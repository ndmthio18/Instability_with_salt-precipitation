// SPDX-License-Identifier: GPL-3.0-or-later
//
// P2-D: isolated C++ port of the minimum crust functionality required to
// couple the Prat (2026) porous-salt growth law to the P1-R2 conservative
// Richards-NaCl solver.
//
// ============================== STATUS ======================================
// THIS FILE HAS NEVER BEEN COMPILED OR EXECUTED.
//
// No C++ toolchain exists in the environment this was written in: g++,
// clang++, cl, cmake, make and ninja are all absent, and the pinned
// DuMuX/Dune stack recorded in
// model/p1_r2_conservative_richards/provenance/source_manifest.json is not
// installed. It is therefore a transcription of the verified Python model in
// model/p2_a_prat_surface/surface/graded_crust.py and its dependencies, not a
// tested artefact. Treat every number it would produce as unverified until
// verify_vectors.cc passes against reference_vectors.txt on a machine that
// can build it. That check is the gate; see cpp/README.md.
// ============================================================================
//
// Design constraints taken from the task and from the repository rules:
//  * header-only and dependency-free (only the C++17 standard library), so it
//    can be dropped into the DuMuX application without touching its build;
//  * no scipy/numpy analogue and no linear-algebra dependency: the crust
//    transport problem is tridiagonal and is solved with the same Thomas
//    algorithm the Python model uses;
//  * every constitutive closure is Prat (2026) as transcribed in
//    docs/PRAT_PHYSICS_AUDIT.md, including equation (17)
//    a_v = 36 (1 - eps) / d_sub and the tortuosity closure
//    D_s,eff = sqrt(eps) D_s;
//  * the soil-side fluxes are the INVENTORY route, so a coupled solver
//    conserves NaCl exactly; the transport route is also returned so it can
//    be audited independently.
//
// Interface contract with P1-R2
// -----------------------------
// P1-R2's two equations are total liquid mass (conti0EqIdx) and dissolved
// NaCl mass (conti1EqIdx), both per bulk volume, and its Neumann values are
// positive when mass LEAVES the domain. `InterfaceFluxes::totalMass` and
// `::saltMass` are exactly those two numbers for the top boundary, in
// kg/(m^2 s). Nothing else in the solver needs to change.

#ifndef PRAT_CRUST_HH
#define PRAT_CRUST_HH

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace prat {

// ---------------------------------------------------------------- errors ---

//! The model left its declared validity envelope. Never a silent fallback.
class DomainExit : public std::runtime_error
{
public:
    explicit DomainExit(const std::string& what) : std::runtime_error(what) {}
};

// ------------------------------------------------------------ constants ----
// Provenance: DuMuX Constants::R, Air::molarMass, and the frozen Case 36877
// input file (Problem.MolarMassH2O, Problem.MolarMassNaCl).

inline constexpr double R_GAS = 8.314472;      // J mol^-1 K^-1
inline constexpr double M_H2O = 0.01801528;    // kg mol^-1
inline constexpr double M_NACL = 0.05844;      // kg mol^-1
inline constexpr double M_AIR = 0.02896;       // kg mol^-1
//! Largest molality the NaCl activity correlation is declared valid at.
inline constexpr double MOLALITY_LIMIT = 8.0;

// --------------------------------------------------------------- thermo ----

//! Pure-water saturation pressure [Pa]; IAPWS-IF97 Region 4, as in DuMuX
//! dumux/material/components/iapws/region4.hh.
inline double saturationPressure(double temperature)
{
    if (!std::isfinite(temperature) || temperature <= 0.0)
        throw DomainExit("saturationPressure requires a positive finite temperature");
    static constexpr double n[10] = {
        0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
        0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
        -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
        0.65017534844798e3};
    const double sigma = temperature + n[8] / (temperature - n[9]);
    const double a = (sigma + n[0]) * sigma + n[1];
    const double b = (n[2] * sigma + n[3]) * sigma + n[4];
    const double c = (n[5] * sigma + n[6]) * sigma + n[7];
    double tmp = 2.0 * c / (std::sqrt(b * b - 4.0 * a * c) - b);
    tmp *= tmp;
    tmp *= tmp;
    return 1.0e6 * tmp;
}

//! Binary H2O-air gas diffusion coefficient [m^2/s]; DuMuX H2O_Air::gasDiffCoeff.
inline double vapourDiffusionCoefficient(double temperature, double pressure)
{
    if (!std::isfinite(temperature) || temperature <= 0.0)
        throw DomainExit("vapour diffusion requires a positive finite temperature");
    if (!std::isfinite(pressure) || pressure <= 0.0)
        throw DomainExit("vapour diffusion requires a positive finite pressure");
    return 2.13e-5 * (1.0e5 / pressure) * std::pow(temperature / 273.15, 1.8);
}

//! Mass fraction of water in an ideal H2O/air gas at vapour mole fraction y.
inline double gasWaterMassFraction(double yWater)
{
    if (!std::isfinite(yWater) || yWater < 0.0 || yWater > 1.0)
        throw DomainExit("gas water mole fraction must lie in [0, 1]");
    return M_H2O * yWater / (M_H2O * yWater + M_AIR * (1.0 - yWater));
}

inline double idealGasDensity(double temperature, double pressure, double yWater)
{
    const double molarMass = yWater * M_H2O + (1.0 - yWater) * M_AIR;
    return pressure * molarMass / (R_GAS * temperature);
}

//! w_s from x_s; identical to the P1-R2 Problem::massFraction closure.
inline double massFractionFromMoleFraction(double xSalt)
{
    const double mixture = (1.0 - xSalt) * M_H2O + xSalt * M_NACL;
    return xSalt * M_NACL / mixture;
}

inline double moleFractionFromMassFraction(double wSalt)
{
    if (wSalt < 0.0 || wSalt >= 1.0)
        throw DomainExit("salt mass fraction must lie in [0, 1)");
    const double nSalt = wSalt / M_NACL;
    const double nWater = (1.0 - wSalt) / M_H2O;
    return nSalt / (nSalt + nWater);
}

inline double molalityFromMassFraction(double wSalt)
{
    if (!std::isfinite(wSalt) || wSalt < 0.0 || wSalt >= 1.0)
        throw DomainExit("salt mass fraction must lie in [0, 1)");
    return wSalt / (M_NACL * (1.0 - wSalt));
}

//! The frozen P1 Richards-salt liquid density law, kept bit for bit.
inline double brineDensity(double wSalt,
                           double pureWaterDensity = 998.203,
                           double linear = 0.670804,
                           double quadratic = 0.373854)
{
    return pureWaterDensity + 1000.0 * wSalt * (linear + quadratic * wSalt);
}

// ------------------------------------------------------- water activity ----
// Pure-NaCl Pitzer osmotic coefficient, Rose et al. (2008) Eqs. A15-A16 with
// Pitzer and Mayorga (1973) parameters at 298.15 K, reused unchanged from the
// heterogeneity project's water_activity.hh. Its own caveat is inherited: use
// at 293.15 K is a stated near-room-temperature approximation, not a fit.

inline double osmoticCoefficient(double molality)
{
    if (!std::isfinite(molality) || molality < 0.0 || molality > MOLALITY_LIMIT)
        throw DomainExit("NaCl water activity requires molality in [0, 8]");
    const double r = std::sqrt(molality);
    return 1.0 - 0.3915 * r / (1.0 + 1.2 * r)
           + molality * (0.0765 + 0.2664 * std::exp(-2.0 * r))
           + molality * molality * 0.00127;
}

inline double waterActivity(double molality)
{
    return std::exp(-2.0 * molality * M_H2O * osmoticCoefficient(molality));
}

inline double waterActivityFromMassFraction(double wSalt)
{
    return waterActivity(molalityFromMassFraction(wSalt));
}

// ------------------------------------------------------------- gas film ----

//! Stagnant gas-film evaporation. No calibration, no normalisation of local
//! weights: the total rate is free to change with salinity.
struct GasFilm
{
    double temperature = 293.15;          //!< [K]
    double relativeHumidity = 0.4;        //!< [-]
    double thickness = 0.0083;            //!< delta [m]
    double atmosphericPressure = 1.0e5;   //!< [Pa]
    //! <= 0 means "use the documented DuMuX H2O-air correlation".
    double vapourDiffusion = -1.0;

    double diffusivity() const
    {
        return vapourDiffusion > 0.0
                   ? vapourDiffusion
                   : vapourDiffusionCoefficient(temperature, atmosphericPressure);
    }

    double saturationPressureValue() const { return saturationPressure(temperature); }

    double ambientVapourMoleFraction() const
    {
        return relativeHumidity * saturationPressureValue() / atmosphericPressure;
    }

    double referenceGasDensity() const
    {
        return idealGasDensity(temperature, atmosphericPressure,
                               ambientVapourMoleFraction());
    }

    //! rho_ref D_v / delta [kg m^-2 s^-1 per unit mass-fraction difference].
    double transferCoefficient() const
    {
        return referenceGasDensity() * diffusivity() / thickness;
    }

    //! Water-component mass flux [kg m^-2 s^-1]; SIGNED, negative is condensation.
    double fluxFromActivity(double activity) const
    {
        if (!std::isfinite(activity) || activity < 0.0 || activity > 1.0)
            throw DomainExit("water activity must lie in [0, 1]");
        const double ys = activity * saturationPressureValue() / atmosphericPressure;
        if (!(ys >= 0.0 && ys < 1.0))
            throw DomainExit("invalid film interfacial vapour state");
        const double driving = gasWaterMassFraction(ys)
                               - gasWaterMassFraction(ambientVapourMoleFraction());
        return transferCoefficient() * driving;
    }

    double fluxFromMassFraction(double wSalt) const
    {
        return fluxFromActivity(waterActivityFromMassFraction(wSalt));
    }
};

// ------------------------------------------------------ crust parameters ----

enum class AreaClosure { Constant, PratEquation17 };
enum class DiffusivityClosure { Constant, SqrtPorosity };

struct CrustParameters
{
    double equilibriumMassFraction = 0.0;   //!< C_eq
    double depositionPorosity = 0.43;       //!< eps_0, Prat's base value
    double crystalDensity = 2165.0;         //!< rho_cr [kg/m3]
    double specificArea = 0.0;              //!< a_v0 [1/m], used by Constant
    double reactionCoefficient = 2.3e-3;    //!< k_r [m/s], Prat's ref. 23
    double effectiveSaltDiffusivity = 0.0;  //!< D_ref [m2/s], used by Constant
    double beta = 2.0;                      //!< Prat appendix cubic construction
    double solutionDensity = 0.0;           //!< rho; <= 0 means the brine law
    double vapourDensity = 0.0;             //!< rho_v [kg/m3]
    double crustSaturation = 1.0;           //!< S_c
    double cloggingPorosity = 0.01;
    double connectivityPorosity = 0.02;
    double substrateSize = 25.0e-6;         //!< d_sub [m], Prat equation (17)
    double molecularDiffusivity = 1.3e-9;   //!< D_s [m2/s]
    AreaClosure areaClosure = AreaClosure::PratEquation17;
    DiffusivityClosure diffusivityClosure = DiffusivityClosure::SqrtPorosity;
    //! Reinstates the C ~ C_eq top linearisation of Prat equation (8).
    //! Only for reproducing equation (28); production runs keep it false.
    bool linearisedTopClosure = false;

    double rho() const
    {
        return solutionDensity > 0.0 ? solutionDensity
                                     : brineDensity(equilibriumMassFraction);
    }

    //! Prat equation (17): a_v = 36 (1 - eps) / d_sub.
    double area(double porosity) const
    {
        const double eps = std::max(porosity, 0.0);
        if (areaClosure == AreaClosure::PratEquation17)
            return 36.0 * (1.0 - eps) / substrateSize;
        return specificArea;
    }

    //! D_s,eff = D_s / tau with tau ~ eps^(-1/2).
    double diffusivity(double porosity) const
    {
        const double eps = std::max(porosity, 0.0);
        if (diffusivityClosure == DiffusivityClosure::SqrtPorosity)
            return std::sqrt(eps) * molecularDiffusivity;
        return effectiveSaltDiffusivity;
    }

    //! delta_r = sqrt( eps D / ((1 - C_eq) a_v k_r) ) = h / sqrt(Da).
    double reactionLength(double porosity) const
    {
        const double eps = porosity > 0.0 ? porosity : depositionPorosity;
        return std::sqrt(eps * diffusivity(eps)
                         / ((1.0 - equilibriumMassFraction) * area(eps)
                            * reactionCoefficient));
    }

    double damkohlerStar(double porosity) const
    {
        const double eps = porosity > 0.0 ? porosity : depositionPorosity;
        return beta * beta * (1.0 - eps) * (1.0 - eps)
               * (1.0 - equilibriumMassFraction) * reactionCoefficient
               / (eps * diffusivity(eps) * area(eps));
    }

    void check() const
    {
        if (!(equilibriumMassFraction > 0.0 && equilibriumMassFraction < 1.0))
            throw DomainExit("C_eq must lie strictly in (0, 1)");
        if (!(depositionPorosity > 0.0 && depositionPorosity < 1.0))
            throw DomainExit("deposition porosity must lie strictly in (0, 1)");
        if (!(crystalDensity > 0.0 && reactionCoefficient > 0.0 && beta > 0.0))
            throw DomainExit("crystal density, k_r and beta must be positive");
        if (!(cloggingPorosity >= 0.0 && cloggingPorosity < depositionPorosity))
            throw DomainExit("clogging porosity must lie in [0, eps_0)");
    }
};

// ----------------------------------------------------------- crust state ----

struct CrustState
{
    double topCellSize = 0.0;
    std::vector<double> thicknesses;
    std::vector<double> solid;     //!< solid volume per projected area [m]
    double height = 0.0;
    double time = 0.0;
    bool active = false;
    double solidSaltMass = 0.0;    //!< [kg/m2]
    double dissolvedSalt = 0.0;    //!< [kg/m2]
    double water = 0.0;            //!< [kg/m2]
    double cumulativeInternalReaction = 0.0;
    double cumulativeTopSkeleton = 0.0;
    double seedSalt = 0.0;
    double seedWater = 0.0;
    bool sealed = false;
    bool disconnected = false;
    bool stopped = false;
    std::string stopReason;

    double porosity(std::size_t index) const
    {
        const double t = thicknesses[index];
        return t <= 0.0 ? 0.0 : 1.0 - solid[index] / t;
    }

    std::vector<double> centres() const
    {
        std::vector<double> out(thicknesses.size());
        double lower = 0.0;
        for (std::size_t i = 0; i < thicknesses.size(); ++i)
        {
            out[i] = lower + 0.5 * thicknesses[i];
            lower += thicknesses[i];
        }
        return out;
    }

    double solidVolume() const
    {
        double sum = 0.0;
        for (const double value : solid) sum += value;
        return sum;
    }

    double minimumPorosity() const
    {
        double value = 1.0;
        for (std::size_t i = 0; i < thicknesses.size(); ++i)
            value = std::min(value, porosity(i));
        return value;
    }
};

// ------------------------------------------------------- tridiagonal solve --

//! Thomas algorithm, matching the Python `_thomas` term for term.
inline std::vector<double> thomas(const std::vector<double>& lower,
                                  const std::vector<double>& diag,
                                  const std::vector<double>& upper,
                                  const std::vector<double>& rhs)
{
    const std::size_t n = diag.size();
    std::vector<double> c(n, 0.0), d(n, 0.0), x(n, 0.0);
    double beta = diag[0];
    if (beta == 0.0) throw DomainExit("singular crust transport system");
    c[0] = upper[0] / beta;
    d[0] = rhs[0] / beta;
    for (std::size_t i = 1; i < n; ++i)
    {
        beta = diag[i] - lower[i] * c[i - 1];
        if (beta == 0.0) throw DomainExit("singular crust transport system");
        c[i] = (i < n - 1) ? upper[i] / beta : 0.0;
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / beta;
    }
    x[n - 1] = d[n - 1];
    for (std::size_t i = n - 1; i-- > 0;)
        x[i] = d[i] - c[i] * x[i + 1];
    return x;
}

// ----------------------------------------------------------- graded grid ----

//! Cell thicknesses, bottom to top, tiling [0, height] exactly.
inline std::vector<double> gradedSeedGrid(double height, double topCell,
                                          double window = 0.0, double ratio = 1.0,
                                          double maximumCell = 0.0)
{
    if (height <= 0.0) throw DomainExit("seed height must be positive");
    if (topCell <= 0.0) throw DomainExit("top cell size must be positive");
    if (ratio < 1.0) throw DomainExit("grading ratio must be at least 1");
    topCell = std::min(topCell, height);
    const double cap = maximumCell > 0.0 ? maximumCell : height;
    std::vector<double> fromTop;
    double covered = 0.0, size = topCell;
    while (covered < height - 1.0e-15 * height)
    {
        if (covered >= window) size = std::min(size * ratio, cap);
        const double take = std::min(size, height - covered);
        fromTop.push_back(take);
        covered += take;
        if (fromTop.size() > 10000000u)
            throw DomainExit("graded grid exceeded the cell cap");
    }
    // Absorb a sliver bottom cell into its neighbour. The pop must happen
    // first; indexing after the pop would merge into the wrong cell.
    if (fromTop.size() > 1 && fromTop.back() < 0.5 * fromTop[fromTop.size() - 2])
    {
        const double sliver = fromTop.back();
        fromTop.pop_back();
        fromTop.back() += sliver;
    }
    std::vector<double> out(fromTop.rbegin(), fromTop.rend());
    return out;
}

// ------------------------------------------------------ profile solution ----

struct ProfileSolution
{
    double height = 0.0;
    double evaporationFlux = 0.0;
    double baseMassFraction = 0.0;
    double topSupersaturation = 0.0;
    double topMassFraction = 0.0;
    double growthVelocity = 0.0;
    double internalReaction = 0.0;
    double topSkeletonRate = 0.0;
    double dissolvedSalt = 0.0;
    double water = 0.0;
    double solutionMassFluxTop = 0.0;
    double solutionMassFluxBase = 0.0;
    double baseGradient = 0.0;
    double saltFluxBaseTransport = 0.0;   //!< independent transport route
    double minimumPorosity = 0.0;
    double topPorosity = 0.0;
    double reactionLength = 0.0;
    bool growthArrested = false;
    std::vector<double> profile;          //!< C*(z) at cell centres
    std::vector<double> porosityProfile;
    std::vector<double> centres;
    std::vector<double> thicknesses;
};

namespace detail {

//! Bracketed bisection for the scalar top closure; no Newton tuning needed.
template<class Residual>
double solveTopSupersaturation(Residual residual, double cEq, double flux)
{
    if (flux <= 0.0) return 0.0;
    double low = 0.0;
    if (residual(low) >= 0.0) return low;
    double high = std::max(1.0e-12, 1.0e-3 * cEq);
    bool bracketed = false;
    for (int i = 0; i < 200; ++i)
    {
        if (residual(high) >= 0.0) { bracketed = true; break; }
        high *= 2.0;
        if (cEq + high >= 1.0)
            throw DomainExit("the crust-top closure has no root below C_top = 1");
    }
    if (!bracketed) throw DomainExit("failed to bracket the crust-top closure");
    for (int i = 0; i < 200; ++i)
    {
        const double mid = 0.5 * (low + high);
        if (residual(mid) < 0.0) low = mid; else high = mid;
        if (high - low <= 1.0e-16 * std::max(high, 1.0e-30)) break;
    }
    return 0.5 * (low + high);
}

} // namespace detail

//! Quasi-steady crust transport with an exact top closure.
//!
//! `baseMassFraction < 0` means "use C_eq", i.e. Prat's halite-equilibrated
//! base. Any other value is the soil/crust interface concentration, which the
//! coupled problem solves for rather than prescribes.
inline ProfileSolution solveProfile(const CrustState& state,
                                    const CrustParameters& params,
                                    double evaporationFlux,
                                    double baseMassFraction = -1.0)
{
    if (!state.active || state.height <= 0.0)
        throw DomainExit("the crust transport problem requires an active crust");
    const std::size_t n = state.thicknesses.size();
    if (n == 0) throw DomainExit("no occupied crust cells");

    const double rho = params.rho();
    const double sat = params.crustSaturation;
    const double cEq = params.equilibriumMassFraction;
    const double one = 1.0 - cEq;

    const std::vector<double> z = state.centres();
    std::vector<double> eps(n), cond(n), lam(n);
    for (std::size_t k = 0; k < n; ++k)
    {
        eps[k] = state.porosity(k);
        cond[k] = eps[k] * params.diffusivity(eps[k]);
        lam[k] = one * params.area(eps[k]) * params.reactionCoefficient;
    }

    std::vector<double> lower(n, 0.0), diag(n, 0.0), upper(n, 0.0), rhs(n, 0.0);
    for (std::size_t k = 0; k < n; ++k)
    {
        if (k == 0)
        {
            const double distance = z[0];
            const double coefficient = distance > 0.0 ? cond[0] / distance : 0.0;
            diag[k] -= coefficient;                       // phi(0) = 0
        }
        else
        {
            const double distance = z[k] - z[k - 1];
            const double coefficient = 0.5 * (cond[k] + cond[k - 1]) / distance;
            diag[k] -= coefficient;
            lower[k] += coefficient;
        }
        if (k == n - 1)
        {
            const double distance = state.height - z[k];
            const double coefficient = distance > 0.0 ? cond[k] / distance : 0.0;
            diag[k] -= coefficient;
            rhs[k] -= coefficient;                        // phi(h) = 1
        }
        else
        {
            const double distance = z[k + 1] - z[k];
            const double coefficient = 0.5 * (cond[k] + cond[k + 1]) / distance;
            diag[k] -= coefficient;
            upper[k] += coefficient;
        }
        diag[k] -= lam[k] * state.thicknesses[k];
    }

    const std::vector<double> phi = thomas(lower, diag, upper, rhs);

    const double cBase = baseMassFraction < 0.0 ? cEq : baseMassFraction;
    const double uBase = cBase - cEq;
    std::vector<double> psi(n, 0.0);
    if (uBase != 0.0)
    {
        std::vector<double> rhsPsi(n, 0.0);
        const double distance0 = z[0];
        rhsPsi[0] -= distance0 > 0.0 ? cond[0] / distance0 : 0.0;
        psi = thomas(lower, diag, upper, rhsPsi);
    }

    const double topDistance = state.height - z[n - 1];
    if (topDistance <= 0.0) throw DomainExit("degenerate top cell in the crust grid");
    const double epsTop = eps[n - 1];
    const double aCoefficient = rho * cond[n - 1] * (1.0 - phi[n - 1]) / topDistance;
    const double bCoefficient = rho * cond[n - 1] * (0.0 - psi[n - 1]) / topDistance;

    const double kKinetic =
        params.beta * params.reactionCoefficient * rho / params.crystalDensity;
    const double skeleton = params.crystalDensity * (1.0 - params.depositionPorosity);

    auto residual = [&](double s) {
        const double vs = kKinetic * s;
        const double cTop = params.linearisedTopClosure ? cEq : cEq + s;
        return aCoefficient * s + bCoefficient * uBase - cTop * evaporationFlux
               + skeleton * (1.0 - cTop) * vs + cTop * params.vapourDensity * vs;
    };

    ProfileSolution out;
    out.growthArrested = residual(0.0) >= 0.0;
    const double sTop = detail::solveTopSupersaturation(residual, cEq, evaporationFlux);
    const double vs = kKinetic * sTop;

    out.profile.resize(n);
    double internal = 0.0, dissolved = 0.0, water = 0.0;
    for (std::size_t k = 0; k < n; ++k)
    {
        out.profile[k] = sTop * phi[k] + uBase * psi[k];
        internal += params.area(eps[k]) * rho * params.reactionCoefficient
                    * out.profile[k] * state.thicknesses[k];
        const double pore = eps[k] * sat * rho * state.thicknesses[k];
        dissolved += pore * (cEq + out.profile[k]);
        water += pore * (1.0 - cEq - out.profile[k]);
    }

    // Moving-interface solution flux. The new pore volume is created at the
    // DEPOSITION porosity, and precipitation frees pore volume that refills
    // with solution, so d(rho U)/dz = -R (1 - rho S / rho_cr).
    const double rhoUTop =
        (rho * params.depositionPorosity * sat + skeleton - params.vapourDensity) * vs
        + evaporationFlux;
    const double release = 1.0 - rho * sat / params.crystalDensity;
    const double rhoUBase = rhoUTop + internal * release;

    const double baseGradient = z[0] > 0.0 ? (out.profile[0] - uBase) / z[0] : 0.0;
    const double diffusiveBase =
        -rho * eps[0] * sat * params.diffusivity(eps[0]) * baseGradient;

    out.height = state.height;
    out.evaporationFlux = evaporationFlux;
    out.baseMassFraction = cBase;
    out.topSupersaturation = sTop;
    out.topMassFraction = cEq + sTop;
    out.growthVelocity = vs;
    out.internalReaction = internal;
    out.topSkeletonRate = skeleton * vs;
    out.dissolvedSalt = dissolved;
    out.water = water;
    out.solutionMassFluxTop = rhoUTop;
    out.solutionMassFluxBase = rhoUBase;
    out.baseGradient = baseGradient;
    out.saltFluxBaseTransport = cBase * rhoUBase + diffusiveBase;
    out.minimumPorosity = *std::min_element(eps.begin(), eps.end());
    out.topPorosity = epsTop;
    out.reactionLength = params.reactionLength(epsTop);
    out.porosityProfile = eps;
    out.centres = z;
    out.thicknesses = state.thicknesses;
    return out;
}

// ------------------------------------------------------------- activation ---

//! Create the seed crust at the declared deposition porosity.
inline void activate(CrustState& state, const CrustParameters& params,
                     double seedHeight, double evaporationFlux,
                     double window = 0.0, double ratio = 1.0,
                     double maximumCell = 0.0, double baseMassFraction = -1.0)
{
    if (state.active) return;
    if (seedHeight <= 0.0)
        throw DomainExit("seed height must be positive; h = 0 is not a crust");
    state.thicknesses =
        gradedSeedGrid(seedHeight, state.topCellSize, window, ratio, maximumCell);
    state.solid.resize(state.thicknesses.size());
    state.height = 0.0;
    for (std::size_t i = 0; i < state.thicknesses.size(); ++i)
    {
        state.solid[i] = (1.0 - params.depositionPorosity) * state.thicknesses[i];
        state.height += state.thicknesses[i];
    }
    state.active = true;
    const ProfileSolution solution =
        solveProfile(state, params, evaporationFlux, baseMassFraction);
    state.solidSaltMass = params.crystalDensity * state.solidVolume();
    state.dissolvedSalt = solution.dissolvedSalt;
    state.water = solution.water;
    state.seedSalt = state.solidSaltMass + state.dissolvedSalt;
    state.seedWater = state.water;
}

// ------------------------------------------------------------ time step ----

struct StepRecord
{
    double time = 0.0;
    double height = 0.0;
    double growthVelocity = 0.0;
    double topSupersaturation = 0.0;
    double minimumPorosity = 0.0;
    double internalReaction = 0.0;
    //! INVENTORY route: d(M_solid + M_dissolved)/dt. This is what the soil
    //! must be given for the coupled pair to conserve NaCl exactly.
    double saltFluxSoil = 0.0;
    //! j - rho_v v_s + dM_w/dt.
    double waterFluxSoil = 0.0;
    //! Independent transport-route reconstruction of the same salt flux.
    double saltFluxBaseTransport = 0.0;
    bool internalReactionClamped = false;
    bool growthArrested = false;
    std::size_t cells = 0;
};

//! Advance porosity (Prat equation 2, pointwise in the material coordinate z),
//! the moving top and the inventories by dt.
inline StepRecord step(CrustState& state, const CrustParameters& params,
                       double evaporationFlux, double dt,
                       double baseMassFraction = -1.0,
                       std::size_t maximumCells = 2000000u)
{
    if (dt <= 0.0 || !std::isfinite(dt))
        throw DomainExit("time step must be positive and finite");
    if (!state.active) throw DomainExit("cannot step an inactive crust");
    if (state.stopped)
        throw DomainExit("crust stopped: " + state.stopReason
                         + "; the saturated capillary-supply closure is "
                           "disabled rather than extrapolated");

    const ProfileSolution before =
        solveProfile(state, params, evaporationFlux, baseMassFraction);
    const double solidBefore = params.crystalDensity * state.solidVolume();
    const double dissolvedBefore = before.dissolvedSalt;
    const double waterBefore = before.water;
    const double rho = params.rho();

    StepRecord record;
    record.internalReactionClamped = false;

    // 1. internal reaction, cell by cell
    double internalVolume = 0.0;
    for (std::size_t k = 0; k < state.thicknesses.size(); ++k)
    {
        const double rate = params.area(before.porosityProfile[k]) * rho
                            * params.reactionCoefficient * before.profile[k];
        double added = rate * state.thicknesses[k] * dt / params.crystalDensity;
        const double available = state.thicknesses[k] - state.solid[k];
        if (added > available)
        {
            added = std::max(0.0, available);
            record.internalReactionClamped = true;
        }
        state.solid[k] += added;
        internalVolume += added;
    }

    // 2. moving top: new material at the declared deposition porosity
    const double growth = before.growthVelocity * dt;
    double remaining = growth;
    const double skeletonFraction = 1.0 - params.depositionPorosity;
    const double room = state.topCellSize - state.thicknesses.back();
    if (room > 0.0 && remaining > 0.0)
    {
        const double take = std::min(room, remaining);
        state.thicknesses.back() += take;
        state.solid.back() += skeletonFraction * take;
        remaining -= take;
    }
    while (remaining > 0.0)
    {
        if (state.thicknesses.size() >= maximumCells)
            throw DomainExit("crust grid exceeded the declared cell cap");
        const double take = std::min(state.topCellSize, remaining);
        state.thicknesses.push_back(take);
        state.solid.push_back(skeletonFraction * take);
        remaining -= take;
    }
    state.height += growth;
    state.time += dt;

    const ProfileSolution after =
        solveProfile(state, params, evaporationFlux, baseMassFraction);
    const double solidAfter = params.crystalDensity * state.solidVolume();
    state.solidSaltMass = solidAfter;
    state.dissolvedSalt = after.dissolvedSalt;
    state.water = after.water;
    state.cumulativeInternalReaction += internalVolume * params.crystalDensity;
    state.cumulativeTopSkeleton += skeletonFraction * params.crystalDensity * growth;

    record.saltFluxSoil =
        ((solidAfter - solidBefore) + (after.dissolvedSalt - dissolvedBefore)) / dt;
    record.waterFluxSoil = evaporationFlux
                           - params.vapourDensity * before.growthVelocity
                           + (after.water - waterBefore) / dt;
    record.saltFluxBaseTransport = before.saltFluxBaseTransport;
    record.time = state.time;
    record.height = state.height;
    record.growthVelocity = before.growthVelocity;
    record.topSupersaturation = before.topSupersaturation;
    record.minimumPorosity = after.minimumPorosity;
    record.internalReaction = before.internalReaction;
    record.growthArrested = before.growthArrested;
    record.cells = state.thicknesses.size();

    // Connectivity is the PRIMARY stop: beyond it the saturated capillary
    // supply the reduction assumes is void. Clogging is recorded too.
    if (record.minimumPorosity <= params.connectivityPorosity && !state.disconnected)
    {
        state.disconnected = true;
        state.stopped = true;
        state.stopReason = "hydraulic_connectivity_loss";
    }
    if (record.minimumPorosity <= params.cloggingPorosity && !state.sealed)
    {
        state.sealed = true;
        state.stopped = true;
        state.stopReason = "geometric_clogging";
    }
    return record;
}

// -------------------------------------------------- DuMuX-facing coupling ---

//! The two Neumann values P1-R2 needs at its top boundary, plus diagnostics.
//! Sign convention: POSITIVE means mass LEAVES the soil, matching
//! `Problem::neumann` in the P1-R2 application.
struct InterfaceFluxes
{
    double totalMass = 0.0;        //!< conti0EqIdx [kg/(m2 s)]
    double saltMass = 0.0;         //!< conti1EqIdx [kg/(m2 s)]
    double waterMass = 0.0;
    double evaporation = 0.0;
    double growthVelocity = 0.0;
    double interfaceMassFraction = 0.0;
    double crustTopMassFraction = 0.0;
    double saltFluxTransportRoute = 0.0;   //!< independent audit of saltMass
    double minimumPorosity = 0.0;
    double height = 0.0;
};

//! Close the gas-film fixed point: j = film(C_top(j)).
//!
//! The mismatch is decreasing in j, is positive at j = 0 and negative at the
//! free-water flux, so the bracket is guaranteed. A crust-top concentration
//! beyond the activity correlation's range is reported as "this flux is too
//! large" rather than extrapolated.
inline double gasFilmFixedPoint(const CrustState& state,
                                const CrustParameters& params,
                                const GasFilm& film,
                                double interfaceMassFraction)
{
    auto mismatch = [&](double flux) {
        const ProfileSolution s =
            solveProfile(state, params, flux, interfaceMassFraction);
        const double wTop = s.topMassFraction;
        if (wTop >= 1.0 || molalityFromMassFraction(wTop) > MOLALITY_LIMIT)
            return flux > 0.0 ? -flux : -1.0;
        return film.fluxFromMassFraction(wTop) - flux;
    };
    double low = 0.0;
    double high = film.fluxFromMassFraction(0.0);
    for (int i = 0; i < 200; ++i)
    {
        const double mid = 0.5 * (low + high);
        if (mismatch(mid) > 0.0) low = mid; else high = mid;
        if (high - low <= 1.0e-14 * std::max(high, 1.0e-30)) break;
    }
    const double flux = 0.5 * (low + high);
    const ProfileSolution s = solveProfile(state, params, flux, interfaceMassFraction);
    if (molalityFromMassFraction(s.topMassFraction) > MOLALITY_LIMIT)
        throw DomainExit("the crust-top molality at the gas-film fixed point "
                         "exceeds the declared Pitzer limit");
    return flux;
}

//! Crust-side component fluxes for a trial interface concentration.
//!
//! `state` is advanced on a COPY, so this can be called repeatedly inside the
//! interface root solve without committing anything.
inline InterfaceFluxes crustResponse(const CrustState& state,
                                     const CrustParameters& params,
                                     const GasFilm& film,
                                     double interfaceMassFraction,
                                     double dt)
{
    const double flux = gasFilmFixedPoint(state, params, film, interfaceMassFraction);
    const ProfileSolution solution =
        solveProfile(state, params, flux, interfaceMassFraction);
    CrustState trial = state;
    const StepRecord record =
        step(trial, params, flux, dt, interfaceMassFraction);

    InterfaceFluxes out;
    out.saltMass = record.saltFluxSoil;
    out.waterMass = record.waterFluxSoil;
    out.totalMass = out.saltMass + out.waterMass;
    out.evaporation = flux;
    out.growthVelocity = solution.growthVelocity;
    out.interfaceMassFraction = interfaceMassFraction;
    out.crustTopMassFraction = solution.topMassFraction;
    out.saltFluxTransportRoute = solution.saltFluxBaseTransport;
    out.minimumPorosity = record.minimumPorosity;
    out.height = record.height;
    return out;
}

//! Solve for the interface concentration that makes the soil-side and
//! crust-side NaCl fluxes agree, then commit the crust step.
//!
//! `soilSaltFlux(totalMassFlux, interfaceMassFraction)` must return the soil's
//! own reconstruction of the top-face salt flux, using the same upwinding and
//! mass-averaged Fickian closure as its interior faces. Neither the interface
//! concentration nor the salt flux is prescribed: the first is solved for and
//! the second follows.
//!
//! `undersaturationLimit` is the declared validity envelope: the Prat
//! reduction assumes a halite-equilibrated base, so a root outside
//! C_eq (1 +/- limit) is refused rather than extrapolated.
inline InterfaceFluxes solveInterface(
    CrustState& state, const CrustParameters& params, const GasFilm& film,
    const std::function<double(double, double)>& soilSaltFlux, double dt,
    double undersaturationLimit = 0.02, double tolerance = 1.0e-9,
    int maximumIterations = 40)
{
    const double cEq = params.equilibriumMassFraction;
    auto residualAt = [&](double w) {
        const InterfaceFluxes f = crustResponse(state, params, film, w, dt);
        return std::make_pair(f, soilSaltFlux(f.totalMass, w) - f.saltMass);
    };

    const double maxMassFraction =
        MOLALITY_LIMIT * M_NACL / (1.0 + MOLALITY_LIMIT * M_NACL);
    double low = 0.0, high = 0.0, gLow = 0.0, gHigh = 0.0;
    bool bracketed = false;
    for (const double offset : {1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, undersaturationLimit})
    {
        const double o = std::min(offset, undersaturationLimit);
        const double lo = cEq * (1.0 - o);
        const double hi = std::min(cEq * (1.0 + o), 0.99 * maxMassFraction);
        try
        {
            const auto a = residualAt(lo);
            const auto b = residualAt(hi);
            if (a.second * b.second <= 0.0)
            {
                low = lo; high = hi; gLow = a.second; gHigh = b.second;
                bracketed = true;
                break;
            }
        }
        catch (const DomainExit&) { continue; }
    }
    if (!bracketed)
        throw DomainExit("no compatible interface concentration inside the "
                         "declared envelope; the reduced Prat closure is not "
                         "defensible here");

    // Illinois: keeps the bracket, converges superlinearly.
    double w = 0.5 * (low + high);
    InterfaceFluxes fluxes = crustResponse(state, params, film, w, dt);
    for (int i = 0; i < maximumIterations; ++i)
    {
        if (gHigh == gLow) break;
        double mid = high - gHigh * (high - low) / (gHigh - gLow);
        if (!(mid > std::min(low, high) && mid < std::max(low, high)))
            mid = 0.5 * (low + high);
        const auto result = residualAt(mid);
        fluxes = result.first;
        w = mid;
        const double g = result.second;
        if (std::abs(g) <= tolerance * std::max(std::abs(fluxes.saltMass), 1.0e-300))
            break;
        if (g * gLow < 0.0) { high = mid; gHigh = g; gLow *= 0.5; }
        else { low = mid; gLow = g; gHigh *= 0.5; }
        if (std::abs(high - low) <= 1.0e-14 * std::max(std::abs(high), 1.0e-30)) break;
    }

    if (w < cEq * (1.0 - 0.999 * undersaturationLimit))
        throw DomainExit("interface undersaturated: the Prat reduction's "
                         "halite-equilibrated base is not defensible");

    // Commit the crust with the converged interface value.
    step(state, params, fluxes.evaporation, dt, w);
    return fluxes;
}

} // namespace prat

#endif // PRAT_CRUST_HH
