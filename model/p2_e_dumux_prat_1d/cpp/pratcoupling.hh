// SPDX-License-Identifier: GPL-3.0-or-later
// Transactional DuMuX-facing wrapper around the frozen P2-D Prat model.
#ifndef PRAT_COUPLING_HH
#define PRAT_COUPLING_HH

#include "../../p2_d_coupling/cpp/pratcrust.hh"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <utility>

namespace prat::coupling {

struct RootDiagnostics
{
    bool converged = false;
    int iterations = 0;
    double residual = std::numeric_limits<double>::quiet_NaN();
    double scale = std::numeric_limits<double>::quiet_NaN();
    double lower = std::numeric_limits<double>::quiet_NaN();
    double upper = std::numeric_limits<double>::quiet_NaN();
};

struct SolverOptions
{
    double undersaturationLimit = 0.02;
    double relativeTolerance = 1.0e-9;
    double absoluteFluxTolerance = 1.0e-14;
    double subrootRelativeTolerance = 1.0e-10;
    double subrootAbsoluteTolerance = 1.0e-14;
    int maximumIterations = 40;
};

struct InterfaceTrial
{
    CrustState baseState;
    CrustState proposedState;
    InterfaceFluxes fluxes;
    RootDiagnostics interfaceRoot;
    double gasFilmResidual = std::numeric_limits<double>::quiet_NaN();
    double topClosureResidual = std::numeric_limits<double>::quiet_NaN();
    bool committed = false;
    bool discarded = false;
};

inline void enforceSubrootConvergence(const CrustState& state,
                                      const CrustParameters& params,
                                      const GasFilm& film,
                                      const InterfaceFluxes& fluxes,
                                      const SolverOptions& options,
                                      double& gasResidual,
                                      double& topResidual);

// Box/CVFE places the boundary degree of freedom at the soil/crust interface.
// Its mass fraction is therefore the interface unknown of the global Newton
// system; introducing a second scalar root would duplicate that unknown.
inline InterfaceTrial evaluateBoundaryState(
    const CrustState& state,
    const CrustParameters& params,
    const GasFilm& film,
    double interfaceMassFraction,
    double dt,
    const SolverOptions& options = {})
{
    if (dt <= 0.0 || !std::isfinite(dt))
        throw DomainExit("trial time step must be positive and finite");
    const double cEq = params.equilibriumMassFraction;
    if (interfaceMassFraction
        < cEq * (1.0 - options.undersaturationLimit))
        throw DomainExit("CVFE boundary state is outside the audited Prat "
                         "undersaturation envelope");

    InterfaceTrial trial;
    trial.baseState = state;
    trial.fluxes = crustResponse(
        state, params, film, interfaceMassFraction, dt);
    enforceSubrootConvergence(state, params, film, trial.fluxes, options,
                              trial.gasFilmResidual, trial.topClosureResidual);
    trial.proposedState = state;
    const StepRecord record = step(
        trial.proposedState, params, trial.fluxes.evaporation, dt,
        interfaceMassFraction);
    const double scale = std::max({std::abs(record.saltFluxSoil),
                                   std::abs(trial.fluxes.saltMass),
                                   options.absoluteFluxTolerance});
    if (std::abs(record.saltFluxSoil - trial.fluxes.saltMass)
        > options.absoluteFluxTolerance + options.relativeTolerance * scale)
        throw DomainExit("CVFE boundary trial replay changed inventory flux");
    // The interface equation is part of the global DuMuX Newton residual.
    trial.interfaceRoot = {true, 0, 0.0, scale,
                           interfaceMassFraction, interfaceMassFraction};
    return trial;
}

inline bool exactlyEqual(const CrustState& a, const CrustState& b)
{
    return a.topCellSize == b.topCellSize
        && a.thicknesses == b.thicknesses
        && a.solid == b.solid
        && a.height == b.height
        && a.time == b.time
        && a.active == b.active
        && a.solidSaltMass == b.solidSaltMass
        && a.dissolvedSalt == b.dissolvedSalt
        && a.water == b.water
        && a.cumulativeInternalReaction == b.cumulativeInternalReaction
        && a.cumulativeTopSkeleton == b.cumulativeTopSkeleton
        && a.seedSalt == b.seedSalt
        && a.seedWater == b.seedWater
        && a.sealed == b.sealed
        && a.disconnected == b.disconnected
        && a.stopped == b.stopped
        && a.stopReason == b.stopReason;
}

inline double topClosureResidual(const CrustState& state,
                                 const CrustParameters& params,
                                 const ProfileSolution& solution)
{
    const std::size_t top = state.thicknesses.size() - 1;
    const double distance = state.height - solution.centres[top];
    const double conductive = params.rho() * solution.porosityProfile[top]
                            * params.diffusivity(solution.porosityProfile[top])
                            * (solution.topSupersaturation - solution.profile[top])
                            / distance;
    const double cTop = solution.topMassFraction;
    const double skeleton = params.crystalDensity * (1.0 - params.depositionPorosity);
    return conductive - cTop * solution.evaporationFlux
         + skeleton * (1.0 - cTop) * solution.growthVelocity
         + cTop * params.vapourDensity * solution.growthVelocity;
}

inline void enforceSubrootConvergence(const CrustState& state,
                                      const CrustParameters& params,
                                      const GasFilm& film,
                                      const InterfaceFluxes& fluxes,
                                      const SolverOptions& options,
                                      double& gasResidual,
                                      double& topResidual)
{
    const ProfileSolution solution = solveProfile(
        state, params, fluxes.evaporation, fluxes.interfaceMassFraction);
    gasResidual = film.fluxFromMassFraction(solution.topMassFraction)
                - fluxes.evaporation;
    const double gasScale = std::max({std::abs(fluxes.evaporation),
                                      std::abs(film.fluxFromMassFraction(
                                          solution.topMassFraction)),
                                      options.subrootAbsoluteTolerance});
    if (std::abs(gasResidual) > options.subrootAbsoluteTolerance
                               + options.subrootRelativeTolerance * gasScale)
        throw DomainExit("gas-film fixed point returned without convergence");

    topResidual = topClosureResidual(state, params, solution);
    const double topScale = std::max({
        std::abs(solution.saltFluxBaseTransport),
        std::abs(solution.evaporationFlux),
        std::abs(params.crystalDensity * (1.0 - params.depositionPorosity)
                 * solution.growthVelocity),
        options.subrootAbsoluteTolerance});
    if (std::abs(topResidual) > options.subrootAbsoluteTolerance
                               + options.subrootRelativeTolerance * topScale)
        throw DomainExit("crust-top closure returned without convergence");
}

inline InterfaceTrial evaluate(
    const CrustState& state,
    const CrustParameters& params,
    const GasFilm& film,
    const std::function<double(double, double)>& soilSaltFlux,
    double dt,
    const SolverOptions& options = {})
{
    if (dt <= 0.0 || !std::isfinite(dt))
        throw DomainExit("trial time step must be positive and finite");
    if (options.maximumIterations <= 0)
        throw DomainExit("interface root iteration limit must be positive");
    if (!(options.undersaturationLimit > 0.0))
        throw DomainExit("interface validity envelope must be positive");

    InterfaceTrial trial;
    trial.baseState = state;

    auto residualAt = [&](double w) {
        const InterfaceFluxes f = crustResponse(state, params, film, w, dt);
        return std::make_pair(f, soilSaltFlux(f.totalMass, w) - f.saltMass);
    };

    const double cEq = params.equilibriumMassFraction;
    const double maxMassFraction =
        MOLALITY_LIMIT * M_NACL / (1.0 + MOLALITY_LIMIT * M_NACL);
    double low = 0.0, high = 0.0, gLow = 0.0, gHigh = 0.0;
    bool bracketed = false;
    for (const double offset : {1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2,
                                options.undersaturationLimit})
    {
        const double o = std::min(offset, options.undersaturationLimit);
        const double lo = cEq * (1.0 - o);
        const double hi = std::min(cEq * (1.0 + o), 0.99 * maxMassFraction);
        try
        {
            const auto a = residualAt(lo);
            const auto b = residualAt(hi);
            if (a.second == 0.0 || b.second == 0.0
                || std::signbit(a.second) != std::signbit(b.second))
            {
                low = lo;
                high = hi;
                gLow = a.second;
                gHigh = b.second;
                bracketed = true;
                break;
            }
        }
        catch (const DomainExit&) {}
    }
    if (!bracketed)
        throw DomainExit("no compatible interface concentration inside the "
                         "declared Prat validity envelope");

    InterfaceFluxes acceptedFluxes;
    double acceptedResidual = std::numeric_limits<double>::quiet_NaN();
    bool converged = false;
    int iterations = 0;
    for (; iterations < options.maximumIterations; ++iterations)
    {
        double w = 0.5 * (low + high);
        if (gHigh != gLow)
        {
            const double candidate = high - gHigh * (high - low) / (gHigh - gLow);
            if (candidate > low && candidate < high)
                w = candidate;
        }
        const auto result = residualAt(w);
        acceptedFluxes = result.first;
        acceptedResidual = result.second;
        const double scale = std::max({std::abs(soilSaltFlux(
                                           acceptedFluxes.totalMass, w)),
                                       std::abs(acceptedFluxes.saltMass),
                                       options.absoluteFluxTolerance});
        if (std::abs(acceptedResidual) <= options.absoluteFluxTolerance
                                        + options.relativeTolerance * scale)
        {
            converged = true;
            ++iterations;
            break;
        }
        if (std::signbit(acceptedResidual) != std::signbit(gLow))
        {
            high = w;
            gHigh = acceptedResidual;
            gLow *= 0.5;
        }
        else
        {
            low = w;
            gLow = acceptedResidual;
            gHigh *= 0.5;
        }
    }

    trial.interfaceRoot = {converged, iterations, acceptedResidual,
        std::max({std::abs(acceptedFluxes.saltMass),
                  std::abs(acceptedResidual + acceptedFluxes.saltMass),
                  options.absoluteFluxTolerance}), low, high};
    if (!converged)
        throw DomainExit("interface flux root did not satisfy the declared "
                         "absolute and relative residual tolerances");
    if (acceptedFluxes.interfaceMassFraction
        < cEq * (1.0 - 0.999 * options.undersaturationLimit))
        throw DomainExit("interface undersaturated outside the audited Prat "
                         "validity envelope");

    enforceSubrootConvergence(state, params, film, acceptedFluxes, options,
                              trial.gasFilmResidual, trial.topClosureResidual);
    trial.fluxes = acceptedFluxes;
    trial.proposedState = state;
    const StepRecord committedRecord = step(
        trial.proposedState, params, acceptedFluxes.evaporation, dt,
        acceptedFluxes.interfaceMassFraction);
    const double inventoryScale = std::max({std::abs(committedRecord.saltFluxSoil),
                                            std::abs(acceptedFluxes.saltMass),
                                            options.absoluteFluxTolerance});
    if (std::abs(committedRecord.saltFluxSoil - acceptedFluxes.saltMass)
        > options.absoluteFluxTolerance + options.relativeTolerance * inventoryScale)
        throw DomainExit("trial replay changed the crust inventory flux");
    return trial;
}

inline void commit(CrustState& state, InterfaceTrial& trial)
{
    if (trial.committed) throw DomainExit("interface trial was already committed");
    if (trial.discarded) throw DomainExit("discarded interface trial cannot be committed");
    if (!trial.interfaceRoot.converged)
        throw DomainExit("unconverged interface trial cannot be committed");
    if (!exactlyEqual(state, trial.baseState))
        throw DomainExit("stale interface trial: persistent crust state changed");
    state = trial.proposedState;
    trial.committed = true;
}

inline void rollback(CrustState& state, InterfaceTrial& trial)
{
    if (trial.committed) throw DomainExit("committed interface trial cannot be rolled back");
    if (!exactlyEqual(state, trial.baseState))
        throw DomainExit("rollback found a modified persistent crust state");
    trial.discarded = true;
}

} // namespace prat::coupling

#endif
