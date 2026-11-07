// SPDX-License-Identifier: GPL-3.0-or-later
//
// Standalone checker for the P2-D C++ crust port.
//
// It recomputes every entry of reference_vectors.txt from pratcrust.hh and
// reports the largest relative difference. The vectors were produced by the
// Python model that P2-A, P2-B and P2-C verified, so a pass here is what
// transfers that verification to the C++ module. Until it is run on a machine
// with a compiler, the header is UNVERIFIED.
//
//   g++ -std=c++17 -O2 -o verify_vectors verify_vectors.cc
//   ./verify_vectors reference_vectors.txt
//
// Exit status 0 means every vector matched within the tolerance below.

#include "pratcrust.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

// Declared vector configuration, mirrored from
// verification/generate_cpp_vectors.py. If these drift apart the comparison
// is meaningless, so they are kept adjacent and both are documented.
constexpr double SEED_HEIGHT = 1.0e-5;
constexpr double CELLS_PER_DELTA = 8.0;
constexpr double WINDOW_LENGTHS = 40.0;
constexpr double RATIO = 1.25;
constexpr double MAXIMUM_CELL_LENGTHS = 10.0;
constexpr double DT = 5.0;
constexpr int STEPS = 10;

//! Relative tolerance. Both sides evaluate the same closed-form expressions
//! in double precision but not necessarily in the same operation order, and
//! both contain bisections with a 1e-16 relative stopping rule, so exact
//! equality is not expected. 1e-11 is far tighter than any physical claim.
constexpr double TOLERANCE = 1.0e-11;

prat::CrustParameters declaredParameters()
{
    prat::CrustParameters p;
    p.equilibriumMassFraction = prat::massFractionFromMoleFraction(0.1);
    p.depositionPorosity = 0.43;
    p.crystalDensity = 2165.0;
    p.reactionCoefficient = 2.3e-3;
    p.beta = 2.0;
    p.solutionDensity = prat::brineDensity(p.equilibriumMassFraction);
    p.vapourDensity = 0.0;
    p.crustSaturation = 1.0;
    p.cloggingPorosity = 0.01;
    p.connectivityPorosity = 0.02;
    p.substrateSize = 25.0e-6;
    p.molecularDiffusivity = 1.3e-9;
    p.areaClosure = prat::AreaClosure::PratEquation17;
    p.diffusivityClosure = prat::DiffusivityClosure::SqrtPorosity;
    p.check();
    return p;
}

prat::GasFilm declaredFilm()
{
    prat::GasFilm f;
    f.temperature = 293.15;
    f.relativeHumidity = 0.4;
    f.thickness = 0.0083;
    f.atmosphericPressure = 1.0e5;
    return f;
}

void add(std::map<std::string, double>& out, const std::string& name, double value)
{
    out[name] = value;
}

std::map<std::string, double> compute()
{
    std::map<std::string, double> v;
    const prat::CrustParameters params = declaredParameters();
    const prat::GasFilm film = declaredFilm();
    const double cEq = params.equilibriumMassFraction;

    add(v, "constant_C_eq", cEq);
    add(v, "constant_rho", params.rho());
    add(v, "constant_a_v0", params.area(params.depositionPorosity));
    add(v, "thermo_psat_293_15", prat::saturationPressure(293.15));
    add(v, "thermo_dv_293_15_1e5", prat::vapourDiffusionCoefficient(293.15, 1.0e5));
    add(v, "thermo_gaswmf_0_02", prat::gasWaterMassFraction(0.02));
    add(v, "thermo_w_from_x_0_1", prat::massFractionFromMoleFraction(0.1));
    add(v, "thermo_x_from_w_Ceq", prat::moleFractionFromMassFraction(cEq));
    add(v, "thermo_molality_Ceq", prat::molalityFromMassFraction(cEq));
    add(v, "thermo_brine_density_Ceq", prat::brineDensity(cEq));
    add(v, "activity_osmotic_1", prat::osmoticCoefficient(1.0));
    add(v, "activity_osmotic_6", prat::osmoticCoefficient(6.0));
    add(v, "activity_from_w_Ceq", prat::waterActivityFromMassFraction(cEq));
    add(v, "film_transfer_coefficient", film.transferCoefficient());
    add(v, "film_flux_fresh", film.fluxFromMassFraction(0.0));
    add(v, "film_flux_Ceq", film.fluxFromMassFraction(cEq));
    add(v, "closure_area_0_43", params.area(0.43));
    add(v, "closure_area_0_25", params.area(0.25));
    add(v, "closure_diffusivity_0_43", params.diffusivity(0.43));
    add(v, "closure_diffusivity_0_25", params.diffusivity(0.25));
    add(v, "closure_reaction_length_0_43", params.reactionLength(0.43));
    add(v, "closure_damkohler_star_0_43", params.damkohlerStar(0.43));

    const double delta = params.reactionLength(params.depositionPorosity);
    prat::CrustState state;
    state.topCellSize = std::min(delta / CELLS_PER_DELTA, SEED_HEIGHT / 8.0);
    const double flux = film.fluxFromMassFraction(cEq);
    prat::activate(state, params, SEED_HEIGHT, flux, WINDOW_LENGTHS * delta, RATIO,
                   MAXIMUM_CELL_LENGTHS * delta);

    add(v, "seed_cells", static_cast<double>(state.thicknesses.size()));
    add(v, "seed_height", state.height);
    add(v, "seed_top_cell", state.topCellSize);
    add(v, "seed_solid_salt", state.solidSaltMass);
    add(v, "seed_dissolved_salt", state.dissolvedSalt);
    add(v, "seed_water", state.water);

    const std::vector<std::pair<std::string, double>> bases = {
        {"Ceq", -1.0}, {"0_99Ceq", 0.99 * cEq}, {"1_01Ceq", 1.01 * cEq}};
    for (const auto& entry : bases)
    {
        const prat::ProfileSolution s =
            prat::solveProfile(state, params, flux, entry.second);
        const std::string prefix = "profile_" + entry.first + "_";
        add(v, prefix + "s_top", s.topSupersaturation);
        add(v, prefix + "v_s", s.growthVelocity);
        add(v, prefix + "internal_reaction", s.internalReaction);
        add(v, prefix + "dissolved_salt", s.dissolvedSalt);
        add(v, prefix + "water", s.water);
        add(v, prefix + "rho_u_base", s.solutionMassFluxBase);
        add(v, prefix + "salt_flux_transport", s.saltFluxBaseTransport);
        add(v, prefix + "base_gradient", s.baseGradient);
    }

    prat::StepRecord record;
    for (int i = 0; i < STEPS; ++i)
        record = prat::step(state, params, flux, DT);

    add(v, "step10_time", state.time);
    add(v, "step10_height", state.height);
    add(v, "step10_cells", static_cast<double>(state.thicknesses.size()));
    add(v, "step10_minimum_porosity", record.minimumPorosity);
    add(v, "step10_solid_salt", state.solidSaltMass);
    add(v, "step10_dissolved_salt", state.dissolvedSalt);
    add(v, "step10_water", state.water);
    add(v, "step10_salt_flux_soil", record.saltFluxSoil);
    add(v, "step10_water_flux_soil", record.waterFluxSoil);
    add(v, "step10_salt_flux_transport", record.saltFluxBaseTransport);
    add(v, "step10_growth_velocity", record.growthVelocity);
    add(v, "step10_internal_reaction", record.internalReaction);
    return v;
}

} // namespace

int main(int argc, char** argv)
{
    const std::string path = argc > 1 ? argv[1] : "reference_vectors.txt";
    std::ifstream file(path);
    if (!file)
    {
        std::cerr << "cannot open " << path << "\n";
        return 2;
    }

    std::map<std::string, double> expected;
    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream stream(line);
        std::string name;
        double value = 0.0;
        if (!(stream >> name >> value)) continue;
        expected[name] = value;
    }

    std::map<std::string, double> actual;
    try
    {
        actual = compute();
    }
    catch (const std::exception& error)
    {
        std::cerr << "computation threw: " << error.what() << "\n";
        return 2;
    }

    int failures = 0, missing = 0;
    double worst = 0.0;
    std::string worstName;
    std::cout << std::scientific << std::setprecision(6);
    for (const auto& entry : expected)
    {
        const auto found = actual.find(entry.first);
        if (found == actual.end())
        {
            std::cout << "MISSING  " << entry.first << "\n";
            ++missing;
            continue;
        }
        const double denominator = std::max(std::abs(entry.second), 1.0e-300);
        const double relative = std::abs(found->second - entry.second) / denominator;
        if (relative > worst) { worst = relative; worstName = entry.first; }
        if (!(relative <= TOLERANCE))
        {
            std::cout << "FAIL     " << entry.first << "  expected "
                      << entry.second << "  got " << found->second
                      << "  rel " << relative << "\n";
            ++failures;
        }
    }
    for (const auto& entry : actual)
        if (expected.find(entry.first) == expected.end())
            std::cout << "EXTRA    " << entry.first << "\n";

    std::cout << "\nchecked " << expected.size() << " vectors, " << failures
              << " failures, " << missing << " missing\n"
              << "worst relative difference " << worst << " at " << worstName
              << " (tolerance " << TOLERANCE << ")\n";
    return (failures == 0 && missing == 0) ? 0 : 1;
}
