// SPDX-License-Identifier: GPL-3.0-or-later
#include "pratcoupling.hh"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

prat::CrustParameters parameters()
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
    p.check();
    return p;
}

prat::GasFilm film()
{
    prat::GasFilm f;
    f.temperature = 293.15;
    f.relativeHumidity = 0.4;
    f.thickness = 0.0083;
    f.atmosphericPressure = 1.0e5;
    return f;
}

prat::CrustState seeded(const prat::CrustParameters& p, const prat::GasFilm& f)
{
    prat::CrustState state;
    const double delta = p.reactionLength(p.depositionPorosity);
    state.topCellSize = std::min(delta / 8.0, 1.0e-5 / 8.0);
    prat::activate(state, p, 1.0e-5,
                   f.fluxFromMassFraction(p.equilibriumMassFraction),
                   40.0 * delta, 1.25, 10.0 * delta);
    return state;
}

auto compatibleSoilFlux(const prat::CrustState& state,
                        const prat::CrustParameters& p,
                        const prat::GasFilm& f,
                        double dt)
{
    const double target = p.equilibriumMassFraction;
    const double salt = prat::crustResponse(state, p, f, target, dt).saltMass;
    return [target, salt](double, double w) {
        return salt + 0.1 * (w - target);
    };
}

template<class Callable>
void requireThrows(Callable&& call, const std::string& message)
{
    bool threw = false;
    try { call(); }
    catch (const prat::DomainExit&) { threw = true; }
    require(threw, message);
}

} // namespace

int main()
{
    try
    {
        const auto p = parameters();
        const auto f = film();
        const double dt = 5.0;
        const auto initial = seeded(p, f);
        const auto soilFlux = compatibleSoilFlux(initial, p, f, dt);

        // Repeated residual/Jacobian evaluations are pure and bit-repeatable.
        auto a = prat::coupling::evaluate(initial, p, f, soilFlux, dt);
        auto b = prat::coupling::evaluate(initial, p, f, soilFlux, dt);
        require(prat::coupling::exactlyEqual(initial, a.baseState),
                "evaluation changed the persistent state");
        require(prat::coupling::exactlyEqual(a.proposedState, b.proposedState),
                "repeat evaluation produced a different proposed state");
        require(a.fluxes.saltMass == b.fluxes.saltMass
                && a.fluxes.totalMass == b.fluxes.totalMass,
                "repeat evaluation produced different flux bits");
        require(a.interfaceRoot.converged,
                "interface root was not marked converged");

        // A rejected Newton attempt is represented by discarding its trial.
        auto persistent = initial;
        prat::coupling::rollback(persistent, a);
        require(prat::coupling::exactlyEqual(persistent, initial),
                "rollback changed the persistent state");
        requireThrows([&] { prat::coupling::commit(persistent, a); },
                      "a discarded trial was accepted");

        // Newton time-step retry: the failed large-dt trial never reaches state.
        auto failedAttempt = prat::coupling::evaluate(
            persistent, p, f, compatibleSoilFlux(persistent, p, f, dt), dt);
        prat::coupling::rollback(persistent, failedAttempt);
        const double retryDt = 0.5 * dt;
        auto retry = prat::coupling::evaluate(
            persistent, p, f, compatibleSoilFlux(persistent, p, f, retryDt), retryDt);
        const auto expectedRetry = retry.proposedState;
        prat::coupling::commit(persistent, retry);
        require(prat::coupling::exactlyEqual(persistent, expectedRetry),
                "accepted retry did not commit exactly once");
        requireThrows([&] { prat::coupling::commit(persistent, retry); },
                      "the same trial committed twice");

        // A trial becomes stale when another accepted step advances the state.
        auto staleBase = initial;
        auto stale = prat::coupling::evaluate(
            staleBase, p, f, compatibleSoilFlux(staleBase, p, f, dt), dt);
        auto winner = prat::coupling::evaluate(
            staleBase, p, f, compatibleSoilFlux(staleBase, p, f, retryDt), retryDt);
        prat::coupling::commit(staleBase, winner);
        requireThrows([&] { prat::coupling::commit(staleBase, stale); },
                      "a stale trial was committed");

        // Forced root failure must leave the caller-owned state untouched.
        prat::coupling::SolverOptions forced;
        forced.maximumIterations = 1;
        forced.relativeTolerance = 0.0;
        forced.absoluteFluxTolerance = 0.0;
        auto forcedState = initial;
        requireThrows([&] {
            prat::coupling::evaluate(forcedState, p, f, soilFlux, dt, forced);
        }, "forced non-convergence was accepted");
        require(prat::coupling::exactlyEqual(forcedState, initial),
                "forced failure changed persistent state");

        std::cout << "PASS: repeat evaluation, rollback, Newton retry, stale "
                     "trial, single commit, and forced failure\n";
        std::cout << "interface iterations=" << b.interfaceRoot.iterations
                  << " residual=" << b.interfaceRoot.residual
                  << " gas residual=" << b.gasFilmResidual
                  << " top residual=" << b.topClosureResidual << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
