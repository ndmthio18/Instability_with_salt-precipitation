// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef P2F_ANALYTICAL_GROWTH_HH
#define P2F_ANALYTICAL_GROWTH_HH

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace p2f {

constexpr double MWater = 0.01801528;
constexpr double MNaCl = 0.05844;

struct DomainExit : public std::runtime_error
{ using std::runtime_error::runtime_error; };

inline double massFractionFromMoleFraction(double x)
{
    if (!(x >= 0.0 && x < 1.0)) throw DomainExit("invalid mole fraction");
    return x*MNaCl / ((1.0-x)*MWater + x*MNaCl);
}

inline double brineDensity(double w)
{
    if (!(w >= 0.0 && w < 1.0)) throw DomainExit("invalid mass fraction");
    return 998.203 + 1000.0*w*(0.670804 + 0.373854*w);
}

struct Parameters
{
    double equilibriumMassFraction = massFractionFromMoleFraction(0.1);
    double porosity = 0.43;
    double crystalDensity = 2165.0;
    double reactionCoefficient = 2.3e-3;
    double beta = 2.0;
    double molecularDiffusivity = 1.3e-9;
    double solutionDensity = 0.0;
    double vapourDensity = 0.0;
    double minimumInterfaceSaturation = 0.99;
    double interfaceRelativeTolerance = 0.02;
    double maximumSolutionPeclet = 0.1;
    double maximumMovingPeclet = 0.1;

    void check() const
    {
        if (!(equilibriumMassFraction > 0.0 && equilibriumMassFraction < 1.0))
            throw DomainExit("C_eq must lie in (0,1)");
        if (!(porosity > 0.0 && porosity < 1.0))
            throw DomainExit("porosity must lie in (0,1)");
        if (!(crystalDensity > 0.0 && reactionCoefficient > 0.0
              && beta > 0.0 && molecularDiffusivity > 0.0))
            throw DomainExit("positive physical parameters required");
    }

    double rho() const
    { return solutionDensity > 0.0 ? solutionDensity : brineDensity(equilibriumMassFraction); }

    double effectiveDiffusivity() const
    { return std::sqrt(porosity)*molecularDiffusivity; }

    double heightCoefficient() const
    {
        return beta*(1.0-porosity)*(1.0-equilibriumMassFraction)*reactionCoefficient
               /(porosity*effectiveDiffusivity());
    }

    double damkohlerTop(double height) const
    {
        if (!(height > 0.0 && std::isfinite(height)))
            throw DomainExit("active crust requires finite h > 0");
        return heightCoefficient()*height;
    }
};

inline double growthPrefactor(const Parameters& p, double evaporation)
{
    if (!(evaporation >= 0.0 && std::isfinite(evaporation)))
        throw DomainExit("evaporation must be finite and nonnegative");
    const double c = p.equilibriumMassFraction;
    return c*evaporation/(p.crystalDensity*(1.0-p.porosity)*(1.0-c));
}

inline double growthVelocity(const Parameters& p, double height, double evaporation)
{
    const double daTop = p.damkohlerTop(height);
    return growthPrefactor(p, evaporation)*daTop/(1.0+daTop);
}

inline double largeDaTopVelocity(const Parameters& p, double evaporation)
{ return growthPrefactor(p, evaporation); }

inline double solutionMassFlux(const Parameters& p, double evaporation, double velocity)
{
    return evaporation + (p.rho()*p.porosity
        + p.crystalDensity*(1.0-p.porosity) - p.vapourDensity)*velocity;
}

inline double solutionPeclet(const Parameters& p, double height,
                             double evaporation, double velocity)
{
    return solutionMassFlux(p, evaporation, velocity)*height
           /(p.rho()*p.porosity*p.effectiveDiffusivity());
}

inline double movingPeclet(const Parameters& p, double height, double velocity)
{ return velocity*height/p.effectiveDiffusivity(); }

inline double topMassFraction(const Parameters& p, double velocity)
{
    return p.equilibriumMassFraction
        + p.crystalDensity*velocity/(p.beta*p.rho()*p.reactionCoefficient);
}

struct State
{
    double height = 0.0;
    double solidSalt = 0.0;
    double dissolvedSalt = 0.0;
    double water = 0.0;
    double time = 0.0;
    double cumulativeEvaporation = 0.0;
    bool active = false;
    std::uint64_t revision = 0;

    double totalSalt() const { return solidSalt + dissolvedSalt; }
};

inline State activate(const Parameters& p, double seedHeight)
{
    if (!(seedHeight > 0.0 && std::isfinite(seedHeight)))
        throw DomainExit("seed height must be positive and finite");
    State s;
    s.height = seedHeight;
    s.solidSalt = p.crystalDensity*(1.0-p.porosity)*seedHeight;
    s.dissolvedSalt = p.rho()*p.porosity*p.equilibriumMassFraction*seedHeight;
    s.water = p.rho()*p.porosity*(1.0-p.equilibriumMassFraction)*seedHeight;
    s.active = true;
    return s;
}

struct Trial
{
    std::uint64_t baseRevision = 0;
    State proposed;
    double evaporation = 0.0;
    double velocity = 0.0;
    double saltFlux = 0.0;
    double waterFlux = 0.0;
    double totalFlux = 0.0;
    double topMassFraction = 0.0;
    double solutionPe = 0.0;
    double movingPe = 0.0;
    bool committed = false;
    bool discarded = false;
};

inline Trial evaluate(const State& state, const Parameters& p,
                      double evaporation, double dt)
{
    if (!state.active) throw DomainExit("cannot grow inactive crust");
    if (!(dt > 0.0 && std::isfinite(dt))) throw DomainExit("invalid dt");
    Trial t;
    t.baseRevision = state.revision;
    t.evaporation = evaporation;
    t.velocity = growthVelocity(p, state.height, evaporation);
    t.proposed = state;
    t.proposed.height += t.velocity*dt;
    t.proposed.time += dt;
    t.proposed.cumulativeEvaporation += evaporation*dt;
    t.proposed.revision += 1;
    t.proposed.solidSalt = p.crystalDensity*(1.0-p.porosity)*t.proposed.height;
    t.proposed.dissolvedSalt = p.rho()*p.porosity*p.equilibriumMassFraction*t.proposed.height;
    t.proposed.water = p.rho()*p.porosity*(1.0-p.equilibriumMassFraction)*t.proposed.height;
    t.saltFlux = (t.proposed.totalSalt()-state.totalSalt())/dt;
    const double waterStorage = (t.proposed.water-state.water)/dt;
    t.waterFlux = evaporation - p.vapourDensity*t.velocity + waterStorage;
    t.totalFlux = t.saltFlux + t.waterFlux;
    t.topMassFraction = p2f::topMassFraction(p, t.velocity);
    t.solutionPe = solutionPeclet(p, state.height, evaporation, t.velocity);
    t.movingPe = movingPeclet(p, state.height, t.velocity);
    return t;
}

inline void commit(State& state, Trial& trial)
{
    if (trial.committed || trial.discarded) throw DomainExit("trial not committable");
    if (state.revision != trial.baseRevision) throw DomainExit("stale trial");
    state = trial.proposed;
    trial.committed = true;
}

inline void rollback(const State& state, Trial& trial)
{
    if (trial.committed) throw DomainExit("committed trial cannot roll back");
    if (state.revision != trial.baseRevision) throw DomainExit("state changed before rollback");
    trial.discarded = true;
}

struct Validity
{
    double interfaceRelativeError = 0.0;
    double interfaceSaturation = 0.0;
    double solutionPe = 0.0;
    double movingPe = 0.0;
    bool valid = false;
};

inline Validity validity(const Parameters& p, const Trial& t,
                         double interfaceMassFraction,
                         double interfaceSaturation)
{
    Validity v;
    v.interfaceRelativeError = std::abs(
        interfaceMassFraction/p.equilibriumMassFraction - 1.0);
    v.interfaceSaturation = interfaceSaturation;
    v.solutionPe = std::abs(t.solutionPe);
    v.movingPe = std::abs(t.movingPe);
    v.valid = v.interfaceRelativeError <= p.interfaceRelativeTolerance
           && interfaceSaturation >= p.minimumInterfaceSaturation
           && v.solutionPe <= p.maximumSolutionPeclet
           && v.movingPe <= p.maximumMovingPeclet;
    return v;
}

} // namespace p2f

#endif
