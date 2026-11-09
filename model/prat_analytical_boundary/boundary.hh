#pragma once
#include "analyticalgrowth.hh"
#include "progressive.hh"
#include <functional>

// Reaction-free reduction of Prat interface balances, NOT paper equation 28.
namespace pratboundary {
using Parameters=p2f::Parameters;
using State=p2f::State;
using Trial=p2f::Trial;
inline State activate(const Parameters& p,double solid,double salt,double water) {
    if(!(solid>0 && salt>=0 && water>0)) throw p2f::DomainExit("invalid activation inventory");
    State s; s.height=solid/(p.crystalDensity*(1-p.porosity));
    s.solidSalt=solid; s.dissolvedSalt=salt; s.water=water; s.active=true;
    const double c=salt/(salt+water),rho=p2f::brineDensity(c);
    if(std::abs(c-p.equilibriumMassFraction)>p.interfaceRelativeTolerance*p.equilibriumMassFraction)
        throw p2f::DomainExit("activation composition outside equilibrium envelope");
    if(std::abs((salt+water)/rho-p.porosity*s.height)>1e-10*p.porosity*s.height)
        throw p2f::DomainExit("activation pore volume inconsistent");
    return s;
}
inline Trial evaluate(const State& s,const Parameters& p,double dt,const std::function<double(double)>& film) {
    if(!s.active || !(dt>0)) throw p2f::DomainExit("invalid active trial");
    const double G=p2f::growthPrefactor(p,1)*p.damkohlerTop(s.height)/(1+p.damkohlerTop(s.height));
    double lo=0,hi=film(p.equilibriumMassFraction);
    if(!(std::isfinite(hi) && hi>=0)) throw p2f::DomainExit("condensation outside scope");
    for(int i=0;i<64;++i) {
        const double mid=.5*(lo+hi),c=p2f::topMassFraction(p,G*mid);
        if(mid>film(c)) hi=mid; else lo=mid;
    }
    const double j=.5*(lo+hi),v=G*j,dh=v*dt;
    Trial t; t.baseRevision=s.revision; t.proposed=s; t.proposed.revision++;
    t.proposed.height+=dh; t.proposed.time+=dt;
    const double solid=p.crystalDensity*(1-p.porosity)*dh;
    const double salt=p.rho()*p.porosity*p.equilibriumMassFraction*dh;
    const double water=p.rho()*p.porosity*(1-p.equilibriumMassFraction)*dh;
    // Increment only: never reset the actual activation pore composition.
    t.proposed.solidSalt+=solid; t.proposed.dissolvedSalt+=salt; t.proposed.water+=water;
    t.proposed.cumulativeEvaporation+=j*dt;
    t.evaporation=j; t.velocity=v; t.saltFlux=(solid+salt)/dt;
    t.waterFlux=j+water/dt; t.totalFlux=t.saltFlux+t.waterFlux;
    t.topMassFraction=p2f::topMassFraction(p,v);
    if(std::abs(t.topMassFraction-p.equilibriumMassFraction)>p.interfaceRelativeTolerance*p.equilibriumMassFraction)
        throw p2f::DomainExit("crust top outside equilibrium envelope");
    t.solutionPe=p2f::solutionPeclet(p,t.proposed.height,j,v);
    t.movingPe=p2f::movingPeclet(p,t.proposed.height,v);
    if(t.solutionPe>p.maximumSolutionPeclet || t.movingPe>p.maximumMovingPeclet)
        throw p2f::DomainExit("analytical Peclet validity exceeded");
    return t;
}
inline void validateSupply(const Parameters& p,double c,double sw) {
    if(sw<p.minimumInterfaceSaturation) throw p2f::DomainExit("hydraulic supply screen failed");
    if(std::abs(c-p.equilibriumMassFraction)>p.interfaceRelativeTolerance*p.equilibriumMassFraction)
        throw p2f::DomainExit("soil interface outside equilibrium envelope");
}
}
