#pragma once
#include "analyticalgrowth.hh"
#include "progressive.hh"
#include <functional>
#include <limits>
#include <string>

// Reaction-free reduction of Prat interface balances, NOT paper equation 28.
namespace pratboundary {
using Parameters=p2f::Parameters;
using State=p2f::State;

enum class GrowthLaw { M1, M2, M3 };

inline GrowthLaw parseGrowthLaw(const std::string& value) {
    if(value=="M1" || value=="m1") return GrowthLaw::M1;
    if(value=="M2" || value=="m2") return GrowthLaw::M2;
    if(value=="M3" || value=="m3") return GrowthLaw::M3;
    throw p2f::DomainExit("unknown analytical growth law: "+value);
}

inline const char* growthLawName(GrowthLaw law) {
    switch(law) {
        case GrowthLaw::M1: return "M1";
        case GrowthLaw::M2: return "M2";
        case GrowthLaw::M3: return "M3";
    }
    throw p2f::DomainExit("unreachable analytical growth law");
}

struct GrowthDiagnostics {
    double damkohler=0;
    double damkohlerStar=0;
    double damkohlerTop=0;
    double specificArea=0;
    double efficiencyM1=0;
    double efficiencyM2=0;
    double efficiencyM3=0;
    double selectedEfficiency=0;
    double m3RelativeError=0;
    bool m3Valid=false;
};

inline GrowthDiagnostics growthDiagnostics(const Parameters& p,double height,double substrateGrainSize,GrowthLaw law) {
    if(!(height>0 && std::isfinite(height))) throw p2f::DomainExit("active crust requires finite h > 0");
    if(!(substrateGrainSize>0 && std::isfinite(substrateGrainSize))) throw p2f::DomainExit("invalid substrate grain size");
    GrowthDiagnostics d;
    const double e=p.porosity,c=p.equilibriumMassFraction,kr=p.reactionCoefficient;
    const double deff=p.effectiveDiffusivity();
    d.specificArea=36*(1-e)/substrateGrainSize;
    d.damkohler=height*height*(1-c)*d.specificArea*kr/(e*deff);
    d.damkohlerStar=p.beta*p.beta*(1-e)*(1-e)*(1-c)*kr/(e*deff*d.specificArea);
    d.damkohlerTop=p.damkohlerTop(height);
    const double rootStar=std::sqrt(d.damkohlerStar);
    const double finiteFactor=rootStar*std::tanh(std::sqrt(d.damkohler));
    d.efficiencyM1=d.damkohlerTop/(1+d.damkohlerTop);
    d.efficiencyM2=finiteFactor/(1+finiteFactor);
    d.efficiencyM3=rootStar/(1+rootStar);
    d.selectedEfficiency=law==GrowthLaw::M1?d.efficiencyM1:(law==GrowthLaw::M2?d.efficiencyM2:d.efficiencyM3);
    d.m3RelativeError=(d.efficiencyM3-d.efficiencyM2)/d.efficiencyM2;
    d.m3Valid=d.damkohler>=10.0; // Prat (2026): use tanh(sqrt(Da)) ~= 1 for Da > 10.
    return d;
}

struct Trial : p2f::Trial { GrowthDiagnostics growth; GrowthLaw law=GrowthLaw::M1; };
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
inline Trial evaluate(const State& s,const Parameters& p,double dt,const std::function<double(double)>& film,
                      GrowthLaw law=GrowthLaw::M1,double substrateGrainSize=25e-6) {
    if(!s.active || !(dt>0)) throw p2f::DomainExit("invalid active trial");
    const auto diagnostics=growthDiagnostics(p,s.height,substrateGrainSize,law);
    const double G=p2f::growthPrefactor(p,1)*diagnostics.selectedEfficiency;
    double lo=0,hi=film(p.equilibriumMassFraction);
    if(!(std::isfinite(hi) && hi>=0)) throw p2f::DomainExit("condensation outside scope");
    for(int i=0;i<64;++i) {
        const double mid=.5*(lo+hi),c=p2f::topMassFraction(p,G*mid);
        if(mid>film(c)) hi=mid; else lo=mid;
    }
    const double j=.5*(lo+hi),v=G*j,dh=v*dt;
    Trial t; t.growth=diagnostics; t.law=law; t.baseRevision=s.revision; t.proposed=s; t.proposed.revision++;
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
