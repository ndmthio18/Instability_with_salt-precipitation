#include "boundary.hh"
#include <iostream>
#include <iomanip>
void require(bool ok,const char* name) { if(!ok) throw std::runtime_error(name); }
int main() {
    pratboundary::Parameters p; p.reactionCoefficient=.001; p.check();
    p2h::PrecursorInventory pore;
    const double h0=3.5e-6,solid=p.crystalDensity*(1-p.porosity)*h0;
    for(int i=0;i<100;++i) p2h::acceptProgressiveIncrement(pore,
        p2h::progressiveIncrement(solid/100,1,p.porosity,p.crystalDensity,p.equilibriumMassFraction,p.rho()));
    auto s=pratboundary::activate(p,solid,pore.poreSalt,pore.poreWater);
    require(s.solidSalt==solid && s.water==pore.poreWater && s.dissolvedSalt==pore.poreSalt,"activation continuity");
    auto a=pratboundary::evaluate(s,p,.5,[](double){return 1e-5;},pratboundary::GrowthLaw::M1,25e-6);
    auto b=pratboundary::evaluate(s,p,.5,[](double){return 1e-5;});
    require(a.proposed.height==b.proposed.height && s.height==h0,"repeat evaluation purity");
    require(std::abs((a.proposed.totalSalt()-s.totalSalt())-.5*a.saltFlux)<1e-17,"salt flux inventory");
    require(std::abs((a.proposed.water-s.water)+.5*a.evaporation-.5*a.waterFlux)<1e-17,"water flux inventory");
    p2f::rollback(s,a); bool refused=false;
    try { p2f::commit(s,a); } catch(const p2f::DomainExit&) { refused=true; }
    require(refused && s.height==h0,"discarded transaction"); p2f::commit(s,b);
    refused=false; try { p2f::commit(s,b); } catch(const p2f::DomainExit&) { refused=true; }
    require(refused,"double commit");
    refused=false; try { pratboundary::validateSupply(p,p.equilibriumMassFraction,.98); } catch(const p2f::DomainExit&) {refused=true;}
    require(refused,"insufficient supply");
    refused=false; try { pratboundary::validateSupply(p,p.equilibriumMassFraction*1.021,1); } catch(const p2f::DomainExit&) {refused=true;}
    require(refused,"concentration screen");
    const auto d=pratboundary::growthDiagnostics(p,h0,25e-6,pratboundary::GrowthLaw::M1);
    require(d.damkohler>10 && d.m3Valid,"large Da validity at activation");
    require(d.efficiencyM1>d.efficiencyM3 && d.efficiencyM3>=d.efficiencyM2,"reference efficiency ordering");
    const auto thin=pratboundary::growthDiagnostics(p,1e-7,25e-6,pratboundary::GrowthLaw::M3);
    require(!thin.m3Valid,"M3 invalid below Da threshold");
    for(auto law:{pratboundary::GrowthLaw::M1,pratboundary::GrowthLaw::M2,pratboundary::GrowthLaw::M3}) {
        auto q=pratboundary::activate(p,solid,pore.poreSalt,pore.poreWater);
        auto t=pratboundary::evaluate(q,p,.5,[](double){return 1e-5;},law,25e-6);
        require(std::abs((t.proposed.totalSalt()-q.totalSalt())-.5*t.saltFlux)<1e-17,"law salt flux inventory");
        require(std::abs((t.proposed.water-q.water)+.5*t.evaporation-.5*t.waterFlux)<1e-17,"law water flux inventory");
    }
    auto zero=pratboundary::evaluate(s,p,1,[](double){return 0.;});
    require(zero.velocity==0 && zero.saltFlux==0 && zero.waterFlux==0,"zero evaporation");
    const double tiny=1e-14,huge=1.,k=p.heightCoefficient();
    require(std::abs(p2f::growthVelocity(p,tiny,1e-5)/(p2f::growthPrefactor(p,1e-5)*k*tiny)-1)<1e-6,"small Da");
    require(std::abs(p2f::growthVelocity(p,huge,1e-5)/p2f::growthPrefactor(p,1e-5)-1)<1e-6,"large Da");
    std::cout<<std::setprecision(17)<<"dt,h_final,implicit_integral_error,da,da_star,da_top,efficiency_m1,efficiency_m2,efficiency_m3\n";
    for(double dt:{1.,.5,.25,.125}) {
        auto q=pratboundary::activate(p,solid,pore.poreSalt,pore.poreWater);
        for(int i=0;i<int(1000/dt);++i) { auto t=pratboundary::evaluate(q,p,dt,[](double){return 1e-5;}); p2f::commit(q,t); }
        const double exactIntegral=q.height-h0+std::log(q.height/h0)/k;
        std::cout<<dt<<','<<q.height<<','<<exactIntegral-p2f::growthPrefactor(p,1e-5)*1000<<','
            <<d.damkohler<<','<<d.damkohlerStar<<','<<d.damkohlerTop<<','<<d.efficiencyM1<<','<<d.efficiencyM2<<','<<d.efficiencyM3<<'\n';
    }
    std::cerr<<"PASS M1 legacy path, M2/M3 formulae, Da>10 flag, filling, fluxes, rollback, supply screens\n";
}
