#include "bdf.hh"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    p2h::Parameters p;
    p.ceq=0.2649408697676898;
    p.bottomC=p.ceq-1e-5;
    p.reaction=1e-3;
    p.crystal=2165;
    p.splitReactionLedgers=true;
    p2h::ALE model(16,2,p);
    auto before=model.initial(3.5e-6,p.ceq);
    p2h::BDF solver(model,before,0.05,1e-9,true);
    const auto after=solver.advance(0.5);
    const auto& salt=solver.integratedSaltFaces();
    const auto& reaction=solver.integratedReaction();
    if(salt.size()!=17 || reaction.size()!=16)
        throw std::runtime_error("quadrature size mismatch");
    double worst=0;
    for(int i=0;i<model.n;++i) {
        const double q0=p.rho*before[i]*before[model.n+i]*before[2*model.n]*model.widths[i];
        const double q1=p.rho*after[i]*after[model.n+i]*after[2*model.n]*model.widths[i];
        const double expected=salt[i]-salt[i+1]-reaction[i];
        const double scale=std::max({std::abs(q1-q0),std::abs(salt[i]),std::abs(salt[i+1]),
                                     std::abs(reaction[i]),1e-30});
        worst=std::max(worst,std::abs((q1-q0-expected)/scale));
    }
    if(!(worst<2e-7))
        throw std::runtime_error("quadrature local salt balance failed");
    p2h::BDF plain(model,before,0.05);
    plain.advance(0.1);
    if(!plain.integratedSaltFaces().empty() || !plain.integratedReaction().empty())
        throw std::runtime_error("diagnostic-disabled solver produced quadrature");
    p2h::Parameters limited=p;
    limited.suppressDissolution=true;
    limited.bottomC=limited.ceq;
    p2h::ALE limitedModel(16,2,limited);
    const auto limitedUnder=limitedModel.flux(limitedModel.initial(3.5e-6,limited.ceq-1e-4));
    if(!std::all_of(limitedUnder.reaction.begin(),limitedUnder.reaction.end(),
                    [](double value){return value==0;}))
        throw std::runtime_error("no-dissolution reaction clamp failed");
    const auto limitedOver=limitedModel.flux(limitedModel.initial(3.5e-6,limited.ceq+1e-4));
    if(!std::all_of(limitedOver.reaction.begin(),limitedOver.reaction.end(),
                    [](double value){return value>0;}))
        throw std::runtime_error("no-dissolution clamp removed precipitation");
    std::cout<<"PASS accepted-step quadrature local salt balance worst_relative="<<worst<<'\n';
    std::cout<<"PASS saturated ALE boundary adds no phase source and dissolution is suppressed\n";
}
