#include "ale.hh"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition,const char* message) {
    if(!condition) throw std::runtime_error(message);
}
bool close(double a,double b,double tolerance=2e-13) {
    return std::abs(a-b)<=tolerance*std::max({1.,std::abs(a),std::abs(b)});
}
}

int main() {
    try {
        p2h::Parameters p;
        require(p.reaction==2.3e-3,"historical reaction default changed");
        require(p.avMultiplier==1.,"historical area default changed");
        p.reaction=1e-3;
        p2h::ALE reference(64,2,p);
        auto y=reference.initial(7e-6,p.ceq+2e-4);
        for(int i=0;i<reference.n;++i) {
            y[i]=.25+.15*double(i)/(reference.n-1);
            y[reference.n+i]=p.ceq+1e-5+2e-4*double(i)/(reference.n-1);
        }
        y[2*reference.n]=7e-6;
        auto state=[&](double multiplier) {
            auto q=p; q.avMultiplier=multiplier;
            p2h::ALE model(64,2,q);
            return model.flux(y);
        };
        const auto low=state(.1),base=state(1),high=state(10);
        for(int i=0;i<reference.n;++i) {
            require(close(low.reaction[i]*10,base.reaction[i]),"0.1av reaction scaling failed");
            require(close(high.reaction[i],base.reaction[i]*10),"10av reaction scaling failed");
        }
        require(close(low.topC,base.topC) && close(high.topC,base.topC),"AvMultiplier changed direct top concentration law");
        require(close(low.speed,base.speed) && close(high.speed,base.speed),"AvMultiplier changed direct top kinetic law");
        auto undersaturated=y;
        std::fill(undersaturated.begin()+reference.n,undersaturated.begin()+2*reference.n,p.ceq-1e-5);
        const auto negative=reference.flux(undersaturated);
        require(std::all_of(negative.reaction.begin(),negative.reaction.end(),[](double value){return value<0;}),
                "internal dissolution sign handling failed");
        auto splitParameters=p; splitParameters.splitReactionLedgers=true;
        p2h::ALE split(64,2,splitParameters);
        auto splitState=y;
        splitState.resize(split.size(),0.);
        const auto derivative=split.rhs(splitState);
        require(close(derivative[split.positiveInv()]-derivative[split.dissolutionInv()],derivative[split.inv()+1]),
                "positive/dissolution ledger rates do not close to net reaction");
        std::cout<<"PASS av scaling and top-law independence\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';
        return 1;
    }
}
