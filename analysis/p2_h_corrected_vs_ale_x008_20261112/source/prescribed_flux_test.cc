#include "ale.hh"
#include <algorithm>
#include <cmath>
#include <iostream>

int main() {
    try {
        p2h::Parameters p;
        p.ceq=0.264940869768;
        p.bottomC=p.ceq;
        p.initialHeight=3.5e-6;
        p2h::ALE legacy(16,3,p), forced(16,3,p);
        const auto y=legacy.initial();
        const double t=123.0, back=1.7e-6, evap=1.8e-5;
        forced.bottomBackDiffusion=[&](double time){return back+0.0*time;};
        forced.evaporationHistory=[&](double time){return evap+0.0*time;};
        const auto fl=legacy.flux(y,t), ff=forced.flux(y,t);
        const double expected=p.ceq*ff.liquid[0]-back;
        if(std::abs(ff.salt[0]-expected)>1e-18) throw std::runtime_error("Neumann salt flux mismatch");
        if(std::abs(ff.evaporation-evap)>1e-18) throw std::runtime_error("evaporation history mismatch");
        if(std::abs(fl.salt[0]-p.ceq*fl.liquid[0])>1e-18) throw std::runtime_error("legacy path changed");
        const auto rhs=forced.rhs(y,t);
        if(!std::all_of(rhs.begin(),rhs.end(),[](double x){return std::isfinite(x);}))
            throw std::runtime_error("non-finite forced RHS");
        auto dissolved=y; std::fill(dissolved.begin(),dissolved.begin()+forced.n,1.0);
        forced.p.allowCompleteDissolution=true;
        forced.validate(dissolved);
        const auto fd=forced.flux(dissolved,t);
        if(!std::isfinite(fd.salt[0])) throw std::runtime_error("epsilon=1 flux invalid");
        std::cout<<"PASS prescribed-flux operator\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
