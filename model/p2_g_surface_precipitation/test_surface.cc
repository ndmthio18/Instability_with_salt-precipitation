#include "surface.hh"
#include <iostream>
int main() {
    auto check=[](bool ok){if(!ok)throw std::runtime_error("test failed");};
    check(p2g::precipitate(.1,1,.25)==0);
    const double p=p2g::precipitate(.5,1,.25);
    check(std::abs(p-1./6)<1e-15);
    check(std::abs((.5-p)/(1.5-p)-.25)<1e-15);
    check(p2g::precipitate(.5-p,1,.25)<1e-15);
    check(p2g::velocity(1e-6,0,.265)==0);
    check(p2g::velocity(1e-6,1e-5,.265)>0);
    check(p2g::velocity(1e-6,1e-5,.265)<p2g::velocity(1e-6,1e-5,.265,false));
    std::cout<<"PASS 7 precipitation/kinetics checks\n";
}
