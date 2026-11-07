#include "ale.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=2) return 2;
    std::ifstream in(argv[1]); int n; double j; in>>n>>j;
    p2h::ALE m(n,3); m.p.evaporation=j;
    auto y=m.initial(); for(double& v:y) in>>v;
    if(!in) return 3;
    auto f=m.flux(y); auto dy=m.rhs(y);
    std::cout<<std::setprecision(17)<<f.topC<<' '<<f.speed<<' ';
    for(const auto& a:{f.liquid,f.salt,f.mesh,dy}) for(double v:a) std::cout<<v<<' ';
    std::cout<<'\n';
}
