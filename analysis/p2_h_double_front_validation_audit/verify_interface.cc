#include "two_front.hh"
#include <iostream>
#include <iomanip>
using namespace audit28;
int main(){
 p2h::Parameters p;p.ceq=.264940869768;p.rho=1202.168577;p.crystal=2165;p.birth=.43;
 p.diffusion=1.3e-9;p.grain=25e-6;p.reaction=.001;p.beta=2;p.initialHeight=3.5e-6;p.evaporation=1.78e-5;
 Model m(32,8,p,1,1e-9);auto initial=m.initial();auto f0=m.flux(initial);double gc=p.rho*std::pow(p.birth,1.5)*p.diffusion/(p.initialHeight*m.centers[0]);
 double gl=p.rho*p.diffusion/(m.gap0*(1-m.lc.back())),A=p.rho*p.reaction*(1-p.birth);
 double expected=-p.rho/p.crystal*p.reaction*(gc+gl)/(gc+gl+A*(1-p.ceq));
 for(double e:{-1e-9,-1e-10,1e-10,1e-9}){
  std::vector<Dual> y(m.size());for(int i=0;i<m.size();++i)y[i]=initial[i];
  for(int i=0;i<m.n;++i)y[m.n+i]=Dual(p.rho*(f0.vol[i]-initial[i])*e,p.rho*(f0.vol[i]-initial[i]));
  for(int i=0;i<m.l;++i)y[2*m.n+i]=Dual(p.rho*f0.vl[i]*e,p.rho*f0.vl[i]);
  auto f=m.flux(y);long double target=e<0?expected:0;
  long double error=std::abs(f.vb.d-target)/(std::abs(expected));
  std::cout<<std::setprecision(17)<<e<<" derivative="<<(double)f.vb.d<<" expected="<<(double)target<<" relative_error="<<(double)error<<'\n';
  if(error>1e-8)return 1;
 }
 if(f0.vb!=0)return 2;
 std::cout<<"PASS one-sided interface derivatives on both sides of saturation; derivative at kink is not unique\n";
}
