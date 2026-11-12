#include "two_bdf.hh"
#include "bdf.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace twofront;
void require(bool ok,const char* text){if(!ok)throw std::runtime_error(text);}
p2h::Parameters parameters(){p2h::Parameters p;p.ceq=.264940869768;p.rho=1202.168577;p.crystal=2165;
 p.birth=.43;p.diffusion=1.3e-9;p.grain=25e-6;p.reaction=.001;p.beta=2;p.initialHeight=3.5e-6;p.evaporation=1.78e-5;
 p.allowCompleteDissolution=true;p.splitReactionLedgers=true;p.prescribedDiffusionUsesEquilibriumInflow=false;return p;}
double budgets(const Model& m,const Vec& y,const Inventory& a){auto b=m.inventory(y);int k=m.ledger();
 double rs=b.solid+b.crustSalt+b.liquidSalt-a.solid-a.crustSalt-a.liquidSalt-y[k+SOIL_SALT];
 double rw=b.water-a.water-y[k+SOIL_MASS]+y[k+SOIL_SALT]+y[k+EVAP];
 double rc=b.solid-a.solid-y[k+TOP]+y[k+DB]-y[k+RP]+y[k+RD];
 double rl=b.liquidSalt-a.liquidSalt-y[k+SOIL_SALT]+y[k+SMINUS];
 double rp=b.crustSalt-a.crustSalt-y[k+SPLUS]+y[k+TOP]+y[k+RP]-y[k+RD];
 return std::max({std::abs(rs),std::abs(rw),std::abs(rc),std::abs(rl),std::abs(rp)});}
int main(){try{
 std::cout<<std::setprecision(17);
 auto p=parameters();Model m(32,0,p,0,0,false);auto y=m.initial();auto a=m.inventory(y);
 p2h::ALE old(32,2,p);old.bottomBackDiffusion=[](double){return 0.;};auto yo=old.initial();
 p2h::BDF prior(old,yo,1,1e-9,true);Solver current(m,y,1e-9,1);
 double recovery=0,budget=0;
 for(double t:{.01,.1,1.,10.,100.}){y=current.advance(t);yo=prior.advance(t);auto f=m.flux(y,t);auto fo=old.flux(yo,t);
   std::cout<<"RECOVERY_SAMPLE t="<<t<<" new="<<f.vt<<" old="<<fo.speed<<" rel="<<std::abs(f.vt/fo.speed-1)<<std::endl;
   recovery=std::max(recovery,std::abs(f.vt/fo.speed-1));budget=std::max(budget,budgets(m,y,a));
   require(y[m.bottomIndex()]==0,"fixed bottom moved");}
 require(recovery<2e-6,"fixed-bottom recovery mismatch");require(budget<2e-9,"recovery budget");
 std::cout<<"RECOVERY "<<recovery<<" BUDGET "<<budget<<'\n';
 // Independent constant-undersaturation asymptote, made accurate by large D.
 auto pd=p;pd.beta=0;pd.avMultiplier=0;pd.evaporation=0;pd.diffusion=1e-6;
 Model dissolve(32,8,pd,1,1e-9);dissolve.bottomDirichlet=true;dissolve.bathC=.26;dissolve.referenceC=.26;
 auto yd=dissolve.initial(.26);auto ad=dissolve.inventory(yd);Solver sd(dissolve,yd,1e-9,.001);
 yd=sd.advance(.1);double exact=pd.rho/pd.crystal*pd.reaction*(pd.ceq-.26)*.1;
 double recession=yd[dissolve.bottomIndex()]-dissolve.gap0;
 double err=std::abs(recession/exact-1);require(err<1e-4,"pure dissolution asymptote");
 require(std::abs(yd[dissolve.topIndex()]-(pd.initialHeight+dissolve.gap0))<1e-14,"top moved in dissolution test");
 require(budgets(dissolve,yd,ad)<2e-9,"dissolution budget");
 double massloss=ad.solid-dissolve.inventory(yd).solid;
 require(std::abs(massloss/(pd.crystal*(1-pd.birth)*recession)-1)<1e-8,"independent crystal loss");
 std::cout<<"DISSOLUTION "<<err<<" BUDGET "<<budgets(dissolve,yd,ad)<<'\n';
 double speeds[2],net[2];
 for(int k=0;k<2;++k){Model both(48,12,p,k?2:1,1e-9);both.Bhistory=[&](double){return .6*p.ceq*p.evaporation;};
   auto yy=both.initial();auto aa=both.inventory(yy);Solver ss(both,yy,1e-9,1);
   for(double t:{.01,.1,1.,10.,100.}){yy=ss.advance(t);auto f=both.flux(yy,t);auto rhs=both.rhs(yy,t);
     require(std::abs(rhs[both.topIndex()]-rhs[both.bottomIndex()]-(f.vt-f.vb))<1e-20,"geometry rate");
     require(budgets(both,yy,aa)<2e-9,"two-region budget");require(std::abs(f.waterJump)<1e-12,"water jump");}
   auto f=both.flux(yy,100);speeds[k]=f.vb;net[k]=f.vt-f.vb;
   require(f.vt>0 && f.vb>0,"simultaneous motion missing");
   std::cout<<"BOTH kappa="<<both.kappa<<" vt="<<f.vt<<" vb="<<f.vb<<" budget="<<budgets(both,yy,aa)<<'\n';}
 require(std::abs(speeds[1]/speeds[0]-1)>.05,"bottom closure perturbation hidden");
 require(std::abs(net[1]-net[0])>1e-11,"net comparison not independent");
 Model remap(16,4,p,1,1e-9);auto yr=remap.initial();
 // Manufactured exactly depleted prefix; charge its full liquid inventory.
 for(int i=0;i<2;++i){yr[i]=0;yr[remap.n+i]=0;}
 yr[8]=0;yr[remap.n+8]=0;
 auto before=remap.inventory(yr);auto converted=convertZeroPrefix(remap,yr);auto after=converted.first.inventory(converted.second);
 require(converted.first.n==14 && converted.first.l==6,"connected prefix selection");
 require(converted.second[6]==0,"isolated zero cell was deleted");
 require(std::abs(before.solid-after.solid)<1e-16,"remap crystal");
 require(std::abs(before.water-after.water)<1e-16,"remap water");
 require(std::abs(before.crustSalt+before.liquidSalt-after.crustSalt-after.liquidSalt)<1e-16,"remap salt");
 require(converted.second[converted.first.topIndex()]==yr[remap.topIndex()],"remap moved top");
 std::cout<<"REMAP conservative connected-prefix transfer PASS\n";
 std::cout<<"PASS P2 transport/recovery/dissolution/independence gates\n";return 0;
 }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
