#include "two_bdf.hh"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <memory>
using namespace twofront;
int main(int argc,char** argv){try{
 if(argc!=11)throw std::runtime_error("id nc nl kappa gap rtol maxstep u end forcing");
 std::string id=argv[1];int nc=std::stoi(argv[2]),nl=std::stoi(argv[3]);
 double kb=std::stod(argv[4]),gap=std::stod(argv[5]),rtol=std::stod(argv[6]),step=std::stod(argv[7]);
 double u=std::stod(argv[8]),end=std::stod(argv[9]);
 p2h::Parameters p;p.ceq=.264940869768;p.rho=1202.168577;p.crystal=2165;p.birth=.43;
 p.diffusion=1.3e-9;p.grain=25e-6;p.reaction=.001;p.beta=2;p.initialHeight=3.5e-6;p.evaporation=1.78e-5;
 Vec times,js,bs;
 if(u<0){std::ifstream input(argv[10]);std::string line;std::getline(input,line);
   while(std::getline(input,line)){std::stringstream ss(line);std::string x;Vec row;
     while(std::getline(ss,x,','))row.push_back(std::stod(x));times.push_back(row.at(0));js.push_back(row.at(2));bs.push_back(row.at(3));}
   if(times.empty())throw std::runtime_error("missing forcing");end=std::min(end,times.back());}
 auto interpolation=[&](double t,const Vec& v){auto it=std::upper_bound(times.begin(),times.end(),t);
   if(it==times.begin())return v.front();if(it==times.end())return v.back();auto k=it-times.begin();double a=(t-times[k-1])/(times[k]-times[k-1]);return (1-a)*v[k-1]+a*v[k];};
 Model m(nc,nl,p,kb,gap);m.jhistory=[&](double t){return u<0?interpolation(t,js):p.evaporation;};
 m.Bhistory=[&](double t){return u<0?interpolation(t,bs):u*p.ceq*p.evaporation;};
 Vec y=m.initial();auto initial=m.inventory(y);double initialb=gap;
 auto solver=std::make_unique<Solver>(m,y,rtol,step);
 std::filesystem::create_directories("raw");std::ofstream out("raw/"+id+"_history.csv"),pro("raw/"+id+"_profiles.csv");
 out<<std::setprecision(17)<<"t,zt,zb,hc,vt,vb,vnet,G,Db,Rplus,Dminus,solid,crust_salt,liquid_salt,water,epsmin,epsmax,phib,Cb,Ct,Csoil,j,B,Fplus,Fminus,Jsoil,msoil,Splus,Sminus,core_phi,core_std,core_volume,salt_res,water_res,solid_res,liquid_salt_res,crust_salt_res,liquid_water_res,crust_water_res,contact_res,water_jump";
 for(auto name:{"TOP","DB","RP","RD","SOIL_SALT","SOIL_MASS","EVAP","SPLUS","SMINUS","MPLUS","MMINUS"})out<<",Q_"<<name;
 out<<'\n';pro<<std::setprecision(17)<<"t,region,cell,zlo,zhi,phi,C,R\n";
 double last=0,previous=0;int transfers=0;std::string status="end";
 auto write=[&](double t,bool profile){auto f=m.flux(y,t);auto a=m.inventory(y);int k=m.ledger();
  double h=y[m.topIndex()]-y[m.bottomIndex()],rp=0,rd=0;
  double ell=std::sqrt(std::pow(p.birth,1.5)*p.diffusion/((1-p.ceq)*(36*(1-p.birth)/p.grain)*p.reaction));
  double cv=0,cphi=0,cphi2=0;
  for(int i=0;i<m.n;++i){rp+=std::max(0.,f.R[i])*f.vol[i];rd+=std::max(0.,-f.R[i])*f.vol[i];
   double v=std::max(0.,std::min(h*m.faces[i+1],h-3*ell)-std::max(h*m.faces[i],3*ell));
   cv+=v;cphi+=v*(1-f.eps[i]);cphi2+=v*std::pow(1-f.eps[i],2);}
  double sr=a.solid+a.crustSalt+a.liquidSalt-initial.solid-initial.crustSalt-initial.liquidSalt-y[k+SOIL_SALT];
  double wr=a.water-initial.water-y[k+SOIL_MASS]+y[k+SOIL_SALT]+y[k+EVAP];
  double cr=a.solid-initial.solid-y[k+TOP]+y[k+DB]-y[k+RP]+y[k+RD];
  double lr=a.liquidSalt-initial.liquidSalt-y[k+SOIL_SALT]+y[k+SMINUS];
  double pr=a.crustSalt-initial.crustSalt-y[k+SPLUS]+y[k+TOP]+y[k+RP]-y[k+RD];
  double lw=p.rho*(y[m.bottomIndex()]-initialb)-a.liquidSalt+initial.liquidSalt-y[k+SOIL_MASS]+y[k+SOIL_SALT]+y[k+MMINUS]-y[k+SMINUS];
  double pw=wr-lw;
  double emin=*std::min_element(f.eps.begin(),f.eps.end()),emax=*std::max_element(f.eps.begin(),f.eps.end());
  double Fp=(f.m[0]-p.rho*f.eps[0]*f.vb)*f.Cb-f.S[0],Fm=(f.ml-p.rho*f.vb)*f.Cb-(f.S[0]-f.Db);
  out<<t<<','<<y[m.topIndex()]<<','<<y[m.bottomIndex()]<<','<<h<<','<<f.vt<<','<<f.vb<<','<<f.vt-f.vb<<','<<f.S[m.n]<<','<<f.Db<<','<<rp<<','<<rd<<','<<a.solid<<','<<a.crustSalt<<','<<a.liquidSalt<<','<<a.water<<','<<emin<<','<<emax<<','<<1-f.eps[0]<<','<<f.Cb<<','<<f.Ct<<','<<f.Csoil<<','<<f.j<<','<<f.B<<','<<Fp<<','<<Fm<<','<<f.Sl[0]<<','<<f.ml<<','<<f.S[0]<<','<<f.S[0]-f.Db<<',';
  if(cv>0)out<<cphi/cv<<','<<std::sqrt(std::max(0.,cphi2/cv-std::pow(cphi/cv,2)));else out<<"nan,nan";
  out<<','<<cv<<','<<sr<<','<<wr<<','<<cr<<','<<lr<<','<<pr<<','<<lw<<','<<pw<<','<<f.contactResidual<<','<<f.waterJump;
  for(int i=0;i<NLEDGER;++i)out<<','<<y[k+i];out<<'\n';out.flush();
  if(profile || emin<=.05 || 1-f.eps[0]<1e-8){
   for(int i=0;i<m.l;++i)pro<<t<<",liquid,"<<i<<','<<y[m.bottomIndex()]*m.lf[i]<<','<<y[m.bottomIndex()]*m.lf[i+1]<<",0,"<<f.cl[i]<<",0\n";
   for(int i=0;i<m.n;++i)pro<<t<<",crust,"<<i<<','<<y[m.bottomIndex()]+h*m.faces[i]<<','<<y[m.bottomIndex()]+h*m.faces[i+1]<<','<<1-f.eps[i]<<','<<f.cc[i]<<','<<f.R[i]<<'\n';pro.flush();}
  if(std::max({std::abs(sr),std::abs(wr),std::abs(cr),std::abs(lr),std::abs(pr),std::abs(lw),std::abs(pw)})>2e-9)return std::string("budget_failed");
  if(emin<=.05)return std::string("pore_filling_screen");
  if(h<=1e-9)return std::string("vanishing_active_crust");
  if(1-f.eps[0]<1e-8 && y[0]>0)return std::string("intact_trace_lost");
  return std::string("end");};
 write(0,true);Vec schedule={.001,.01,.1,1,2,5,10,30,60};
 for(double t=120;t<=end;t+=60)schedule.push_back(t);schedule.push_back(end);
 std::sort(schedule.begin(),schedule.end());schedule.erase(std::unique(schedule.begin(),schedule.end()),schedule.end());
 for(double t:schedule){if(t>end)break;try{y=solver->advance(t);last=t;
   if(y[0]==0){solver.reset();auto c=convertZeroPrefix(m,y);m=std::move(c.first);y=std::move(c.second);
      solver=std::make_unique<Solver>(m,y,rtol,step,t);++transfers;}
   status=write(t,t<=60 || std::fmod(t,3600)==0 || t==end);if(status!="end")break;previous=t;
  }catch(const std::exception& e){status=std::string("solver_stop:")+e.what();break;}}
 std::ofstream report("raw/"+id+"_status.txt");report<<std::setprecision(17)<<"id="<<id<<"\nstatus="<<status<<"\nlast_s="<<last<<"\nprevious_s="<<previous<<"\ntarget_s="<<end<<"\ntransfers="<<transfers<<'\n';
 std::cout<<id<<' '<<status<<" last="<<last<<'\n';return status=="end"?0:2;
 }catch(const std::exception& e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}}
