#include "bdf.hh"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <filesystem>

int main(int argc,char** argv) {
 try {
  if(argc!=11) throw std::runtime_error("id n power rtol maxstep u beta av jfactor forcing_csv");
  const std::string id=argv[1]; int n=std::stoi(argv[2]); double power=std::stod(argv[3]);
  const bool denseEvents=id.find("_event")!=std::string::npos;
  const double rtol=std::stod(argv[4]),maxstep=std::stod(argv[5]),u=std::stod(argv[6]);
  p2h::Parameters p; p.ceq=.264940869768;p.rho=1202.168577;p.crystal=2165;p.birth=.43;
  p.initialHeight=3.5e-6;p.reaction=1e-3;p.beta=std::stod(argv[7]);p.avMultiplier=std::stod(argv[8]);
  p.diffusion=1.3e-9;p.grain=25e-6;p.evaporation=1.78e-5*std::stod(argv[9]);
  p.splitReactionLedgers=true;p.allowCompleteDissolution=true;p.prescribedDiffusionUsesEquilibriumInflow=false;
  const double clockScale=std::max(1.0,std::stod(argv[9]));
  std::vector<double> times,js,bs; double end=180000/clockScale;
  if(u<0) { std::ifstream file(argv[10]);std::string line;std::getline(file,line);
   while(std::getline(file,line)){std::stringstream ss(line);std::string s;std::vector<double> v;
    while(std::getline(ss,s,','))v.push_back(std::stod(s));
    times.push_back(v.at(0));js.push_back(v.at(2));bs.push_back(v.at(3));}
   if(times.empty())throw std::runtime_error("missing forcing");end=times.back(); }
  auto interp=[&](double t,const std::vector<double>&v){auto it=std::upper_bound(times.begin(),times.end(),t);
   if(it==times.begin())return v.front();if(it==times.end())return v.back();
   auto k=it-times.begin();double w=(t-times[k-1])/(times[k]-times[k-1]);return v[k-1]*(1-w)+v[k]*w;};
  p2h::ALE m(n,power,p);m.evaporationHistory=[&](double t){return u<0?interp(t,js):p.evaporation;};
  m.bottomBackDiffusion=[&](double t){return u<0?interp(t,bs):u*p.ceq*p.evaporation;};
  auto y=m.initial(),initial=y;auto inv0=m.inventory(y);p2h::BDF solver(m,y,maxstep,rtol,true);
  std::filesystem::create_directories("raw");
  solver.rejectedState=[&](const p2h::Vec& state,double time){
   std::ofstream file("raw/"+id+"_rejected.csv");
   file<<std::setprecision(17)<<"t,cell,eps,C\n";
   for(int i=0;i<n;i++)file<<time<<','<<i<<','<<state[i]<<','<<state[n+i]<<'\n';
  };
  std::ofstream hist("raw/"+id+"_history.csv"),prof("raw/"+id+"_profiles.csv");
  hist<<std::setprecision(17)<<"t,h,v,ct,cb,j,B,Fb,Jb,advb,diffb,G,Rplus,Dminus,solid,dissolved,water,epsmin,epsmax,epsmean,Qb,Qtop,Qplus,Qminus,Qliquid,Qevap,salt_res,water_res,diss_res,local_res\n";
  prof<<std::setprecision(17)<<"t,cell,xi,width,h,eps,C,R,storage,advection,diffusion,mesh,reaction\n";
  double previous=0,last=0;std::string status="end";
  auto write=[&](double t,bool profile){auto f=m.flux(y,t);auto iv=m.inventory(y);int k=m.inv();
   double pos=0,neg=0,mean=0,local=0;
   for(int i=0;i<n;i++){pos+=std::max(0.,f.reaction[i])*f.volume[i];neg+=std::max(0.,-f.reaction[i])*f.volume[i];mean+=y[i]*m.widths[i];
    if(t>0){double q=p.rho*y[i]*y[n+i]*f.volume[i],q0=p.rho*initial[i]*initial[n+i]*p.initialHeight*m.widths[i];
     local=std::max(local,std::abs(q-q0-solver.integratedSaltFaces()[i]+solver.integratedSaltFaces()[i+1]+solver.integratedReaction()[i]));}}
   double sr=iv.solid+iv.salt-inv0.solid-inv0.salt-y[k+2],wr=iv.water-inv0.water-y[k+3]+y[k+2]+y[k+4];
   double dr=iv.salt-inv0.salt-y[k+2]+y[k]+y[m.positiveInv()]-y[m.dissolutionInv()];
   double emin=*std::min_element(y.begin(),y.begin()+n),emax=*std::max_element(y.begin(),y.begin()+n);
   hist<<t<<','<<y[2*n]<<','<<f.speed<<','<<f.topC<<','<<f.bottomC<<','<<f.evaporation<<','<<m.bottomBackDiffusion(t)<<','<<f.liquid[0]<<','<<f.salt[0]<<','<<f.bottomAdvectiveSalt<<','<<f.bottomDiffusiveSalt<<','<<f.salt[n]<<','<<pos<<','<<neg<<','<<iv.solid<<','<<iv.salt<<','<<iv.water<<','<<emin<<','<<emax<<','<<mean<<','<<y[k+2]<<','<<y[k]<<','<<y[m.positiveInv()]<<','<<y[m.dissolutionInv()]<<','<<y[k+3]<<','<<y[k+4]<<','<<sr<<','<<wr<<','<<dr<<','<<local<<'\n';
   if(std::max({std::abs(sr),std::abs(wr),std::abs(dr),local})>2e-9)throw std::runtime_error("budget_failed");
   if(profile || emax>=.999 || emin<=.05){
    p2h::Vec adv(n+1),diff(n+1),mesh(n+1);adv[0]=f.bottomAdvectiveSalt;diff[0]=f.bottomDiffusiveSalt;mesh[0]=0;
    for(int a=1;a<n;a++){double rel=f.liquid[a]+p.rho*f.mesh[a],c=rel>=0?y[n+a-1]:y[n+a];
     adv[a]=f.liquid[a]*c;mesh[a]=p.rho*f.mesh[a]*c;diff[a]=f.salt[a]-adv[a]-mesh[a];}
    adv[n]=f.liquid[n]*f.topC;mesh[n]=p.rho*f.mesh[n]*f.topC;diff[n]=f.salt[n]-adv[n]-mesh[n];
    for(int i=0;i<n;i++){double vol=f.volume[i],a=(adv[i]-adv[i+1])/vol,d=(diff[i]-diff[i+1])/vol,g=(mesh[i]-mesh[i+1])/vol,r=-f.reaction[i];
     prof<<t<<','<<i<<','<<m.centers[i]<<','<<m.widths[i]<<','<<y[2*n]<<','<<y[i]<<','<<y[n+i]<<','<<f.reaction[i]<<','<<a+d+g+r<<','<<a<<','<<d<<','<<g<<','<<r<<'\n';}}
   hist.flush();prof.flush();
   if(emax>=.999)return std::string("solid_framework_limit_eps_0.999");
   if(emin<=.05)return std::string("pore_filling_limit_eps_0.05");return std::string("end");
  };
  write(0,true);
  std::vector<double> schedule={.001,.01,.1,1,2,5,10,20,30,60};
  for(double& t:schedule)t/=clockScale;
  if(denseEvents)for(double t=1/clockScale;t<=600/clockScale;t+=1/clockScale)schedule.push_back(t);
  for(double t=120/clockScale;t<=end;t+=60/clockScale)schedule.push_back(t);
  schedule.push_back(end);
  std::sort(schedule.begin(),schedule.end());schedule.erase(std::unique(schedule.begin(),schedule.end()),schedule.end());
  for(double t:schedule){if(t>end)break;try{y=solver.advance(t);last=t;
    status=write(t,t<=60/clockScale || std::abs(t*clockScale/600-std::round(t*clockScale/600))<1e-7 || t==end);
    if(status!="end")break;previous=t;
   }catch(const std::exception&e){status=std::string("solver_stop:")+e.what();break;}}
  std::ofstream s("raw/"+id+"_status.txt");s<<std::setprecision(17)<<"id="<<id<<"\nstatus="<<status<<"\nprevious_s="<<previous<<"\nlast_s="<<last<<"\ntarget_s="<<end<<"\nsteps="<<solver.steps()<<'\n';
  std::cout<<id<<' '<<status<<" at "<<last<<'\n';return status.find("solver_stop")==0?2:0;
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
