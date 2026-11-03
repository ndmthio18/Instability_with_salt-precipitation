// Isolated WP8 v2: equation-scale diagnostics and budget-bounded convergence.
#pragma once
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <dumux/common/parameters.hh>

struct Wp8LinearInfo { double minPivot=1e300,maxPivot=0,growth=0,etaNorm=0,etaComponent=0; };
inline std::vector<double> wp8SolveV2(const std::vector<double>& original,const std::vector<double>& rhs,Wp8LinearInfo& info)
{
    auto a=original,b=rhs;const auto n=b.size();double before=0,after=0;
    for(std::size_t i=0;i<n;++i){double scale=0;for(std::size_t j=0;j<n;++j)scale=std::max(scale,std::abs(a[i*n+j]));if(!(scale>0))throw std::runtime_error("zero Jacobian row");for(std::size_t j=0;j<n;++j){a[i*n+j]/=scale;before=std::max(before,std::abs(a[i*n+j]));}b[i]/=scale;}
    for(std::size_t k=0;k<n;++k){auto pivot=k;for(auto i=k+1;i<n;++i)if(std::abs(a[i*n+k])>std::abs(a[pivot*n+k]))pivot=i;
        double piv=std::abs(a[pivot*n+k]);if(piv<1e-30)throw std::runtime_error("singular reduced Jacobian");info.minPivot=std::min(info.minPivot,piv);info.maxPivot=std::max(info.maxPivot,piv);
        if(pivot!=k){for(std::size_t j=k;j<n;++j)std::swap(a[k*n+j],a[pivot*n+j]);std::swap(b[k],b[pivot]);}
        for(auto i=k+1;i<n;++i){double q=a[i*n+k]/a[k*n+k];a[i*n+k]=0;for(auto j=k+1;j<n;++j){a[i*n+j]-=q*a[k*n+j];after=std::max(after,std::abs(a[i*n+j]));}b[i]-=q*b[k];}}
    for(std::size_t i=n;i-->0;){for(auto j=i+1;j<n;++j)b[i]-=a[i*n+j]*b[j];b[i]/=a[i*n+i];if(!std::isfinite(b[i]))throw std::runtime_error("nonfinite correction");}
    double rn=0,an=0,bn=0,dn=0;
    for(auto v:b)dn=std::max(dn,std::abs(v));
    for(std::size_t i=0;i<n;++i){long double e=rhs[i],den=std::abs(rhs[i]),ar=0;for(std::size_t j=0;j<n;++j){e-=static_cast<long double>(original[i*n+j])*b[j];den+=std::abs(static_cast<long double>(original[i*n+j])*b[j]);ar+=std::abs(original[i*n+j]);}info.etaComponent=std::max(info.etaComponent,static_cast<double>(std::abs(e)/std::max(den,1e-300L)));rn=std::max(rn,static_cast<double>(std::abs(e)));an=std::max(an,static_cast<double>(ar));bn=std::max(bn,std::abs(rhs[i]));}
    info.etaNorm=rn/std::max(an*dn+bn,1e-300);info.growth=std::max(before,after)/before;return b;
}

template<class Solution,class Geometry>
void wp8DiagnosticReplay(Solution& x,const Geometry& geometry)
{
    auto file=Dumux::getParam<std::string>("WP8.ReplayProfile","");if(file.empty())return;
    if(!Dumux::getParam<bool>("WP8.AllowDiagnosticReplay",false))throw std::runtime_error("replay is diagnostic only");
    std::ifstream f(file);if(!f)throw std::runtime_error("replay file missing");double z,p,c;std::map<double,std::pair<double,double>> values;
    while(f>>z>>p>>c){if(!std::isfinite(z+p+c))throw std::runtime_error("nonfinite replay");values.emplace(z,std::make_pair(p,c));}
    if(values.size()!=201)throw std::runtime_error("replay levels must match 201 original heights");
    for(const auto& v:vertices(geometry.gridView())){double y=v.geometry().center()[1];auto it=values.lower_bound(y-1e-12);if(it==values.end()||std::abs(it->first-y)>1e-12)throw std::runtime_error("replay height mismatch; no interpolation");auto i=geometry.vertexMapper().index(v);x[i][0]=it->second.first;x[i][1]=it->second.second;}
    std::cout<<"WP8_DIAGNOSTIC_REPLAY approximate saved state; no restart-equivalence claim\n";
}

template<class Assembler,class Solution,class Geometry>
void paper2ConstrainedStep(Assembler& assembler,Solution& x,const Geometry& geometry,double dt)
{
    if(geometry.gridView().comm().size()!=1)throw std::runtime_error("one rank required");
    std::map<double,std::size_t> levels;std::vector<std::size_t> group(x.size());for(const auto& v:vertices(geometry.gridView()))levels.emplace(v.geometry().center()[1],0);
    std::size_t l=0;std::vector<double> z;for(auto& kv:levels){kv.second=l++;z.push_back(kv.first);}for(const auto& v:vertices(geometry.gridView()))group[geometry.vertexMapper().index(v)]=levels.at(v.geometry().center()[1]);
    std::vector<std::size_t> first(l,x.size()),count(l,0);for(std::size_t i=0;i<x.size();++i){++count[group[i]];if(first[group[i]]==x.size())first[group[i]]=i;else for(int c=0;c<2;++c)if(x[i][c]!=x[first[group[i]]][c])throw std::runtime_error("nonuniform input rejected");}for(auto c:count)if(c!=2)throw std::runtime_error("one horizontal cell required");
    std::size_t n=2*l;std::vector<double> u(n);for(std::size_t j=0;j<l;++j)for(int c=0;c<2;++c)u[2*j+c]=x[first[j]][c];
    const auto criterion=Dumux::getParam<std::string>("WP8.Criterion","legacy");if(criterion!="legacy"&&criterion!="scaled")throw std::runtime_error("unknown criterion");
    const double factor=Dumux::getParam<double>("WP8.RoundoffFactor",32.0);if(factor<=0||factor>64)throw std::runtime_error("roundoff factor outside tested bounds");
    const int maxiter=Dumux::getParam<int>("WP8.MaxIterations",20);if(maxiter<0||maxiter>20)throw std::runtime_error("iteration cap outside bounded protocol");
    auto cuts=Dumux::getParam<std::vector<double>>("Transport.CutCoordinates");cuts.push_back((z[0]+z[1])/2);double W=geometry.bBoxMax()[0]-geometry.bBoxMin()[0];
    const double saltAllocation=1e-12/4800,waterAllocation=1e-10/4800;static int step=0;++step;
    std::ofstream diag("convergence_v2.csv",step==1?std::ios::out:std::ios::app);diag<<std::setprecision(17);if(step==1)diag<<"step,iteration,water_max,salt_max,relative_shift,scaled_nonlinear,linear_eta_norm,linear_eta_component,min_pivot,max_pivot,pivot_growth,regional_salt_error_m,regional_water_error_m,accepted\n";
    auto residualAt=[&](const std::vector<double>& v){auto trial=x;for(std::size_t i=0;i<x.size();++i)for(int c=0;c<2;++c)trial[i][c]=v[2*group[i]+c];assembler.gridVariables().update(trial);assembler.assembleJacobianAndResidual(trial);std::vector<double> out(n,0);for(std::size_t i=0;i<x.size();++i)for(int c=0;c<2;++c)out[2*group[i]+c]+=assembler.residual()[i][c];return out;};
    double shift=1;
    for(int iter=0;iter<=maxiter;++iter){
        for(std::size_t i=0;i<x.size();++i)for(int c=0;c<2;++c)x[i][c]=u[2*group[i]+c];assembler.gridVariables().update(x);assembler.assembleJacobianAndResidual(x);
        auto& f=assembler.residual();auto& jac=assembler.jacobian();std::vector<double> r(n,0),a(n*n,0);
        for(std::size_t i=0;i<x.size();++i){for(int c=0;c<2;++c)r[2*group[i]+c]+=f[i][c];for(auto it=jac[i].begin();it!=jac[i].end();++it)for(int c=0;c<2;++c)for(int d=0;d<2;++d)a[(2*group[i]+c)*n+2*group[it.index()]+d]+=(*it)[c][d];}
        double water=0,salt=0,scaled=0,regionalS=0,regionalW=0;
        for(std::size_t i=0;i<n;++i){if(!std::isfinite(r[i]))throw std::runtime_error("nonfinite residual");if(i%2)salt=std::max(salt,std::abs(r[i]));else water=std::max(water,std::abs(r[i]));long double scale=0;for(std::size_t j=0;j<n;++j)scale+=std::abs(static_cast<long double>(a[i*n+j]))*std::max(std::abs(u[j]),j%2?.0036:1.0);scaled=std::max(scaled,static_cast<double>(std::abs(r[i])/std::max(scale,1e-300L)));}
        for(double cut:cuts){long double rs=0,rw=0;for(std::size_t j=1;j<l;++j)if(z[j]>=cut){rw+=r[2*j];rs+=r[2*j+1];}regionalS=std::max(regionalS,static_cast<double>(std::abs(rs)*dt/W));regionalW=std::max(regionalW,static_cast<double>(std::abs(rw)*dt/W));}
        Wp8LinearInfo info;auto du=wp8SolveV2(a,r,info);
        bool accept=criterion=="legacy"?(water<=2e-18&&salt<=2e-20&&shift<=1e-10):(shift<=1e-10&&scaled<=factor*std::numeric_limits<double>::epsilon()&&info.etaComponent<=1e-12&&regionalS<=saltAllocation&&regionalW<=waterAllocation);
        diag<<step<<","<<iter<<","<<water<<","<<salt<<","<<shift<<","<<scaled<<","<<info.etaNorm<<","<<info.etaComponent<<","<<info.minPivot<<","<<info.maxPivot<<","<<info.growth<<","<<regionalS<<","<<regionalW<<","<<accept<<"\n";diag.flush();
        if(Dumux::getParam<bool>("WP8.Probe",false)&&step==1&&iter==2){
            std::ofstream mat("jacobian.csv"),state("state.csv"),fd("fd_columns.csv"),ulp("ulp_probes.csv");mat<<std::setprecision(17);state<<std::setprecision(17);fd<<std::setprecision(17);ulp<<std::setprecision(17);
            for(std::size_t i=0;i<n;++i){for(std::size_t j=0;j<n;++j)mat<<a[i*n+j]<<(j+1<n?',':'\n');state<<u[i]<<","<<r[i]<<","<<du[i]<<"\n";}
            fd<<"column,row,h,fd,assembled\n";
            for(auto level:std::vector<std::size_t>{1,l/2,l-2,l-1})for(int comp=0;comp<2;++comp)for(double rel:{1e-7,1e-6,1e-5}){auto col=2*level+comp;double h=rel*std::max(std::abs(u[col]),comp?.0036:1.0);auto up=u,dn=u;up[col]+=h;dn[col]-=h;auto fp=residualAt(up),fm=residualAt(dn);for(std::size_t i=0;i<n;++i)fd<<col<<","<<i<<","<<h<<","<<(fp[i]-fm[i])/(2*h)<<","<<a[i*n+col]<<"\n";}
            ulp<<"component,direction,row,residual,baseline\n";
            for(int comp=0;comp<2;++comp)for(int sign:{-1,1}){auto v=u;for(std::size_t j=1;j<l;++j)v[2*j+comp]=std::nextafter(v[2*j+comp],sign>0?std::numeric_limits<double>::infinity():-std::numeric_limits<double>::infinity());auto fp=residualAt(v);for(std::size_t i=0;i<n;++i)ulp<<comp<<","<<sign<<","<<i<<","<<fp[i]<<","<<r[i]<<"\n";}
            assembler.gridVariables().update(x);assembler.assembleJacobianAndResidual(x);
        }
        std::cout<<"WP8_V2 step="<<step<<" iteration="<<iter<<" accepted="<<accept<<" scaled="<<scaled<<" salt_budget="<<regionalS<<" water_budget="<<regionalW<<"\n";
        if(accept){std::cout<<"WP8_CONSTRAINT iterations="<<iter<<" residual="<<std::max(water,salt)<<" shift="<<shift<<"\n";return;}
        if(iter==maxiter)throw std::runtime_error("WP8 v2 rejected Newton step; no acceptance");shift=0;for(std::size_t j=0;j<n;++j){shift=std::max(shift,std::abs(du[j])/std::max(std::abs(u[j]),j%2?.0036:1.0));u[j]-=du[j];}
    }
}
