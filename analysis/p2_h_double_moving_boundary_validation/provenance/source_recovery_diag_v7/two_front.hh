#pragma once
#include "ale.hh"
#include <numeric>

namespace twofront {
using p2h::Vec;
enum Ledger { TOP, DB, RP, RD, SOIL_SALT, SOIL_MASS, EVAP, SPLUS, SMINUS, MPLUS, MMINUS, NLEDGER };
struct Flux {
    Vec eps,cc,cl,ec,el,vol,vl,R,m,H,S,Sl,U,Ul;
    double vt=0,vb=0,Db=0,Cb=0,Ct=0,Csoil=0,j=0,B=0,ml=0;
    double contactResidual=0,waterJump=0;
};
struct Inventory { double solid=0,crustSalt=0,liquidSalt=0,water=0,pore=0; };
class Model {
public:
    p2h::Parameters p;
    int n,l;
    double kappa=1,gap0=1e-9,referenceC=0;
    bool bottomDirichlet=false;
    double bathC=.26;
    Vec faces,centers,widths,lf,lc,lw;
    std::function<double(double)> jhistory,Bhistory;
    p2h::ALE original;
    Model(int nc,int nl,p2h::Parameters parameters,double kb=1,double seed=1e-9,bool symmetric=true)
      :p(parameters),n(nc),l(nl),kappa(kb),gap0(seed),original(nc,2,parameters) {
        referenceC=p.ceq;
        if(n<4 || l<0 || kappa<0 || (l==0 && (gap0!=0 || kappa!=0)))
            throw std::runtime_error("invalid two-front configuration");
        faces=original.faces;
        if(symmetric)for(int i=0;i<=n;++i)faces[i]=.5*(1-std::cos(std::acos(-1.)*i/n));
        for(int i=0;i<n;++i){widths.push_back(faces[i+1]-faces[i]);centers.push_back(.5*(faces[i]+faces[i+1]));}
        original.faces=faces;original.centers=centers;original.widths=widths;
        for(int i=0;i<=l && l>0;++i)lf.push_back(double(i)/l);
        for(int i=0;i<l;++i){lw.push_back(1./l);lc.push_back((i+.5)/l);}
    }
    int topIndex() const{return 2*n+l;}
    int bottomIndex() const{return topIndex()+1;}
    int ledger() const{return topIndex()+2;}
    int size() const{return ledger()+NLEDGER;}
    Vec initial(double concentration=-1) const {
        const double c=concentration<0?p.ceq:concentration;
        Vec y(size(),0);y[topIndex()]=gap0+p.initialHeight;y[bottomIndex()]=gap0;
        for(int i=0;i<n;++i){double v=p.initialHeight*widths[i];y[i]=(1-p.birth)*v;y[n+i]=p.rho*p.birth*v*(c-referenceC);}
        for(int i=0;i<l;++i)y[2*n+i]=p.rho*gap0*lw[i]*(c-referenceC);
        return y;
    }
    Flux flux(const Vec& y,double time=0) const {
        const double b=y[bottomIndex()],h=y[topIndex()]-b;
        if(!(h>0) || (l>0 && !(b>0)))throw std::runtime_error("nonpositive region length");
        Flux f;f.eps.resize(n);f.cc.resize(n);f.cl.resize(l);f.ec.resize(n);f.el.resize(l);f.vol.resize(n);f.vl.resize(l);f.R.resize(n);
        f.m.resize(n+1);f.H.resize(n+1);f.S.resize(n+1);f.Sl.resize(l+1);f.U.resize(n+1);f.Ul.resize(l+1);
        f.j=jhistory?jhistory(time):p.evaporation;f.B=Bhistory?Bhistory(time):0;
        for(int i=0;i<n;++i){
            f.vol[i]=h*widths[i];
            // Domain extension for unaccepted Newton iterates only.
            const double solid=std::clamp(y[i],0.,f.vol[i]*(1-1e-12));
            f.eps[i]=1-solid/f.vol[i];f.ec[i]=y[n+i]/(p.rho*(f.vol[i]-solid));f.cc[i]=referenceC+f.ec[i];
            f.R[i]=p.rho*36*(1-f.eps[i])/p.grain*p.avMultiplier*p.reaction*(f.ec[i]+(referenceC-p.ceq));
        }
        for(int i=0;i<l;++i){f.vl[i]=b*lw[i];f.el[i]=y[2*n+i]/(p.rho*f.vl[i]);f.cl[i]=referenceC+f.el[i];}
        if(p.beta==0){
            if(f.j!=0)throw std::runtime_error("disabled top test requires zero evaporation");
            f.Ct=f.cc.back();f.vt=0;
        }else{
            auto result=original.top(f.cc.back(),f.eps.back(),std::sqrt(f.eps.back())*p.diffusion,h,f.j);
            f.Ct=result.first;f.vt=result.second;
        }
        f.m[n]=f.j+(p.rho*p.birth+p.crystal*(1-p.birth))*f.vt;
        for(int i=n-1;i>=0;--i)f.m[i]=f.m[i+1]+(1-p.rho/p.crystal)*f.R[i]*f.vol[i];
        const double gc=p.rho*std::pow(f.eps[0],1.5)*p.diffusion/(h*centers[0]);
        const double gl=l?p.rho*p.diffusion/(b*(1-lc.back())):0;
        const double phi=1-f.eps[0],A=p.rho*p.reaction*kappa*phi;
        const double rhs=gc*f.ec[0]+(l?gl*f.el.back():-f.B);
        const double deficit=p.ceq-referenceC,q=1-referenceC;
        const double M=gc+gl+A*(deficit+q),forced=rhs+A*deficit*q;
        const double excess=rhs/(gc+gl)<deficit?2*forced/(M+std::sqrt(M*M-4*A*forced)):rhs/(gc+gl);
        f.Cb=referenceC+excess;
        f.vb=p.rho/p.crystal*p.reaction*kappa*std::max(deficit-excess,0.);
        f.Db=p.crystal*phi*f.vb;
        f.ml=f.m[0]-(1-p.rho/p.crystal)*f.Db;
        f.H[0]=-f.vb*f.eps[0];f.H[n]=-f.vt*p.birth;
        const double Fplus=gc*(f.ec[0]-excess);
        f.U[0]=(f.m[0]-p.rho*f.eps[0]*f.vb)*excess-Fplus;
        f.S[0]=(f.m[0]-p.rho*f.eps[0]*f.vb)*referenceC+f.U[0];
        // One shared jump flux prevents subtractive roundoff at a vanishing seed.
        const double Sminus=f.S[0]-f.Db;
        const double Uminus=f.U[0]-(1-referenceC)*f.Db;
        const double Fminus=(f.ml-p.rho*f.vb)*excess-Uminus;
        f.contactResidual=l?Fminus-gl*(excess-f.el.back()):Fminus-f.B;
        f.waterJump=(f.m[0]-p.rho*f.eps[0]*f.vb-f.S[0])-(f.ml-p.rho*f.vb-Sminus);
        if(l){
            const double soilg=p.rho*p.diffusion/(b*lc[0]);
            const double es=bottomDirichlet?bathC-referenceC:f.el[0]-f.B/soilg;
            f.Csoil=referenceC+es;
            const double soilF=bottomDirichlet?soilg*(f.el[0]-(bathC-referenceC)):f.B;
            f.Ul[0]=f.ml*es-soilF;f.Ul[l]=Uminus;
            f.Sl[0]=f.ml*referenceC+f.Ul[0];f.Sl[l]=Sminus;
            for(int i=1;i<l;++i){
                double relative=f.ml-p.rho*lf[i]*f.vb;
                f.Ul[i]=relative*(relative>=0?f.el[i-1]:f.el[i])-p.rho*p.diffusion*(f.el[i]-f.el[i-1])/(b*(lc[i]-lc[i-1]));
                f.Sl[i]=relative*referenceC+f.Ul[i];
            }
        }else{f.Csoil=f.Cb;f.Sl[0]=f.S[0];f.Ul[0]=f.U[0];}
        for(int i=1;i<n;++i){
            double w=f.vb+faces[i]*(f.vt-f.vb);
            f.H[i]=-w*(w>=0?f.eps[i]:f.eps[i-1]);
            double relative=f.m[i]+p.rho*f.H[i];
            double resistance=h*((faces[i]-centers[i-1])/(std::pow(f.eps[i-1],1.5)*p.diffusion)
                                +(centers[i]-faces[i])/(std::pow(f.eps[i],1.5)*p.diffusion));
            f.U[i]=relative*(relative>=0?f.ec[i-1]:f.ec[i])-p.rho*(f.ec[i]-f.ec[i-1])/resistance;
            f.S[i]=relative*referenceC+f.U[i];
        }
        f.S[n]=p.crystal*(1-p.birth)*f.vt;
        f.U[n]=(1-referenceC)*f.S[n]-referenceC*f.j;
        return f;
    }
    Vec rhs(const Vec& y,double time) const {
        auto f=flux(y,time);Vec dy(size(),0);
        for(int i=0;i<n;++i){
            double dE=f.H[i]-f.H[i+1]-f.R[i]*f.vol[i]/p.crystal;
            dy[i]=(f.vt-f.vb)*widths[i]-dE;
            dy[n+i]=f.U[i]-f.U[i+1]-(1-referenceC)*f.R[i]*f.vol[i];
            dy[ledger()+RP]+=std::max(0.,f.R[i])*f.vol[i];
            dy[ledger()+RD]+=std::max(0.,-f.R[i])*f.vol[i];
        }
        for(int i=0;i<l;++i)dy[2*n+i]=f.Ul[i]-f.Ul[i+1];
        dy[topIndex()]=f.vt;dy[bottomIndex()]=f.vb;
        dy[ledger()+TOP]=f.S[n];dy[ledger()+DB]=f.Db;
        dy[ledger()+SOIL_SALT]=f.Sl[0];dy[ledger()+SOIL_MASS]=f.ml;dy[ledger()+EVAP]=f.j;
        dy[ledger()+SPLUS]=f.S[0];dy[ledger()+SMINUS]=f.S[0]-f.Db;
        dy[ledger()+MPLUS]=f.m[0]-p.rho*f.eps[0]*f.vb;dy[ledger()+MMINUS]=f.ml-p.rho*f.vb;
        return dy;
    }
    Inventory inventory(const Vec& y) const {
        Inventory a;double h=y[topIndex()]-y[bottomIndex()];
        for(int i=0;i<n;++i){a.solid+=p.crystal*y[i];a.crustSalt+=y[n+i]+p.rho*referenceC*(h*widths[i]-y[i]);a.pore+=h*widths[i]-y[i];}
        a.pore+=y[bottomIndex()];
        for(int i=0;i<l;++i)a.liquidSalt+=y[2*n+i]+p.rho*referenceC*y[bottomIndex()]*lw[i];
        a.water=p.rho*a.pore-a.crustSalt-a.liquidSalt;return a;
    }
    void validate(const Vec& y,double time=0) const {
        auto f=flux(y,time);double h=y[topIndex()]-y[bottomIndex()];
        for(int i=0;i<n;++i)if(y[i]<0 || y[i]>h*widths[i]*.98 || f.cc[i]<0 || f.cc[i]>.95)
            throw std::runtime_error("crust accepted-state bound");
        for(double c:f.cl)if(c<0 || c>.95)throw std::runtime_error("liquid accepted-state bound");
        if(f.Csoil<0 || f.Csoil>.95 || f.Cb<0 || f.Cb>.95)throw std::runtime_error("boundary concentration bound");
    }
};
// Transfer a connected exactly liquid-only prefix without any phase change.
// Preserve physical cell faces, all inventories and accumulated ledgers.
inline std::pair<Model,Vec> convertZeroPrefix(const Model& m,const Vec& y){
    int cut=0;while(cut<m.n && y[cut]==0)++cut;
    if(!cut)return {m,y};
    if(m.n-cut<4)throw std::runtime_error("remaining crust too small for active-grid conversion");
    double oldb=y[m.bottomIndex()],oldt=y[m.topIndex()],h=oldt-oldb;
    double newb=oldb+h*m.faces[cut],newh=oldt-newb;
    Model a(m.n-cut,m.l+cut,m.p,m.kappa,newb,false);a.bottomDirichlet=m.bottomDirichlet;a.bathC=m.bathC;a.referenceC=m.referenceC;
    a.jhistory=m.jhistory;a.Bhistory=m.Bhistory;
    for(int i=0;i<=a.n;++i)a.faces[i]=(oldb+h*m.faces[i+cut]-newb)/newh;
    for(int i=0;i<a.n;++i){a.widths[i]=a.faces[i+1]-a.faces[i];a.centers[i]=.5*(a.faces[i]+a.faces[i+1]);}
    for(int i=0;i<=m.l;++i)a.lf[i]=m.l?oldb*m.lf[i]/newb:0;
    for(int i=1;i<=cut;++i)a.lf[m.l+i]=(oldb+h*m.faces[i])/newb;
    for(int i=0;i<a.l;++i){a.lw[i]=a.lf[i+1]-a.lf[i];a.lc[i]=.5*(a.lf[i+1]+a.lf[i]);}
    a.original.faces=a.faces;a.original.widths=a.widths;a.original.centers=a.centers;
    Vec z(a.size(),0);z[a.topIndex()]=oldt;z[a.bottomIndex()]=newb;
    for(int i=0;i<a.n;++i){z[i]=y[i+cut];z[a.n+i]=y[m.n+i+cut];}
    for(int i=0;i<m.l;++i)z[2*a.n+i]=y[2*m.n+i];
    for(int i=0;i<cut;++i)z[2*a.n+m.l+i]=y[m.n+i];
    for(int i=0;i<NLEDGER;++i)z[a.ledger()+i]=y[m.ledger()+i];
    return {a,z};
}
}
