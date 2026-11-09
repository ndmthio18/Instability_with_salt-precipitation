#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace p2h {
using Vec = std::vector<double>;
struct Parameters {
    double ceq=.264, rho=1200, crystal=2170, birth=.43, diffusion=1.3e-9;
    double reaction=.0023, avMultiplier=1, beta=2, grain=25e-6, initialHeight=.001;
    double evaporation=1200*.01/86400, bottomC=.264;
    bool splitReactionLedgers=false, suppressDissolution=false;
};
struct Fluxes {
    Vec volume, reaction, liquid, salt, mesh;
    double topC, speed, evaporation;
};

// Exact M3 spatial operator. Primitive variables reproduce the pinned M3
// tolerance convention. No explicit time update is implemented in this class.
class ALE {
public:
    Parameters p;
    int n;
    Vec faces, centers, widths;
    std::function<double(double)> gasFilm;
    explicit ALE(int cells=640, double power=3, Parameters params={}) : p(params), n(cells) {
        if(n<4 || power<1) throw std::runtime_error("invalid ALE mesh");
        for(int i=0;i<=n;++i) faces.push_back(1-std::pow(1-double(i)/n,power));
        for(int i=0;i<n;++i) { centers.push_back((faces[i]+faces[i+1])/2); widths.push_back(faces[i+1]-faces[i]); }
    }
    int size() const { return 2*n+6+(p.splitReactionLedgers?2:0); }
    int inv() const { return 2*n+1; }
    int positiveInv() const { return 2*n+6; }
    int dissolutionInv() const { return 2*n+7; }
    double specificArea(double epsilon) const {
        return p.avMultiplier*36*(1-epsilon)/p.grain;
    }
    double reactionRate(double epsilon,double concentration) const {
        const double rate=p.rho*specificArea(epsilon)*p.reaction*(concentration-p.ceq);
        return p.suppressDissolution?std::max(0.,rate):rate;
    }
    Vec initial(double h=-1,double concentration=-1) const {
        Vec y(size(),0);
        std::fill(y.begin(),y.begin()+n,p.birth);
        std::fill(y.begin()+n,y.begin()+2*n,concentration<0?p.ceq:concentration);
        y[2*n]=h<0?p.initialHeight:h;
        return y;
    }
    std::pair<double,double> top(double c,double eps,double deff,double h,double j) const {
        const double g=p.rho*eps*deff/(h*(1-centers.back()));
        const double a=p.beta*p.rho*(1-p.birth)*p.reaction;
        const double b=j-a*(1-p.ceq)-g, d=j*p.ceq+g*(c-p.ceq);
        const double discriminant=b*b-4*a*d;
        if(!(discriminant>=0)) throw std::runtime_error("top discriminant invalid");
        const double q=-.5*(b+std::copysign(std::sqrt(discriminant),b));
        const double candidates[2]={d/q,q/a};
        bool found=false; double excess=0;
        for(double e:candidates) if(std::isfinite(e) && p.ceq+e>=0 && p.ceq+e<=1)
            if(!found || std::abs(e)<std::abs(excess)) { excess=e; found=true; }
        if(!found) throw std::runtime_error("no physical top root");
        return {p.ceq+excess,p.beta*p.rho/p.crystal*p.reaction*excess};
    }
    Fluxes flux(const Vec& y) const {
        if(int(y.size())!=size()) throw std::runtime_error("ALE state size");
        const double h=y[2*n];
        if(!(h>0)) throw std::runtime_error("nonpositive ALE height");
        Vec d(n), eface(n+1);
        Fluxes f{Vec(n),Vec(n),Vec(n+1),Vec(n+1),Vec(n+1),0,0,p.evaporation};
        for(int i=0;i<n;++i) {
            if(!(y[i]>0 && y[i]<1)) throw std::runtime_error("invalid porosity");
            d[i]=std::sqrt(y[i])*p.diffusion;
            f.volume[i]=h*widths[i];
            f.reaction[i]=reactionRate(y[i],y[n+i]);
        }
        auto t=top(y[2*n-1],y[n-1],d[n-1],h,f.evaporation);
        if(gasFilm) {
            bool done=false;
            for(int k=0;k<40;++k) {
                const double next=gasFilm(t.first);
                if(!(next>=0)) throw std::runtime_error("condensation outside scope");
                const double difference=std::abs(next-f.evaporation);
                f.evaporation=next; t=top(y[2*n-1],y[n-1],d[n-1],h,next);
                if(difference<1e-14) { done=true; break; }
            }
            if(!done) throw std::runtime_error("top film fixed point failed");
        }
        f.topC=t.first; f.speed=t.second;
        f.liquid[n]=(p.rho*p.birth+p.crystal*(1-p.birth))*f.speed+f.evaporation;
        for(int i=n-1;i>=0;--i)
            f.liquid[i]=f.liquid[i+1]+(1-p.rho/p.crystal)*f.reaction[i]*f.volume[i];
        eface[0]=y[0]; eface[n]=p.birth;
        f.mesh[0]=0; f.mesh[n]=-f.speed*p.birth;
        for(int i=1;i<n;++i) {
            eface[i]=f.speed>=0?y[i]:y[i-1];
            f.mesh[i]=-faces[i]*f.speed*eface[i];
            const double resistance=h*((faces[i]-centers[i-1])/(y[i-1]*d[i-1])+
                                          (centers[i]-faces[i])/(y[i]*d[i]));
            const double relative=f.liquid[i]+p.rho*f.mesh[i];
            f.salt[i]=relative*(relative>=0?y[n+i-1]:y[n+i])-p.rho*(y[n+i]-y[n+i-1])/resistance;
        }
        f.salt[0]=p.bottomC*f.liquid[0]-p.rho*y[0]*d[0]*(y[n]-p.bottomC)/(h*centers[0]);
        f.salt[n]=p.crystal*(1-p.birth)*f.speed;
        return f;
    }
    Vec rhs(const Vec& y) const {
        const auto f=flux(y); Vec v(size(),0);
        double internal=0,positive=0,dissolution=0;
        for(int i=0;i<n;++i) {
            const double e=f.mesh[i]-f.mesh[i+1]-f.reaction[i]*f.volume[i]/p.crystal;
            const double q=f.salt[i]-f.salt[i+1]-f.reaction[i]*f.volume[i];
            v[i]=(e-y[i]*widths[i]*f.speed)/f.volume[i];
            v[n+i]=(q-p.rho*y[n+i]*e)/(p.rho*y[i]*f.volume[i]);
            internal+=f.reaction[i]*f.volume[i];
            positive+=std::max(0.,f.reaction[i])*f.volume[i];
            dissolution+=std::max(0.,-f.reaction[i])*f.volume[i];
        }
        v[2*n]=f.speed;
        v[inv()]=f.salt[n]; v[inv()+1]=internal; v[inv()+2]=f.salt[0];
        v[inv()+3]=f.liquid[0]; v[inv()+4]=f.evaporation;
        if(p.splitReactionLedgers) {
            v[positiveInv()]=positive;
            v[dissolutionInv()]=dissolution;
        }
        return v;
    }
    struct Inventory { double salt=0,water=0,solid=0,pore=0; };
    Inventory inventory(const Vec& y) const {
        Inventory m;
        for(int i=0;i<n;++i) {
            const double volume=y[2*n]*widths[i], pore=y[i]*volume;
            m.pore+=pore; m.salt+=p.rho*pore*y[n+i];
            m.water+=p.rho*pore*(1-y[n+i]); m.solid+=p.crystal*(volume-pore);
        }
        return m;
    }
    void validate(const Vec& y) const {
        for(int i=0;i<n;++i)
            if(!(y[i]>.02 && y[i]<.95 && y[n+i]>=0 && y[n+i]<.95))
                throw std::runtime_error("ALE accepted state outside M3 physical bounds");
    }
};
}
