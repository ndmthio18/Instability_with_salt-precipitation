// WP8 equal-height constrained Newton. Original FV residual/constitutive laws retained.
// Serial one-cell extrusion only. Failed Newton steps throw; no long-run admission here.
#pragma once
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

inline std::vector<double> wp8DenseSolve(std::vector<double> a, std::vector<double> b)
{
    const auto n=b.size();
    for (std::size_t i=0;i<n;++i) {
        double scale=0;for(std::size_t j=0;j<n;++j)scale=std::max(scale,std::abs(a[i*n+j]));
        if (!(scale>0))throw std::runtime_error("WP8 zero Jacobian row");
        for(std::size_t j=0;j<n;++j)a[i*n+j]/=scale;b[i]/=scale;
    }
    for(std::size_t k=0;k<n;++k) {
        auto pivot=k;for(auto i=k+1;i<n;++i)if(std::abs(a[i*n+k])>std::abs(a[pivot*n+k]))pivot=i;
        if(std::abs(a[pivot*n+k])<1e-30)throw std::runtime_error("WP8 singular reduced Jacobian");
        if(pivot!=k){for(std::size_t j=k;j<n;++j)std::swap(a[k*n+j],a[pivot*n+j]);std::swap(b[k],b[pivot]);}
        for(auto i=k+1;i<n;++i){const double q=a[i*n+k]/a[k*n+k];a[i*n+k]=0;for(auto j=k+1;j<n;++j)a[i*n+j]-=q*a[k*n+j];b[i]-=q*b[k];}
    }
    for(std::size_t ii=n;ii-->0;){for(auto j=ii+1;j<n;++j)b[ii]-=a[ii*n+j]*b[j];b[ii]/=a[ii*n+ii];if(!std::isfinite(b[ii]))throw std::runtime_error("WP8 nonfinite correction");}
    return b;
}

template<class Assembler, class Solution, class Geometry>
void paper2ConstrainedStep(Assembler& assembler, Solution& x, const Geometry& geometry)
{
    if(geometry.gridView().comm().size()!=1)throw std::runtime_error("WP8 constraint requires one MPI rank");
    std::map<double,std::size_t> levels;std::vector<std::size_t> group(x.size());
    for(const auto& v:vertices(geometry.gridView()))levels.emplace(v.geometry().center()[1],0);
    std::size_t l=0;for(auto& kv:levels)kv.second=l++;
    for(const auto& v:vertices(geometry.gridView()))group[geometry.vertexMapper().index(v)]=levels.at(v.geometry().center()[1]);
    std::vector<std::size_t> first(l,x.size()),count(l,0);
    for(std::size_t i=0;i<x.size();++i){++count[group[i]];if(first[group[i]]==x.size())first[group[i]]=i;else for(int c=0;c<2;++c)if(x[i][c]!=x[first[group[i]]][c])throw std::runtime_error("WP8 nonuniform input; post-step projection forbidden");}
    for(auto n:count)if(n!=2)throw std::runtime_error("WP8 requires exactly one horizontal cell");
    const std::size_t n=2*l;std::vector<double> u(n);
    for(std::size_t j=0;j<l;++j)for(int c=0;c<2;++c)u[2*j+c]=x[first[j]][c];
    double lastShift=1;
    for(int iter=0;iter<=20;++iter){
        // Every accepted/Newton state is an exact lift; no horizontal averaging.
        for(std::size_t i=0;i<x.size();++i)for(int c=0;c<2;++c)x[i][c]=u[2*group[i]+c];
        assembler.gridVariables().update(x);assembler.assembleJacobianAndResidual(x);
        const auto& f=assembler.residual();const auto& jac=assembler.jacobian();
        std::vector<double> r(n,0),a(n*n,0);
        for(std::size_t i=0;i<x.size();++i){
            for(int c=0;c<2;++c)r[2*group[i]+c]+=f[i][c];
            for(auto it=jac[i].begin();it!=jac[i].end();++it)for(int c=0;c<2;++c)for(int d=0;d<2;++d)a[(2*group[i]+c)*n+2*group[it.index()]+d]+=(*it)[c][d];
        }
        // Native FVAssembler retains its Dirichlet row residuals. Two identical
        // copies sum to 2*(u-u_D); this scaling leaves the constraint unchanged.
        double norm=0,waterNorm=0,saltNorm=0;for(std::size_t j=0;j<n;++j){const auto v=std::abs(r[j]);if(!std::isfinite(v))throw std::runtime_error("WP8 nonfinite residual");norm=std::max(norm,v);if(j%2)saltNorm=std::max(saltNorm,v);else waterNorm=std::max(waterNorm,v);}
        // r01/r02 exposed a water residual rounding floor near 1e-18 with
        // relative shifts below 2e-15. These are internal algebraic tolerances;
        // unchanged independent salt/water budgets and 1e-14 uniformity gate
        // remain mandatory. Never use one raw residual cutoff for both scales.
        std::cout<<"WP8_ITER "<<iter<<" water_residual="<<std::scientific<<waterNorm<<" salt_residual="<<saltNorm<<" shift="<<lastShift<<"\n";
        if(waterNorm<=2e-18 && saltNorm<=2e-20 && lastShift<=1e-10){std::cout<<"WP8_CONSTRAINT iterations="<<iter<<" residual="<<norm<<" shift="<<lastShift<<"\n";return;}
        if(iter==20)throw std::runtime_error("WP8 reduced Newton failed; no accepted step");
        auto du=wp8DenseSolve(a,r);lastShift=0;
        for(std::size_t j=0;j<n;++j){lastShift=std::max(lastShift,std::abs(du[j])/std::max(std::abs(u[j]),j%2 ? .0036 : 1.0));u[j]-=du[j];}
    }
}
