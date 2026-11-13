#pragma once
#include "two_front.hh"
#include <iostream>
namespace audit28 {
// Alternative physical initial condition: dc/dt=0 at prescribed initial solid geometry.
// This does not replace the uniform initial condition or remove its recorded transient.
inline Vec consistentInitial(const Model& m,Vec y){
 const int N=m.n+m.l;auto f=m.flux(y);Vec mass(N);std::vector<long double> e(N,0);
 for(int i=0;i<m.n;++i)mass[i]=m.p.rho*(f.vol[i]-y[i]);
 for(int i=0;i<m.l;++i)mass[m.n+i]=m.p.rho*f.vl[i];
 auto residual=[&](const std::vector<Dual>& ec){std::vector<Dual> state(m.size());
  for(int i=0;i<m.size();++i)state[i]=y[i];
  for(int i=0;i<N;++i)state[m.n+i]=mass[i]*ec[i];
  auto d=m.rhs(state,0);std::vector<Dual> r(N);
  for(int i=0;i<N;++i){Dual de=i<m.n?(d[m.topIndex()]-d[m.bottomIndex()])*m.widths[i]-d[i]:d[m.bottomIndex()]*m.lw[i-m.n];
   r[i]=d[m.n+i]-m.p.rho*ec[i]*de;}
  return r;};
 for(int it=0;it<12;++it){std::vector<Dual> ec(N);for(int i=0;i<N;++i)ec[i]=e[i];auto r=residual(ec);long double norm=0;
  for(auto z:r)norm=std::max(norm,std::abs(z.v));
  if(norm<1e-15){for(int i=0;i<N;++i)y[m.n+i]=mass[i]*e[i];m.validate(y);std::cout<<"CONSISTENT_INIT iterations="<<it<<" max_residual="<<(double)norm<<'\n';return y;}
  std::vector<std::vector<long double>> A(N,std::vector<long double>(N+1));
  for(int j=0;j<N;++j){ec[j].d=1;auto dr=residual(ec);ec[j].d=0;for(int i=0;i<N;++i)A[i][j]=dr[i].d;}
  for(int i=0;i<N;++i)A[i][N]=-r[i].v;
  for(int k=0;k<N;++k){int pivot=k;for(int i=k+1;i<N;++i)if(std::abs(A[i][k])>std::abs(A[pivot][k]))pivot=i;
   std::swap(A[k],A[pivot]);if(std::abs(A[k][k])<1e-30)throw std::runtime_error("initialization singular Jacobian");
   for(int i=k+1;i<N;++i){long double a=A[i][k]/A[k][k];for(int j=k+1;j<=N;++j)A[i][j]-=a*A[k][j];}}
  std::vector<long double> delta(N);for(int i=N-1;i>=0;--i){long double b=A[i][N];for(int j=i+1;j<N;++j)b-=A[i][j]*delta[j];delta[i]=b/A[i][i];e[i]+=delta[i];}
 }
 throw std::runtime_error("consistent initialization failed");
}
}
