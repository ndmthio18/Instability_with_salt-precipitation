#pragma once
#include <cmath>
namespace audit28 {
struct Dual {
    double v=0,d=0;
    Dual()=default; Dual(double value,double derivative=0):v(value),d(derivative){}
    Dual& operator+=(Dual b){v+=b.v;d+=b.d;return *this;}
};
inline Dual operator+(Dual a,Dual b){return {a.v+b.v,a.d+b.d};}
inline Dual operator-(Dual a,Dual b){return {a.v-b.v,a.d-b.d};}
inline Dual operator-(Dual a){return {-a.v,-a.d};}
inline Dual operator*(Dual a,Dual b){return {a.v*b.v,a.d*b.v+a.v*b.d};}
inline Dual operator/(Dual a,Dual b){return {a.v/b.v,(a.d-(a.v/b.v)*b.d)/b.v};}
inline double value(double a){return a;}
inline double value(Dual a){return a.v;}
inline double root(double a){return std::sqrt(a);}
inline Dual root(Dual a){double r=std::sqrt(a.v);return {r,a.d/(2*r)};}
inline double power(double a,double b){return std::pow(a,b);}
inline Dual power(Dual a,double b){return {std::pow(a.v,b),b*std::pow(a.v,b-1)*a.d};}
template<class T>T positive(T a){return value(a)>0?a:T(0);}
template<class T>T bounded(T a,T upper){return value(a)<0?T(0):(value(a)>value(upper)?upper:a);}
}
