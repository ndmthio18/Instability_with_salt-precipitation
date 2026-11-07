#pragma once
#include <cmath>
#include <stdexcept>
#include <algorithm>
namespace p2g {
inline double precipitate(double salt, double water, double ceq) {
    if (!(salt>=0 && water>0 && ceq>0 && ceq<1))
        throw std::runtime_error("invalid precipitation inventory");
    return std::max(0.0, salt-ceq/(1-ceq)*water);
}
inline double velocity(double h,double j,double ceq, bool finiteReaction=true) {
    if (!(h>0 && j>=0)) throw std::runtime_error("invalid growth reference");
    const double eps=.43, kr=2.3e-3, beta=2, D=std::sqrt(eps)*1.3e-9;
    const double av=36*(1-eps)/25e-6;
    const double dat=beta*(1-eps)*(1-ceq)*h*kr/(eps*D);
    const double da=h*h*(1-ceq)*av*kr/(eps*D);
    const double factor=finiteReaction ? dat*std::tanh(std::sqrt(da))/std::sqrt(da) : dat;
    return ceq*j/(2165*(1-eps)*(1-ceq))*factor/(1+factor);
}
}
