#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace p2h {
struct Donor {
    std::size_t id;
    double depth, deepest, volume, porosity, saturation, density, concentration;
    double available=0, extracted=0;
};
struct Allocation {
    std::vector<Donor> donors;
    double required=0, available=0, salt=0, water=0, usedDepth=0;
    bool sufficient=false;
};
// Pure trial allocation: no state mutation. Equal-depth donors share a partial
// layer proportionally, preserving a horizontally homogeneous one-cell column.
inline Allocation allocate(std::vector<Donor> donors,double required,double maxDepth,double minSw) {
    if(!(required>0 && maxDepth>0 && minSw>0 && minSw<1))
        throw std::runtime_error("invalid donor controls");
    std::sort(donors.begin(),donors.end(),[](const Donor& a,const Donor& b){
        return a.depth==b.depth ? a.id<b.id : a.depth<b.depth;
    });
    Allocation result; result.required=required;
    for(auto d:donors) {
        if(d.deepest>maxDepth+1e-14) continue;
        // A 1e-13 upward roundoff guard keeps reconstructed Sw >= minSw.
        d.available=d.density*d.porosity*d.volume*std::max(0.,d.saturation-(minSw+1e-13));
        result.available+=d.available; result.donors.push_back(d);
    }
    double remaining=required;
    for(std::size_t a=0;a<result.donors.size() && remaining>0;) {
        std::size_t b=a+1;
        while(b<result.donors.size() && std::abs(result.donors[b].depth-result.donors[a].depth)<1e-14) ++b;
        double capacity=0;
        for(auto i=a;i<b;++i) capacity+=result.donors[i].available;
        const double take=std::min(capacity,remaining);
        if(capacity>0) for(auto i=a;i<b;++i) {
            auto& d=result.donors[i]; d.extracted=take*d.available/capacity;
            result.salt+=d.extracted*d.concentration;
            result.water+=d.extracted*(1-d.concentration);
            if(d.extracted>0) result.usedDepth=std::max(result.usedDepth,d.deepest);
        }
        remaining-=take; a=b;
    }
    // Never pad a deficit with invented liquid; this only tolerates arithmetic ULPs.
    result.sufficient=remaining<=required*1e-14;
    return result;
}
}
