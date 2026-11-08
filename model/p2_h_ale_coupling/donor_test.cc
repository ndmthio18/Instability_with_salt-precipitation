#include "donor.hh"
#include <cassert>
#include <iostream>
int main() {
    using namespace p2h;
    std::vector<Donor> d={{3,.003,.004,2e-6,.4,.999,1200,.25},
        {1,0,.001,1e-6,.4,.999,1190,.2},{2,0,.001,1e-6,.4,.999,1190,.2},
        {4,.008,.009,4e-6,.4,.98,1180,.1}};
    auto a=allocate(d,1e-5,.005,.99);
    assert(a.sufficient && a.donors.size()==3);
    assert(a.donors[0].extracted==a.donors[1].extracted);
    assert(a.donors[2].extracted>0 && a.donors[2].extracted<a.donors[2].available);
    assert(std::abs(a.salt+a.water-1e-5)<1e-19);
    assert(a.usedDepth==.004);
    auto b=allocate(d,1e-5,.002,.99); assert(!b.sufficient && b.donors.size()==2);
    auto c=allocate(d,1e-5,.010,.99); assert(c.donors.back().available==0);
    for(const auto& q:a.donors) assert(q.saturation-q.extracted/(q.density*q.porosity*q.volume)>=.99);
    // Reordering input and equal-depth IDs cannot create a lateral bias.
    std::reverse(d.begin(),d.end()); auto e=allocate(d,1e-5,.005,.99);
    assert(e.salt==a.salt && e.water==a.water);
    assert(d.front().extracted==0); // caller state is immutable even on refusal
    std::cout<<"PASS: graded volumes, depth cutoff, partial layer, mixing, dry donor, immutable trials\n";
}
