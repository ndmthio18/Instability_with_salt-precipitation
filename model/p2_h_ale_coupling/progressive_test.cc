#include "progressive.hh"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

namespace {
bool close(double a,double b,double relative=2e-14) {
    return std::abs(a-b)<=relative*std::max({1.0,std::abs(a),std::abs(b)});
}
}

int main() {
    constexpr double j=2.7e-6,dt=.25,eps=.43,rhoc=2165,c=.264940870,rho=1200;
    const auto one=p2h::progressiveIncrement(j,dt,eps,rhoc,c,rho);
    assert(close(one.solid,j*dt));
    assert(close(one.poreVolume,one.solid*eps/(rhoc*(1-eps))));
    assert(close(one.poreSalt+one.poreWater,rho*one.poreVolume));
    assert(close(one.poreSalt/(one.poreSalt+one.poreWater),c));

    p2h::PrecursorInventory many;
    for(int i=0;i<20;++i) p2h::acceptProgressiveIncrement(many,one);
    const auto all=p2h::progressiveIncrement(j,20*dt,eps,rhoc,c,rho);
    assert(close(many.poreSalt,all.poreSalt));
    assert(close(many.poreWater,all.poreWater));
    assert(close(many.requiredPoreVolume,all.poreVolume));
    assert(close(many.filledPoreVolume,all.poreVolume));
    assert(close(p2h::precursorVolumeResidual(many),0));

    // Rejected hydraulic trial: work on a copy and discard it.  No inventory
    // compensation, duplicate transfer, or waiting-time fill is introduced.
    const auto accepted=many;
    auto rejected=many;
    p2h::acceptProgressiveIncrement(rejected,one);
    rejected=accepted;
    assert(close(rejected.poreSalt,many.poreSalt));
    assert(close(rejected.filledPoreVolume,many.filledPoreVolume));

    // The kinematic demand is independent of an arbitrary top BOX volume.
    for(double topBoxVolume: {1e-9,1e-12,1e-15}) {
        (void)topBoxVolume;
        const auto value=p2h::progressiveIncrement(j,dt,eps,rhoc,c,rho);
        assert(close(value.poreVolume,one.poreVolume));
    }

    // Activation is a one-time conservative transfer from the accepted ledger.
    const double transferredSalt=many.poreSalt;
    const double transferredWater=many.poreWater;
    assert(close(transferredSalt+transferredWater,rho*many.filledPoreVolume));
    assert(close(p2h::precursorConcentration(many),c));
    std::cout<<"progressive precursor tests passed\n";
}
