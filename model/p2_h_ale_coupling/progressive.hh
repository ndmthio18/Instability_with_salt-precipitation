#pragma once
#include <cmath>
#include <stdexcept>

namespace p2h {

// Per-unit-soil-surface inventories.  Volumes therefore have units m3/m2=m.
// This helper contains no soil-control-volume quantity and no time-step-based
// availability rule: its continuum limit is set by the physical precipitation
// flux supplied by the surface complementarity condition.
struct PrecursorInventory {
    double poreSalt=0;          // kg/m2
    double poreWater=0;         // kg/m2
    double requiredPoreVolume=0;// m3/m2
    double filledPoreVolume=0;  // m3/m2
};

struct ProgressiveIncrement {
    double solid=0;             // kg/m2
    double poreSalt=0;          // kg/m2
    double poreWater=0;         // kg/m2
    double poreVolume=0;        // m3/m2
};

inline ProgressiveIncrement progressiveIncrement(double precipitationFlux,
                                                  double dt,
                                                  double birthPorosity,
                                                  double crystalDensity,
                                                  double liquidConcentration,
                                                  double liquidDensity) {
    if(!(precipitationFlux>=0 && dt>0 && birthPorosity>0 && birthPorosity<1 &&
         crystalDensity>0 && liquidConcentration>=0 && liquidConcentration<1 &&
         liquidDensity>0))
        throw std::runtime_error("invalid progressive precursor input");
    ProgressiveIncrement value;
    value.solid=precipitationFlux*dt;
    value.poreVolume=value.solid*birthPorosity/
                     (crystalDensity*(1-birthPorosity));
    const double liquid=liquidDensity*value.poreVolume;
    value.poreSalt=liquidConcentration*liquid;
    value.poreWater=(1-liquidConcentration)*liquid;
    return value;
}

inline void acceptProgressiveIncrement(PrecursorInventory& inventory,
                                       const ProgressiveIncrement& value) {
    inventory.requiredPoreVolume+=value.poreVolume;
    inventory.filledPoreVolume+=value.poreVolume;
    inventory.poreSalt+=value.poreSalt;
    inventory.poreWater+=value.poreWater;
}

inline double precursorConcentration(const PrecursorInventory& inventory) {
    const double liquid=inventory.poreSalt+inventory.poreWater;
    if(!(liquid>0)) throw std::runtime_error("empty precursor pore liquid");
    return inventory.poreSalt/liquid;
}

inline double precursorVolumeResidual(const PrecursorInventory& inventory) {
    return inventory.requiredPoreVolume-inventory.filledPoreVolume;
}

} // namespace p2h
