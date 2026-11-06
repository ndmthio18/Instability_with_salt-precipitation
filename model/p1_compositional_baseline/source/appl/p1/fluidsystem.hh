#ifndef P1_COMPOSITIONAL_FLUIDSYSTEM_HH
#define P1_COMPOSITIONAL_FLUIDSYSTEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidsystems/brineair.hh>

namespace Dumux {

template<class Scalar>
class P1FluidSystem
: public FluidSystems::BrineAir<Scalar, Components::H2O<Scalar>>
{
    using Parent = FluidSystems::BrineAir<Scalar, Components::H2O<Scalar>>;

public:
    template<class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        const bool liquidSaltPair = phaseIdx == Parent::liquidPhaseIdx
            && ((compIIdx == Parent::H2OIdx && compJIdx == Parent::NaClIdx)
                || (compJIdx == Parent::H2OIdx && compIIdx == Parent::NaClIdx));
        if (liquidSaltPair)
            return getParam<Scalar>("Problem.SaltDiffusivity", 1.5e-9);
        return Parent::binaryDiffusionCoefficient(fluidState, phaseIdx, compIIdx, compJIdx);
    }

    template<class FluidState, class ParameterCache>
    static Scalar binaryDiffusionCoefficient(const FluidState& fluidState,
                                             const ParameterCache&,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
    {
        return binaryDiffusionCoefficient(fluidState, phaseIdx, compIIdx, compJIdx);
    }
};

} // namespace Dumux

#endif
