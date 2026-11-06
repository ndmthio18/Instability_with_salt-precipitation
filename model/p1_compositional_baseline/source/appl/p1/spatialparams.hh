#ifndef P1_COMPOSITIONAL_SPATIALPARAMS_HH
#define P1_COMPOSITIONAL_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class P1SpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar,
                                           P1SpatialParams<GridGeometry, Scalar>>
{
    using ThisType = P1SpatialParams<GridGeometry, Scalar>;
    using Parent = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    explicit P1SpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : Parent(gridGeometry)
    , permeability_(getParam<Scalar>("Problem.Permeability"))
    , porosity_(getParam<Scalar>("Problem.Porosity"))
    , temperature_(getParam<Scalar>("Problem.Temperature", 293.15))
    , gravityMagnitude_(getParam<Scalar>("Problem.GravitationalAcceleration", 9.8))
    , enableGravity_(getParam<bool>("Problem.EnableGravity", true))
    , pcKrSwCurve_("SpatialParams")
    {}

    PermeabilityType permeabilityAtPos(const GlobalPosition&) const
    { return permeability_; }

    Scalar porosityAtPos(const GlobalPosition&) const
    { return porosity_; }

    Scalar temperatureAtPos(const GlobalPosition&) const
    { return temperature_; }

    GlobalPosition gravity(const GlobalPosition&) const
    {
        GlobalPosition g(0.0);
        if (enableGravity_)
            g[GridView::dimensionworld - 1] = -gravityMagnitude_;
        return g;
    }

    auto fluidMatrixInteractionAtPos(const GlobalPosition&) const
    { return makeFluidMatrixInteraction(pcKrSwCurve_); }

    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition&) const
    { return FluidSystem::liquidPhaseIdx; }

private:
    Scalar permeability_;
    Scalar porosity_;
    Scalar temperature_;
    Scalar gravityMagnitude_;
    bool enableGravity_;
    const PcKrSwCurve pcKrSwCurve_;
};

} // namespace Dumux

#endif
