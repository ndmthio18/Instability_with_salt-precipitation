#ifndef P1_COMPOSITIONAL_PROPERTIES_HH
#define P1_COMPOSITIONAL_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/porousmediumflow/2pnc/model.hh>

#include "fluidsystem.hh"
#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct P1Compositional { using InheritsFrom = std::tuple<TwoPNC>; };
struct P1Compositional2D { using InheritsFrom = std::tuple<P1Compositional, CCTpfaModel>; };
struct P1Compositional1D { using InheritsFrom = std::tuple<P1Compositional, CCTpfaModel>; };
}

template<class TypeTag>
struct Grid<TypeTag, TTag::P1Compositional2D>
{
    using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2>>;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::P1Compositional1D>
{
    using type = Dune::YaspGrid<1, Dune::TensorProductCoordinates<double, 1>>;
};

template<class TypeTag>
struct Problem<TypeTag, TTag::P1Compositional>
{
    using type = P1Problem<TypeTag>;
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::P1Compositional>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = P1FluidSystem<Scalar>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::P1Compositional>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = P1SpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::P1Compositional>
{
    static constexpr bool value = true;
};

template<class TypeTag>
struct Formulation<TypeTag, TTag::P1Compositional>
{
    static constexpr auto value = TwoPFormulation::p0s1;
};

} // namespace Dumux::Properties

#endif
