// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux-Course contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXERCISE_RICHARDS_SALT_MODEL_HH
#define DUMUX_EXERCISE_RICHARDS_SALT_MODEL_HH

// In the file `model.hh`, we define the volume variables, the model equations, and
// set all default model properties.

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/volumevariables.hh>
#include <dumux/discretization/method.hh>

// The property tag is simply an empty struct with the name `RichardsSaltModel`
namespace Dumux::Properties::TTag {
//! The nonlinear diffusion model tag that we can specialize properties for
struct RichardsSaltModel {};
} // end namespace Dumux::Properties::TTag

namespace Dumux {

/*!
 * \ingroup RichardsSaltModel
 * \brief Volume averaged quantities required by the nonlinear diffusion model.
 * This contains the quantities which are related to a control volume
 */
template <class Traits>
class RichardsSaltVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)

{
    // Set the primary variables and some derived properties
    priVars_ = elemSol[scv.indexInElement()];
    sw_ = problem.sw(priVars_[Indices::pressureIdx]);
    rhoW_ = problem.rhoW(priVars_[Indices::moleFractionIdx]); 
}   

    Scalar pw() const
    { return priVars_[Indices::pressureIdx]; }

    Scalar moleFraction() const
    { return priVars_[Indices::moleFractionIdx]; }

    Scalar sw() const
    { return sw_; }

    Scalar rhoW() const
    { return rhoW_; }


    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
    Scalar sw_;
    Scalar rhoW_; 
    Scalar vn_; 
};

/*!
 * \ingroup RichardsSaltModel
 * \brief Element-wise calculation of the residual and its derivatives
 *        for a nonlinear diffusion problem.
 */
template<class TypeTag>
class RichardsSaltModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        return normalizedStorage(problem, volVars);
    }

    // Shared by the residual and Paper 2 diagnostics: no molar-density factor.
    static NumEqVector normalizedStorage(const Problem& problem, const VolumeVariables& volVars)
    {
        NumEqVector storage;
        storage[Indices::conti0EqIdx] = problem.porosity() * volVars.sw();
        storage[Indices::conti1EqIdx] = problem.porosity() * volVars.sw() * volVars.moleFraction();
        return storage;
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto f = transportFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        NumEqVector result;
        result[Indices::conti0EqIdx] = f.water;
        result[Indices::conti1EqIdx] = f.total;
        return result;
    }

    struct TransportFlux {
        Scalar water, advective, diffusive, total, upwindMoleFraction;
    };

    // Same gradients, arithmetic mean D/rho, upstream mobility and concentration
    // as the original residual. Values already include SCVF area.
    static TransportFlux transportFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "This local residual is hard-coded to control-volume finite element schemes");

        // Compute ∇c at the integration point of this sub control volume face.
        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradMoleFraction(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradP(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            gradMoleFraction.axpy(
                volVars.moleFraction(),
                fluxVarCache.gradN(scv.indexInElement())
            );

            gradP.axpy(
                volVars.pw(),
                fluxVarCache.gradN(scv.indexInElement())
            );
        }

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto K = problem.permeability();
        const auto insideD = problem.effectiveDiffusionCoefficient(insideVolVars.sw());
        const auto outsideD = problem.effectiveDiffusionCoefficient(outsideVolVars.sw());
        const auto D = 0.5*(insideD + outsideD);

        const auto insideRhoW = problem.rhoW(insideVolVars.moleFraction());
        const auto outsideRhoW = problem.rhoW(outsideVolVars.moleFraction());
        const auto rhoW = 0.5*(insideRhoW + outsideRhoW);

        static const bool enableGravity = getParam<bool>("Problem.EnableGravity");
        GlobalPosition gravity(0.0);
        if (enableGravity)
            gravity[dimWorld-1] = -problem.gravitationalAcceleration();

        auto vn = -vtmv(scvf.unitOuterNormal(), K, gradP - rhoW*gravity);
        const auto upVolVar = (vn > 0) ? insideVolVars : outsideVolVars;
        vn *= problem.krw(upVolVar.sw()) / problem.viscosity();

        // Compute the flux with `vtmv` (vector transposed times matrix times vector) or -n^T D ∇c A.
        // The diffusion coefficient comes from the `problem` (see Part 2 of the example).
        const auto diffusiveDensity = -vtmv(scvf.unitOuterNormal(), D, gradMoleFraction);
        return {scvf.area() * vn,
                scvf.area() * vn * upVolVar.moleFraction(),
                scvf.area() * diffusiveDensity,
                scvf.area() * (vn * upVolVar.moleFraction() + diffusiveDensity),
                upVolVar.moleFraction()};
    }

};
} // end namespace Dumux

// The model properties
namespace Dumux::Properties {

// The type of the local residual is the class defined above
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::RichardsSaltModel>
{ using type = RichardsSaltModelLocalResidual<TypeTag>; };

// The default scalar type is double
// we compute with double precision floating point numbers
template<class TypeTag>
struct Scalar<TypeTag, TTag::RichardsSaltModel>
{ using type = double; };

// The model traits specify some information about our equation system.
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::RichardsSaltModel>
{
    struct type
    {
        struct Indices
        {
            static constexpr int pressureIdx = 0;
            static constexpr int moleFractionIdx = 1;
            static constexpr int conti0EqIdx = 0;
            static constexpr int conti1EqIdx = 1;
        };

        static constexpr int numEq() { return 2; }
    };
};

// The primary variable vector has entries of type `Scalar` and is
// as large as the number of equations
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::RichardsSaltModel>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

// The volume variables that are used in our model
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::RichardsSaltModel>
{
    struct Traits
    {
        using ModelTraits
            = GetPropType<TypeTag, Properties::ModelTraits>;
        using PrimaryVariables
            = GetPropType<TypeTag, Properties::PrimaryVariables>;
    };

    using type = RichardsSaltVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties

#endif
