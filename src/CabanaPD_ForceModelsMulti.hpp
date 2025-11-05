/****************************************************************************
 * Copyright (c) 2022 by Oak Ridge National Laboratory                      *
 * All rights reserved.                                                     *
 *                                                                          *
 * This file is part of CabanaPD. CabanaPD is distributed under a           *
 * BSD 3-clause license. For the licensing terms see the LICENSE file in    *
 * the top-level directory.                                                 *
 *                                                                          *
 * SPDX-License-Identifier: BSD-3-Clause                                    *
 ****************************************************************************/

#ifndef FORCE_MODELS_MULTI_H
#define FORCE_MODELS_MULTI_H

#include <CabanaPD_Output.hpp>

namespace CabanaPD
{

namespace Impl
{
template <std::size_t current, typename FunctorType, typename ParameterPackType,
          typename... ARGS,
          std::enable_if_t<current<ParameterPackType::size - 1, int> = 0>
              KOKKOS_INLINE_FUNCTION auto
                  recursion_functor_for_index_in_pack_with_args(
                      FunctorType const& functor, std::size_t index,
                      ParameterPackType& parameterPack, ARGS... args )
{
    if ( index == current )
    {
        return functor( Cabana::get<current>( parameterPack ), args... );
    }
    else
    {
        return recursion_functor_for_index_in_pack_with_args<current + 1>(
            functor, index, parameterPack, args... );
    }
}

template <std::size_t current, typename FunctorType, typename ParameterPackType,
          typename... ARGS,
          std::enable_if_t<current == ParameterPackType::size - 1, int> = 0>
KOKKOS_INLINE_FUNCTION auto recursion_functor_for_index_in_pack_with_args(
    FunctorType const& functor, std::size_t index,
    ParameterPackType& parameterPack, ARGS... args )
{
    if ( index == current )
    {
        return functor( Cabana::get<current>( parameterPack ), args... );
    }
    else
    {
        Kokkos::abort( "Requested index not contained in ParameterPack" );
        return functor( Cabana::get<current>( parameterPack ), args... );
    }
}

template <typename FunctorType, typename ParameterPackType, typename... ARGS>
KOKKOS_INLINE_FUNCTION auto run_functor_for_index_in_pack_with_args(
    FunctorType const& functor, std::size_t index,
    ParameterPackType& parameterPack, ARGS... args )
{
    return recursion_functor_for_index_in_pack_with_args<0>(
        functor, index, parameterPack, args... );
}

struct IdentityFunctor
{
    template <typename Model, typename... ARGS>
    KOKKOS_FUNCTION auto operator()( Model& model, ARGS... args ) const
    {
        return model( args... );
    }
};

template <typename IfType, typename ElseType, bool Condition>
struct IfElseType
{
    static_assert( true, "We shouuld never end up here" );
};

template <typename IfType, typename ElseType>
struct IfElseType<IfType, ElseType, true>
{
    using type = IfType;
};

template <typename IfType, typename ElseType>
struct IfElseType<IfType, ElseType, false>
{
    using type = ElseType;
};

template <typename ParameterPackType, typename Sequence>
struct CheckTemperatureDependence
{
};

template <typename ParameterPackType, std::size_t... Indices>
struct CheckTemperatureDependence<ParameterPackType,
                                  std::index_sequence<Indices...>>
{
    using type = typename IfElseType<
        TemperatureDependent, TemperatureIndependent,
        std::disjunction_v<std::is_same<
            typename ParameterPackType::value_type<Indices>::thermal_type,
            TemperatureDependent>...>>::type;
};

template <typename ParameterPackType, typename Sequence>
struct CheckFractureModel
{
};

template <typename ParameterPackType, std::size_t... Indices>
struct CheckFractureModel<ParameterPackType, std::index_sequence<Indices...>>
{
    using type = typename IfElseType<
        Fracture, NoFracture,
        std::disjunction_v<std::is_same<
            typename ParameterPackType::value_type<Indices>::fracture_type,
            Fracture>...>>::type;
};

template <typename FractureType, typename ParameterPackType, unsigned Index,
          bool Condition>
struct FirstModelWithFractureTypeImpl
{
};

template <typename FractureType, typename ParameterPackType, unsigned Index>
struct FirstModelWithFractureTypeImpl<FractureType, ParameterPackType, Index,
                                      true>
{
    using type = typename ParameterPackType::value_type<Index>;
};

template <typename FractureType, typename ParameterPackType, unsigned Index>
struct FirstModelWithFractureTypeImpl<FractureType, ParameterPackType, Index,
                                      false>
{
    using type = typename FirstModelWithFractureTypeImpl<
        FractureType, ParameterPackType, Index + 1,
        std::is_same_v<
            typename ParameterPackType::value_type<Index>::fracture_type,
            FractureType>>::type;
};

template <typename FractureType, typename ParameterPackType>
struct FirstModelWithFractureType
{
    using type = typename FirstModelWithFractureTypeImpl<
        FractureType, ParameterPackType, 0,
        std::is_same_v<typename ParameterPackType::value_type<0>::fracture_type,
                       FractureType>>::type;
};

// Wrap multiple models in a single object.
template <typename MaterialType, typename ParameterPackType, typename Sequence>
struct ForceModelsImpl
{
};

template <typename MaterialType, typename ParameterPackType,
          std::size_t... Indices>
struct ForceModelsImpl<MaterialType, ParameterPackType,
                       std::index_sequence<Indices...>>
{
    using material_type = MultiMaterial;

    static_assert(
        (std::conjunction_v<std::is_same<
             typename ParameterPackType::value_type<Indices>::base_model,
             typename ParameterPackType::value_type<0>::base_model>...>),
        "All models need the same base model" );
    using base_model = typename ParameterPackType::value_type<0>::base_model;

    using fracture_type = typename CheckFractureModel<
        ParameterPackType,
        std::make_index_sequence<ParameterPackType::size - 1>>::type;
    using thermal_type = typename CheckTemperatureDependence<
        ParameterPackType,
        std::make_index_sequence<ParameterPackType::size - 1>>::type;
    using model_type =
        typename FirstModelWithFractureType<fracture_type,
                                            ParameterPackType>::type;

    ForceModelsImpl( MaterialType t, ParameterPackType const& m )
        : type( t )
        , models( m )
    {
        setHorizon();
    }

    auto cutoff() const { return delta; }

    void updateBonds( const int num_local, const int max_neighbors )
    {
        ( Cabana::get<Indices>( models ).updateBounds( num_local,
                                                       max_neighbors ),
          ... );
    }

    KOKKOS_INLINE_FUNCTION int getIndex( const int i, const int j ) const
    {
        // TODO
        const int type_i = type( i );
        const int type_j = type( j );
        // FIXME: only for binary.
        if ( type_i == type_j )
            return type_i;
        else
            return 2;
    }

    // Extract model index and hide template indexing.
    template <typename Tag, typename... Args>
    KOKKOS_INLINE_FUNCTION auto operator()( Tag tag, const int i, const int j,
                                            Args... args ) const
    {
        auto t = getIndex( i, j );
        // Call individual model.
        return run_functor_for_index_in_pack_with_args(
            IdentityFunctor{}, t, models, tag, i, j, args... );
    }

    // This is only for LPS force/energy, currently the only cases that require
    // type information. When running models individually, the SingleMaterial
    // tag is used in the model directly; here it is replaced with the
    // MultiMaterial tag instead.
    template <typename Tag, typename... Args>
    KOKKOS_INLINE_FUNCTION auto operator()( Tag tag, SingleMaterial,
                                            const int i, const int j,
                                            Args... args ) const
    {
        const int type_i = type( i );
        const int type_j = type( j );

        auto t = getIndex( i, j );
        MultiMaterial mtag;
        // Call individual model.
        return run_functor_for_index_in_pack_with_args(
            IdentityFunctor{}, t, models, tag, mtag, type_i, type_j, args... );
    }

    auto horizon( const int ) { return delta; }
    auto maxHorizon() { return delta; }

    void update( const MaterialType _type ) { type = _type; }

    double delta;
    MaterialType type;
    ParameterPackType models;

  protected:
    void setHorizon()
    {
        // Enforce equal cutoff for now.
        const double tol = 1e-10;
        delta = Cabana::get<1>( models ).delta;

        if ( ( ( std::abs( Cabana::get<Indices>( models ).delta - delta ) >
                 tol ) ||
               ... ) )
        {
            log_err( std::cout, "Horizon for each model must match for "
                                "multi-material systems." );
        }
    }
};
} // namespace Impl

/******************************************************************************
  Multi-material models
******************************************************************************/
struct AverageTag
{
};
template <typename MaterialType, typename ParameterPackType>
struct ForceModels : CabanaPD::Impl::ForceModelsImpl<
                         MaterialType, ParameterPackType,
                         std::make_index_sequence<ParameterPackType::size - 1>>
{
    ForceModels( MaterialType t, ParameterPackType const& m )
        : CabanaPD::Impl::ForceModelsImpl<
              MaterialType, ParameterPackType,
              std::make_index_sequence<ParameterPackType::size - 1>>( t, m )
    {
    }
};

template <typename ParticleType, typename... ModelTypes>
auto createMultiForceModel( ParticleType particles, ModelTypes... m )
{
    auto type = particles.sliceType();
    return ForceModels( type, Cabana::makeParameterPack( m... ) );
}

template <typename ParticleType, typename ModelType1, typename ModelType2>
auto createMultiForceModel( ParticleType particles, AverageTag, ModelType1 m1,
                            ModelType2 m2 )
{
    ModelType1 m12( m1, m2 );
    return createMultiForceModel( particles, m1, m2, m12 );
}

} // namespace CabanaPD

#endif
