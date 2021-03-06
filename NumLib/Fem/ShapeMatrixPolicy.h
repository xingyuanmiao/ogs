/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef SHAPE_MATRIX_POLICY_H_
#define SHAPE_MATRIX_POLICY_H_

#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"

#ifdef OGS_USE_EIGEN
#include <Eigen/Dense>

namespace detail
{
    /// Forwards the Eigen::Matrix type for general N and M.
    /// There is a partial specialization for M = 1 to store the matrix in
    /// column major storage order.
    template <int N, int M>
    struct EigenMatrixType
    {
        using type = Eigen::Matrix<double, N, M, Eigen::RowMajor>;
    };

    /// Specialization for Nx1 matrices which can be stored only in column major
    /// form in Eigen-3.2.5.
    template <int N>
    struct EigenMatrixType<N, 1>
    {
        using type = Eigen::Matrix<double, N, 1, Eigen::ColMajor>;
    };

    /// Specialization for 0xM matrices. Using fixed size Eigen matrices here
    /// would lead to zero sized arrays, which cause compilation errors on
    /// some compilers.
    template <int M>
    struct EigenMatrixType<0, M>
    {
        using type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    };

    template <>
    struct EigenMatrixType<0, 1>
    {
        using type = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    };

}   // detail

/// An implementation of ShapeMatrixPolicy using fixed size (compile-time) eigen
/// matrices and vectors.
template <typename ShapeFunction, unsigned GlobalDim>
struct EigenFixedShapeMatrixPolicy
{
    template <int N>
    using VectorType = typename ::detail::EigenMatrixType<N, 1>::type;

    template <int N, int M>
    using MatrixType = typename ::detail::EigenMatrixType<N, M>::type;

    using NodalMatrixType = MatrixType<ShapeFunction::NPOINTS, ShapeFunction::NPOINTS>;
    using NodalVectorType = VectorType<ShapeFunction::NPOINTS>;
    using DimNodalMatrixType = MatrixType<ShapeFunction::DIM, ShapeFunction::NPOINTS>;
    using DimMatrixType = MatrixType<ShapeFunction::DIM, ShapeFunction::DIM>;
    using GlobalDimNodalMatrixType = MatrixType<GlobalDim, ShapeFunction::NPOINTS>;
    using GlobalDimMatrixType = MatrixType<GlobalDim, GlobalDim>;

    using ShapeMatrices =
        NumLib::ShapeMatrices<
            NodalVectorType,
            DimNodalMatrixType,
            DimMatrixType,
            GlobalDimNodalMatrixType>;
};

/// An implementation of ShapeMatrixPolicy using dynamic size eigen matrices and
/// vectors.
template <typename ShapeFunction, unsigned GlobalDim>
struct EigenDynamicShapeMatrixPolicy
{
    // Dynamic size local matrices are much slower in allocation than their
    // fixed counterparts.

    template<int N, int M>
    using MatrixType =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    template<int N>
    using VectorType =
        Eigen::Matrix<double, Eigen::Dynamic, 1>;

    using NodalMatrixType = MatrixType<0,0>;
    using NodalVectorType = VectorType<0>;
    using DimNodalMatrixType = MatrixType<0,0>;
    using DimMatrixType = MatrixType<0,0>;
    using GlobalDimNodalMatrixType = MatrixType<0,0>;
    using GlobalDimMatrixType = MatrixType<0,0>;

    using ShapeMatrices =
        NumLib::ShapeMatrices<
            NodalVectorType,
            DimNodalMatrixType,
            DimMatrixType,
            GlobalDimNodalMatrixType>;
};

#ifdef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
template <typename ShapeFunction, unsigned GlobalDim>
using ShapeMatrixPolicyType = EigenDynamicShapeMatrixPolicy<ShapeFunction, GlobalDim>;

const unsigned OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG = 1;
#else
template <typename ShapeFunction, unsigned GlobalDim>
using ShapeMatrixPolicyType = EigenFixedShapeMatrixPolicy<ShapeFunction, GlobalDim>;

const unsigned OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG = 0;
#endif

#endif  // OGS_USE_EIGEN

//static_assert(std::is_class<ShapeMatrixPolicyType<>::value,
        //"ShapeMatrixPolicyType was not defined.");

#endif  // SHAPE_MATRIX_POLICY_H_
