/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PHASEFIELD_FEM_H_
#define PROCESS_LIB_PHASEFIELD_FEM_H_

#include <iostream>
#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "PhaseFieldProcessData.h"
#include "LowerDimShapeTable.h"

namespace ProcessLib
{
namespace PhaseField
{
template <typename BMatricesType, typename ShapeMatrixTypeDisplacement,
          typename ShapeMatricesTypePhaseField, int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : _solid_material(solid_material),
          _material_state_variables(
              _solid_material.createMaterialStateVariables())
    {
    }

#if defined(_MSC_VER) && _MSC_VER < 1900
    // The default generated move-ctor is correctly generated for other
    // compilers.
    explicit IntegrationPointData(IntegrationPointData&& other)
        : _b_matrices(std::move(other._b_matrices)),
          _sigma(std::move(other._sigma)),
          _sigma_prev(std::move(other._sigma_prev)),
          _eps(std::move(other._eps)),
          _eps_prev(std::move(other._eps_prev)),
          _solid_material(other._solid_material),
          _material_state_variables(std::move(other._material_state_variables)),
          _C(std::move(other._C)),
          integration_weight(std::move(other.integration_weight))
    {
    }
#endif  // _MSC_VER

    typename ShapeMatrixTypeDisplacement::NodalRowVectorType _N_u;
    typename BMatricesType::BMatrixType _b_matrices;
    typename BMatricesType::KelvinVectorType _sigma, _sigma_prev;
    typename BMatricesType::KelvinVectorType _eps, _eps_prev;

    typename ShapeMatricesTypePhaseField::NodalRowVectorType _N_d;
    typename ShapeMatricesTypePhaseField::GlobalDimNodalMatrixType _dNdx_d;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& _solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        _material_state_variables;

    typename BMatricesType::KelvinMatrixType _C;
    double integration_weight;

    void pushBackState()
    {
        _eps_prev = _eps;
        _sigma_prev = _sigma;
        _material_state_variables->pushBackState();
    }

    void updateConstitutiveRelation(
        double const t,
        SpatialPosition const& x_position,
        double const dt,
        Eigen::Map<typename ShapeMatrixTypeDisplacement::
                       NodalDOFVectorType const> const& u)
    {
        _eps.noalias() = _b_matrices * u;
        _solid_material.computeConstitutiveRelation(
            t, x_position, dt, _eps_prev, _eps, _sigma_prev, _sigma, _C,
            *_material_state_variables);
    }
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class PhaseFieldLocalAssembler : public ProcessLib::LocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    // Types for displacement.
    // (Higher order elements = ShapeFunctionDisplacement).
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename BMatricesType::NodalForceVectorType;

    // Types for phasefield.
    // (Lower order elements = Order(ShapeFunctionDisplacement) - 1).
    using ShapeFunctionPhaseField =
        typename LowerDim<ShapeFunctionDisplacement>::type;
    using ShapeMatricesTypePhaseField =
        ShapeMatrixPolicyType<ShapeFunctionPhaseField, DisplacementDim>;

    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler const&) = delete;
    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler&&) = delete;

    PhaseFieldLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        PhaseFieldProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);

        auto const shape_matrices_u =
            initShapeMatrices<ShapeFunctionDisplacement, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        auto const shape_matrices_d =
            initShapeMatrices<ShapeFunctionPhaseField, ShapeMatricesTypePhaseField,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices_u[ip].detJ;
            ip_data._b_matrices.resize(
                kelvin_vector_size,
                ShapeFunctionDisplacement::NPOINTS * DisplacementDim);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunctionDisplacement,
                                       ShapeMatricesType>(
                    e, shape_matrices_u[ip].N);
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS>(
                shape_matrices_u[ip].dNdx, ip_data._b_matrices,
                is_axially_symmetric, shape_matrices_u[ip].N, x_coord);

            ip_data._sigma.resize(kelvin_vector_size);
            ip_data._sigma_prev.resize(kelvin_vector_size);
            ip_data._eps.resize(kelvin_vector_size);
            ip_data._eps_prev.resize(kelvin_vector_size);
            ip_data._C.resize(kelvin_vector_size, kelvin_vector_size);

            ip_data._N_u = shape_matrices_u[ip].N;
            ip_data._N_d = shape_matrices_d[ip].N;
            ip_data._dNdx_d = shape_matrices_d[ip].dNdx;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "PhaseFieldLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override
    {
        /*
        auto const local_matrix_size = local_x.size();
        assert(local_matrix_size == phasefield_size + displacement_size);

        auto d =
            Eigen::Map<typename ShapeMatricesTypePhaseField::template VectorType<
                phasefield_size> const>(local_x.data() + phasefield_index,
                                      phasefield_size);

        auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

        auto d_dot =
            Eigen::Map<typename ShapeMatricesTypePhaseField::template VectorType<
                phasefield_size> const>(local_xdot.data() + phasefield_index,
                                      phasefield_size);
        auto u_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
                                      displacement_size);

        auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_rhs =
            MathLib::createZeroedVector<NodalDisplacementVectorType>(
                local_rhs_data, local_matrix_size);

        double const& dt = _process_data.dt;

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N_d = _ip_data[ip]._N_d;
            auto const& dNdx_d = _ip_data[ip]._dNdx_d;

            auto const& B = _ip_data[ip]._b_matrices;
            auto const& sigma = _ip_data[ip]._sigma;
            auto const& eps = _ip_data[ip]._eps;

            auto const& C = _ip_data[ip]._C;

            // auto const [&](member){ return _process_data.member(t,
            // x_position); };
            double const k =
                _process_data.residual_stiffness(t, x_position)[0];
            double const gc =
                _process_data.crack_resistance(t, x_position)[0];
            double const ls =
                _process_data.crack_length_scale(t, x_position)[0];
            double const M =
                _process_data.kinetic_coefficient(t, x_position)[0];
            auto const rho_sr = _process_data.solid_density(t, x_position)[0];
            auto const body_force =
                _process_data.specific_body_force(t, x_position);
            assert(body_force.size() == DisplacementDim);
            auto const b =
                Eigen::Map<typename ShapeMatricesType::template VectorType<
                    DisplacementDim> const>(body_force.data(), DisplacementDim);

            //
            // displacement equation, displacement part
            //
            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u);

            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() += B.transpose() * (d * d + k) * C * B * w;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);

            for (int i = 0; i < DisplacementDim; ++i)
                N_u.template block<1, displacement_size / DisplacementDim>(
                       i, i * displacement_size / DisplacementDim)
                    .noalias() = _ip_data[ip]._N_u;

            local_rhs
                .template block<displacement_size, 1>(displacement_index, 0)
                .noalias() -=
                (B.transpose() * sigma - N_u.transpose() * rho_sr * b) * w;

            //
            // displacement equation, phasefield part
            //
            Kud.noalias() += B.transpose() * 2. * d * C * eps * N_d * w;

            //
            // phasefield equation, displacement part
            //
            Kdu.noalias() += N_d.transpose() * 2. * d * (C * eps).transpose() * B * w;

            //
            // phasefield equation, phasefield part.
            //
            Kdd.noalias() += (dNdx_d.transpose() * 2. * gc * ls * dNdx_d +
                              N_d.transpose() * eps.transpose() * C * eps +
                              N_d.transpose() * gc / 2. / ls * N_d) * w;

            Ddd.noalias() += N_d.transpose() / M * N_d * w;

            local_rhs.template block<_size, 1>(phasefield_index, 0)
                .noalias() += (N_d.transpose() * d_dot / M +
                               dNdx_d.transpose() * 2 * gc * ls * dNdx_d * d +
                               N_d.transpose() * d * eps.transpose() * C * eps -
                               N_d.transpose() * gc / 2. / ls * (1 - d)) * w;

            //
            // phasefield equation, displacement part.
            //
            // Reusing Kud.transpose().
        }
        // displacement equation, phasefield part
        local_Jac
            .template block<displacement_size, phasefield_size>(
                displacement_index, phasefield_index)
            .noalias() -= Kud;

        // phasefield equation, phasefield part.
        local_Jac
            .template block<phasefield_size, phasefield_size>(phasefield_index,
                                                          phasefield_index)
            .noalias() += Kdd + Ddd / dt;

        // phasefield equation, displacement part.
        local_Jac
            .template block<phasefield_size, displacement_size>(
                phasefield_index, displacement_index)
            .noalias() += Kdu;
        */
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

private:
    PhaseFieldProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType,
                             ShapeMatricesTypePhaseField, DisplacementDim>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;

    static const int phasefield_index = 0;
    static const int phasefield_size = ShapeFunctionPhaseField::NPOINTS;
    static const int displacement_index = ShapeFunctionPhaseField::NPOINTS;
    static const int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;
    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerData final
    : public PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                          DisplacementDim>
{
public:
    LocalAssemblerData(LocalAssemblerData const&) = delete;
    LocalAssemblerData(LocalAssemblerData&&) = delete;

    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       PhaseFieldProcessData<DisplacementDim>& process_data)
        : PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                       DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace PhaseField
}  // namespace ProcessLib

#endif  // PROCESS_LIB_PHASEFIELD_FEM_H_
