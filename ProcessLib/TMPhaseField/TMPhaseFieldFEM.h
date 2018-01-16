/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>
#include <Eigen/Eigenvalues>

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropicTMPhaseField.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "TMPhaseFieldProcessData.h"

namespace ProcessLib
{
namespace TMPhaseField
{
template <typename BMatricesType, typename ShapeMatrixType, int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    typename ShapeMatrixType::NodalRowVectorType N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType dNdx;
    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;
    typename BMatricesType::KelvinVectorType eps_m, eps_m_prev;

    typename BMatricesType::KelvinVectorType sigma_tensile, sigma_compressive,
        sigma_real_prev, sigma_real;
    double strain_energy_tensile;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables> material_state_variables;

    typename BMatricesType::KelvinMatrixType C_tensile, C_compressive;
    double integration_weight;
    double history_variable;
    double history_variable_prev;

    void pushBackState()
    {
        if (history_variable_prev < history_variable)
        {
            history_variable_prev = history_variable;
        }
        eps_m_prev = eps_m;
        eps_prev = eps;
        sigma_prev = sigma;
        sigma_real_prev = sigma_real;
        material_state_variables->pushBackState();
    }

    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants =
        MaterialLib::SolidModels::Invariants<kelvin_vector_size>;

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(
        double const t,
        SpatialPosition const& x_position,
        double const dt,
        DisplacementVectorType const& u,
        double const alpha,
        double const delta_T,
        double const degradation)
    {
        eps_m.noalias() = eps - alpha * delta_T * Invariants::identity2;
        solid_material.integrateStress(
            t, x_position, dt, eps_m_prev, eps_m, sigma_prev,
            *material_state_variables);

        static_cast<MaterialLib::Solids::TMPhaseFieldExtension<DisplacementDim>&>(
            solid_material)
            .calculateDegradedStress(t, x_position, eps_m,
                             strain_energy_tensile, sigma_tensile,
                             sigma_compressive, C_tensile, C_compressive,
                             sigma_real, degradation);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class TMPhaseFieldLocalAssembler : public TMPhaseFieldLocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;

    // Types for displacement.
    // (Higher order elements = ShapeFunction).
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using RhsVector = typename ShapeMatricesType::template VectorType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim>;
    using JacobianMatrix = typename ShapeMatricesType::template MatrixType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim,
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim>;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
            ShapeMatricesType, ShapeFunction::NPOINTS,
            2 * ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim, DisplacementDim>;

    using NodalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

    TMPhaseFieldLocalAssembler(TMPhaseFieldLocalAssembler const&) = delete;
    TMPhaseFieldLocalAssembler(TMPhaseFieldLocalAssembler&&) = delete;

    TMPhaseFieldLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        TMPhaseFieldProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            ip_data.integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ;

            ip_data.sigma.setZero(kelvin_vector_size);
            ip_data.sigma_prev.resize(kelvin_vector_size);
            ip_data.eps.setZero(kelvin_vector_size);
            ip_data.eps_prev.resize(kelvin_vector_size);
            ip_data.eps_m.setZero(kelvin_vector_size);
            ip_data.eps_m_prev.resize(kelvin_vector_size);
            ip_data.C_tensile.setZero(kelvin_vector_size, kelvin_vector_size);
            ip_data.C_compressive.setZero(kelvin_vector_size,
                                         kelvin_vector_size);
            ip_data.sigma_tensile.setZero(kelvin_vector_size);
            ip_data.sigma_compressive.setZero(kelvin_vector_size);
            ip_data.history_variable =
                process_data.history_field(0, x_position)[0];
            ip_data.history_variable_prev =
                process_data.history_field(0, x_position)[0];
            ip_data.sigma_real.resize(kelvin_vector_size);

            ip_data.N = shape_matrices[ip].N;
            ip_data.dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "TMPhaseFieldLocalAssembler: assembly without jacobian is not "
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
        auto const local_matrix_size = local_x.size();
        assert(local_matrix_size == temperature_size + phasefield_size + displacement_size);

        auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                    temperature_size);

        auto d = Eigen::Map<typename ShapeMatricesType::template VectorType<
            phasefield_size> const>(local_x.data() + phasefield_index,
                                    phasefield_size);

        auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

        auto T_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                    temperature_size);

        auto d_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + phasefield_index,
                                    phasefield_size);

        auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_rhs = MathLib::createZeroedVector<RhsVector>(
            local_rhs_data, local_matrix_size);

        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        phasefield_size>
            Kud;
        Kud.setZero(displacement_size, phasefield_size);

        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        temperature_size>
            KuT;
        KuT.setZero(displacement_size, temperature_size);

        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        displacement_size>
            Kdu;
        Kdu.setZero(phasefield_size, displacement_size);

        typename ShapeMatricesType::template MatrixType<phasefield_size,
                                                        temperature_size>
            KdT;
        KdT.setZero(phasefield_size, temperature_size);

        typename ShapeMatricesType::template MatrixType<temperature_size,
                                                        phasefield_size>
            KTd;
        KTd.setZero(temperature_size, phasefield_size);

        typename ShapeMatricesType::NodalMatrixType KTT;
        KTT.setZero(temperature_size, temperature_size);

        typename ShapeMatricesType::NodalMatrixType DTT;
        DTT.setZero(temperature_size, temperature_size);

        typename ShapeMatricesType::NodalMatrixType Kdd;
        Kdd.setZero(phasefield_size, phasefield_size);

        typename ShapeMatricesType::NodalMatrixType Ddd;
        Ddd.setZero(phasefield_size, phasefield_size);

        double const& dt = _process_data.dt;

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& dNdx = _ip_data[ip].dNdx;
            auto const& N = _ip_data[ip].N;

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
            auto const& B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);

            auto& eps_m = _ip_data[ip].eps_m;
            auto& eps = _ip_data[ip].eps;

            auto const& C_tensile = _ip_data[ip].C_tensile;
            auto const& C_compressive = _ip_data[ip].C_compressive;

            auto const& strain_energy_tensile =
                _ip_data[ip].strain_energy_tensile;
            auto const& sigma_tensile = _ip_data[ip].sigma_tensile;

            auto& history_variable = _ip_data[ip].history_variable;
            auto& history_variable_prev = _ip_data[ip].history_variable_prev;

            auto const& sigma_real = _ip_data[ip].sigma_real;

            // auto const [&](member){ return _process_data.member(t,
            // x_position); };
            auto const k = _process_data.residual_stiffness(t, x_position)[0];
            auto const gc = _process_data.crack_resistance(t, x_position)[0];
            auto const ls = _process_data.crack_length_scale(t, x_position)[0];
            auto const M = _process_data.kinetic_coefficient(t, x_position)[0];
            auto rho_sr = _process_data.solid_density(t, x_position)[0];
            auto const alpha = _process_data.linear_thermal_expansion_coefficient(t, x_position)[0];
            double const c = _process_data.specific_heat_capacity(t, x_position)[0];
            auto const lambda = _process_data.thermal_conductivity(t, x_position)[0];
            auto const lambda_res = _process_data.residual_thermal_conductivity(t, x_position)[0];
            double const T0 = _process_data.reference_temperature;
            auto const& b = _process_data.specific_body_force;

            double const T_ip = N.dot(T);
            double const delta_T = T_ip - T0;
            // calculate real density
            double const rho_s = rho_sr * (1 - 3 * alpha * delta_T);
            // calculate thermally induced strain

            // Kdd_1 defines one term which both used in Kdd and local_rhs for phase field
            typename ShapeMatricesType::NodalMatrixType const Kdd_1 = dNdx.transpose() * 2 * gc * ls * dNdx;

            //
            // displacement equation, displacement part
            //
            auto const d_prev = d - d_dot*dt;

            double const d_ip = N.dot(d);
            double const d_ip_prev = N.dot(d_prev);
            double const degradation = d_ip * d_ip * (1 - k) + k;
            double const degradation_prev = d_ip_prev * d_ip_prev * (1 - k) + k;
            eps.noalias() = B * u;
            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u, alpha, delta_T, degradation_prev);

            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() +=
                B.transpose() * (degradation_prev * C_tensile + C_compressive) * B * w;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);

            for (int i = 0; i < DisplacementDim; ++i)
                N_u.template block<1, displacement_size / DisplacementDim>(
                     i, i * displacement_size / DisplacementDim)
                    .noalias() = N;

            using Invariants =
                MaterialLib::SolidModels::Invariants<kelvin_vector_size>;

            local_rhs
                .template block<displacement_size, 1>(displacement_index, 0)
                .noalias() -=
                (B.transpose() * sigma_real - N_u.transpose() * rho_s * b) * w;

            //
            // displacement equation, temperature part
            //
            KuT.noalias() += B.transpose() * (degradation_prev * C_tensile + C_compressive) *
                                             alpha * Invariants::identity2 * N * w;

            //
            // displacement equation, phasefield part
            //
            Kud.noalias() += B.transpose() * 2 * d_ip * sigma_tensile * N * w;

            if (history_variable_prev < strain_energy_tensile)
            {
                // INFO("History variable %g:", history_variable);
                // INFO("History variable previous %g:", history_variable_prev);
                history_variable = strain_energy_tensile;
                Kdu.noalias() += N.transpose() * 2 * d_ip * sigma_tensile.transpose() * B * w;
                KdT.noalias() += N.transpose() * 2 * d_ip * sigma_tensile.transpose() *
                        alpha * Invariants::identity2 * N * w;
            }
            else
            {
                history_variable = history_variable_prev;
            }

            //
            // phasefield equation, phasefield part
            //

            Kdd.noalias() += (Kdd_1 +
                              N.transpose() * 2 * history_variable * N +
                              N.transpose() * 0.5 * gc / ls * N) *
                             w;

            Ddd.noalias() += N.transpose() / M * N * w;

            //
            // temperature equation, temperature part;
            // temperature equation, phasefield part;
            // phasefield equation, temperature part
            //

            double const epsm_trace = Invariants::trace(eps_m);
            if (epsm_trace >= 0)
            {
                // KTT.noalias() += dNdx.transpose() * (d_ip*d_ip*(1-lambda_res) + lambda_res) *
                //                 lambda * dNdx * w;
                KTT.noalias() += dNdx.transpose() * (d_ip*d_ip*lambda + (1 - d_ip)*(1 - d_ip)*lambda_res) *
                                 dNdx * w;
                KTd.noalias() += dNdx.transpose() * 2 * d_ip * lambda *
                                 dNdx * T * N * w;
            }
            else
            {
                KTT.noalias() += dNdx.transpose() * lambda * dNdx * w;
            }

            DTT.noalias() += N.transpose() * rho_s * c * N * w;

            double const d_dot_ip = N.dot(d_dot);

            local_rhs.template block<phasefield_size, 1>(phasefield_index, 0)
               .noalias() -=
               (N.transpose() * d_dot_ip / M +
                Kdd_1 * d +
                N.transpose() * d_ip * 2 * history_variable -
                N.transpose() * 0.5 * gc / ls * (1 - d_ip)) *
               w;

        }
        // temperature equation, temperature part
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += KTT + DTT / dt;
        // temperature equation, phasefield part
        local_Jac
            .template block<temperature_size, phasefield_size>(
                temperature_index, phasefield_index)
            .noalias() += KTd;
        // displacement equation, temperature part
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= KuT;
        // displacement equation, phasefield part
        local_Jac
            .template block<displacement_size, phasefield_size>(
                displacement_index, phasefield_index)
            .noalias() += Kud;
        // phasefield equation, phasefield part.
        local_Jac
            .template block<phasefield_size, phasefield_size>(
                phasefield_index, phasefield_index)
            .noalias() += Kdd + Ddd / dt;
        // phasefield equation, displacement part.
        local_Jac
            .template block<phasefield_size, displacement_size>(
                phasefield_index, displacement_index)
            .noalias() += Kdu;
        // phasefield equation, temperature part
        local_Jac
            .template block<phasefield_size, temperature_size>(
                phasefield_index, temperature_index)
            .noalias() -= KdT;

        local_rhs.template block<temperature_size, 1>(temperature_index, 0)
           .noalias() -= KTT * T + DTT * T_dot;

        local_rhs.template block<temperature_size, 1>(temperature_index, 0)
           .noalias() -= KTd * d;

        local_rhs.template block<phasefield_size, 1>(phasefield_index, 0)
           .noalias() += KdT * T;

        local_rhs.template block<displacement_size, 1>(displacement_index, 0)
           .noalias() += KuT * T;

    }

//template <typename ShapeFunction,typename IntegrationMethod,
//          int DisplacementDim>
//std::vector<double> const&
//TMPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
//DisplacementDim>::
//    getIntPtHeatFlux(const double t,
//                     GlobalVector const& current_solution,
//                     NumLib::LocalToGlobalIndexMap const& dof_table,
//                     std::vector<double>& cache) const
//{
//    auto const num_intpts = _ip_data.size();

//    auto const indices = NumLib::getIndices(_element.getID(), dof_table);
//    assert(!indices.empty());
//    auto const local_x = current_solution.get(indices);

//    cache.clear();
//    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
//            double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
//                cache, DisplacementDim, num_intpts);

//    SpatialPosition pos;
//    pos.setElementID(_element.getID());

//    auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
//            temperature_size> const>(local_x.data() + temperature_index,
//                                     temperature_size);

//    auto d = Eigen::Map<typename ShapeMatricesType::template VectorType<
//            phasefield_size> const>(local_x.data() + phasefield_index,
//                                    phasefield_size);

//    unsigned const n_integration_points =
//            _integration_method.getNumberOfPoints();

//    SpatialPosition x_position;
//    x_position.setElementID(_element.getID());
//    for (unsigned ip = 0; ip < n_integration_points; ip++)
//    {
//        x_position.setIntegrationPoint(ip);
//        auto const& dNdx = _ip_data[ip].dNdx;
//        auto const& N = _ip_data[ip].N;
//        double const d_ip = N.dot(d);
//        auto const lambda = _process_data.thermal_conductivity(t, x_position)[0];
//        auto const lambda_res = _process_data.residual_thermal_conductivity(t, x_position)[0];

//        // Compute the heat flux
//        cache_matrix.col(ip).noalias() =
//                -(d_ip*d_ip*lambda + (1 - d_ip)*(1 - d_ip)*lambda_res) * dNdx * T;
//    }

//    return cache;
//}
    /*void computeSecondaryVariableConcrete(
            const double t, std::vector<double> const& local_x) override
        {
            auto const local_matrix_size = local_x.size();
            assert(local_matrix_size == temperature_size + phasefield_size + displacement_size);

            auto d = Eigen::Map<typename ShapeMatricesType::template VectorType<
                phasefield_size> const>(local_x.data() + phasefield_index,
                                        phasefield_size);

            unsigned const n_integration_points =
                _integration_method.getNumberOfPoints();

            SpatialPosition x_position;
            x_position.setElementID(_element.getID());
            const auto local_x_vec =
                MathLib::toVector<NodalVectorType>(local_x, local_matrix_size);

            for (unsigned ip = 0; ip < n_integration_points; ip++)
            {
                x_position.setIntegrationPoint(ip);
                auto const& dNdx = _ip_data[ip].dNdx;
                auto const& N = _ip_data[ip].N;
                double const d_ip = N.dot(d);
                auto const lambda = _process_data.thermal_conductivity(t, x_position)[0];
                auto const lambda_res = _process_data.residual_thermal_conductivity(t, x_position)[0];
                // heat flux only computed for output.
                GlobalDimVectorType const heat_flux = -(d_ip*d_ip*lambda + (1 - d_ip)*(1 - d_ip)*lambda_res) * dNdx * local_x_vec;

                for (unsigned dm = 0; dm < DisplacementDim; ++dm)
                {
                    _heat_fluxes[dm][ip] = heat_flux[dm];
                }
            }
        }*/

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

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSigmaXX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 0);
    }

    std::vector<double> const& getIntPtSigmaYY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 1);
    }

    std::vector<double> const& getIntPtSigmaZZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 2);
    }

    std::vector<double> const& getIntPtSigmaXY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 3);
    }

    std::vector<double> const& getIntPtSigmaYZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 4);
    }

    std::vector<double> const& getIntPtSigmaXZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 5);
    }

    std::vector<double> const& getIntPtEpsilonXX(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 0);
    }

    std::vector<double> const& getIntPtEpsilonYY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 1);
    }

    std::vector<double> const& getIntPtEpsilonZZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 2);
    }

    std::vector<double> const& getIntPtEpsilonXY(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 3);
    }

    std::vector<double> const& getIntPtEpsilonYZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 4);
    }

    std::vector<double> const& getIntPtEpsilonXZ(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 5);
    }

    std::vector<double> const& getIntPtHeatFlux(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override
    {
        auto const num_intpts = _ip_data.size();

        auto const indices = NumLib::getIndices(_element.getID(), dof_table);
        assert(!indices.empty());
        auto const local_x = current_solution.get(indices);

        cache.clear();
        auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, num_intpts);

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                    temperature_size);

        auto d = Eigen::Map<typename ShapeMatricesType::template VectorType<
            phasefield_size> const>(local_x.data() + phasefield_index,
                                    phasefield_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition x_position;
        x_position.setElementID(_element.getID());
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& dNdx = _ip_data[ip].dNdx;
            auto const& N = _ip_data[ip].N;
            double const d_ip = N.dot(d);
            auto const lambda = _process_data.thermal_conductivity(t, x_position)[0];
            auto const lambda_res = _process_data.residual_thermal_conductivity(t, x_position)[0];

            // Compute the heat flux
            cache_matrix.col(ip).noalias() =
                -(d_ip*d_ip*lambda + (1 - d_ip)*(1 - d_ip)*lambda_res) * dNdx * T;
        }

        return cache;
    }

private:
    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                             std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.sigma_real[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.sigma_real[component] / std::sqrt(2));
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        std::vector<double>& cache, std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            cache.push_back(ip_data.eps[component]);
        }

        return cache;
    }

    TMPhaseFieldProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;

    static const int temperature_index = 0;
    static const int temperature_size = ShapeFunction::NPOINTS;
    static const int phasefield_index = ShapeFunction::NPOINTS;
    static const int phasefield_size = ShapeFunction::NPOINTS;
    static const int displacement_index = 2 * ShapeFunction::NPOINTS;
    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
};

/*template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerData final
    : public TMPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                      DisplacementDim>
{
public:
    LocalAssemblerData(LocalAssemblerData const&) = delete;
    LocalAssemblerData(LocalAssemblerData&&) = delete;

    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       TMPhaseFieldProcessData<DisplacementDim>& process_data)
        : TMPhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};*/

}  // namespace TMPhaseField
}  // namespace ProcessLib
