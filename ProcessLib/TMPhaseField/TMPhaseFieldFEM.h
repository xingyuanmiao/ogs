/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TMPHASEFIELD_FEM_H_
#define PROCESS_LIB_TMPHASEFIELD_FEM_H_

#include <iostream>
#include <memory>
#include <vector>
#include <Eigen/Eigenvalues>

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropicTMPhaseField.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
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
          _eps_m(std::move(other._eps_m)),
          _eps_m_prev(std::move(other._eps_m_prev)),
          _solid_material(other._solid_material),
          _material_state_variables(std::move(other._material_state_variables)),
          _C(std::move(other._C)),
          integration_weight(std::move(other.integration_weight)),
          strain_energy_tensile(std::move(other.strain_energy_tensile)),
          _sigma_tensile(std::move(other._sigma_tensile)),
          _sigma_compressive(std::move(other._sigma_compressive))
    {
    }
#endif  // _MSC_VER

    typename ShapeMatrixType::NodalRowVectorType _N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType _dNdx;
    typename BMatricesType::BMatrixType _b_matrices;
    typename BMatricesType::KelvinVectorType _sigma, _sigma_prev;
    typename BMatricesType::KelvinVectorType _eps, _eps_prev;
    typename BMatricesType::KelvinVectorType _eps_m, _eps_m_prev;

    typename BMatricesType::KelvinVectorType _sigma_tensile, _sigma_compressive;
    double strain_energy_tensile;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& _solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        _material_state_variables;

    typename BMatricesType::KelvinMatrixType _C;
    double integration_weight;

    void pushBackState()
    {
        _eps_m_prev = _eps_m;
        _eps_prev = _eps;
        _sigma_prev = _sigma;
        _material_state_variables->pushBackState();
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
        double& delta_T)
    {
        _eps.noalias() = _b_matrices * u;
        _eps_m.noalias() = _eps - alpha * delta_T * Invariants::identity2;
        _solid_material.computeConstitutiveRelation(
            t, x_position, dt, _eps_m_prev, _eps_m, _sigma_prev, _sigma, _C,
            *_material_state_variables);

        static_cast<MaterialLib::Solids::TMPhaseFieldExtension<DisplacementDim>&>(
            _solid_material)
            .specialFunction(t, x_position, _eps,
                             strain_energy_tensile, _sigma_tensile,
                             _sigma_compressive);
    }
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType> N;
};

struct PhaseFieldLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
    virtual std::vector<double> const& getIntPtSigmaXX(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaZZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYZ(
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class PhaseFieldLocalAssembler : public PhaseFieldLocalAssemblerInterface
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;

    // Types for displacement.
    // (Higher order elements = ShapeFunction).
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using RhsVector = typename ShapeMatricesType::template VectorType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim>;
    using JacobianMatrix = typename ShapeMatricesType::template MatrixType<
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim,
        ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim>;

    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler const&) = delete;
    PhaseFieldLocalAssembler(PhaseFieldLocalAssembler&&) = delete;

    PhaseFieldLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        TMPhaseFieldProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, DisplacementDim>(
                e, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // displacement (subscript u)
            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].detJ;
            ip_data._b_matrices.resize(
                kelvin_vector_size, ShapeFunction::NPOINTS * DisplacementDim);

            auto const x_coord =
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    e, shape_matrices[ip].N);
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunction::NPOINTS>(
                shape_matrices[ip].dNdx, ip_data._b_matrices,
                is_axially_symmetric, shape_matrices[ip].N, x_coord);

            ip_data._sigma.resize(kelvin_vector_size);
            ip_data._sigma_prev.resize(kelvin_vector_size);
            ip_data._eps.resize(kelvin_vector_size);
            ip_data._eps_prev.resize(kelvin_vector_size);
            ip_data._eps_m.resize(kelvin_vector_size);
            ip_data._eps_m_prev.resize(kelvin_vector_size);
            ip_data._C.resize(kelvin_vector_size, kelvin_vector_size);
            ip_data._sigma_tensile.resize(kelvin_vector_size);
            ip_data._sigma_compressive.resize(kelvin_vector_size);
            _ip_data[ip].strain_energy_tensile;

            ip_data._N = shape_matrices[ip].N;
            ip_data._dNdx = shape_matrices[ip].dNdx;

            _secondary_data.N[ip] = shape_matrices[ip].N;
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

        Eigen::MatrixXd local_Jac_numerical = Eigen::MatrixXd::Zero(
            local_matrix_size, local_matrix_size);

        // auto local_K = MathLib::createZeroedMatrix<JacobianMatrix>(
        //     local_K_data, local_matrix_size, local_matrix_size);

        // auto local_M = MathLib::createZeroedMatrix<JacobianMatrix>(
        //     local_M_data, local_matrix_size, local_matrix_size);

        Eigen::MatrixXd local_b_p = Eigen::VectorXd::Zero(local_matrix_size);

        Eigen::MatrixXd local_b_m = Eigen::VectorXd::Zero(local_matrix_size);

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

            auto const& dNdx = _ip_data[ip]._dNdx;
            auto const& N = _ip_data[ip]._N;

            auto const& B = _ip_data[ip]._b_matrices;
            auto const& sigma = _ip_data[ip]._sigma;
            auto const& eps = _ip_data[ip]._eps;

            auto const& C = _ip_data[ip]._C;

            auto const& strain_energy_tensile = _ip_data[ip].strain_energy_tensile;
            auto const& sigma_tensile = _ip_data[ip]._sigma_tensile;
            auto const& sigma_compressive = _ip_data[ip]._sigma_compressive;

            // auto const [&](member){ return _process_data.member(t,
            // x_position); };
            double const k = _process_data.residual_stiffness;
            double const gc = _process_data.crack_resistance;
            double const ls = _process_data.crack_length_scale;
            double const M = _process_data.kinetic_coefficient;
            double const gamma = _process_data.penalty_constant;
            auto rho_sr = _process_data.solid_density(t, x_position)[0];
            double const alpha = _process_data.linear_thermal_expansion_coefficient;
            double const c = _process_data.specific_heat_capacity;
            auto const lambda = _process_data.thermal_conductivity(t, x_position)[0];
            double const T0 = _process_data.reference_temperature;
            auto const& b = _process_data.specific_body_force;

            double T_ip = N.dot(T);
            double delta_T = T_ip - T0;
            // calculate real density
            rho_sr = rho_sr * (1 - 3 * alpha * delta_T);
            // calculate thermally induced strain

            // Kdd_1 defines one term which both used in Kdd and local_rhs for phase field
            typename ShapeMatricesType::NodalMatrixType const Kdd_1 = dNdx.transpose() * 2 * gc * ls * dNdx;

            //
            // displacement equation, displacement part
            //
            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u, alpha, delta_T);

            double const d_ip = N.dot(d);
            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() +=
                B.transpose() * (d_ip*d_ip + k) * C * B * w;

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
                (B.transpose() * ((d_ip*d_ip + k) * sigma_tensile + sigma_compressive)
                                  - N_u.transpose() * rho_sr * b) * w;
            local_rhs
                .template block<displacement_size, 1>(displacement_index, 0)
                .noalias() -=
                B.transpose() * C * alpha * T0 * Invariants::identity2 * w;

            // local_rhs
            //    .template block<displacement_size, 1>(displacement_index, 0)
            //    .noalias() -=
            //    B.transpose() * (d_ip*d_ip + k) * (C * alpha * T0 * Invariants::identity2) * w;
            //
            // displacement equation, temperature part
            //
            KuT.noalias() += B.transpose() * (d_ip*d_ip + k) * (C * alpha * Invariants::identity2) * N * w;

            //
            // displacement equation, phasefield part
            //
            Kud.noalias() += (B.transpose() * 2 * d_ip * sigma_tensile * N) * w;

            double const scalar = (Invariants::identity2).transpose() * C * Invariants::identity2;
            //
            // phasefield equation, phasefield part
            //
            Kdd.noalias() += (Kdd_1 +
                              N.transpose() * 2 * strain_energy_tensile * N +
                              N.transpose() * 0.5 * gc / ls * N) *
                             w;

            Ddd.noalias() += N.transpose() / M * N * w;

            //
            // temperature equation, temperature part;
            // temperature equation, phasefield part;
            // phasefield equation, temperature part
            //

            double const eps_trace = Invariants::trace(eps);
            if (eps_trace >= 0)
            {
                KTT.noalias() += dNdx.transpose() * (d_ip*d_ip + k) *
                        lambda * dNdx * w;
                KTd.noalias() += dNdx.transpose() * 2 * d_ip * lambda * T_ip * dNdx * w;
                KdT.noalias() += N.transpose() * 2 * d_ip * scalar * alpha * alpha * T_ip * N * w;
            }
            else
            {
                KTT.noalias() += dNdx.transpose() * lambda * dNdx * w;
            }

            DTT.noalias() += N.transpose() * rho_sr * c * N * w;

            double const d_dot_ip = N.dot(d_dot);

            local_rhs.template block<phasefield_size, 1>(phasefield_index, 0)
               .noalias() -=
               (N.transpose() * d_dot_ip / M +
                Kdd_1 * d +
                N.transpose() * d_ip * 2 * strain_energy_tensile -
                N.transpose() * 0.5 * gc / ls * (1 - d_ip)) *
               w;
            local_rhs.template block<phasefield_size, 1>(phasefield_index, 0)
               .noalias() -=
               N.transpose() * 0.5 * alpha * T0 * scalar * alpha * T0 * w;

            double const T_dot_ip = N.dot(T_dot);

            local_rhs.template block<temperature_size, 1>(temperature_index, 0)
               .noalias() -= KTT * T +
               N.transpose() * rho_sr * c * T_dot_ip * w;

            // calculate numerical Jac
            double num_p = 1e-8;
            Eigen::VectorXd num_vec = Eigen::VectorXd::Zero(local_matrix_size);
            std::vector<double> local_perturbed = local_x;
            for (Eigen::MatrixXd::Index i = 0; i < local_matrix_size; i++)
            {
                num_vec[i] = num_p;
                local_perturbed[i] += num_vec[i];
                auto T_p = Eigen::Map<typename ShapeMatricesType::template VectorType<
                    temperature_size> const>(local_x.data() + temperature_index,
                                            temperature_size);
                auto d_p = Eigen::Map<typename ShapeMatricesType::template VectorType<
                    phasefield_size> const>(local_perturbed.data() + phasefield_index,
                                            phasefield_size);
                auto u_p = Eigen::Map<typename ShapeMatricesType::template VectorType<
                    displacement_size> const>(local_perturbed.data() + displacement_index,
                                              displacement_size);
                double const d_ip_p = N.dot(d_p);
                double const T_ip_p = N.dot(T_p);
                _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u_p, alpha, delta_T);
                local_b_p
                    .template block<displacement_size, 1>(displacement_index, 0)
                    .noalias() =
                    (B.transpose() * ((d_ip_p*d_ip_p + k) * sigma_tensile + sigma_compressive)
                                      - N_u.transpose() * rho_sr * b) * w;
                local_b_p
                    .template block<displacement_size, 1>(displacement_index, 0)
                    .noalias() =
                    B.transpose() * C * alpha * T0 * Invariants::identity2 * w;
                local_b_p
                   .template block<phasefield_size, 1>(phasefield_index, 0)
                   .noalias() =
                   (N.transpose() * d_ip_p / dt / M +
                    Kdd_1 * d_p +
                    N.transpose() * d_ip_p * 2 * strain_energy_tensile -
                    N.transpose() * 0.5 * gc / ls * (1 - d_ip_p)) *
                   w;
                local_b_p.template block<phasefield_size, 1>(phasefield_index, 0)
                   .noalias() =
                   N.transpose() * 0.5 * alpha * T0 * scalar * alpha * T0 * w;
                local_b_p
                   .template block<temperature_size, 1>(temperature_index, 0)
                   .noalias() = KTT * T_p + N.transpose() * rho_sr * c * T_ip_p / dt * w;

                local_perturbed[i] = local_x[i] - num_vec[i];
                auto T_m = Eigen::Map<typename ShapeMatricesType::template VectorType<
                    temperature_size> const>(local_x.data() + temperature_index,
                                            temperature_size);
                auto d_m = Eigen::Map<typename ShapeMatricesType::template VectorType<
                    phasefield_size> const>(local_perturbed.data() + phasefield_index,
                                            phasefield_size);
                auto u_m = Eigen::Map<typename ShapeMatricesType::template VectorType<
                    displacement_size> const>(local_perturbed.data() + displacement_index,
                                              displacement_size);
                double const d_ip_m = N.dot(d_m);
                double const T_ip_m = N.dot(T_m);
                _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u_m, alpha, delta_T);
                local_b_m
                    .template block<displacement_size, 1>(displacement_index, 0)
                    .noalias() =
                    (B.transpose() * ((d_ip_m*d_ip_m + k) * sigma_tensile + sigma_compressive)
                                      - N_u.transpose() * rho_sr * b) * w;
                local_b_m
                    .template block<displacement_size, 1>(displacement_index, 0)
                    .noalias() =
                    B.transpose() * C * alpha * T0 * Invariants::identity2 * w;
                local_b_m
                   .template block<phasefield_size, 1>(phasefield_index, 0)
                   .noalias() =
                   (N.transpose() * d_ip_m / dt / M +
                    Kdd_1 * d_m +
                    N.transpose() * d_ip_m * 2 * strain_energy_tensile -
                    N.transpose() * 0.5 * gc / ls * (1 - d_ip_m)) *
                   w;
                local_b_m.template block<phasefield_size, 1>(phasefield_index, 0)
                   .noalias() =
                   N.transpose() * 0.5 * alpha * T0 * scalar * alpha * T0 * w;
                local_b_m
                   .template block<temperature_size, 1>(temperature_index, 0)
                   .noalias() = KTT * T_m + N.transpose() * rho_sr * c * T_ip_m / dt * w;
                local_perturbed[i] = local_x[i];
                local_Jac_numerical.col(i).noalias() += (local_b_p - local_b_m) / (2.0 * num_vec[i]);

            }

            // local_rhs.template block<phasefield_size, 1>(phasefield_index, 0)
            //    .noalias() -= N.transpose() * d_ip *
            //         alpha * T0 * scalar * alpha * T0;

            // local_rhs.template block<phasefield_size, 1>(phasefield_index, 0)
            //     .noalias() -=
            //     (Kdd_1 * d +
            //      N.transpose() * d_ip * eps.transpose() * C * eps -
            //      N.transpose() * 0.5 * gc / ls * (1 - d_ip)) *
            //     w;

            // Additional penalty term if inside the damaged region
            // if (damaged_region)
            // {
            //     Kdd.template block<phasefield_size, phasefield_size>(phasefield_index, 0)
            //         .noalias() += N.transpose() / gamma * N * w;
            //     local_rhs
            //         .template block<phasefield_size, 1>(phasefield_index, 0)
            //         .noalias() -= N.transpose() * d_ip / gamma * w;
                // damping
                // local_rhs
                //     .template block<displacement_size, 1>(displacement_index, 0)
                //     .noalias() = 0.5 * local_rhs
                //         .template block<displacement_size, 1>(displacement_index, 0);
                // local_rhs
                //     .template block<phasefield_size, 1>(phasefield_index, 0)
                //     .noalias() = 0.5 * local_rhs
                //         .template block<phasefield_size, 1>(phasefield_index, 0);
            // }

            // phasefield equation, displacement part.
            //
            // Reusing Kud.transpose().
        }
        // temperature equation, temperature part
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += KTT + DTT / dt;
        // displacement equation, phasefield part
        local_Jac
            .template block<displacement_size, phasefield_size>(
                displacement_index, phasefield_index)
            .noalias() += Kud;

        // displacement equation, temperature part
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= KuT;

        // phasefield equation, phasefield part.
        local_Jac
            .template block<phasefield_size, phasefield_size>(phasefield_index,
                                                              phasefield_index)
            .noalias() += Kdd + Ddd / dt;

        // phasefield equation, displacement part.
        local_Jac
            .template block<phasefield_size, displacement_size>(
                phasefield_index, displacement_index)
            .noalias() += Kud.transpose();

        // phasefield equation, temperature part.
        local_Jac
            .template block<phasefield_size, temperature_size>(
                phasefield_index, temperature_index)
            .noalias() -= KdT;

        // temperature equation, phasefield part.
        local_Jac
            .template block<temperature_size, phasefield_size>(
                temperature_index, phasefield_index)
            .noalias() += KTd;

        // temperature equation, temperature part.
        local_Jac
            .template block<temperature_size, temperature_size>(temperature_index,
                                                              temperature_index)
            .noalias() -= KTT + DTT / dt;

        // compare analytical Jac to numerical Jac
        for (Eigen::MatrixXd::Index i = 0; i < local_matrix_size; i++)
        {
            for (Eigen::MatrixXd::Index j = 0; j < local_matrix_size; j++)
            {
                double Jac_analytical = local_Jac(i,j);
                double Jac_numerical = local_Jac_numerical(i,j);
                if (Jac_analytical != 0 && Jac_numerical != 0)
                {
                    double relative_deviation =
                            std::abs(Jac_numerical - Jac_analytical) / Jac_analytical;
                    if (relative_deviation > 1e-4)
                    {
                        OGS_FATAL("Deviation larger than the tolerance.");
                    }
                }
                else if (Jac_analytical = 0 && Jac_numerical > 1e-2)
                {
                    OGS_FATAL("Numerical Jacobian element does not equal zero.");
                }
            }
        }
        // Eigen::EigenSolver<JacobianMatrix> eigensolver(local_Jac);
        // std::cout << "eigenvalues" << eigensolver.eigenvalues() << "\n";

    }

    void preTimestepConcrete(std::vector<double> const& local_x,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        // Update damaged region.
        auto const d =
            Eigen::Map<typename ShapeMatricesType::template VectorType<
                phasefield_size> const>(local_x.data() + phasefield_index,
                                        phasefield_size);
        double const crtol = _process_data.critical_tolerance;
        damaged_region = d.squaredNorm() < crtol*crtol;

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
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 0);
    }

    std::vector<double> const& getIntPtSigmaYY(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 1);
    }

    std::vector<double> const& getIntPtSigmaZZ(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 2);
    }

    std::vector<double> const& getIntPtSigmaXY(
        std::vector<double>& cache) const override
    {
        return getIntPtSigma(cache, 3);
    }

    std::vector<double> const& getIntPtSigmaXZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 4);
    }

    std::vector<double> const& getIntPtSigmaYZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtSigma(cache, 5);
    }

private:
    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                             std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data) {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data._sigma[component]);
            else    // mixed xy, yz, xz components
                cache.push_back(ip_data._sigma[component] / std::sqrt(2));
        }

        return cache;
    }

    TMPhaseFieldProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>>
        _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;

    /// CR_{l-1} := { x \in \Omega s.t. d < CRTOL }. Damaged region indicator.
    /// true means "is damaged". Updated after a timestep.
    bool damaged_region;

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
                       TMPhaseFieldProcessData<DisplacementDim>& process_data)
        : PhaseFieldLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace TMPhaseField
}  // namespace ProcessLib

#endif  // PROCESS_LIB_TMPHASEFIELD_FEM_H_
