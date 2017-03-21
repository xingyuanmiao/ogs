/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_THERMOMECHANICS_FEM_H_
#define PROCESS_LIB_THERMOMECHANICS_FEM_H_

#include <iostream>
#include <memory>
#include <vector>
#include <Eigen/Eigenvalues>

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
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

#include "ThermoMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanics
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
    {
    }
#endif  // _MSC_VER

    typename ShapeMatrixType::NodalRowVectorType _N;
    typename ShapeMatrixType::GlobalDimNodalMatrixType _dNdx;
    typename BMatricesType::BMatrixType _b_matrices;
    typename BMatricesType::KelvinVectorType _sigma, _sigma_prev;
    typename BMatricesType::KelvinVectorType _eps, _eps_prev;
    typename BMatricesType::KelvinVectorType _eps_m, _eps_m_prev;

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

    }
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType> N;
};

struct ThermoMechanicsLocalAssemblerInterface
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

    virtual std::vector<double> const& getIntPtEpsilonXX(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonZZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXY(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXZ(
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYZ(
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class ThermoMechanicsLocalAssembler : public ThermoMechanicsLocalAssemblerInterface
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

    ThermoMechanicsLocalAssembler(ThermoMechanicsLocalAssembler const&) = delete;
    ThermoMechanicsLocalAssembler(ThermoMechanicsLocalAssembler&&) = delete;

    ThermoMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        ThermoMechanicsProcessData<DisplacementDim>& process_data)
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
            "ThermoMechanicsLocalAssembler: assembly without jacobian is not "
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
        assert(local_matrix_size == temperature_size + displacement_size);

        auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_x.data() + temperature_index,
                                    temperature_size);

        auto u = Eigen::Map<typename ShapeMatricesType::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

        auto T_dot = Eigen::Map<typename ShapeMatricesType::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                    temperature_size);

        auto local_Jac = MathLib::createZeroedMatrix<JacobianMatrix>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_rhs = MathLib::createZeroedVector<RhsVector>(
            local_rhs_data, local_matrix_size);

        typename ShapeMatricesType::template MatrixType<displacement_size,
                                                        temperature_size>
            KuT;
        KuT.setZero(displacement_size, temperature_size);

        typename ShapeMatricesType::NodalMatrixType KTT;
        KTT.setZero(temperature_size, temperature_size);

        typename ShapeMatricesType::NodalMatrixType DTT;
        DTT.setZero(temperature_size, temperature_size);

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

            // auto const [&](member){ return _process_data.member(t,
            // x_position); };
            auto rho_sr = _process_data.solid_density(t, x_position)[0];
            double const alpha = _process_data.linear_thermal_expansion_coefficient;
            double const c = _process_data.specific_heat_capacity;
            auto const lambda = _process_data.thermal_conductivity(t, x_position)[0];
            double const T0 = _process_data.reference_temperature;
            auto const& b = _process_data.specific_body_force;

            double T_ip = N.dot(T);
            double delta_T = T_ip - T0;
            // calculate real density
            auto rho_s = rho_sr * (1 - 3 * alpha * delta_T);
            // calculate thermally induced strain

            //
            // displacement equation, displacement part
            //
            _ip_data[ip].updateConstitutiveRelation(t, x_position, dt, u, alpha, delta_T);

            local_Jac
                .template block<displacement_size, displacement_size>(
                    displacement_index, displacement_index)
                .noalias() +=
                B.transpose() * C * B * w;

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
                (B.transpose() * (sigma + C * alpha * delta_T * Invariants::identity2)
                 - N_u.transpose() * rho_s * b) * w;
            local_rhs
                .template block<displacement_size, 1>(displacement_index, 0)
                .noalias() -=
                B.transpose() * C * alpha * T0 * Invariants::identity2 * w;

            //
            // displacement equation, temperature part
            //
            KuT.noalias() += B.transpose() * C * alpha * Invariants::identity2 * N * w;

            //
            // temperature equation, temperature part;
            //

            KTT.noalias() += dNdx.transpose() * lambda * dNdx * w;

            DTT.noalias() += N.transpose() * rho_s * c * N * w;

        }
        // temperature equation, temperature part
        local_Jac
            .template block<temperature_size, temperature_size>(
                temperature_index, temperature_index)
            .noalias() += KTT + DTT / dt;

        // displacement equation, temperature part
        local_Jac
            .template block<displacement_size, temperature_size>(
                displacement_index, temperature_index)
            .noalias() -= KuT;

        local_rhs.template block<temperature_size, 1>(temperature_index, 0)
           .noalias() -= KTT * T + DTT * T_dot;

        local_rhs.template block<displacement_size, 1>(displacement_index, 0)
           .noalias() += KuT * T;

        // Eigen::EigenSolver<JacobianMatrix> eigensolver(local_Jac);
        // std::cout << "eigenvalues" << eigensolver.eigenvalues() << "\n";

    }

    void preTimestepConcrete(std::vector<double> const& local_x,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        // Update damaged region.

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

    std::vector<double> const& getIntPtEpsilonXX(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 0);
    }

    std::vector<double> const& getIntPtEpsilonYY(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 1);
    }

    std::vector<double> const& getIntPtEpsilonZZ(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 2);
    }

    std::vector<double> const& getIntPtEpsilonXY(
        std::vector<double>& cache) const override
    {
        return getIntPtEpsilon(cache, 3);
    }

    std::vector<double> const& getIntPtEpsilonXZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 4);
    }

    std::vector<double> const& getIntPtEpsilonYZ(
        std::vector<double>& cache) const override
    {
        assert(DisplacementDim == 3);
        return getIntPtEpsilon(cache, 5);
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
    std::vector<double> const& getIntPtEpsilon(
        std::vector<double>& cache, std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data) {
            cache.push_back(ip_data._eps[component]);
        }

        return cache;
    }

    ThermoMechanicsProcessData<DisplacementDim>& _process_data;

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
    static const int displacement_index = ShapeFunction::NPOINTS;
    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim, int DisplacementDim>
class LocalAssemblerData final
    : public ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                      DisplacementDim>
{
public:
    LocalAssemblerData(LocalAssemblerData const&) = delete;
    LocalAssemblerData(LocalAssemblerData&&) = delete;

    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       ThermoMechanicsProcessData<DisplacementDim>& process_data)
        : ThermoMechanicsLocalAssembler<ShapeFunction, IntegrationMethod,
                                   DisplacementDim>(
              e, local_matrix_size, is_axially_symmetric, integration_order,
              process_data)
    {
    }
};

}  // namespace ThermoMechanics
}  // namespace ProcessLib

#endif  // PROCESS_LIB_THERMOMECHANICS_FEM_H_
