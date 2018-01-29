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
    typename ShapeMatrixType::GlobalDimVectorType heatflux, heatflux_prev;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

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
        heatflux_prev = heatflux;
        eps_m_prev = eps_m;
        eps_prev = eps;
        sigma_prev = sigma;
        sigma_real_prev = sigma_real;
        material_state_variables->pushBackState();
    }

    static const int kelvin_vector_size =
        KelvinVectorDimensions<DisplacementDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<kelvin_vector_size>;

    template <typename DisplacementVectorType>
    void updateConstitutiveRelation(double const t,
                                    SpatialPosition const& x_position,
                                    double const dt,
                                    DisplacementVectorType const& u,
                                    double const alpha,
                                    double const delta_T,
                                    double const degradation)
    {
        eps_m.noalias() = eps - alpha * delta_T * Invariants::identity2;
        solid_material.integrateStress(t, x_position, dt, eps_m_prev, eps_m,
                                       sigma_prev, *material_state_variables);

        static_cast<
            MaterialLib::Solids::TMPhaseFieldExtension<DisplacementDim>&>(
            solid_material)
            .calculateDegradedStress(t, x_position, eps_m,
                                     strain_energy_tensile, sigma_tensile,
                                     sigma_compressive, C_tensile,
                                     C_compressive, sigma_real, degradation);
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
        2 * ShapeFunction::NPOINTS + ShapeFunction::NPOINTS * DisplacementDim,
        DisplacementDim>;

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
        TMPhaseFieldProcessData<DisplacementDim>& process_data,
        int const mechanics_related_process_id,
        int const phase_field_process_id,
        int const heat_conduction_process_id)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric),
          _mechanics_related_process_id(mechanics_related_process_id),
          _phase_field_process_id(phase_field_process_id),
          _heat_conduction_process_id(heat_conduction_process_id)
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
            ip_data.sigma_prev.setZero(kelvin_vector_size);
            ip_data.eps.setZero(kelvin_vector_size);
            ip_data.eps_prev.setZero(kelvin_vector_size);
            ip_data.eps_m.setZero(kelvin_vector_size);
            ip_data.eps_m_prev.setZero(kelvin_vector_size);
            ip_data.C_tensile.setZero(kelvin_vector_size, kelvin_vector_size);
            ip_data.C_compressive.setZero(kelvin_vector_size,
                                          kelvin_vector_size);
            ip_data.sigma_tensile.setZero(kelvin_vector_size);
            ip_data.sigma_compressive.setZero(kelvin_vector_size);
            ip_data.heatflux.setZero();
            ip_data.heatflux_prev.setZero();
            ip_data.history_variable =
                process_data.history_field(0, x_position)[0];
            ip_data.history_variable_prev =
                process_data.history_field(0, x_position)[0];
            ip_data.sigma_real.setZero(kelvin_vector_size);

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
                              std::vector<double>& local_Jac_data) override;

    void assembleWithJacobianForStaggeredScheme(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions) override;

    /*void computeSecondaryVariableConcrete(
            const double t, std::vector<double> const& local_x) override
        {
            auto const local_matrix_size = local_x.size();
            assert(local_matrix_size == temperature_size + phasefield_size +
       displacement_size);

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
                auto const lambda = _process_data.thermal_conductivity(t,
       x_position)[0];
                auto const lambda_res =
       _process_data.residual_thermal_conductivity(t, x_position)[0];
                // heat flux only computed for output.
                GlobalDimVectorType const heat_flux = -(d_ip*d_ip*lambda + (1 -
       d_ip)*(1 - d_ip)*lambda_res) * dNdx * local_x_vec;

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
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
//        auto const kelvin_vector_size =
//            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, DisplacementDim, num_intpts);

        // TODO make a general implementation for converting KelvinVectors
        // back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& heatflux = _ip_data[ip].heatflux;

            for (typename KelvinVectorType::Index component = 0;
                 component < DisplacementDim;
                 ++component)
            {  // x, y, z components
                cache_mat(component, ip) = heatflux[component];
            }
        }

        return cache;
    }


//    std::vector<double> const& getIntPtHeatFlux(
//        const double t,
//        GlobalVector const& current_solution,
//        NumLib::LocalToGlobalIndexMap const& dof_table,
//        std::vector<double>& cache) const override
//    {
//        auto const num_intpts = _ip_data.size();

//        auto const indices = NumLib::getIndices(_element.getID(), dof_table);
//        assert(!indices.empty());
//        auto const local_x = current_solution.get(indices);

//        cache.clear();
//        auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
//            double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
//            cache, DisplacementDim, num_intpts);

//        SpatialPosition pos;
//        pos.setElementID(_element.getID());

//        auto T = Eigen::Map<typename ShapeMatricesType::template VectorType<
//            temperature_size> const>(local_x.data() + temperature_index,
//                                     temperature_size);

//        auto d = Eigen::Map<typename ShapeMatricesType::template VectorType<
//            phasefield_size> const>(local_x.data() + phasefield_index,
//                                    phasefield_size);

//        unsigned const n_integration_points =
//            _integration_method.getNumberOfPoints();

//        SpatialPosition x_position;
//        x_position.setElementID(_element.getID());
//        for (unsigned ip = 0; ip < n_integration_points; ip++)
//        {
//            x_position.setIntegrationPoint(ip);
//            auto const& dNdx = _ip_data[ip].dNdx;
//            auto const& N = _ip_data[ip].N;
//            double const d_ip = N.dot(d);
//            auto const lambda =
//                _process_data.thermal_conductivity(t, x_position)[0];
//            auto const lambda_res =
//                _process_data.residual_thermal_conductivity(t, x_position)[0];

//            // Compute the heat flux
//            cache_matrix.col(ip).noalias() =
//                -(d_ip * d_ip * lambda + (1 - d_ip) * (1 - d_ip) * lambda_res) *
//                dNdx * T;
//            std::cout << "phasefield" << d << std::endl;
//            std::cout << "temperature" << T << std::endl;
//        }

//        return cache;
//    }

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

    void assembleWithJacobianForDeformationEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    void assembleWithJacobianForHeatConductionEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

    void assembleWithJacobianForPhaseFieldEquations(
        double const t, std::vector<double> const& local_xdot,
        const double dxdot_dx, const double dx_dx,
        std::vector<double>& local_M_data, std::vector<double>& local_K_data,
        std::vector<double>& local_b_data, std::vector<double>& local_Jac_data,
        LocalCoupledSolutions const& local_coupled_solutions);

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

    /// ID of the processes that contains mechanical process.
    int const _mechanics_related_process_id;

    /// ID of phase field process.
    int const _phase_field_process_id;

    /// ID of heat conduction process.
    int const _heat_conduction_process_id;
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

#include "TMPhaseFieldFEM-impl.h"
