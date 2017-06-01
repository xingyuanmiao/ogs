/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <vector>


#include "HCProcessData.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace HC
{
template < typename NodalRowVectorType, typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType const& N_,
                         GlobalDimNodalMatrixType const& dNdx_,
                         double const& integration_weight_)
        : N(N_), dNdx(dNdx_), integration_weight(integration_weight_)
    {}

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;
};

const unsigned NUM_NODAL_DOF = 2;

class HCLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public HCLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType = typename ShapeMatricesType::template MatrixType<
        NUM_NODAL_DOF * ShapeFunction::NPOINTS,
        NUM_NODAL_DOF * ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<NUM_NODAL_DOF *
                                                        ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       HCProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _darcy_velocities(
              GlobalDim,
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        (void)local_matrix_size;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);
        }
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<LocalVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const num_nodes = ShapeFunction::NPOINTS;
        auto p_nodal_values =
            Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);

        auto const & b = _process_data.specific_body_force.head(GlobalDim);

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);

            auto const fluid_reference_density =
                _process_data.fluid_reference_density(t, pos)[0];

            auto const density_solid = _process_data.density_solid(t, pos)[0];
            // \todo the argument to getValue() has to be changed for non
            // constant storage model
            auto const specific_storage =
                _process_data.porous_media_properties.getSpecificStorage(t, pos)
                    .getValue(0.0);

            auto const thermal_conductivity_solid =
                _process_data.thermal_conductivity_solid(t, pos)[0];
            auto const thermal_conductivity_fluid =
                _process_data.thermal_conductivity_fluid(t, pos)[0];

            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First T, then P!
            NumLib::shapeFunctionInterpolate(local_x, N, T_int_pt, p_int_pt);

            // \todo the first argument has to be changed for non constant
            // porosity model
            auto const porosity =
                _process_data.porous_media_properties.getPorosity(t, pos)
                    .getValue(0.0, T_int_pt);

            double const thermal_conductivity =
                thermal_conductivity_solid * (1 - porosity) +
                 thermal_conductivity_fluid * porosity;

            auto const specific_heat_capacity_solid =
                _process_data.specific_heat_capacity_solid(t, pos)[0];
            auto const specific_heat_capacity_fluid =
                _process_data.specific_heat_capacity_fluid(t, pos)[0];

            auto const thermal_dispersivity_longitudinal =
                _process_data.thermal_dispersivity_longitudinal(t, pos)[0];
            auto const thermal_dispersivity_transversal =
                _process_data.thermal_dispersivity_transversal(t, pos)[0];

            auto Ktt = local_K.template block<num_nodes, num_nodes>(0, 0);
            auto Mtt = local_M.template block<num_nodes, num_nodes>(0, 0);
            auto Kpp = local_K.template block<num_nodes, num_nodes>(num_nodes,
                                                                    num_nodes);
            auto Mpp = local_M.template block<num_nodes, num_nodes>(num_nodes,
                                                                    num_nodes);
            auto Bp = local_b.template block<num_nodes, 1>(num_nodes, 0);

            // Use the fluid density model to compute the density
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::T)] = T_int_pt;
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;
            auto const density_water_T =
                _process_data.fluid_density->getValue(vars);

            auto const& K =
                _process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos);
            // Use the viscosity model to compute the viscosity
            auto const mu = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);

            GlobalDimMatrixType const K_over_mu = K / mu;

            GlobalDimVectorType const velocity =
                -perm_over_visc * (dNdx * p_nodal_values - density_water_T * b);

            double const velocity_magnitude = velocity.norm();
            GlobalDimMatrixType const& I(
                GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));
            GlobalDimMatrixType thermal_dispersivity =
                fluid_reference_density * specific_heat_capacity_fluid *
                (thermal_dispersivity_transversal * velocity_magnitude *
                     I +
                 (thermal_dispersivity_longitudinal -
                  thermal_dispersivity_transversal) /
                     velocity_magnitude * velocity * velocity.transpose());

            auto const hydrodynamic_thermodispersion =
                thermal_conductivity * I + thermal_dispersivity;

            double const heat_capacity =
                density_solid * specific_heat_capacity_solid * (1 - porosity) +
                fluid_reference_density * specific_heat_capacity_fluid * porosity;

            // matrix assembly
            Ktt.noalias() +=
                (dNdx.transpose() * hydrodynamic_thermodispersion * dNdx +
                 N.transpose() * velocity.transpose() * dNdx *
                     fluid_reference_density * specific_heat_capacity_fluid) *
                w;
            Kpp.noalias() += w * dNdx.transpose() * perm_over_visc * dNdx;
            Mtt.noalias() += w * N.transpose() * heat_capacity * N;
            Mpp.noalias() += w * N.transpose() * specific_storage * N;
            Bp += w * density_water_T * dNdx.transpose() * perm_over_visc * b;
            /* with Oberbeck-Boussing assumption density difference only exists
             * in buoyancy effects */

            // velocity computed for output.
            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = velocity[d];
            }
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 0);
        return _darcy_velocities[0];
    }

    std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 1);
        return _darcy_velocities[1];
    }

    std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 2);
        return _darcy_velocities[2];
    }

private:
    MeshLib::Element const& _element;
    HCProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<
            IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>>>
        _ip_data;
    std::vector<std::vector<double>> _darcy_velocities;
};

}  // namespace HC
}  // namespace ProcessLib
