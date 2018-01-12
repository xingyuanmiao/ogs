/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "BaseLib/Functional.h"
#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"
#include "ProcessLib/Process.h"

#include "TMPhaseFieldFEM.h"
#include "TMPhaseFieldProcessData.h"
#include "TMPhaseFieldProcess.h"

namespace ProcessLib
{
namespace TMPhaseField
{
template <int DisplacementDim>
TMPhaseFieldProcess<DisplacementDim>::TMPhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    TMPhaseFieldProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
}

template <int DisplacementDim>
bool TMPhaseFieldProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_local_to_global_index_map_single_component);
    const int process_id = 0;  // Only one process in the monolithic scheme.
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
    std::generate_n(
        std::back_inserter(all_mesh_subsets),
        getProcessVariables(process_id)[1].get().getNumberOfComponents() +
            1,
        [&]() {
            return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
        });

    std::vector<int> const vec_n_components{1, DisplacementDim};
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION);
    assert(_local_to_global_index_map);
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, TMPhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table,
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
        _process_data);

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &TMPhaseFieldLocalAssemblerInterface::getIntPtSigmaXZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &TMPhaseFieldLocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtEpsilonXX));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtEpsilonYY));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtEpsilonZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtEpsilonXY));
    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &TMPhaseFieldLocalAssemblerInterface::getIntPtEpsilonYZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &TMPhaseFieldLocalAssemblerInterface::getIntPtEpsilonXZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "heat_flux",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &TMPhaseFieldLocalAssemblerInterface::getIntPtHeatFlux));
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble the equation for TMPhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleJacobian TMPhaseFieldProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_table, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, _coupled_solutions);
}
template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int /*process_id*/)
{
    DBUG("PreTimestep TMPhaseFieldProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &TMPhaseFieldLocalAssemblerInterface::preTimestep, _local_assemblers,
        *_local_to_global_index_map, x, t, dt);
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, int const process_id)
{
    DBUG("PostTimestep TMPhaseFieldProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &TMPhaseFieldLocalAssemblerInterface::postTimestep, _local_assemblers,
        *_local_to_global_index_map, x);
}

}  // namespace TMPhaseField
}  // namespace ProcessLib
