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
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _process_data(std::move(process_data))
{
}

template <int DisplacementDim>
bool TMPhaseFieldProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
MathLib::MatrixSpecifications
TMPhaseFieldProcess<DisplacementDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the TM process (temperature-deformation) in the staggered
    // scheme.
    if (_use_monolithic_scheme || process_id == 1)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and phase field process.
    auto const& l = *_local_to_global_index_map_single_component;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_single_component};
}

template <int DisplacementDim>
NumLib::LocalToGlobalIndexMap const&
TMPhaseFieldProcess<DisplacementDim>::getDOFTable(const int process_id) const
{
    // If monolithic scheme is used or the equations of temperature-deformation is solved in
    // the staggered scheme.
    if (_use_monolithic_scheme || process_id == 1)
    {
        return *_local_to_global_index_map;
    }

    // For the equation of phasefield
    return *_local_to_global_index_map_single_component;
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

    if (_use_monolithic_scheme)
    {
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        all_mesh_subsets.emplace_back(_mesh_subset_all_nodes.get());
        all_mesh_subsets.emplace_back(_mesh_subset_all_nodes.get());

        const int monolithic_process_id = 0;  // Only one process in the monolithic scheme.
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(monolithic_process_id)[2].get().getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{1, 1, DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
    else
    {
        // For displacement equation.
        const int process_id = 0;
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(process_id)[0].get().getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{DisplacementDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);

        // For phase field equation.
        _sparsity_pattern_with_single_component =
            NumLib::computeSparsityPattern(
                *_local_to_global_index_map_single_component, _mesh);
    }
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, TMPhaseFieldLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order,
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
void TMPhaseFieldProcess<DisplacementDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_pf = 0;
        initializeProcessBoundaryCondition(*_local_to_global_index_map,
                                           process_id_of_pf);
        return;
    }

    // for the equations of temperature-deformation.
    const int process_id_of_u = 0;
    initializeProcessBoundaryCondition(*_local_to_global_index_map,
                                       process_id_of_u);
    // Staggered scheme:
    // for the phase field
    const int process_id_of_pf = 1;
    initializeProcessBoundaryCondition(
        *_local_to_global_index_map_single_component, process_id_of_pf);

}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble the equations for TMPhaseFieldProcess.");

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
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG(
            "AssembleJacobian TMPhaseFieldProcess for the monolithic scheme.");
        std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
            dof_table = {std::ref(*_local_to_global_index_map)};
        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, dof_table, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
            Jac, _coupled_solutions);
        return;
    }

    // For the staggered scheme
    if (_coupled_solutions->process_id == 1)
    {
        DBUG(
            "Assemble the Jacobian equations of phase field in "
            "TMPhaseFieldProcess for the staggered scheme.");
    }
    else
    {
        DBUG(
            "Assemble the Jacobian equations of temperature-deformation in "
            "TMPhaseFieldProcess for the staggered scheme.");
    }

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables = {std::ref(*_local_to_global_index_map),
                      std::ref(*_local_to_global_index_map_single_component)};
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
        Jac, _coupled_solutions);
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep TMPhaseFieldProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &TMPhaseFieldLocalAssemblerInterface::preTimestep, _local_assemblers,
        getDOFTable(process_id), x, t, dt);
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x, int const process_id)
{
    DBUG("PostTimestep TMPhaseFieldProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &TMPhaseFieldLocalAssemblerInterface::postTimestep, _local_assemblers,
        getDOFTable(process_id), x);
}

template <int DisplacementDim>
void TMPhaseFieldProcess<DisplacementDim>::postNonLinearSolverProcess(
    GlobalVector const& x, const double t, const int process_id)
{
    if (!_use_monolithic_scheme && process_id == 0)
    {
        return;
    }

    DBUG("PostNonLinearSolver TMPhaseFieldProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        getDOFTable(process_id), x, t, _use_monolithic_scheme);
}

}  // namespace TMPhaseField
}  // namespace ProcessLib
