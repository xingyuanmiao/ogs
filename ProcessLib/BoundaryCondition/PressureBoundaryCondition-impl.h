/**
 * \file
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/MeshSearch/NodeSearch.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "PressureBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
namespace PressureBoundaryCondition
{
template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
PressureBoundaryCondition<LocalAssemblerImplementation>::
    PressureBoundaryCondition(
        bool const is_axially_symmetric, unsigned const integration_order,
        unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, unsigned const global_dim,
        std::vector<MeshLib::Element*>&& elements,
        Parameter<double> const& pressure)
    : _elements(std::move(elements)),
      _integration_order(integration_order),
      _pressure(pressure)
{
    std::vector<MeshLib::Node*> const nodes =
        MeshLib::getUniqueNodes(_elements);
    DBUG("Found %d nodes for Natural BCs for the variable %d", nodes.size(),
         variable_id);

    // Assume that the mesh subsets are equal for all components of the
    // variable.
    auto const& mesh_subsets = dof_table_bulk.getMeshSubsets(variable_id, 0);

    // TODO extend the node intersection to all parts of mesh_subsets, i.e.
    // to each of the MeshSubset in the mesh_subsets.
    _mesh_subset_all_nodes.reset(
        mesh_subsets.getMeshSubset(0).getIntersectionByNodes(nodes));
    std::unique_ptr<MeshLib::MeshSubsets> all_mesh_subsets{
        new MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}};

    // Create component ids vector for the current variable.
    auto const& number_of_components =
        dof_table_bulk.getNumberOfVariableComponents(variable_id);
    std::vector<int> component_ids(number_of_components);
    std::iota(std::begin(component_ids), std::end(component_ids), 0);

    // Create local DOF table from intersected mesh subsets for the given
    // variable and component ids.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, component_ids, std::move(all_mesh_subsets), _elements));

    createLocalAssemblers<LocalAssemblerImplementation>(
        global_dim, _elements, *_dof_table_boundary, shapefunction_order,
        _local_assemblers, is_axially_symmetric, _integration_order, _pressure);
}

template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
PressureBoundaryCondition<
    LocalAssemblerImplementation>::~PressureBoundaryCondition()
{
    for (auto e : _elements)
        delete e;
}

template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
void PressureBoundaryCondition<LocalAssemblerImplementation>::applyNaturalBC(
    const double t, const GlobalVector& x, GlobalMatrix& K, GlobalVector& b)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &PressureBoundaryConditionLocalAssemblerInterface::assemble,
        _local_assemblers, *_dof_table_boundary, t, x, K, b);
}

std::unique_ptr<
    PressureBoundaryCondition<PressureBoundaryConditionLocalAssembler>>
createPressureBoundaryCondition(
    BaseLib::ConfigTree const& config,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    bool is_axially_symmetric, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters)
{
    DBUG("Constructing PressureBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "Pressure");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Pressure__parameter}
    auto const parameter_name =
        config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s", parameter_name.c_str());

    auto const& pressure = findParameter<double>(parameter_name, parameters, 1);
    return std::unique_ptr<
        PressureBoundaryCondition<PressureBoundaryConditionLocalAssembler>>{
        new PressureBoundaryCondition<PressureBoundaryConditionLocalAssembler>{
            is_axially_symmetric, integration_order, shapefunction_order,
            dof_table, variable_id, global_dim, std::move(elements), pressure}};
}

}  // namespace PressureBoundaryCondition
}  // ProcessLib
