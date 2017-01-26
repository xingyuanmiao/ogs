/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreatePhaseFieldProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/CreateLinearElasticIsotropicPhaseField.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"

#include "PhaseFieldProcess.h"
#include "PhaseFieldProcessData.h"

namespace ProcessLib
{
namespace PhaseField
{
template <int DisplacementDim>
class PhaseFieldProcess;

extern template class PhaseFieldProcess<2>;
extern template class PhaseFieldProcess<3>;

template <int DisplacementDim>
std::unique_ptr<Process> createPhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "PHASE_FIELD");
    DBUG("Create PhaseFieldProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__PHASE_FIELD__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__PHASE_FIELD__process_variables__pressure}
         "phasefield",
        //! \ogs_file_param_special{prj__processes__process__PHASE_FIELD__process_variables__displacement}
         "displacement"});

    DBUG("Associate displacement with process variable \'%s\'.",
         process_variables[1].get().getName().c_str());

    if (process_variables[1].get().getNumberOfComponents() !=
        DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            process_variables[1].get().getName().c_str(),
            process_variables[1].get().getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate phase field with process variable \'%s\'.",
         process_variables[0].get().getName().c_str());
    if (process_variables[0].get().getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phase field process variable '%s' is not a scalar variable but has "
            "%d components.",
            process_variables[0].get().getName().c_str(),
            process_variables[0].get().getNumberOfComponents(),
            DisplacementDim);
    }


    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{process__PHASE_FIELD_constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__PHASE_FIELD__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::PhaseFieldExtension<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropicPhaseField")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropicPhaseField<DisplacementDim>(
                parameters, constitutive_relation_config);
    }
    else
    {
        OGS_FATAL(
            "Cannot construct constitutive relation of given type \'%s\'.",
            type.c_str());
    }

    // Residual stiffness
    const double residual_stiffness =
        //! \ogs_file_param_special{process__PHASE_FIELD_residual_stiffness}
        config.getConfigParameter<double>("residual_stiffness");
    // Crack resistance
    const double crack_resistance =
        //! \ogs_file_param_special{process__PHASE_FIELD_crack_resistance}
        config.getConfigParameter<double>("crack_resistance");
    // Crack length scale
    const double crack_length_scale =
        //! \ogs_file_param_special{process__PHASE_FIELD_crack_length_scale}
        config.getConfigParameter<double>("crack_length_scale");
    // Kinetic coefficient
    const double kinetic_coefficient =
        //! \ogs_file_param_special{process__PHASE_FIELD_kinetic_coefficient}
        config.getConfigParameter<double>("kinetic_coefficient");
    // Penalty constant
    const double penalty_constant =
        //! \ogs_file_param_special{process__PHASE_FIELD_penalty_constant}
        config.getConfigParameter<double>("penalty_constant");
    // Critical tolerance
    const double critical_tolerance =
        //! \ogs_file_param_special{process__PHASE_FIELD_critical_tolerance}
        config.getConfigParameter<double>("critical_tolerance");

    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__PHASE_FIELD_solid_density}
        "solid_density",
        parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.",
         solid_density.name.c_str());

    // Specific body force
    Eigen::Matrix<double, DisplacementDim, 1> specific_body_force;
    {
        std::vector<double> const b =
            //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__specific_body_force}
            config.getConfigParameter<std::vector<double>>(
                "specific_body_force");
        if (specific_body_force.size() != DisplacementDim)
            OGS_FATAL(
                "The size of the specific body force vector does not match the "
                "displacement dimension. Vector size is %d, displacement "
                "dimension is %d",
                specific_body_force.size(), DisplacementDim);

        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    PhaseFieldProcessData<DisplacementDim> process_data{
        std::move(material), residual_stiffness, crack_resistance,
        crack_length_scale, kinetic_coefficient, penalty_constant,
        critical_tolerance, solid_density, specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"PhaseField_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<PhaseFieldProcess<DisplacementDim>>{
        new PhaseFieldProcess<DisplacementDim>{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data),
            std::move(secondary_variables), std::move(named_function_caller)}};
}

template std::unique_ptr<Process> createPhaseFieldProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createPhaseFieldProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace PhaseField
}  // namespace ProcessLib
