/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTMPhaseFieldProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/CreateLinearElasticIsotropicTMPhaseField.h"
#include "ProcessLib/Utils/ParseSecondaryVariables.h"

#include "TMPhaseFieldProcess.h"
#include "TMPhaseFieldProcessData.h"

namespace ProcessLib
{
namespace TMPhaseField
{
template <int DisplacementDim>
class TMPhaseFieldProcess;

extern template class TMPhaseFieldProcess<2>;
extern template class TMPhaseFieldProcess<3>;

template <int DisplacementDim>
std::unique_ptr<Process> createTMPhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__type}
    config.checkConfigParameter("type", "TMPHASE_FIELD");
    DBUG("Create TMPhaseFieldProcess.");

    // Process variable.

    //! \ogs_file_param{prj__processes__process__TMPHASE_FIELD__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__TMPHASE_FIELD__process_variables__temperature}
         "temperature",
         //! \ogs_file_param_special{prj__processes__process__TMPHASE_FIELD__process_variables__phasefield}
         "phasefield",
         //! \ogs_file_param_special{prj__processes__process__TMPHASE_FIELD__process_variables__displacement}
         "displacement"});

    DBUG("Associate displacement with process variable \'%s\'.",
         process_variables[2].get().getName().c_str());

    if (process_variables[2].get().getNumberOfComponents() !=
        DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the displacement dimension: got %d, expected %d",
            process_variables[2].get().getName().c_str(),
            process_variables[2].get().getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate phase field with process variable \'%s\'.",
         process_variables[1].get().getName().c_str());
    if (process_variables[1].get().getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Phase field process variable '%s' is not a scalar variable but has "
            "%d components.",
            process_variables[1].get().getName().c_str(),
            process_variables[1].get().getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate temperature with process variable \'%s\'.",
         process_variables[0].get().getName().c_str());
    if (process_variables[0].get().getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Temperature process variable '%s' is not a scalar variable but has "
            "%d components.",
            process_variables[0].get().getName().c_str(),
            process_variables[0].get().getNumberOfComponents(),
            DisplacementDim);
    }


    // Constitutive relation.
    // read type;
    auto const constitutive_relation_config =
        //! \ogs_file_param{process__TMPHASE_FIELD_constitutive_relation}
        config.getConfigSubtree("constitutive_relation");

    auto const type =
        //! \ogs_file_param{prj__processes__process__TMPHASE_FIELD__constitutive_relation__type}
        constitutive_relation_config.peekConfigParameter<std::string>("type");

    std::unique_ptr<MaterialLib::Solids::TMPhaseFieldExtension<DisplacementDim>>
        material = nullptr;
    if (type == "LinearElasticIsotropicTMPhaseField")
    {
        material =
            MaterialLib::Solids::createLinearElasticIsotropicTMPhaseField<DisplacementDim>(
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
        //! \ogs_file_param_special{process__TMPHASE_FIELD_residual_stiffness}
        config.getConfigParameter<double>("residual_stiffness");
    // Crack resistance
    const double crack_resistance =
        //! \ogs_file_param_special{process__TMPHASE_FIELD_crack_resistance}
        config.getConfigParameter<double>("crack_resistance");
    // Crack length scale
    const double crack_length_scale =
        //! \ogs_file_param_special{process__TMPHASE_FIELD_crack_length_scale}
        config.getConfigParameter<double>("crack_length_scale");
    // Kinetic coefficient
    const double kinetic_coefficient =
        //! \ogs_file_param_special{process__TMPHASE_FIELD_kinetic_coefficient}
        config.getConfigParameter<double>("kinetic_coefficient");
    // Penalty constant
    const double penalty_constant =
        //! \ogs_file_param_special{process__TMPHASE_FIELD_penalty_constant}
        config.getConfigParameter<double>("penalty_constant");
    // Critical tolerance
    const double critical_tolerance =
        //! \ogs_file_param_special{process__TMPHASE_FIELD_critical_tolerance}
        config.getConfigParameter<double>("critical_tolerance");

    // Solid density
    auto& solid_density = findParameter<double>(
        config,
        //! \ogs_file_param_special{process__TMPHASE_FIELD_solid_density}
        "solid_density",
        parameters, 1);
    DBUG("Use \'%s\' as solid density parameter.",
         solid_density.name.c_str());

    // Linear thermal expansion coefficient
    const double linear_thermal_expansion_coefficient =
        //! \ogs_file_param_special{process__TMPHASE_FIELD_linear_thermal_expansion_coefficient}
        config.getConfigParameter<double>("linear_thermal_expansion_coefficient");
    // Specific heat capacity
    const double specific_heat_capacity =
        //! \ogs_file_param_special{process__TMPHASE_FIELD_specific_heat_capacity}
        config.getConfigParameter<double>("specific_heat_capacity");
    // Thermal conductivity
    auto& thermal_conductivity = findParameter<double>(
            config,
            //! \ogs_file_param_special{process__TMPHASE_FIELD_thermal_conductivity}
            "thermal_conductivity", parameters, 1);
        DBUG("Use \'%s\' as thermal conductivity parameter.",
             thermal_conductivity.name.c_str());
    // Reference temperature
    const double reference_temperature =
        //! \ogs_file_param_special{process__TMPHASE_FIELD_reference_temperature}
        config.getConfigParameter<double>("reference_temperature");

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

    TMPhaseFieldProcessData<DisplacementDim> process_data{
        std::move(material), residual_stiffness, crack_resistance,
        crack_length_scale, kinetic_coefficient, penalty_constant,
        critical_tolerance, solid_density, linear_thermal_expansion_coefficient,
        specific_heat_capacity, thermal_conductivity, reference_temperature,
        specific_body_force};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"PhaseField_temperature_displacement"});

    ProcessLib::parseSecondaryVariables(config, secondary_variables,
                                        named_function_caller);

    return std::unique_ptr<TMPhaseFieldProcess<DisplacementDim>>{
        new TMPhaseFieldProcess<DisplacementDim>{
            mesh, std::move(jacobian_assembler), parameters, integration_order,
            std::move(process_variables), std::move(process_data),
            std::move(secondary_variables), std::move(named_function_caller)}};
}

template std::unique_ptr<Process> createTMPhaseFieldProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createTMPhaseFieldProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace TMPhaseField
}  // namespace ProcessLib
