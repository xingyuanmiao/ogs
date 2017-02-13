/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "LinearElasticIsotropicTMPhaseField.h"
#include "ProcessLib/Utils/ProcessUtils.h"  // required for findParameter

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<LinearElasticIsotropicTMPhaseField<DisplacementDim>>
createLinearElasticIsotropicTMPhaseField(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{process__SMALL_DEFORMATION__constitutive_relation__type}
    config.checkConfigParameter("type", "LinearElasticIsotropicTMPhaseField");
    DBUG("Create LinearElasticIsotropicTMPhaseField material");

    // Youngs modulus
    auto& youngs_modulus = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{process__SMALL_DEFORMATION__constitutive_relation__LinearElasticIsotropic__youngs_modulus}
        config, "youngs_modulus", parameters, 1);

    DBUG("Use '%s' as youngs_modulus parameter.", youngs_modulus.name.c_str());

    // Poissons ratio
    auto& poissons_ratio = ProcessLib::findParameter<double>(
        //! \ogs_file_param_special{process__SMALL_DEFORMATION__constitutive_relation__LinearElasticIsotropic__poissons_ratio}
        config, "poissons_ratio", parameters, 1);

    DBUG("Use '%s' as poissons_ratio parameter.", poissons_ratio.name.c_str());

    typename LinearElasticIsotropicTMPhaseField<
        DisplacementDim>::MaterialProperties mp{youngs_modulus, poissons_ratio};

    return std::unique_ptr<LinearElasticIsotropicTMPhaseField<DisplacementDim>>{
        new LinearElasticIsotropicTMPhaseField<DisplacementDim>{mp}};
}

}  // namespace Solids
}  // namespace MaterialLib
