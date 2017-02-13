/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_CREATETMPHASEFIELDPROCESS_H_
#define PROCESS_LIB_CREATETMPHASEFIELDPROCESS_H_

#include "ProcessLib/Process.h"

namespace ProcessLib
{
namespace TMPhaseField
{
template <int DisplacementDim>
std::unique_ptr<Process> createTMPhaseFieldProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace TMPhaseField
}  // namespace ProcessLib

#endif  // PROCESS_LIB_CREATETMPHASEFIELDPROCESS_H_
