/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TMPhaseFieldProcess-fwd.h"
#include "TMPhaseFieldProcess.h"

namespace ProcessLib
{
namespace TMPhaseField
{

template class TMPhaseFieldProcess<2>;
template class TMPhaseFieldProcess<3>;

}   // namespace TMPhaseField
}   // namespace ProcessLib
