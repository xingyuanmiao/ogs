/**
 * \file
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{

template <int DisplacementDim>
struct TMPhaseFieldExtension : public MechanicsBase<DisplacementDim>
{
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    virtual bool specialFunction(double const t,
                                 ProcessLib::SpatialPosition const& x,
                                 KelvinVector const& eps_m,
                                 double& strain_energy_tensile,
                                 KelvinVector& sigma_tensile,
                                 KelvinVector& sigma_compressive) const = 0;

    /// Dynamic size Kelvin vector and matrix wrapper for the polymorphic
    /// constitutive relation compute function.
    bool specialFunction(
        double const t,
        ProcessLib::SpatialPosition const& x,
        Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps_m,
        double& strain_energy_tensile,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& sigma_tensile,
        Eigen::Matrix<double, Eigen::Dynamic, 1>& sigma_compressive) const
    {
        // TODO Avoid copies of data:
        // Using MatrixBase<Derived> not possible because template functions
        // cannot be virtual. Maybe there is a workaround for this.  Using
        // Map<Matrix<double, ...>> makes the interface (for the material model
        // implementation) unnecessary difficult.
        KelvinVector const eps_m_{eps_m};
        KelvinVector sigma_tensile_{sigma_tensile};
        KelvinVector sigma_compressive_{sigma_compressive};

        bool const result = specialFunction(t,
                                            x,
                                            eps_m_,
                                            strain_energy_tensile,
                                            sigma_tensile_,
                                            sigma_compressive_);

        sigma_tensile = sigma_tensile_;
        sigma_compressive = sigma_compressive_;
        return result;
    }

    virtual ~TMPhaseFieldExtension() = default;
};

}  // namespace Solids
}  // namespace MaterialLib
