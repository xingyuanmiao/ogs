/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_TMPHASEFIELD_TMPHASEFIELDPROCESSDATA_H_
#define PROCESSLIB_TMPHASEFIELD_TMPHASEFIELDPROCESSDATA_H_

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace TMPhaseField
{
template <int DisplacementDim>
struct TMPhaseFieldProcessData
{
    TMPhaseFieldProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& residual_stiffness_,
        Parameter<double> const& crack_resistance_,
        Parameter<double> const& crack_length_scale_,
        Parameter<double> const& kinetic_coefficient_,
        Parameter<double> const& solid_density_,
        Parameter<double> const& linear_thermal_expansion_coefficient_,
        Parameter<double> const& specific_heat_capacity_,
        Parameter<double> const& thermal_conductivity_,
        double const reference_temperature_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_)
        : material{std::move(material_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          kinetic_coefficient(kinetic_coefficient_),
          solid_density(solid_density_),
          linear_thermal_expansion_coefficient(linear_thermal_expansion_coefficient_),
          specific_heat_capacity(specific_heat_capacity_),
          thermal_conductivity(thermal_conductivity_),
          reference_temperature(reference_temperature_),
          specific_body_force(specific_body_force_)
    {
    }

    TMPhaseFieldProcessData(TMPhaseFieldProcessData&& other)
        : material{std::move(other.material)},
          residual_stiffness(other.residual_stiffness),
          crack_resistance(other.crack_resistance),
          crack_length_scale(other.crack_length_scale),
          kinetic_coefficient(other.kinetic_coefficient),
          solid_density(other.solid_density),
          linear_thermal_expansion_coefficient(other.linear_thermal_expansion_coefficient),
          specific_heat_capacity(other.specific_heat_capacity),
          thermal_conductivity(other.thermal_conductivity),
          reference_temperature(other.reference_temperature),
          specific_body_force(other.specific_body_force),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    TMPhaseFieldProcessData(TMPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TMPhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TMPhaseFieldProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    Parameter<double> const& residual_stiffness;
    Parameter<double> const& crack_resistance;
    Parameter<double> const& crack_length_scale;
    Parameter<double> const& kinetic_coefficient;
    Parameter<double> const& solid_density;
    Parameter<double> const& linear_thermal_expansion_coefficient;
    Parameter<double> const& specific_heat_capacity;
    Parameter<double> const& thermal_conductivity;
    double const reference_temperature;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt;
    double t;
};

}  // namespace TMPhaseField
}  // namespace ProcessLib

#endif  // PROCESSLIB_TMPHASEFIELD_TMPHASEFIELDPROCESSDATA_H_
