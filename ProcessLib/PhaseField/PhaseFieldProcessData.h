/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_PHASEFIELD_PHASEFIELDPROCESSDATA_H_
#define PROCESSLIB_PHASEFIELD_PHASEFIELDPROCESSDATA_H_

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace PhaseField
{
template <int DisplacementDim>
struct PhaseFieldProcessData
{
    PhaseFieldProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        double const residual_stiffness_,
        double const crack_resistance_,
        double const crack_length_scale_,
        double const kinetic_coefficient_,
        double const penalty_constant_,
        double const critical_tolerance_,
        Parameter<double> const& solid_density_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_)
        : material{std::move(material_)},
          residual_stiffness(residual_stiffness_),
          crack_resistance(crack_resistance_),
          crack_length_scale(crack_length_scale_),
          kinetic_coefficient(kinetic_coefficient_),
          penalty_constant(penalty_constant_),
          critical_tolerance(critical_tolerance_),
          solid_density(solid_density_),
          specific_body_force(specific_body_force_)
    {
    }

    PhaseFieldProcessData(PhaseFieldProcessData&& other)
        : material{std::move(other.material)},
          residual_stiffness(other.residual_stiffness),
          crack_resistance(other.crack_resistance),
          crack_length_scale(other.crack_length_scale),
          kinetic_coefficient(other.kinetic_coefficient),
          penalty_constant(other.penalty_constant),
          critical_tolerance(other.critical_tolerance),
          solid_density(other.solid_density),
          specific_body_force(other.specific_body_force),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    PhaseFieldProcessData(PhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(PhaseFieldProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    double const residual_stiffness;
    double const crack_resistance;
    double const crack_length_scale;
    double const kinetic_coefficient;
    double const penalty_constant;
    double const critical_tolerance;
    Parameter<double> const& solid_density;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt;
    double t;
};

}  // namespace PhaseField
}  // namespace ProcessLib

#endif  // PROCESSLIB_PHASEFIELD_PHASEFIELDPROCESSDATA_H_
