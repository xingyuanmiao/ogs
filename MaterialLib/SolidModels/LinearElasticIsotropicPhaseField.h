/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PhaseFieldExtension.h"
#include "LinearElasticIsotropic.h"
#include "KelvinVector.h"
#include <iostream>

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
class LinearElasticIsotropicPhaseField final
    : public LinearElasticIsotropic<DisplacementDim>,
      public PhaseFieldExtension<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    explicit LinearElasticIsotropicPhaseField(
        typename LinearElasticIsotropic<
            DisplacementDim>::MaterialProperties const& material_properties)
        : LinearElasticIsotropic<DisplacementDim>(material_properties)
    {
    }

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() override
    {
        return LinearElasticIsotropic<
            DisplacementDim>::createMaterialStateVariables();
    }

    bool computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        KelvinVector& sigma,
        KelvinMatrix& C,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
            material_state_variables) override
    {
        return LinearElasticIsotropic<DisplacementDim>::
            computeConstitutiveRelation(t,
                                        x,
                                        dt,
                                        eps_prev,
                                        eps,
                                        sigma_prev,
                                        sigma,
                                        C,
                                        material_state_variables);
    }

    bool specialFunction(double const t,
                         ProcessLib::SpatialPosition const& x,
                         KelvinVector const& eps_prev,
                         KelvinVector const& eps,
                         double& strain_energy_tensile,
                         KelvinVector& sigma_tensile,
                         KelvinVector& sigma_compressive) const override
    {
        using Invariants =
            MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
        // calculation of deviatoric parts
        auto const& P_dev = Invariants::deviatoric_projection;
        // KelvinVector const epsd_i = P_dev * (eps - eps_prev);
        KelvinVector const epsd_curr = P_dev * eps;

        // Hydrostatic part for the stress and the tangent.
        // double const eps_i_trace = Invariants::trace(eps - eps_prev);
        double const eps_curr_trace = Invariants::trace(eps);

        auto const& K =
            LinearElasticIsotropic<DisplacementDim>::_mp.bulk_modulus(t, x);
        auto const& mu =
            LinearElasticIsotropic<DisplacementDim>::_mp.mu(t, x);

        /*strain_energy_tensile =
                K / 8 * (eps_curr_trace + std::abs(eps_curr_trace)) *
                                (eps_curr_trace + std::abs(eps_curr_trace)) +
                                mu * epsd_curr.transpose() * epsd_curr;

        sigma_tensile.noalias() =
                (K / 2 * (eps_curr_trace + std::abs(eps_curr_trace)) *
                                Invariants::identity2 + 2 * mu * epsd_curr).eval();

        sigma_compressive.noalias() =
                (K / 2 * (eps_curr_trace - std::abs(eps_curr_trace)) *
                                Invariants::identity2).eval();*/
        if (eps_curr_trace >= 0)
        {
            strain_energy_tensile =
                 K / 2 * eps_curr_trace * eps_curr_trace +
                 mu * epsd_curr.transpose() * epsd_curr;
            sigma_tensile.noalias() =
                 (K * eps_curr_trace * Invariants::identity2 +
                 2 * mu * epsd_curr).eval();
            sigma_compressive.noalias() = KelvinVector::Zero();
        }
        else
        {
            strain_energy_tensile = mu * epsd_curr.transpose() * epsd_curr;
            sigma_tensile.noalias() = (2 * mu * epsd_curr).eval();
            sigma_compressive.noalias() =
                 (K * eps_curr_trace * Invariants::identity2).eval();
        }
        // KelvinVector sigma = K * eps_curr_trace * Invariants::identity2 +
        //                      2 * mu *epsd_curr;
        // std::cout << "Tensile stress" << sigma_tensile << std::endl;
        // std::cout << "Compressive stress" << sigma_compressive << std::endl;
        // std::cout << "Stress" << sigma << std::endl;
        // std::cout << "strain energy tensile" << strain_energy_tensile << std::endl;
        return true;
    }
};

extern template class LinearElasticIsotropicPhaseField<2>;
extern template class LinearElasticIsotropicPhaseField<3>;

}  // namespace Solids
}  // namespace MaterialLib
