/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TMPhaseFieldExtension.h"
#include "LinearElasticIsotropic.h"
#include "KelvinVector.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
class LinearElasticIsotropicTMPhaseField final
    : public LinearElasticIsotropic<DisplacementDim>,
      public TMPhaseFieldExtension<DisplacementDim>
{
public:
    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<DisplacementDim>::value;
    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    explicit LinearElasticIsotropicTMPhaseField(
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
                         KelvinVector const& eps,
                         KelvinVector const& eps_m,
                         double& strain_energy_tensile,
                         KelvinVector& sigma_tensile,
                         KelvinVector& sigma_compressive,
                         KelvinMatrix& C_tensile,
                         KelvinMatrix& C_compressive,
                         KelvinVector& sigma_real,
                         double const degradation) const override
    {
        using Invariants =
            MaterialLib::SolidModels::Invariants<KelvinVectorSize>;
        // thermally induced strain
        // KelvinVector const eps_thermal = alpha * delta_T * Invariants::identity2;

        // calculation of deviatoric parts
        auto const& P_dev = Invariants::deviatoric_projection;
        KelvinVector const epsd_curr = P_dev * eps;
        KelvinVector const epsdm_curr = P_dev * eps_m;

        // Hydrostatic part for the stress and the tangent.
        double const eps_curr_trace = Invariants::trace(eps);
        double const epsm_curr_trace = Invariants::trace(eps_m);

        auto const& K =
            LinearElasticIsotropic<DisplacementDim>::_mp.bulk_modulus(t, x);
        auto const& mu =
            LinearElasticIsotropic<DisplacementDim>::_mp.mu(t, x);

        C_tensile = KelvinMatrix::Zero();
        C_compressive = KelvinMatrix::Zero();

        if (epsm_curr_trace >= 0)
        {
            strain_energy_tensile =
                 K / 2 * epsm_curr_trace * epsm_curr_trace +
                 mu * epsdm_curr.transpose() * epsdm_curr;
            sigma_tensile.noalias() =
                 (K * epsm_curr_trace * Invariants::identity2 +
                 2 * mu * epsdm_curr).eval();
            sigma_compressive.noalias() = KelvinVector::Zero();
            C_tensile.template topLeftCorner<3, 3>().setConstant(K);
            C_tensile.noalias() += 2 * mu * P_dev * KelvinMatrix::Identity();
        }
        else
        {
            strain_energy_tensile = mu * epsdm_curr.transpose() * epsdm_curr;
            sigma_tensile.noalias() = (2 * mu * epsdm_curr).eval();
            sigma_compressive.noalias() =
                 (K * epsm_curr_trace * Invariants::identity2).eval();
            C_tensile.noalias() = 2 * mu * P_dev * KelvinMatrix::Identity();
            C_compressive.template topLeftCorner<3, 3>().setConstant(K);
        }
        sigma_real.noalias() = degradation * sigma_tensile + sigma_compressive;
        return true;
    }
};

extern template class LinearElasticIsotropicTMPhaseField<2>;
extern template class LinearElasticIsotropicTMPhaseField<3>;

}  // namespace Solids
}  // namespace MaterialLib
