/*!
 * \file zx_mixture.cpp
 * \brief Two-phase mixture helpers implementation.
 * \author Colin Macritchie (Ripple Group, LLC)
 */

#include "zx/zx_mixture.h"
#include <algorithm>

extern "C"
{

  /**
   * @brief Darcy drag coefficient β = μ / max(k,k_min).
   * @param mu Dynamic viscosity [Pa·s].
   * @param k Permeability [m^2].
   * @param k_min Minimum permeability clamp [m^2].
   * @return β in [Pa·s/m^2] (0 if clamped denominator is non-positive).
   */
  float ZX_CALL zx_beta_darcy(float mu, float k, float k_min)
  {
    const float kk = std::max(k, k_min);
    if (kk <= 0.0F)
    {
      return 0.0F;
    }
    return mu / kk;
  }

  /**
   * @brief Porosity update proxy with clamp.
   * @param phi Current porosity [0,1].
   * @param solid_divergence Divergence of solid velocity (1/s).
   * @param dt Time step [s].
   * @param phi_min Minimum porosity clamp.
   * @param phi_max Maximum porosity clamp.
   * @return Next porosity after update and clamp.
   */
  float ZX_CALL zx_update_porosity(float phi, float solid_divergence, float dt, float phi_min,
                                   float phi_max)
  {
    if (!(phi_min <= phi_max))
    {
      std::swap(phi_min, phi_max);
    }
    float phi_next = phi - (dt * ((1.0F - phi) * solid_divergence));
    phi_next       = std::clamp(phi_next, phi_min, phi_max);
    return phi_next;
  }

  /**
   * @brief Bog memory beta via Darcy with clamps.
   * @param mu Dynamic viscosity [Pa·s].
   * @param permeability Permeability [m^2].
   * @param k_min Minimum permeability clamp [m^2].
   * @param beta_min Minimum beta clamp.
   * @param beta_max Maximum beta clamp.
   * @return Clamped Darcy beta.
   */
  float ZX_CALL zx_bog_beta_darcy(float mu, float permeability, float k_min, float beta_min,
                                  float beta_max)
  {
    const float k_eff = std::max(permeability, k_min);
    float beta        = (k_eff > 0.0F) ? (mu / k_eff) : 0.0F;
    if (!(beta_min <= beta_max))
    {
      std::swap(beta_min, beta_max);
    }
    beta = std::clamp(beta, beta_min, beta_max);
    return beta;
  }

  /**
   * @brief Bog memory beta via HB-P effective viscosity.
   * @param gamma_dot Shear rate proxy.
   * @param hbp HB-P parameters.
   * @param permeability Permeability [m^2].
   * @param k_min Minimum permeability clamp [m^2].
   * @param mu_min Minimum viscosity clamp.
   * @param mu_max Maximum viscosity clamp.
   * @param beta_min Minimum beta clamp.
   * @param beta_max Maximum beta clamp.
   * @param policy Viscosity clamp policy.
   * @param softness_k Softness factor for clamp.
   * @return Clamped Darcy beta computed from effective viscosity.
   */
  float ZX_CALL zx_bog_beta_hbp(float gamma_dot, const zx_hbp_params* hbp, float permeability,
                                float k_min, float mu_min, float mu_max, float beta_min,
                                float beta_max, zx_mu_clamp_policy policy, float softness_k)
  {
    float mu_eff = zx_hbp_mu_eff_policy(gamma_dot, hbp, mu_min, mu_max, policy, softness_k);
    return zx_bog_beta_darcy(mu_eff, permeability, k_min, beta_min, beta_max);
  }

  /**
   * @brief Isotropic effective pressure p_eff = p_total - alpha * p_pore.
   * @param p_total Total (Cauchy) pressure.
   * @param p_pore Pore pressure.
   * @param alpha_biot Biot coefficient [0,1].
   * @return Effective pressure.
   */
  float ZX_CALL zx_effective_pressure(float p_total, float p_pore, float alpha_biot)
  {
    const float a = std::clamp(alpha_biot, 0.0F, 1.0F);
    return p_total - (a * p_pore);
  }

  /**
   * @brief Symmetric drag forces on solid and fluid phases for a single point.
   * @param beta Drag coefficient.
   * @param vfx Fluid velocity components.
   * @param vfy Fluid velocity components.
   * @param vfz Fluid velocity components.
   * @param vsx Solid velocity components.
   * @param vsy Solid velocity components.
   * @param vsz Solid velocity components.
   * @param volume Cell/point volume.
   * @param dt Time step.
   * @param[out] fs_x Solid drag x-component.
   * @param[out] fs_y Solid drag y-component.
   * @param[out] fs_z Solid drag z-component.
   * @param[out] ff_x Fluid drag x-component.
   * @param[out] ff_y Fluid drag y-component.
   * @param[out] ff_z Fluid drag z-component.
   */
  void ZX_CALL zx_mixture_drag_force(float beta, float vfx, float vfy, float vfz, float vsx,
                                     float vsy, float vsz, float volume, float dt, float* fs_x,
                                     float* fs_y, float* fs_z, float* ff_x, float* ff_y,
                                     float* ff_z)
  {
    if ((fs_x == nullptr) || (fs_y == nullptr) || (fs_z == nullptr) || (ff_x == nullptr) ||
        (ff_y == nullptr) || (ff_z == nullptr))
    {
      return;
    }
    const float vx    = (vfx - vsx);
    const float vy    = (vfy - vsy);
    const float vz    = (vfz - vsz);
    const float coeff = std::max(0.0F, beta) * std::max(0.0F, volume) * std::max(0.0F, dt);
    const float fx    = coeff * vx;
    const float fy    = coeff * vy;
    const float fz    = coeff * vz;
    *fs_x             = +fx;
    *fs_y             = +fy;
    *fs_z             = +fz;
    *ff_x             = -fx;
    *ff_y             = -fy;
    *ff_z             = -fz;
  }

}  // extern "C"
