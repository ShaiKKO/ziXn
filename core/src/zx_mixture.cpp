/*!
 * \file zx_mixture.cpp
 * \brief Two-phase mixture helpers implementation.
 * \author Colin Macritchie (Ripple Group, LLC)
 */

#include "zx/zx_mixture.h"
#include <algorithm>

extern "C"
{

  // NOLINTBEGIN(bugprone-easily-swappable-parameters)
  float ZX_CALL zx_beta_darcy(float mu, float k, float k_min)
  // NOLINTEND(bugprone-easily-swappable-parameters)
  {
    const float kk = std::max(k, k_min);
    if (kk <= 0.0F)
    {
      return 0.0F;
    }
    return mu / kk;
  }

  // NOLINTBEGIN(bugprone-easily-swappable-parameters)
  float ZX_CALL zx_update_porosity(float phi, float solid_divergence, float dt, float phi_min,
                                   float phi_max)
  // NOLINTEND(bugprone-easily-swappable-parameters)
  {
    if (!(phi_min <= phi_max))
    {
      std::swap(phi_min, phi_max);
    }
    float phi_next = phi - (dt * (1.0F - phi) * solid_divergence);
    phi_next       = std::min(std::max(phi_next, phi_min), phi_max);
    return phi_next;
  }

  // NOLINTBEGIN(bugprone-easily-swappable-parameters)
  float ZX_CALL zx_bog_beta_darcy(float mu, float permeability, float k_min, float beta_min,
                                  float beta_max)
  // NOLINTEND(bugprone-easily-swappable-parameters)
  {
    float k_eff = std::max(permeability, k_min);
    float beta  = (k_eff > 0.0F) ? (mu / k_eff) : 0.0F;
    if (!(beta_min <= beta_max))
    {
      std::swap(beta_min, beta_max);
    }
    beta = std::min(std::max(beta, beta_min), beta_max);
    return beta;
  }

  // NOLINTBEGIN(bugprone-easily-swappable-parameters)
  float ZX_CALL zx_bog_beta_hbp(float gamma_dot, const zx_hbp_params* hbp, float permeability,
                                float k_min, float mu_min, float mu_max, float beta_min,
                                float beta_max, zx_mu_clamp_policy policy, float softness_k)
  // NOLINTEND(bugprone-easily-swappable-parameters)
  {
    float mu_eff = zx_hbp_mu_eff_policy(gamma_dot, hbp, mu_min, mu_max, policy, softness_k);
    return zx_bog_beta_darcy(mu_eff, permeability, k_min, beta_min, beta_max);
  }

  // NOLINTBEGIN(bugprone-easily-swappable-parameters)
  float ZX_CALL zx_effective_pressure(float p_total, float p_pore, float alpha_biot)
  // NOLINTEND(bugprone-easily-swappable-parameters)
  {
    const float a = std::min(std::max(alpha_biot, 0.0F), 1.0F);
    return p_total - (a * p_pore);
  }

  // NOLINTBEGIN(bugprone-easily-swappable-parameters)
  void ZX_CALL zx_mixture_drag_force(float beta, float vfx, float vfy, float vfz, float vsx,
                                     float vsy, float vsz, float volume, float dt, float* Fs_x,
                                     float* Fs_y, float* Fs_z, float* Ff_x, float* Ff_y,
                                     float* Ff_z)
  // NOLINTEND(bugprone-easily-swappable-parameters)
  {
    if ((Fs_x == nullptr) || (Fs_y == nullptr) || (Fs_z == nullptr) || (Ff_x == nullptr) ||
        (Ff_y == nullptr) || (Ff_z == nullptr))
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
    *Fs_x             = +fx;
    *Fs_y             = +fy;
    *Fs_z             = +fz;
    *Ff_x             = -fx;
    *Ff_y             = -fy;
    *Ff_z             = -fz;
  }

}  // extern "C"
