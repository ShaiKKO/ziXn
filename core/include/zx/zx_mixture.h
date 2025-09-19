/*!
\file zx_mixture.h
\brief Two-phase (soil-water) mixture helpers: effective stress, drag, porosity.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_MIXTURE_H
#define ZX_MIXTURE_H

#include <stdint.h>

#include "zx_abi.h"
#include "zx_hbp.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct zx_mixture_params
  {
    /* Biot coefficient (0..1) linking pore pressure to total stress */
    float alpha_biot;
    /* Dynamic viscosity of pore fluid [Pa·s] */
    float mu_fluid;
    /* Isotropic permeability [m^2] */
    float permeability;
    /* Fluid density [kg/m^3] */
    float rho_fluid;
    /* Porosity clamps */
    float phi_min;
    float phi_max;
  } zx_mixture_params;

  /** \brief Darcy drag coefficient β = μ / k, with k clamped to k_min. */
  ZX_API float ZX_CALL zx_beta_darcy(float mu, float k, float k_min);

  /** \brief Porosity update using mass conservation proxy.
   * phi_{t+dt} = clamp(phi - dt * (1-phi)*div_vs).
   */
  ZX_API float ZX_CALL zx_update_porosity(float phi, float solid_divergence, float dt,
                                          float phi_min, float phi_max);

  /** \brief Isotropic effective pressure: p_eff = p_total - alpha * p_pore. */
  ZX_API float ZX_CALL zx_effective_pressure(float p_total, float p_pore, float alpha_biot);

  /** \brief Symmetric drag forces on solid and fluid phases for a single point. */
  ZX_API void ZX_CALL zx_mixture_drag_force(float beta, float vfx, float vfy, float vfz, float vsx,
                                            float vsy, float vsz, float volume, float dt,
                                            /* out */ float* Fs_x, float* Fs_y, float* Fs_z,
                                            /* out */ float* Ff_x, float* Ff_y, float* Ff_z);

  /** \brief Expose bog memory knob β ≈ μ/k with explicit clamps.
   * Units: beta [Pa·s/m^2] when mu in [Pa·s], k in [m^2].
   */
  ZX_API float ZX_CALL zx_bog_beta_darcy(float mu, float permeability, float k_min, float beta_min,
                                         float beta_max);

  /** \brief Bog memory beta using HB-P effective viscosity at shear-rate gamma_dot.
   * Computes mu_eff via zx_hbp_mu_eff_policy, then returns clamped beta = mu_eff / max(k,k_min),
   * further clamped into [beta_min, beta_max].
   */
  ZX_API float ZX_CALL zx_bog_beta_hbp(float gamma_dot, const zx_hbp_params* hbp,
                                       float permeability, float k_min, float mu_min, float mu_max,
                                       float beta_min, float beta_max, zx_mu_clamp_policy policy,
                                       float softness_k);

#ifdef __cplusplus
}
#endif

#endif /* ZX_MIXTURE_H */
