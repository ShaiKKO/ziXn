/*!
 * \file zx_mixture.cpp
 * \brief Two-phase mixture helpers implementation.
 * \author Colin Macritchie (Ripple Group, LLC)
 */

#include "zx/zx_mixture.h"
#include <algorithm>

extern "C" {

ZX_API float ZX_CALL zx_beta_darcy(float mu, float k, float k_min)
{
    const float kk = std::max(k, k_min);
    if (kk <= 0.0f) return 0.0f;
    return mu / kk;
}

ZX_API float ZX_CALL zx_update_porosity(float phi,
                                        float solid_divergence,
                                        float dt,
                                        float phi_min,
                                        float phi_max)
{
    if (!(phi_min <= phi_max)) std::swap(phi_min, phi_max);
    float phi_next = phi - dt * (1.0f - phi) * solid_divergence;
    phi_next = std::min(std::max(phi_next, phi_min), phi_max);
    return phi_next;
}

ZX_API float ZX_CALL zx_bog_beta_darcy(float mu,
                                       float permeability,
                                       float k_min,
                                       float beta_min,
                                       float beta_max)
{
    float k_eff = std::max(permeability, k_min);
    float beta = (k_eff > 0.0f) ? (mu / k_eff) : 0.0f;
    if (!(beta_min <= beta_max)) std::swap(beta_min, beta_max);
    beta = std::min(std::max(beta, beta_min), beta_max);
    return beta;
}

ZX_API float ZX_CALL zx_bog_beta_hbp(float gamma_dot,
                                     const zx_hbp_params* hbp,
                                     float permeability,
                                     float k_min,
                                     float mu_min,
                                     float mu_max,
                                     float beta_min,
                                     float beta_max,
                                     zx_mu_clamp_policy policy,
                                     float softness_k)
{
    float mu_eff = zx_hbp_mu_eff_policy(gamma_dot, hbp, mu_min, mu_max, policy, softness_k);
    return zx_bog_beta_darcy(mu_eff, permeability, k_min, beta_min, beta_max);
}

ZX_API float ZX_CALL zx_effective_pressure(float p_total,
                                           float p_pore,
                                           float alpha_biot)
{
    const float a = std::min(std::max(alpha_biot, 0.0f), 1.0f);
    return p_total - a * p_pore;
}

ZX_API void ZX_CALL zx_mixture_drag_force(float beta,
                                          float vfx, float vfy, float vfz,
                                          float vsx, float vsy, float vsz,
                                          float volume,
                                          float dt,
                                          float* Fs_x, float* Fs_y, float* Fs_z,
                                          float* Ff_x, float* Ff_y, float* Ff_z)
{
    if (!Fs_x || !Fs_y || !Fs_z || !Ff_x || !Ff_y || !Ff_z) return;
    const float vx = (vfx - vsx);
    const float vy = (vfy - vsy);
    const float vz = (vfz - vsz);
    const float coeff = std::max(0.0f, beta) * std::max(0.0f, volume) * std::max(0.0f, dt);
    const float Fx = coeff * vx;
    const float Fy = coeff * vy;
    const float Fz = coeff * vz;
    *Fs_x = +Fx; *Fs_y = +Fy; *Fs_z = +Fz;
    *Ff_x = -Fx; *Ff_y = -Fy; *Ff_z = -Fz;
}

} // extern "C"


