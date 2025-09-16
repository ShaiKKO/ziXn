/*!
\file zx_hbp.h
\brief Herschel–Bulkley–Papanastasiou viscous update on coarse grid.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_HBP_H
#define ZX_HBP_H

#include <stdint.h>

#include "zx_abi.h"
#include "zx_tiles.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct zx_hbp_params {
    /* Base (solvent) viscosity [Pa·s] */
    float mu0;
    /* Consistency index [Pa·s^n] */
    float K;
    /* Flow behavior index (n<1 shear-thinning, n>1 shear-thickening) [-] */
    float n;
    /* Yield stress [Pa] */
    float tau_y;
    /* Papanastasiou regularization parameter m [s] */
    float m;
} zx_hbp_params;

typedef enum zx_mu_clamp_policy {
    ZX_MU_CLAMP_NONE = 0,
    ZX_MU_CLAMP_HARD = 1,
    ZX_MU_CLAMP_SMOOTH_TANH = 2
} zx_mu_clamp_policy;

/*!
 * \brief Compute regularized effective viscosity for HB-P model.
 *
 * Uses: mu_eff = mu0 + K * max(gamma_dot,eps)^(n-1)
 *              + tau_y * f(gamma_dot), with f = (1-exp(-m*gamma_dot))/max(gamma_dot,eps)
 * and the correct limit f(0) = m. Clamps into [mu_min, mu_max].
 */
ZX_API float ZX_CALL zx_hbp_mu_eff(float gamma_dot,
                                   const zx_hbp_params* params,
                                   float mu_min,
                                   float mu_max);

/* Variant with explicit clamp policy and softness parameter (k>0 for smoother clamp). */
ZX_API float ZX_CALL zx_hbp_mu_eff_policy(float gamma_dot,
                                          const zx_hbp_params* params,
                                          float mu_min,
                                          float mu_max,
                                          zx_mu_clamp_policy policy,
                                          float softness_k);

/* Validate and sanitize parameters (returns 0 if unchanged, 1 if adjusted, <0 on error). */
ZX_API int ZX_CALL zx_hbp_validate_params(zx_hbp_params* inout);

/* CFL-like stability upper bound for explicit diffusion: dt_max = (h^2 * rho_min) / (6 * mu_max). */
ZX_API float ZX_CALL zx_hbp_dt_stable_upper_bound(float h, float mu_max, float rho_min);

/*!
 * \brief Perform an explicit viscous diffusion step on grid node velocities.
 *
 * Updates node momenta in-place using a 6-point Laplacian stencil over a tile-sparse grid.
 * The Laplacian operates on velocities v = mom/mass. Density is lumped per node
 * as rho = mass / h^3. Missing neighbors use zero-Neumann (copy-center) boundary.
 *
 * Stable for dt <= h^2/(6*nu_max) with nu_max = mu_max / rho_min.
 */
ZX_API void ZX_CALL zx_hbp_update_coarse_grid(zx_tile* pool,
                                              uint32_t tile_count,
                                              float h,
                                              float dt,
                                              const zx_hbp_params* params,
                                              float mu_min,
                                              float mu_max);

#ifdef __cplusplus
}
#endif

#endif /* ZX_HBP_H */


