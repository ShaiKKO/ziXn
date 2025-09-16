/*!
\file zx_constitutive_ref.h
\brief CPU reference constitutive models: Elastic trial and DP/MC(+cap) return mapping.
\author Colin Macritchie (Ripple Group, LLC)
*/

#ifndef ZX_CONSTITUTIVE_REF_H
#define ZX_CONSTITUTIVE_REF_H

#include <stdint.h>
#include "zx_abi.h"

typedef struct zx_elastic_params {
    float young_E;
    float poisson_nu;
} zx_elastic_params;

typedef struct zx_mc_params {
    float friction_deg;  /* φ in degrees */
    float cohesion_kpa;  /* cohesion in kPa */
} zx_mc_params;

typedef struct zx_dp_params {
    float alpha;   /* slope in I1-J2 space */
    float k;       /* yield size */
} zx_dp_params;

typedef struct zx_cap_params {
    int   enabled;      /* 0/1 */
    float i1_min;       /* minimum I1 (compression cap), negative value */
} zx_cap_params;

/* Map Mohr–Coulomb (φ, c) to Drucker–Prager (alpha, k) using triaxial compression fit. */
ZX_API void ZX_CALL zx_mc_to_dp(const zx_mc_params* mc, zx_dp_params* out_dp);

/* Compute Lamé from (E, ν). */
ZX_API void ZX_CALL zx_elastic_lame(const zx_elastic_params* ep, float* lambda, float* mu);

/* Compute invariants and deviatoric part of stress. sigma is 3x3 row-major. */
ZX_API void ZX_CALL zx_stress_invariants(const float sigma[9], float* I1, float* J2);

/* Return map onto DP (+optional compression cap).
 * Inputs: elastic params, DP params, cap params, trial Cauchy stress (row-major).
 * Outputs: projected stress (row-major) and plastic multiplier Δγ.
 */
ZX_API void ZX_CALL zx_dp_return_map(
    const zx_elastic_params* ep,
    const zx_dp_params* dp,
    const zx_cap_params* cap,
    const float sigma_trial[9],
    float sigma_out[9],
    float* delta_gamma_out);

/* Mohr–Coulomb return mapping in principal stress space with non-associated flow.
 * Inputs: elastic params (for scaling), MC params, dilatancy_deg (ψ), trial Cauchy stress.
 * State: eps_v_pl (plastic volumetric strain) updated if non-null.
 */
ZX_API void ZX_CALL zx_mc_return_map(
    const zx_elastic_params* ep,
    const zx_mc_params* mc,
    float dilatancy_deg,
    const float sigma_trial[9],
    float sigma_out[9],
    float* delta_gamma_out,
    float* eps_v_pl /* in/out, may be null */);

/* NorSand parameters and state */
typedef struct zx_norsand_params {
    float M;             /* critical stress ratio */
    float lambda_cs;     /* CSL slope in e-ln p space */
    float kappa;         /* swelling slope */
    float p_ref;         /* reference pressure (Pa) */
    float e_ref;         /* void ratio at p_ref on CSL */
    float n_exp;         /* exponent for e_cs(p) = e_ref - lambda_cs * ln((p+p_ref)/p_ref)^n */
    float dilatancy_scale; /* scale for dilatancy vs state parameter */
} zx_norsand_params;

typedef struct zx_norsand_state {
    float void_ratio_e;  /* current void ratio */
} zx_norsand_state;

/* Critical state void ratio as a function of mean pressure p (compression negative by convention). */
ZX_API float ZX_CALL zx_norsand_e_cs(float p_mean, const zx_norsand_params* ns);

/* State parameter ψ = e - e_cs(p). */
ZX_API float ZX_CALL zx_norsand_state_parameter(const zx_norsand_state* st, float p_mean, const zx_norsand_params* ns);

/* NorSand return map (simplified) that adjusts stress toward critical stress ratio M and
 * updates volumetric plastic strain via a dilatancy law proportional to ψ.
 */
ZX_API void ZX_CALL zx_norsand_return_map(
    const zx_elastic_params* ep,
    const zx_norsand_params* ns,
    zx_norsand_state* state,
    const float sigma_trial[9],
    float sigma_out[9],
    float* delta_gamma_out,
    float* eps_v_pl /* in/out, may be null */);

#endif /* ZX_CONSTITUTIVE_REF_H */


