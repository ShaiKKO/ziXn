/*!
\file zx_integration.h
\brief Integration scenes (dam-break over porous bed) metrics API.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
*/

#ifndef ZX_INTEGRATION_H
#define ZX_INTEGRATION_H

#include <stdint.h>

#include "zx_abi.h"
#include "zx_tiles.h"
#include "zx_hbp.h"
#include "zx_mixture.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct zx_dambreak_params {
    uint32_t tiles;   /* number of tiles in pool (>=1) */
    float h;          /* grid spacing [m] */
    float dt;         /* timestep [s] */
    uint32_t steps;   /* number of steps */
    uint32_t bed_k;   /* bed thickness in node layers (0..B) */
    float init_head_u;/* initial x-velocity in reservoir [m/s] */
    float gamma_dot;  /* shear-rate estimate for mu_eff [1/s] */
    zx_hbp_params hbp;
    float permeability; /* [m^2] */
    float k_min;        /* [m^2] */
    float mu_min; float mu_max; /* for mu_eff clamp */
    float beta_min; float beta_max; /* for beta clamp */
    zx_mu_clamp_policy policy; float softness_k;
} zx_dambreak_params;

typedef struct zx_dambreak_metrics {
    float front_x;      /* estimated front advance [m] */
    float kinetic_j;    /* total kinetic energy [J], unit density */
} zx_dambreak_metrics;

/* Compute dam-break metrics on a single-tile coarse grid with HBP diffusion and Darcy drag in bed. */
ZX_API zx_dambreak_metrics ZX_CALL zx_integration_dambreak_run(const zx_dambreak_params* params);

typedef struct zx_bogging_params {
    float h; float dt; uint32_t steps;
    uint32_t wheel_radius_nodes; /* wheel radius in node units */
    float wheel_pull_u;          /* horizontal pull speed [m/s] */
    float wheel_push_w;          /* downward push speed [m/s] */
    zx_hbp_params hbp;           /* viscous diffusion */
    float mu_min; float mu_max;  /* viscosity clamp */
    float permeability; float k_min; /* Darcy */
    float beta_min; float beta_max;  /* beta clamp */
    float gamma_dot;                   /* HB-P shear rate estimate */
    zx_mu_clamp_policy policy; float softness_k;
} zx_bogging_params;

typedef struct zx_bogging_metrics {
    float sink_depth_m;  /* average sink depth under wheel [m] */
    float drag_N;        /* effective resisting force in x [N] (unit density proxy) */
} zx_bogging_metrics;

ZX_API zx_bogging_metrics ZX_CALL zx_integration_wheel_bogging_run(const zx_bogging_params* params);

typedef struct zx_puddle_params {
    float h; float dt; uint32_t steps;
    float init_head_u;         /* initial puddle push [m/s] */
    zx_hbp_params hbp; float mu_min; float mu_max;
    float permeability; float k_min; float beta_min; float beta_max;
    float gamma_dot; zx_mu_clamp_policy policy; float softness_k;
} zx_puddle_params;

typedef struct zx_puddle_metrics {
    float creep_dist_x; /* puddle creep distance [m] */
} zx_puddle_metrics;

ZX_API zx_puddle_metrics ZX_CALL zx_integration_puddle_creep_run(const zx_puddle_params* params);

#ifdef __cplusplus
}
#endif

#endif /* ZX_INTEGRATION_H */


