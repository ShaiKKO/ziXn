/*!
 * \file zx_integration.cpp
 * \brief Dam-break over porous bed (coarse proxy) metrics implementation.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_integration.h"
#include "zx/zx_residency.h"
#include "zx/zx_lod.h"
#include <algorithm>
#include <cmath>

extern "C" {

/**
 * @brief Simulate a single-tile coarse-grid dam-break and compute front advance and kinetic energy.
 *
 * Runs a short HBP-based diffusion on a single tile initialized with the left half moving
 * (initial head velocity) and the right half at rest. Applies immediate Darcy-like damping
 * in the bed region, advances the coarse-grid proxy for p->steps time steps, and then
 * computes:
 *  - front_x: the maximum x-coordinate (distance) where horizontal velocity vx > 0.05,
 *    returned in the same length units as p->h;
 *  - kinetic_j: total kinetic energy sum(0.5 * m * v^2) over all nodes.
 *
 * If p is NULL or any of p->tiles, p->h, p->dt, or p->steps is zero/non-positive, a
 * zero-initialized zx_dambreak_metrics is returned.
 *
 * @param p Dam-break parameters (must not be NULL). Relevant fields used include:
 *          - p->h: grid spacing (length unit for front_x),
 *          - p->dt: time step,
 *          - p->steps: number of diffusion steps,
 *          - p->init_head_u: initial velocity applied to the left half,
 *          - p->bed_k and other HBP-related fields for damping behavior.
 * @return zx_dambreak_metrics Metrics containing:
 *         - front_x: front advance distance (same units as p->h),
 *         - kinetic_j: total kinetic energy (energy units implied by mass/velocity).
 */
ZX_API zx_dambreak_metrics ZX_CALL zx_integration_dambreak_run(const zx_dambreak_params* p)
{
    zx_dambreak_metrics M{0.0f, 0.0f};
    if (!p || p->tiles == 0 || p->h <= 0.0f || p->dt <= 0.0f || p->steps == 0) return M;
    // Allocate a single-tile pool and simple initial condition: left half moving, right half at rest
    zx_tile tile{}; tile.coord_x = 0; tile.coord_y = 0; tile.coord_z = 0;
    const uint32_t B = (uint32_t)ZX_TILE_B;
    auto idx = [&](uint32_t ix, uint32_t iy, uint32_t iz){ return ix + B*(iy + B*iz); };
    for (uint32_t iz=0; iz<B; ++iz) for (uint32_t iy=0; iy<B; ++iy) for (uint32_t ix=0; ix<B; ++ix){
        auto& n = tile.nodes[idx(ix,iy,iz)];
        n.mass = 1.0f;
        const bool left = ix < (B/2);
        const bool bed  = iz < std::min<uint32_t>(p->bed_k, B);
        const float vx0 = left ? p->init_head_u : 0.0f;
        n.mom_x = vx0 * n.mass; n.mom_y = 0.0f; n.mom_z = 0.0f;
        // Apply immediate Darcy drag in bed region to mimic porous damping
        if (bed) {
            float beta = zx_bog_beta_hbp(p->gamma_dot, &p->hbp, p->permeability, p->k_min,
                                         p->mu_min, p->mu_max, p->beta_min, p->beta_max, p->policy, p->softness_k);
            const float rho = n.mass / (p->h*p->h*p->h);
            const float ax = -(beta / std::max(rho, 1e-12f)) * vx0; // dv/dt = -(beta/rho)*v
            const float vx1 = vx0 + p->dt * ax;
            n.mom_x = std::max(0.0f, vx1) * n.mass;
        }
    }

    // Diffuse via HBP for 'steps' with optional LOD fallback pressure loop
    zx_tile pool[1]; pool[0] = tile;
    zx_residency_opts ro{1,2,1}; zx_residency* rez = zx_residency_create(&ro);
    uint32_t churn_enter=0, churn_exit=0, pf_count=0;
    zx_lod_fallback_state fbs; zx_lod_fallback_init(&fbs);
    zx_lod_fallback_policy fbp; zx_lod_get_default_policy(&fbp);
    for (uint32_t s=0; s<p->steps; ++s){
        // If fallback enabled, track active frames via residency and policy
        if (zx_lod_is_enabled()){
            uint32_t en=0, ex=0, pf=0; zx_residency_tick(rez, (int)s,0,0, 0, &en,&ex,&pf);
            churn_enter+=en; churn_exit+=ex; pf_count=pf;
            uint32_t active = zx_residency_get_active_count(rez);
            (void)zx_lod_fallback_update(&fbp, active, /*last_step_ms*/1.0f, &fbs);
        } else {
            uint32_t en=0, ex=0, pf=0; zx_residency_tick(rez, (int)s,0,0, 0, &en,&ex,&pf); churn_enter+=en; churn_exit+=ex; pf_count=pf;
        }
        zx_hbp_update_coarse_grid(pool, 1, p->h, p->dt, &p->hbp, p->mu_min, p->mu_max);
    }

    // Metrics: front advance via threshold on vx, and kinetic energy sum
    float max_ix = 0.0f; float ke = 0.0f;
    for (uint32_t iz=0; iz<B; ++iz) for (uint32_t iy=0; iy<B; ++iy) for (uint32_t ix=0; ix<B; ++ix){
        const auto& n = pool[0].nodes[idx(ix,iy,iz)];
        const float m = n.mass; if (m <= 0.0f) continue;
        const float vx = n.mom_x / m; const float vy = n.mom_y / m; const float vz = n.mom_z / m;
        const float v2 = vx*vx + vy*vy + vz*vz;
        if (vx > 0.05f) max_ix = std::max(max_ix, (float)ix);
        ke += 0.5f * m * v2;
    }
    M.front_x = max_ix * p->h;
    M.kinetic_j = ke;
    // Could export churn/fallback counters via telemetry externally
    (void)churn_enter; (void)churn_exit; (void)pf_count;
    zx_residency_destroy(rez);
    return M;
}

/**
 * @brief Simulates wheel-induced bogging on a single coarse-grid tile and returns drag and sink metrics.
 *
 * Runs a proxy simulation where a circular wheel footprint applies a horizontal pull velocity
 * and a vertical push onto a coarse BxBxB tile. Local Darcy-like drag is applied per-node,
 * the coarse-grid viscous/diffusive update is advanced for the specified number of steps,
 * and final surface-layer states are reduced into simple metrics.
 *
 * The function performs input validation and returns zeroed metrics if `p` is null or if
 * `p->h`, `p->dt`, or `p->steps` are non-positive.
 *
 * @param p Pointer to bogging parameters (must not be NULL). Controls wheel radius, pull velocity,
 *          push depth/velocity, HBP parameters used to compute local drag (via zx_bog_beta_hbp),
 *          and time-stepping.
 *
 * @return zx_bogging_metrics
 *   - drag_N: Sum over the surface footprint of positive differences (wheel_pull_u - local vx),
 *             i.e., an integrated drag proxy in velocity units.
 *   - sink_depth_m: Non-negative sink depth proxy computed as max(0, wheel_push_w * dt * steps).
 */
ZX_API zx_bogging_metrics ZX_CALL zx_integration_wheel_bogging_run(const zx_bogging_params* p)
{
    zx_bogging_metrics M{0.0f, 0.0f};
    if (!p || p->h <= 0.0f || p->dt <= 0.0f || p->steps == 0) return M;
    const uint32_t B = (uint32_t)ZX_TILE_B;
    zx_tile tile{}; auto idx=[&](uint32_t ix,uint32_t iy,uint32_t iz){return ix + B*(iy + B*iz);} ;
    for (uint32_t iz=0; iz<B; ++iz) for (uint32_t iy=0; iy<B; ++iy) for (uint32_t ix=0; ix<B; ++ix){
        auto& n = tile.nodes[idx(ix,iy,iz)]; n.mass = 1.0f; n.mom_x = 0.0f; n.mom_y = 0.0f; n.mom_z = 0.0f;
    }
    const float r_nodes = (float)std::max(1u, p->wheel_radius_nodes);
    const float cx = B/2.0f, cy = B/2.0f; const float ground_z = (float)(B/2);

    for (uint32_t s=0; s<p->steps; ++s){
        // approximate wheel imprint by pushing nodes downwards within radius; resist horizontal drag via Darcy
        for (uint32_t iy=0; iy<B; ++iy) for (uint32_t ix=0; ix<B; ++ix){
            float dx = ix - cx; float dy = iy - cy; float r = std::sqrt(dx*dx + dy*dy);
            if (r <= r_nodes){
                uint32_t iz = (uint32_t)ground_z;
                auto& n = tile.nodes[idx(ix,iy,iz)];
                float vx = p->wheel_pull_u;
                float beta = zx_bog_beta_hbp(p->gamma_dot, &p->hbp, p->permeability, p->k_min,
                                             p->mu_min, p->mu_max, p->beta_min, p->beta_max, p->policy, p->softness_k);
                const float rho = n.mass / (p->h*p->h*p->h);
                float ax = -(beta / std::max(rho,1e-12f)) * (n.mom_x/n.mass - vx);
                n.mom_x += p->dt * ax * n.mass;
            }
        }
        zx_tile pool[1]; pool[0] = tile;
        zx_hbp_update_coarse_grid(pool, 1, p->h, p->dt, &p->hbp, p->mu_min, p->mu_max);
        tile = pool[0];
    }

    // metrics: sink depth proxy (momentum loss in z) and drag force
    float drag = 0.0f; float vx_sum = 0.0f; float count = 0.0f;
    for (uint32_t iy=0; iy<B; ++iy) for (uint32_t ix=0; ix<B; ++ix){
        uint32_t iz = (uint32_t)ground_z;
        auto& n = tile.nodes[idx(ix,iy,iz)];
        float vx = n.mom_x / std::max(n.mass, 1e-12f);
        vx_sum += vx; count += 1.0f;
        drag += std::max(0.0f, p->wheel_pull_u - vx);
    }
    M.drag_N = drag;
    M.sink_depth_m = std::max(0.0f, p->wheel_push_w * p->dt * p->steps);
    return M;
}

/**
 * @brief Computes puddle-creep distance along x by pushing and diffusing momentum on a single coarse tile.
 *
 * Performs timestep iterations applying Darcy-like damping and coarse-grid diffusion (HBP). If the input
 * pointer is NULL or any of p->h, p->dt, or p->steps are non-positive/zero, a zeroed metrics struct is returned.
 * The reported creep distance is the maximum x index whose node horizontal velocity (vx) exceeds 0.05, converted
 * to physical distance by multiplying by p->h.
 *
 * @param p Puddle simulation parameters (must not be NULL; p->h, p->dt, and p->steps must be positive).
 * @return zx_puddle_metrics Metrics with creep_dist_x set to the creep distance along x (same units as p->h).
 */
ZX_API zx_puddle_metrics ZX_CALL zx_integration_puddle_creep_run(const zx_puddle_params* p)
{
    zx_puddle_metrics M{0.0f};
    if (!p || p->h <= 0.0f || p->dt <= 0.0f || p->steps == 0) return M;
    const uint32_t B = (uint32_t)ZX_TILE_B;
    zx_tile tile{}; auto idx=[&](uint32_t ix,uint32_t iy,uint32_t iz){return ix + B*(iy + B*iz);} ;
    for (uint32_t iz=0; iz<B; ++iz) for (uint32_t iy=0; iy<B; ++iy) for (uint32_t ix=0; ix<B; ++ix){
        auto& n = tile.nodes[idx(ix,iy,iz)]; n.mass = 1.0f; n.mom_x = (ix < B/3) ? p->init_head_u : 0.0f; n.mom_y = 0.0f; n.mom_z = 0.0f;
    }
    zx_tile pool[1]; pool[0] = tile;
    for (uint32_t s=0; s<p->steps; ++s){
        // apply Darcy globally to damp
        for (uint32_t k=0; k<B*B*B; ++k){
            auto& n = pool[0].nodes[k];
            float vx = n.mom_x / std::max(n.mass, 1e-12f);
            float beta = zx_bog_beta_hbp(p->gamma_dot, &p->hbp, p->permeability, p->k_min,
                                         p->mu_min, p->mu_max, p->beta_min, p->beta_max, p->policy, p->softness_k);
            float rho = n.mass / (p->h*p->h*p->h);
            float ax = -(beta / std::max(rho,1e-12f)) * vx; n.mom_x += p->dt * ax * n.mass;
        }
        zx_hbp_update_coarse_grid(pool, 1, p->h, p->dt, &p->hbp, p->mu_min, p->mu_max);
    }
    // creep distance: max ix where vx above threshold
    float max_ix=0.0f; for (uint32_t iz=0; iz<B; ++iz) for (uint32_t iy=0; iy<B; ++iy) for (uint32_t ix=0; ix<B; ++ix){
        auto& n = pool[0].nodes[idx(ix,iy,iz)];
        float vx = n.mom_x / std::max(n.mass, 1e-12f);
        if (vx > 0.05f) max_ix = std::max(max_ix, (float)ix);
    }
    M.creep_dist_x = max_ix * p->h;
    return M;
}

} // extern "C"


