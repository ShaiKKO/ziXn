/**
 * @file zx_integration.cpp
 * @brief Dam-break over porous bed (coarse proxy) metrics implementation.
 * @details Implements coarse proxy integrations for dam-break, wheel bogging, and puddle creep
 *          scenarios using HBP updates and simple drag models. Functions return minimal metrics
 *          used by validation tests. Thread-safe given independent tiles/contexts.
 * @copyright
 *   (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
 *   Licensed for use within the ziXn project under project terms.
 */

#include "zx/zx_integration.h"
#include "zx/zx_lod.h"
#include "zx/zx_residency.h"
#include <algorithm>
#include <cmath>

extern "C"
{

  /**
   * @brief Simulate a single-tile coarse-grid dam-break and compute front advance and kinetic
   * energy.
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
  zx_dambreak_metrics ZX_CALL zx_integration_dambreak_run(const zx_dambreak_params* p)
  {
    zx_dambreak_metrics m{0.0F, 0.0F};
    if ((p == nullptr) || (p->tiles == 0U) || (p->h <= 0.0F) || (p->dt <= 0.0F) || (p->steps == 0U))
    {
      return m;
    }
    constexpr float k_vx_thresh = 0.05F;
    constexpr float k_rho_eps   = 1.0e-12F;
    // Allocate a single-tile pool and simple initial condition: left half moving, right half at
    // rest
    zx_tile tile{};
    tile.coord_x = 0;
    tile.coord_y = 0;
    tile.coord_z = 0;
    const auto b = static_cast<uint32_t>(ZX_TILE_B);
    auto idx = [&](uint32_t ix, uint32_t iy, uint32_t iz) { return ix + (b * (iy + (b * iz))); };
    for (uint32_t iz = 0; iz < b; ++iz)
    {
      for (uint32_t iy = 0; iy < b; ++iy)
      {
        for (uint32_t ix = 0; ix < b; ++ix)
        {
          auto& n         = tile.nodes[idx(ix, iy, iz)];
          n.mass          = 1.0F;
          const bool left = ix < (b / 2U);
          const bool bed  = iz < std::min<uint32_t>(p->bed_k, b);
          const float vx0 = left ? p->init_head_u : 0.0F;
          n.mom_x         = vx0 * n.mass;
          n.mom_y         = 0.0F;
          n.mom_z         = 0.0F;
          // Apply immediate Darcy drag in bed region to mimic porous damping
          if (bed)
          {
            float beta =
                zx_bog_beta_hbp(p->gamma_dot, &p->hbp, p->permeability, p->k_min, p->mu_min,
                                p->mu_max, p->beta_min, p->beta_max, p->policy, p->softness_k);
            const float rho = n.mass / (p->h * p->h * p->h);
            const float ax  = -(beta / std::max(rho, k_rho_eps)) * vx0;  // dv/dt = -(beta/rho)*v
            const float vx1 = vx0 + (p->dt * ax);
            n.mom_x         = std::max(0.0F, vx1) * n.mass;
          }
        }
      }
    }

    // Diffuse via HBP for 'steps' with optional LOD fallback pressure loop
    zx_tile pool[1];
    pool[0] = tile;
    zx_residency_opts ro{1, 2, 1};
    zx_residency* rez    = zx_residency_create(&ro);
    uint32_t churn_enter = 0;
    uint32_t churn_exit  = 0;
    uint32_t pf_count    = 0;
    zx_lod_fallback_state fbs;
    zx_lod_fallback_init(&fbs);
    zx_lod_fallback_policy fbp;
    zx_lod_get_default_policy(&fbp);
    for (uint32_t s = 0; s < p->steps; ++s)
    {
      // If fallback enabled, track active frames via residency and policy
      if (zx_lod_is_enabled() != 0)
      {
        uint32_t en = 0;
        uint32_t ex = 0;
        uint32_t pf = 0;
        zx_residency_tick(rez, static_cast<int>(s), 0, 0, 0, &en, &ex, &pf);
        churn_enter += en;
        churn_exit += ex;
        pf_count        = pf;
        uint32_t active = zx_residency_get_active_count(rez);
        (void) zx_lod_fallback_update(&fbp, active, /*last_step_ms*/ 1.0F, &fbs);
      }
      else
      {
        uint32_t en = 0;
        uint32_t ex = 0;
        uint32_t pf = 0;
        zx_residency_tick(rez, static_cast<int>(s), 0, 0, 0, &en, &ex, &pf);
        churn_enter += en;
        churn_exit += ex;
        pf_count = pf;
      }
      zx_hbp_update_coarse_grid(pool, 1, p->h, p->dt, &p->hbp, p->mu_min, p->mu_max);
    }

    // Metrics: front advance via threshold on vx, and kinetic energy sum
    float max_ix = 0.0F;
    float ke     = 0.0F;
    for (uint32_t iz = 0; iz < b; ++iz)
    {
      for (uint32_t iy = 0; iy < b; ++iy)
      {
        for (uint32_t ix = 0; ix < b; ++ix)
        {
          const auto& n = pool[0].nodes[idx(ix, iy, iz)];
          const float m = n.mass;
          if (m <= 0.0F)
          {
            continue;
          }
          const float vx = n.mom_x / m;
          const float vy = n.mom_y / m;
          const float vz = n.mom_z / m;
          const float v2 = (vx * vx) + (vy * vy) + (vz * vz);
          if (vx > k_vx_thresh)
          {
            max_ix = std::max(max_ix, static_cast<float>(ix));
          }
          ke += 0.5F * m * v2;
        }
      }
    }
    m.front_x   = max_ix * p->h;
    m.kinetic_j = ke;
    // Could export churn/fallback counters via telemetry externally
    (void) churn_enter;
    (void) churn_exit;
    (void) pf_count;
    zx_residency_destroy(rez);
    return m;
  }

  /**
   * @brief Simulates wheel-induced bogging on a single coarse-grid tile and returns drag and sink
   * metrics.
   *
   * Runs a proxy simulation where a circular wheel footprint applies a horizontal pull velocity
   * and a vertical push onto a coarse BxBxB tile. Local Darcy-like drag is applied per-node,
   * the coarse-grid viscous/diffusive update is advanced for the specified number of steps,
   * and final surface-layer states are reduced into simple metrics.
   *
   * The function performs input validation and returns zeroed metrics if `p` is null or if
   * `p->h`, `p->dt`, or `p->steps` are non-positive.
   *
   * @param p Pointer to bogging parameters (must not be NULL). Controls wheel radius, pull
   * velocity, push depth/velocity, HBP parameters used to compute local drag (via zx_bog_beta_hbp),
   *          and time-stepping.
   *
   * @return zx_bogging_metrics
   *   - drag_N: Sum over the surface footprint of positive differences (wheel_pull_u - local vx),
   *             i.e., an integrated drag proxy in velocity units.
   *   - sink_depth_m: Non-negative sink depth proxy computed as max(0, wheel_push_w * dt * steps).
   */
  zx_bogging_metrics ZX_CALL zx_integration_wheel_bogging_run(const zx_bogging_params* p)
  {
    zx_bogging_metrics m{0.0F, 0.0F};
    if ((p == nullptr) || (p->h <= 0.0F) || (p->dt <= 0.0F) || (p->steps == 0U))
    {
      return m;
    }
    const auto b = static_cast<uint32_t>(ZX_TILE_B);
    zx_tile tile{};
    auto idx = [&](uint32_t ix, uint32_t iy, uint32_t iz) { return ix + (b * (iy + (b * iz))); };
    for (uint32_t iz = 0; iz < b; ++iz)
    {
      for (uint32_t iy = 0; iy < b; ++iy)
      {
        for (uint32_t ix = 0; ix < b; ++ix)
        {
          auto& n = tile.nodes[idx(ix, iy, iz)];
          n.mass  = 1.0F;
          n.mom_x = 0.0F;
          n.mom_y = 0.0F;
          n.mom_z = 0.0F;
        }
      }
    }
    const float r_nodes       = static_cast<float>(std::max(1U, p->wheel_radius_nodes));
    const float cx            = static_cast<float>(b) / 2.0F;
    const float cy            = static_cast<float>(b) / 2.0F;
    const auto ground_z       = static_cast<float>(b) / 2.0F;
    constexpr float k_rho_eps = 1.0e-12F;

    for (uint32_t s = 0; s < p->steps; ++s)
    {
      // approximate wheel imprint by pushing nodes downwards within radius; resist horizontal drag
      // via Darcy
      for (uint32_t iy = 0; iy < b; ++iy)
      {
        for (uint32_t ix = 0; ix < b; ++ix)
        {
          const float dx = static_cast<float>(ix) - cx;
          const float dy = static_cast<float>(iy) - cy;
          const float r  = std::sqrt((dx * dx) + (dy * dy));
          if (r <= r_nodes)
          {
            auto iz  = static_cast<uint32_t>(ground_z);
            auto& n  = tile.nodes[idx(ix, iy, iz)];
            float vx = p->wheel_pull_u;
            float beta =
                zx_bog_beta_hbp(p->gamma_dot, &p->hbp, p->permeability, p->k_min, p->mu_min,
                                p->mu_max, p->beta_min, p->beta_max, p->policy, p->softness_k);
            const float rho = n.mass / (p->h * p->h * p->h);
            float ax        = -(beta / std::max(rho, k_rho_eps)) * (n.mom_x / n.mass - vx);
            n.mom_x += p->dt * ax * n.mass;
          }
        }
      }
      zx_tile pool[1];
      pool[0] = tile;
      zx_hbp_update_coarse_grid(pool, 1, p->h, p->dt, &p->hbp, p->mu_min, p->mu_max);
      tile = pool[0];
    }

    // metrics: sink depth proxy (momentum loss in z) and drag force
    float drag = 0.0F;
    for (uint32_t iy = 0; iy < b; ++iy)
    {
      for (uint32_t ix = 0; ix < b; ++ix)
      {
        auto iz  = static_cast<uint32_t>(ground_z);
        auto& n  = tile.nodes[idx(ix, iy, iz)];
        float vx = n.mom_x / std::max(n.mass, k_rho_eps);
        drag += std::max(0.0F, p->wheel_pull_u - vx);
      }
    }
    m.drag_N       = drag;
    m.sink_depth_m = std::max(0.0F, p->wheel_push_w * p->dt * static_cast<float>(p->steps));
    return m;
  }

  /**
   * @brief Computes puddle-creep distance along x by pushing and diffusing momentum on a single
   * coarse tile.
   *
   * Performs timestep iterations applying Darcy-like damping and coarse-grid diffusion (HBP). If
   * the input pointer is NULL or any of p->h, p->dt, or p->steps are non-positive/zero, a zeroed
   * metrics struct is returned. The reported creep distance is the maximum x index whose node
   * horizontal velocity (vx) exceeds 0.05, converted to physical distance by multiplying by p->h.
   *
   * @param p Puddle simulation parameters (must not be NULL; p->h, p->dt, and p->steps must be
   * positive).
   * @return zx_puddle_metrics Metrics with creep_dist_x set to the creep distance along x (same
   * units as p->h).
   */
  zx_puddle_metrics ZX_CALL zx_integration_puddle_creep_run(const zx_puddle_params* p)
  {
    zx_puddle_metrics m{0.0F};
    if ((p == nullptr) || (p->h <= 0.0F) || (p->dt <= 0.0F) || (p->steps == 0U))
    {
      return m;
    }
    constexpr float k_vx_thresh = 0.05F;
    constexpr float k_rho_eps   = 1.0e-12F;
    const auto b                = static_cast<uint32_t>(ZX_TILE_B);
    zx_tile tile{};
    auto idx = [&](uint32_t ix, uint32_t iy, uint32_t iz) { return ix + (b * (iy + (b * iz))); };
    for (uint32_t iz = 0; iz < b; ++iz)
    {
      for (uint32_t iy = 0; iy < b; ++iy)
      {
        for (uint32_t ix = 0; ix < b; ++ix)
        {
          auto& n = tile.nodes[idx(ix, iy, iz)];
          n.mass  = 1.0F;
          n.mom_x = (ix < (b / 3U)) ? p->init_head_u : 0.0F;
          n.mom_y = 0.0F;
          n.mom_z = 0.0F;
        }
      }
    }
    zx_tile pool[1];
    pool[0] = tile;
    for (uint32_t s = 0; s < p->steps; ++s)
    {
      // apply Darcy globally to damp
      for (auto& n : pool[0].nodes)
      {
        float vx   = n.mom_x / std::max(n.mass, k_rho_eps);
        float beta = zx_bog_beta_hbp(p->gamma_dot, &p->hbp, p->permeability, p->k_min, p->mu_min,
                                     p->mu_max, p->beta_min, p->beta_max, p->policy, p->softness_k);
        float rho  = n.mass / (p->h * p->h * p->h);
        float ax   = -(beta / std::max(rho, k_rho_eps)) * vx;
        n.mom_x += p->dt * ax * n.mass;
      }
      zx_hbp_update_coarse_grid(pool, 1, p->h, p->dt, &p->hbp, p->mu_min, p->mu_max);
    }
    // creep distance: max ix where vx above threshold
    float max_ix = 0.0F;
    for (uint32_t iz = 0; iz < b; ++iz)
    {
      for (uint32_t iy = 0; iy < b; ++iy)
      {
        for (uint32_t ix = 0; ix < b; ++ix)
        {
          auto& n  = pool[0].nodes[idx(ix, iy, iz)];
          float vx = n.mom_x / std::max(n.mass, k_rho_eps);
          if (vx > k_vx_thresh)
          {
            max_ix = std::max(max_ix, static_cast<float>(ix));
          }
        }
      }
    }
    m.creep_dist_x = max_ix * p->h;
    return m;
  }

}  // extern "C"
