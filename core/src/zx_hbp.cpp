/**
 * @file zx_hbp.cpp
 * @brief Herschel–Bulkley–Papanastasiou coarse-grid viscous update (CPU reference).
 * @details Implements coarse-grid viscosity update using HBP rheology for integration proxies
 *          (dam-break, bogging, puddle). Stateless wrt global state; functions operate on
 *          provided tile pools. Thread-safety depends on non-overlapping tile buffers.
 * @copyright
 *   (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
 *   Licensed for use within the ziXn project under project terms.
 */

#include "zx/zx_hbp.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace
{

  struct TileKey
  {
    int x, y, z;
    bool operator==(const TileKey& o) const noexcept
    {
      return x == o.x && y == o.y && z == o.z;
    }
  };

  struct TileKeyHash
  {
    size_t operator()(const TileKey& k) const noexcept
    {
      auto mix = [](uint64_t h, uint64_t v)
      {
        h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
      };
      uint64_t h = 1469598103934665603ULL;
      h          = mix(h, static_cast<uint64_t>(static_cast<uint32_t>(k.x)));
      h          = mix(h, static_cast<uint64_t>(static_cast<uint32_t>(k.y)));
      h          = mix(h, static_cast<uint64_t>(static_cast<uint32_t>(k.z)));
      return static_cast<size_t>(h);
    }
  };

  inline uint32_t node_index(uint32_t ix, uint32_t iy, uint32_t iz)
  {
    const auto b = static_cast<uint32_t>(ZX_TILE_B);
    return ix + (b * (iy + (b * iz)));
  }

  struct NeighborRef
  {
    int tile_idx;
    uint32_t local_idx;
  };

  NeighborRef get_neighbor_ref(const zx_tile* pool,
                               const std::unordered_map<TileKey, int, TileKeyHash>& key_to_index,
                               int tile_idx, int ix, int iy, int iz, int dx, int dy, int dz)
  {
    const zx_tile& t = pool[tile_idx];
    int ntx          = t.coord_x;
    int nty          = t.coord_y;
    int ntz          = t.coord_z;

    int nix = ix + dx;
    int niy = iy + dy;
    int niz = iz + dz;

    int carry_x = 0;
    int carry_y = 0;
    int carry_z = 0;
    if (nix < 0)
    {
      nix += ZX_TILE_B;
      carry_x = -1;
    }
    else if (nix >= ZX_TILE_B)
    {
      nix -= ZX_TILE_B;
      carry_x = 1;
    }
    if (niy < 0)
    {
      niy += ZX_TILE_B;
      carry_y = -1;
    }
    else if (niy >= ZX_TILE_B)
    {
      niy -= ZX_TILE_B;
      carry_y = 1;
    }
    if (niz < 0)
    {
      niz += ZX_TILE_B;
      carry_z = -1;
    }
    else if (niz >= ZX_TILE_B)
    {
      niz -= ZX_TILE_B;
      carry_z = 1;
    }

    TileKey nk{ntx + carry_x, nty + carry_y, ntz + carry_z};
    auto it = key_to_index.find(nk);
    if (it == key_to_index.end())
    {
      return NeighborRef{-1, 0};
    }
    return NeighborRef{it->second,
                       node_index(static_cast<uint32_t>(nix), static_cast<uint32_t>(niy),
                                  static_cast<uint32_t>(niz))};
  }

}  // namespace

extern "C"
{

  float ZX_CALL zx_hbp_mu_eff(float gamma_dot, const zx_hbp_params* params, float mu_min,
                              float mu_max)
  {
    if (params == nullptr)
    {
      return std::clamp(0.0F, mu_min, mu_max);
    }
    const float g      = std::max(gamma_dot, 0.0F);
    const float eps    = 1.0e-8F;
    const float g_safe = std::max(g, eps);

    float f_y = 0.0F;
    if ((params->tau_y > 0.0F) && (params->m > 0.0F))
    {
      if (g > eps)
      {
        f_y = (1.0F - std::exp(-params->m * g)) / g;
      }
      else
      {
        f_y = params->m;
      }
    }

    float power_term = 0.0F;
    if (params->K > 0.0F)
    {
      power_term = params->K * std::pow(g_safe, params->n - 1.0F);
    }

    float mu = params->mu0 + power_term + (params->tau_y * f_y);
    if (!(mu_min <= mu_max))
    {
      std::swap(mu_min, mu_max);
    }
    mu = std::min(std::max(mu, mu_min), mu_max);
    return mu;
  }

  float ZX_CALL zx_hbp_mu_eff_policy(float gamma_dot, const zx_hbp_params* params, float mu_min,
                                     float mu_max, zx_mu_clamp_policy policy, float softness_k)
  {
    float mu = zx_hbp_mu_eff(gamma_dot, params, -INFINITY, INFINITY);
    if (!(mu_min <= mu_max))
    {
      std::swap(mu_min, mu_max);
    }
    switch (policy)
    {
    case ZX_MU_CLAMP_NONE:
      break;
    case ZX_MU_CLAMP_HARD:
      mu = std::min(std::max(mu, mu_min), mu_max);
      break;
    case ZX_MU_CLAMP_SMOOTH_TANH:
    {
      const float k = (softness_k > 0.0F) ? softness_k : 1.0F;
      // Smoothly push into range using tanh-based soft clipping near bounds
      // Map to [-1,1] via affine, then invert mapping
      const float mid  = 0.5F * (mu_min + mu_max);
      const float half = 0.5F * (mu_max - mu_min);
      float t          = (mu - mid) / std::max(half, 1.0e-12F);
      t                = std::tanh(k * t);
      mu               = mid + half * t;
      // tiny extra hard clamp to ensure numeric safety
      mu = std::min(std::max(mu, mu_min), mu_max);
      break;
    }
    default:
      mu = std::min(std::max(mu, mu_min), mu_max);
      break;
    }
    return mu;
  }

  int ZX_CALL zx_hbp_validate_params(zx_hbp_params* inout)
  {
    if (inout == nullptr)
    {
      return -1;
    }
    int adjusted = 0;
    if (!(inout->mu0 >= 0.0F))
    {
      inout->mu0 = std::max(inout->mu0, 0.0F);
      adjusted   = 1;
    }
    if (!(inout->K >= 0.0F))
    {
      inout->K = std::max(inout->K, 0.0F);
      adjusted = 1;
    }
    if (!(inout->n > 0.0F))
    {
      inout->n = std::max(inout->n, 1.0e-6F);
      adjusted = 1;
    }
    if (!(inout->tau_y >= 0.0F))
    {
      inout->tau_y = std::max(inout->tau_y, 0.0F);
      adjusted     = 1;
    }
    if (!(inout->m >= 0.0F))
    {
      inout->m = std::max(inout->m, 0.0F);
      adjusted = 1;
    }
    return adjusted;
  }

  float ZX_CALL zx_hbp_dt_stable_upper_bound(float h, float mu_max, float rho_min)
  {
    if ((h <= 0.0F) || (mu_max <= 0.0F) || (rho_min <= 0.0F))
    {
      return 0.0F;
    }
    return (h * h * rho_min) / (6.0F * mu_max);
  }

  void ZX_CALL zx_hbp_update_coarse_grid(zx_tile* pool, uint32_t tile_count, float h, float dt,
                                         const zx_hbp_params* params, float mu_min, float mu_max)
  {
    if ((pool == nullptr) || (tile_count == 0U) || (h <= 0.0F) || (dt <= 0.0F) ||
        (params == nullptr))
    {
      return;
    }

    const auto b                  = static_cast<uint32_t>(ZX_TILE_B);
    const uint32_t nodes_per_tile = b * b * b;
    const uint64_t total_nodes =
        static_cast<uint64_t>(tile_count) * static_cast<uint64_t>(nodes_per_tile);

    std::unordered_map<TileKey, int, TileKeyHash> key_to_index;
    key_to_index.reserve(static_cast<size_t>(tile_count) * 2ULL);
    for (uint32_t t = 0; t < tile_count; ++t)
    {
      key_to_index.insert(
          {TileKey{pool[t].coord_x, pool[t].coord_y, pool[t].coord_z}, static_cast<int>(t)});
    }

    std::vector<float> vx(total_nodes, 0.0F);
    std::vector<float> vy(total_nodes, 0.0F);
    std::vector<float> vz(total_nodes, 0.0F);
    std::vector<float> mass(total_nodes, 0.0F);
    for (uint32_t t = 0; t < tile_count; ++t)
    {
      const zx_tile& tile = pool[t];
      const uint64_t base = static_cast<uint64_t>(t) * nodes_per_tile;
      for (uint32_t k = 0; k < nodes_per_tile; ++k)
      {
        const zx_tile_node& n = tile.nodes[k];
        const float m         = n.mass;
        mass[base + k]        = m;
        if (m > 0.0F)
        {
          vx[base + k] = n.mom_x / m;
          vy[base + k] = n.mom_y / m;
          vz[base + k] = n.mom_z / m;
        }
      }
    }

    std::vector<float> mu(total_nodes, mu_min);
    const float inv2h = 0.5F / h;
    for (uint32_t t = 0; t < tile_count; ++t)
    {
      const uint64_t base = static_cast<uint64_t>(t) * nodes_per_tile;
      for (uint32_t iz = 0; iz < b; ++iz)
      {
        for (uint32_t iy = 0; iy < b; ++iy)
        {
          for (uint32_t ix = 0; ix < b; ++ix)
          {
            const uint32_t k   = node_index(ix, iy, iz);
            const float m_node = mass[base + k];
            if (m_node <= 0.0F)
            {
              mu[base + k] = mu_min;
              continue;
            }

            auto sample_v = [&](int dx, int dy, int dz)
            {
              NeighborRef nr =
                  get_neighbor_ref(pool, key_to_index, static_cast<int>(t), static_cast<int>(ix),
                                   static_cast<int>(iy), static_cast<int>(iz), dx, dy, dz);
              if (nr.tile_idx < 0)
              {
                return std::array<float, 3>{vx[base + k], vy[base + k], vz[base + k]};
              }
              const uint64_t nbase = static_cast<uint64_t>(nr.tile_idx) * nodes_per_tile;
              return std::array<float, 3>{vx[nbase + nr.local_idx], vy[nbase + nr.local_idx],
                                          vz[nbase + nr.local_idx]};
            };

            auto vx_p = sample_v(+1, 0, 0);
            auto vx_m = sample_v(-1, 0, 0);
            auto vy_p = sample_v(0, +1, 0);
            auto vy_m = sample_v(0, -1, 0);
            auto vz_p = sample_v(0, 0, +1);
            auto vz_m = sample_v(0, 0, -1);

            const float dvx_dx = (vx_p[0] - vx_m[0]) * inv2h;
            const float dvx_dy = (vy_p[0] - vy_m[0]) * inv2h;
            const float dvx_dz = (vz_p[0] - vz_m[0]) * inv2h;

            const float dvy_dx = (vx_p[1] - vx_m[1]) * inv2h;
            const float dvy_dy = (vy_p[1] - vy_m[1]) * inv2h;
            const float dvy_dz = (vz_p[1] - vz_m[1]) * inv2h;

            const float dvz_dx = (vx_p[2] - vx_m[2]) * inv2h;
            const float dvz_dy = (vy_p[2] - vy_m[2]) * inv2h;
            const float dvz_dz = (vz_p[2] - vz_m[2]) * inv2h;

            const float sxx = dvx_dx;
            const float syy = dvy_dy;
            const float szz = dvz_dz;
            const float sxy = 0.5F * (dvx_dy + dvy_dx);
            const float sxz = 0.5F * (dvx_dz + dvz_dx);
            const float syz = 0.5F * (dvy_dz + dvz_dy);

            const float frob2 = (sxx * sxx) + (syy * syy) + (szz * szz) +
                                (2.0F * ((sxy * sxy) + (sxz * sxz) + (syz * syz)));
            const float gamma_dot = std::sqrt(std::max(0.0F, 2.0F * frob2));

            mu[base + k] = zx_hbp_mu_eff(gamma_dot, params, mu_min, mu_max);
          }
        }
      }
    }

    const float inv_h2 = 1.0F / (h * h);
    std::vector<float> n_vx(total_nodes, 0.0F);
    std::vector<float> n_vy(total_nodes, 0.0F);
    std::vector<float> n_vz(total_nodes, 0.0F);
    for (uint32_t t = 0; t < tile_count; ++t)
    {
      const uint64_t base = static_cast<uint64_t>(t) * nodes_per_tile;
      for (uint32_t iz = 0; iz < b; ++iz)
      {
        for (uint32_t iy = 0; iy < b; ++iy)
        {
          for (uint32_t ix = 0; ix < b; ++ix)
          {
            const uint32_t k   = node_index(ix, iy, iz);
            const float m_node = mass[base + k];
            if (m_node <= 0.0F)
            {
              n_vx[base + k] = 0.0F;
              n_vy[base + k] = 0.0F;
              n_vz[base + k] = 0.0F;
              continue;
            }

            const float vxc = vx[base + k];
            const float vyc = vy[base + k];
            const float vzc = vz[base + k];
            const float muc = mu[base + k];

            auto accum_axis = [&](int dx, int dy, int dz, float& ax, float& ay, float& az)
            {
              NeighborRef nr =
                  get_neighbor_ref(pool, key_to_index, static_cast<int>(t), static_cast<int>(ix),
                                   static_cast<int>(iy), static_cast<int>(iz), dx, dy, dz);
              float vxn = vxc;
              float vyn = vyc;
              float vzn = vzc;
              float muf = muc;
              if (nr.tile_idx >= 0)
              {
                const uint64_t nbase = static_cast<uint64_t>(nr.tile_idx) * nodes_per_tile;
                vxn                  = vx[nbase + nr.local_idx];
                vyn                  = vy[nbase + nr.local_idx];
                vzn                  = vz[nbase + nr.local_idx];
                muf                  = 0.5F * (muc + mu[nbase + nr.local_idx]);
              }
              ax += muf * (vxn - vxc);
              ay += muf * (vyn - vyc);
              az += muf * (vzn - vzc);
            };

            float ax = 0.0F;
            float ay = 0.0F;
            float az = 0.0F;
            accum_axis(+1, 0, 0, ax, ay, az);
            accum_axis(-1, 0, 0, ax, ay, az);
            accum_axis(0, +1, 0, ax, ay, az);
            accum_axis(0, -1, 0, ax, ay, az);
            accum_axis(0, 0, +1, ax, ay, az);
            accum_axis(0, 0, -1, ax, ay, az);

            const float rho   = m_node / (h * h * h);
            const float scale = (dt * inv_h2) / std::max(rho, 1.0e-12F);
            n_vx[base + k]    = vxc + scale * ax;
            n_vy[base + k]    = vyc + scale * ay;
            n_vz[base + k]    = vzc + scale * az;
          }
        }
      }
    }

    for (uint32_t t = 0; t < tile_count; ++t)
    {
      zx_tile& tile       = pool[t];
      const uint64_t base = static_cast<uint64_t>(t) * nodes_per_tile;
      for (uint32_t k = 0; k < nodes_per_tile; ++k)
      {
        const float m = mass[base + k];
        if (m > 0.0F)
        {
          tile.nodes[k].mom_x = n_vx[base + k] * m;
          tile.nodes[k].mom_y = n_vy[base + k] * m;
          tile.nodes[k].mom_z = n_vz[base + k] * m;
        }
        else
        {
          tile.nodes[k].mom_x = 0.0F;
          tile.nodes[k].mom_y = 0.0F;
          tile.nodes[k].mom_z = 0.0F;
        }
      }
    }
  }

}  // extern "C"
