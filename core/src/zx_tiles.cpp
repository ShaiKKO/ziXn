/*!
 * \file zx_tiles.cpp
 * \brief Allocation and initialization for tiles and particle SoA.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_tiles.h"
#include "zx/zx_tiles_api.h"
#include <cstring>

static void zx_memzero(void* p, size_t n)
{
  if ((p != nullptr) && (n != 0U))
  {
    std::memset(p, 0, n);
  }
}

extern "C"
{

  /* Allocation helpers are exposed as C symbols to ease testing and early wiring. */

  zx_particle_soa ZX_CALL zx_create_particle_so_a(size_t count, void* (*alloc_fn)(size_t, void*),
                                                  void* user)
  {
    zx_particle_soa soa{};
    if (alloc_fn == nullptr || count == 0)
    {
      return soa;
    }
    const size_t n = count;
    soa.pos_x      = static_cast<float*>(alloc_fn(sizeof(float) * n, user));
    zx_memzero(soa.pos_x, sizeof(float) * n);
    soa.pos_y = static_cast<float*>(alloc_fn(sizeof(float) * n, user));
    zx_memzero(soa.pos_y, sizeof(float) * n);
    soa.pos_z = static_cast<float*>(alloc_fn(sizeof(float) * n, user));
    zx_memzero(soa.pos_z, sizeof(float) * n);
    soa.vel_x = static_cast<float*>(alloc_fn(sizeof(float) * n, user));
    zx_memzero(soa.vel_x, sizeof(float) * n);
    soa.vel_y = static_cast<float*>(alloc_fn(sizeof(float) * n, user));
    zx_memzero(soa.vel_y, sizeof(float) * n);
    soa.vel_z = static_cast<float*>(alloc_fn(sizeof(float) * n, user));
    zx_memzero(soa.vel_z, sizeof(float) * n);
    soa.mass = static_cast<float*>(alloc_fn(sizeof(float) * n, user));
    zx_memzero(soa.mass, sizeof(float) * n);
    soa.volume = static_cast<float*>(alloc_fn(sizeof(float) * n, user));
    zx_memzero(soa.volume, sizeof(float) * n);
    constexpr auto k_mat3 = static_cast<size_t>(zx_mat3_size);
    soa.F                 = static_cast<float*>(alloc_fn(sizeof(float) * n * k_mat3, user));
    zx_memzero(soa.F, sizeof(float) * n * k_mat3);
    soa.C = static_cast<float*>(alloc_fn(sizeof(float) * n * k_mat3, user));
    zx_memzero(soa.C, sizeof(float) * n * k_mat3);
    soa.mat_id = static_cast<uint16_t*>(alloc_fn(sizeof(uint16_t) * n, user));
    zx_memzero(soa.mat_id, sizeof(uint16_t) * n);
    soa.flags = static_cast<uint16_t*>(alloc_fn(sizeof(uint16_t) * n, user));
    zx_memzero(soa.flags, sizeof(uint16_t) * n);
    return soa;
  }

  void ZX_CALL zx_destroy_particle_so_a(zx_particle_soa* soa, void (*free_fn)(void*, void*),
                                        void* user)
  {
    if (soa == nullptr || free_fn == nullptr)
    {
      return;
    }
    free_fn(soa->pos_x, user);
    free_fn(soa->pos_y, user);
    free_fn(soa->pos_z, user);
    free_fn(soa->vel_x, user);
    free_fn(soa->vel_y, user);
    free_fn(soa->vel_z, user);
    free_fn(soa->mass, user);
    free_fn(soa->volume, user);
    free_fn(soa->F, user);
    free_fn(soa->C, user);
    free_fn(soa->mat_id, user);
    free_fn(soa->flags, user);
    *soa = zx_particle_soa{};
  }

  zx_tile* ZX_CALL zx_create_tile_pool(uint32_t capacity, void* (*alloc_fn)(size_t, void*),
                                       void* user)
  {
    if (alloc_fn == nullptr || capacity == 0)
    {
      return nullptr;
    }
    auto* pool = static_cast<zx_tile*>(alloc_fn(sizeof(zx_tile) * capacity, user));
    zx_memzero(pool, sizeof(zx_tile) * capacity);
    return pool;
  }

  void ZX_CALL zx_destroy_tile_pool(zx_tile* pool, void (*free_fn)(void*, void*), void* user)
  {
    if (pool == nullptr || free_fn == nullptr)
    {
      return;
    }
    free_fn(pool, user);
  }
}
