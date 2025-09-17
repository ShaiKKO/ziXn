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
  if (p && n)
  {
    std::memset(p, 0, n);
  }
}

extern "C"
{

  /* Allocation helpers are exposed as C symbols to ease testing and early wiring. */

  zx_particle_soa ZX_CALL zxCreateParticleSoA(size_t count, void* (*alloc_fn)(size_t, void*),
                                              void* user)
  {
    zx_particle_soa soa{};
    if (!alloc_fn || count == 0)
      return soa;
    const size_t n = count;
    soa.pos_x      = (float*) alloc_fn(sizeof(float) * n, user);
    zx_memzero(soa.pos_x, sizeof(float) * n);
    soa.pos_y = (float*) alloc_fn(sizeof(float) * n, user);
    zx_memzero(soa.pos_y, sizeof(float) * n);
    soa.pos_z = (float*) alloc_fn(sizeof(float) * n, user);
    zx_memzero(soa.pos_z, sizeof(float) * n);
    soa.vel_x = (float*) alloc_fn(sizeof(float) * n, user);
    zx_memzero(soa.vel_x, sizeof(float) * n);
    soa.vel_y = (float*) alloc_fn(sizeof(float) * n, user);
    zx_memzero(soa.vel_y, sizeof(float) * n);
    soa.vel_z = (float*) alloc_fn(sizeof(float) * n, user);
    zx_memzero(soa.vel_z, sizeof(float) * n);
    soa.mass = (float*) alloc_fn(sizeof(float) * n, user);
    zx_memzero(soa.mass, sizeof(float) * n);
    soa.volume = (float*) alloc_fn(sizeof(float) * n, user);
    zx_memzero(soa.volume, sizeof(float) * n);
    soa.F = (float*) alloc_fn(sizeof(float) * n * 9, user);
    zx_memzero(soa.F, sizeof(float) * n * 9);
    soa.C = (float*) alloc_fn(sizeof(float) * n * 9, user);
    zx_memzero(soa.C, sizeof(float) * n * 9);
    soa.mat_id = (uint16_t*) alloc_fn(sizeof(uint16_t) * n, user);
    zx_memzero(soa.mat_id, sizeof(uint16_t) * n);
    soa.flags = (uint16_t*) alloc_fn(sizeof(uint16_t) * n, user);
    zx_memzero(soa.flags, sizeof(uint16_t) * n);
    return soa;
  }

  void ZX_CALL zxDestroyParticleSoA(zx_particle_soa* soa, void (*free_fn)(void*, void*), void* user)
  {
    if (!soa || !free_fn)
      return;
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

  zx_tile* ZX_CALL zxCreateTilePool(uint32_t capacity, void* (*alloc_fn)(size_t, void*), void* user)
  {
    if (!alloc_fn || capacity == 0)
      return nullptr;
    zx_tile* pool = (zx_tile*) alloc_fn(sizeof(zx_tile) * capacity, user);
    zx_memzero(pool, sizeof(zx_tile) * capacity);
    return pool;
  }

  void ZX_CALL zxDestroyTilePool(zx_tile* pool, void (*free_fn)(void*, void*), void* user)
  {
    if (!pool || !free_fn)
      return;
    free_fn(pool, user);
  }
}
