/**
 * @file zx_tiles.cpp
 * @brief Allocation and initialization for tiles and particle SoA.
 * @details C-ABI allocation helpers for tile pools and particle SoA using host-provided
 *          callbacks. Functions defensively check pointers and sizes; no global state.
 *          Thread-safe when allocators are thread-safe and buffers are disjoint.
 * @copyright
 *   (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
 *   Licensed for use within the ziXn project under project terms.
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

  /**
   * @brief Allocate a particle SoA of length count using a user allocator and zero-initialize.
   * @param count Number of particles.
   * @param alloc_fn User allocator callback: alloc_fn(num_bytes, user) -> pointer.
   * @param user Opaque pointer forwarded to the allocator.
   * @return zx_particle_soa SoA with all pointers set (or zeroed if allocation fails/invalid).
   */
  zx_particle_soa ZX_CALL zx_create_particle_soa(size_t count, void* (*alloc_fn)(size_t, void*),
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

  /**
   * @brief Destroy a particle SoA using a user-provided free function and reset fields.
   * @param soa SoA to destroy; may be null.
   * @param free_fn User free callback: free_fn(ptr, user).
   * @param user Opaque pointer forwarded to the free callback.
   */
  void ZX_CALL zx_destroy_particle_soa(zx_particle_soa* soa, void (*free_fn)(void*, void*),
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

  /**
   * @brief Allocate a tile pool of given capacity using a user allocator and zero-initialize.
   * @param capacity Number of tiles to allocate.
   * @param alloc_fn User allocator callback.
   * @param user Opaque pointer forwarded to the allocator.
   * @return zx_tile* Pointer to zeroed tile array, or nullptr on invalid input.
   */
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

  /**
   * @brief Free a tile pool via user callback (no-op on null inputs).
   * @param pool Pointer to tile pool array.
   * @param free_fn User free callback.
   * @param user Opaque pointer forwarded to the free callback.
   */
  void ZX_CALL zx_destroy_tile_pool(zx_tile* pool, void (*free_fn)(void*, void*), void* user)
  {
    if (pool == nullptr || free_fn == nullptr)
    {
      return;
    }
    free_fn(pool, user);
  }
}
