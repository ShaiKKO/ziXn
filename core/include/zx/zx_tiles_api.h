/*!
\file zx_tiles_api.h
\brief C symbols for tiles/SoA allocation to support early wiring/tests.
\author Colin Macritchie (Ripple Group, LLC)
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.

Windows-first ABI notes:
- All exported functions use ZX_CALL (\__cdecl on Windows) and ZX_API visibility.
*/

#ifndef ZX_TILES_API_H
#define ZX_TILES_API_H

#include "zx_abi.h"
#include "zx_tiles.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /** \brief Create a particle SoA with N elements.
   * @param count Number of particles to allocate
   * @param alloc_fn Host allocator callback (size,user) -> ptr (must not be NULL)
   * @param user User pointer passed to callbacks
   * @return A zx_particle_soa whose arrays are allocated; caller owns and must destroy
   */
  ZX_API zx_particle_soa ZX_CALL zx_create_particle_so_a(size_t count,
                                                         void* (*alloc_fn)(size_t, void*),
                                                         void* user);
  /** \brief Destroy a particle SoA and free its arrays via callback.
   * @param soa SoA to destroy (must not be NULL)
   * @param free_fn Host free callback (ptr,user) (must not be NULL)
   * @param user User pointer passed to callbacks
   */
  ZX_API void ZX_CALL zx_destroy_particle_so_a(zx_particle_soa* soa, void (*free_fn)(void*, void*),
                                               void* user);
  /** \brief Create a tile pool of capacity tiles.
   * @param capacity Number of tiles to allocate
   * @param alloc_fn Host allocator callback
   * @param user User pointer passed to callbacks
   * @return Pointer to tile array; caller owns and must destroy
   */
  ZX_API zx_tile* ZX_CALL zx_create_tile_pool(uint32_t capacity, void* (*alloc_fn)(size_t, void*),
                                              void* user);
  /** \brief Destroy a tile pool via callback.
   * @param pool Tile array pointer
   * @param free_fn Host free callback
   * @param user User pointer passed to callbacks
   */
  ZX_API void ZX_CALL zx_destroy_tile_pool(zx_tile* pool, void (*free_fn)(void*, void*),
                                           void* user);

#ifdef __cplusplus
}
#endif

#endif /* ZX_TILES_API_H */
