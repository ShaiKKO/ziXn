/*!
\file zx_tiles_api.h
\brief C symbols for tiles/SoA allocation to support early wiring/tests.
\author Colin Macritchie (Ripple Group, LLC)
*/

#ifndef ZX_TILES_API_H
#define ZX_TILES_API_H

#include "zx_abi.h"
#include "zx_tiles.h"

#ifdef __cplusplus
extern "C" {
#endif

ZX_API zx_particle_soa ZX_CALL zxCreateParticleSoA(size_t count, void* (*alloc_fn)(size_t, void*), void* user);
ZX_API void ZX_CALL zxDestroyParticleSoA(zx_particle_soa* soa, void (*free_fn)(void*, void*), void* user);
ZX_API zx_tile* ZX_CALL zxCreateTilePool(uint32_t capacity, void* (*alloc_fn)(size_t, void*), void* user);
ZX_API void ZX_CALL zxDestroyTilePool(zx_tile* pool, void (*free_fn)(void*, void*), void* user);

#ifdef __cplusplus
}
#endif

#endif /* ZX_TILES_API_H */


