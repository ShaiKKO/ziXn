/*!
 * \file zx_checksum.cpp
 * \brief Stable 64-bit checksum for grid/node snapshots.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_checksum.h"

static inline uint64_t mix64(uint64_t x){
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x;
}

extern "C" {

ZX_API uint64_t ZX_CALL zx_checksum_tile(const zx_tile* tile){
    if (!tile) return 0ULL;
    uint64_t h = 1469598103934665603ULL;
    h ^= mix64((uint64_t)(uint32_t)tile->coord_x);
    h ^= mix64((uint64_t)(uint32_t)tile->coord_y);
    h ^= mix64((uint64_t)(uint32_t)tile->coord_z);
    for (int i=0;i<ZX_TILE_B*ZX_TILE_B*ZX_TILE_B;++i){
        const zx_tile_node& n = tile->nodes[i];
        const uint64_t* p = reinterpret_cast<const uint64_t*>(&n);
        // consume three floats of momentum and one float of mass as 4x32; fold to 64
        uint64_t a = 0; const uint32_t* u32 = reinterpret_cast<const uint32_t*>(&n);
        a ^= (uint64_t)u32[0] <<  0; // mass
        a ^= (uint64_t)u32[1] << 16; // mom_x (low bits folded)
        a ^= (uint64_t)u32[2] << 32; // mom_y
        a ^= (uint64_t)u32[3] << 48; // mom_z
        h ^= mix64(a + (uint64_t)i);
    }
    return h;
}

}


