/*!
 * \file zx_checksum.cpp
 * \brief Stable 64-bit checksum for grid/node snapshots.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_checksum.h"

namespace
{
  constexpr uint64_t k_mix1             = 0xff51afd7ed558ccdULL;
  constexpr uint64_t k_mix2             = 0xc4ceb9fe1a85ec53ULL;
  constexpr uint64_t k_fnv_offset_basis = 1469598103934665603ULL;
  constexpr int k_shift0                = 0;
  constexpr int k_shift1                = 16;
  constexpr int k_shift2                = 32;
  constexpr int k_shift3                = 48;
}  // namespace

static inline uint64_t mix64(uint64_t x)
{
  x ^= x >> 33;
  x *= k_mix1;
  x ^= x >> 33;
  x *= k_mix2;
  x ^= x >> 33;
  return x;
}

extern "C"
{

  uint64_t ZX_CALL zx_checksum_tile(const zx_tile* tile)
  {
    if (tile == nullptr)
    {
      return 0ULL;
    }
    uint64_t h = k_fnv_offset_basis;
    h ^= mix64(static_cast<uint64_t>(static_cast<uint32_t>(tile->coord_x)));
    h ^= mix64(static_cast<uint64_t>(static_cast<uint32_t>(tile->coord_y)));
    h ^= mix64(static_cast<uint64_t>(static_cast<uint32_t>(tile->coord_z)));
    for (int i = 0; i < ZX_TILE_B * ZX_TILE_B * ZX_TILE_B; ++i)
    {
      const zx_tile_node& n = tile->nodes[i];
      // consume three floats of momentum and one float of mass as 4x32; fold to 64
      uint64_t a = 0ULL;
      // NOLINTBEGIN(cppcoreguidelines-pro-type-reinterpret-cast)
      const auto* u32 = reinterpret_cast<const uint32_t*>(&n);
      // NOLINTEND(cppcoreguidelines-pro-type-reinterpret-cast)
      a ^= static_cast<uint64_t>(u32[0]) << k_shift0;  // mass
      a ^= static_cast<uint64_t>(u32[1]) << k_shift1;  // mom_x (low bits folded)
      a ^= static_cast<uint64_t>(u32[2]) << k_shift2;  // mom_y
      a ^= static_cast<uint64_t>(u32[3]) << k_shift3;  // mom_z
      h ^= mix64(a + static_cast<uint64_t>(i));
    }
    return h;
  }
}
