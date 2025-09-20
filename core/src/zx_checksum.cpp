/**
 * @file zx_checksum.cpp
 * @brief Stable 64-bit checksum for grid/node snapshots.
 * @details Implements a deterministic 64-bit hash over a tile's node payload using
 *          FNV-like mixing and shift-folding of node fields. Thread-safe and pure.
 * @copyright
 *   (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
 *   Licensed for use within the ziXn project under project terms.
 */

#include "zx/zx_checksum.h"
// Local includes for byte operations and fixed-size array view
#include <array>
#include <cstring>

namespace
{
  constexpr uint64_t k_mix1             = 0xff51afd7ed558ccdULL;
  constexpr uint64_t k_mix2             = 0xc4ceb9fe1a85ec53ULL;
  constexpr uint64_t k_fnv_offset_basis = 1469598103934665603ULL;
  constexpr int k_shift0                = 0;
  constexpr int k_shift1                = 16;
  constexpr int k_shift2                = 32;
  constexpr int k_shift3                = 48;
  constexpr int k_fold_shift            = 33;
  constexpr int k_mix_shift             = 33;
}  // namespace

static inline uint64_t mix64(uint64_t x)
{
  x ^= x >> k_mix_shift;
  x *= k_mix1;
  x ^= x >> k_mix_shift;
  x *= k_mix2;
  x ^= x >> k_mix_shift;
  return x;
}

extern "C"
{

  /**
   * @brief Compute a stable 64-bit checksum over a tile's node payload.
   *
   * Mixes tile coordinates and all node records into a deterministic 64-bit hash
   * suitable for parity checks and determinism tests.
   *
   * @param tile Pointer to tile; may be null (returns 0).
   * @return uint64_t Stable checksum value (0 if tile is null).
   */
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
    constexpr int k_tile_node_count = ZX_TILE_B * ZX_TILE_B * ZX_TILE_B;
    for (int i = 0; i < k_tile_node_count; ++i)
    {
      // Access by runtime index is required to walk nodes; struct layout is ABI-fixed.
      const zx_tile_node& n = tile->nodes[i];
      // consume three floats of momentum and one float of mass as 4x32; fold to 64
      uint64_t a = 0ULL;
      // Copy raw bytes into a 4x32 view to avoid array-to-pointer decay and subscript on
      // non-const-expr
      std::array<uint32_t, 4> u32{};
      static_assert(sizeof(u32) == sizeof(n), "size mismatch");
      std::memcpy(u32.data(), &n, sizeof(n));
      a ^= static_cast<uint64_t>(u32[0]) << k_shift0;  // mass
      a ^= static_cast<uint64_t>(u32[1]) << k_shift1;  // mom_x (low bits folded)
      a ^= static_cast<uint64_t>(u32[2]) << k_shift2;  // mom_y
      a ^= static_cast<uint64_t>(u32[3]) << k_shift3;  // mom_z
      // Shift on 64-bit to avoid UB from shifting a 32-bit int by >= 32
      h ^= mix64(a + (static_cast<uint64_t>(i) << k_fold_shift));
    }
    return h;
  }
}
