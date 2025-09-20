/**
 * @file zx_checksum.h
 * @brief Stable 64-bit checksum for small grid/node snapshots.
 * @details Declares C-ABI checksum function to hash a tile's node payload deterministically.
 *          Pure and thread-safe.
 */

#ifndef ZX_CHECKSUM_H
#define ZX_CHECKSUM_H

#include "zx_abi.h"
#include "zx_tiles.h"
#include <cstdint>

#ifdef __cplusplus
extern "C"
{
#endif

  ZX_API uint64_t ZX_CALL zx_checksum_tile(const zx_tile* tile);

#ifdef __cplusplus
}
#endif

#endif /* ZX_CHECKSUM_H */
