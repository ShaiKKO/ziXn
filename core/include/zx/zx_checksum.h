/*!
\file zx_checksum.h
\brief Stable 64-bit checksum for small grid/node snapshots.
\author Colin Macritchie (Ripple Group, LLC)
*/

#ifndef ZX_CHECKSUM_H
#define ZX_CHECKSUM_H

#include <stdint.h>
#include "zx_abi.h"
#include "zx_tiles.h"

#ifdef __cplusplus
extern "C" {
#endif

ZX_API uint64_t ZX_CALL zx_checksum_tile(const zx_tile* tile);

#ifdef __cplusplus
}
#endif

#endif /* ZX_CHECKSUM_H */


