/*!
\file zx_determinism.h
\brief Global determinism controls (for CPU reference and tests).
\author Colin Macritchie (Ripple Group, LLC)
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_DETERMINISM_H
#define ZX_DETERMINISM_H

#include "zx_abi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \brief Enable or disable global determinism. */
ZX_API void ZX_CALL zx_set_determinism(int enable);
/** \brief Get current determinism flag (0/1). */
ZX_API int  ZX_CALL zx_get_determinism(void);
/** \brief Seed the global RNG to a fixed value. */
ZX_API void ZX_CALL zx_seed_rng(uint64_t seed);
/** \brief Retrieve the current RNG seed. */
ZX_API uint64_t ZX_CALL zx_get_rng_seed(void);

#ifdef __cplusplus
}
#endif

#endif /* ZX_DETERMINISM_H */


