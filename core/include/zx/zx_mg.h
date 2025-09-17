/*!
\file zx_mg.h
\brief Geometric Multigrid preconditioner (1D Poisson V-cycle) for PCG.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_MG_H
#define ZX_MG_H

#include <stddef.h>
#include <stdint.h>

#include "zx_abi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct zx_mg_context zx_mg_context; /* opaque */

typedef struct zx_mg_opts {
    uint32_t max_levels;    /* maximum levels including finest */
    uint32_t pre_smooth;    /* Jacobi pre-smoothing iterations */
    uint32_t post_smooth;   /* Jacobi post-smoothing iterations */
    uint32_t coarse_iters;  /* Jacobi iterations on coarsest level */
    float omega;            /* weighted Jacobi relaxation (0,1] */
} zx_mg_opts;

/** \brief Build a 1D Poisson hierarchy for size n (n>=3).
 * @param n Problem size (>=3)
 * @param opts Options (may be NULL for defaults)
 * @return Context pointer or NULL on failure
 */
ZX_API zx_mg_context* ZX_CALL zx_mg_create_poisson1d(size_t n, const zx_mg_opts* opts);

/** \brief Destroy multigrid context. */
ZX_API void ZX_CALL zx_mg_destroy(zx_mg_context* ctx);

/** \brief Preconditioner apply compatible with zx_pcg: z = M^{-1} r. */
ZX_API void ZX_CALL zx_mg_prec_apply(const float* r, float* z, void* user_ctx);

#ifdef __cplusplus
}
#endif

#endif /* ZX_MG_H */


