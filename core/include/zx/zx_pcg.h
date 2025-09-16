/*!
\file zx_pcg.h
\brief Preconditioned Conjugate Gradient (PCG) solver (C-ABI).
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
*/

#ifndef ZX_PCG_H
#define ZX_PCG_H

#include <stddef.h>
#include <stdint.h>

#include "zx_abi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct zx_pcg_opts {
    uint32_t max_iters;
    float tol_abs;
    float tol_rel;
} zx_pcg_opts;

typedef void (*zx_apply_A_fn)(const float* x, float* y, void* user);
typedef void (*zx_apply_prec_fn)(const float* r, float* z, void* user);

/* Returns number of iterations used (<=max_iters). On error returns 0 and leaves x unchanged. */
ZX_API uint32_t ZX_CALL zx_pcg_solve(size_t n,
                                     const float* b,
                                     float* x,
                                     zx_apply_A_fn apply_A,
                                     zx_apply_prec_fn apply_prec,
                                     void* user,
                                     const zx_pcg_opts* opts,
                                     float* out_final_resid);

#ifdef __cplusplus
}
#endif

#endif /* ZX_PCG_H */


