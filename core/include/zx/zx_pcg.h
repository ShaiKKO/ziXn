/*!
\file zx_pcg.h
\brief Preconditioned Conjugate Gradient (PCG) solver (C-ABI).
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
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

/** \brief Solve Ax=b using PCG with optional preconditioner.
 * @param n Dimension of A
 * @param b Right-hand side vector (size n)
 * @param x In: initial guess; Out: solution (size n)
 * @param apply_A Callback to compute y=Ax
 * @param apply_prec Callback to compute z=M^{-1} r (may be NULL for identity)
 * @param user User pointer passed to callbacks
 * @param opts Solver options (must not be NULL)
 * @param out_final_resid Out: final residual norm (may be NULL)
 * @return Iterations used (<= max_iters); 0 on error (x unchanged)
 */
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


