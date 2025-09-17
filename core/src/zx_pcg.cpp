/*!
 * \file zx_pcg.cpp
 * \brief Preconditioned Conjugate Gradient (PCG) solver (CPU reference).
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_pcg.h"
#include <algorithm>
#include <cmath>
#include <vector>

extern "C" {

static inline float dot(size_t n, const float* a, const float* b){
    double s = 0.0; // extended precision accumulator
    for (size_t i=0;i<n;++i) s += static_cast<double>(a[i]) * static_cast<double>(b[i]);
    return static_cast<float>(s);
}

static inline void axpy(size_t n, float alpha, const float* x, float* y){
    for (size_t i=0;i<n;++i) y[i] += alpha * x[i];
}

static inline void scal(size_t n, float alpha, float* x){
    for (size_t i=0;i<n;++i) x[i] *= alpha;
}

/**
 * @brief Copy n float elements from source to destination.
 *
 * Copies exactly n elements from the array pointed to by x to the array
 * pointed to by y.
 *
 * @param n Number of elements to copy.
 * @param x Source pointer (must point to at least n elements).
 * @param y Destination pointer (must point to at least n elements).
 *
 * @note The source and destination ranges must not overlap; overlapping
 * ranges lead to undefined behavior.
 */
static inline void copy(size_t n, const float* x, float* y){
    std::copy(x, x+n, y);
}

/**
 * @brief Solve the symmetric positive-definite linear system A x = b using the
 *        Preconditioned Conjugate Gradient (PCG) method.
 *
 * Runs up to opts->max_iters iterations (at least 1) and stops when the
 * Euclidean norm of the residual r = b - A x is <= max(tol_abs, tol_rel * ||b||).
 * If apply_prec is non-NULL it is used to compute z = M^{-1} r; otherwise the
 * residual is used directly (no preconditioning).
 *
 * @param n Dimension of A (number of rows/columns and vector length).
 * @param b Right-hand side vector (length n).
 * @param x In: initial guess; Out: final solution (length n). On error x is left unchanged.
 * @param apply_A Callback invoked as apply_A(input, output, user) to compute output := A * input.
 * @param apply_prec Optional callback invoked as apply_prec(r, z, user) to compute z := M^{-1} r; may be NULL.
 * @param opts Solver options (must not be NULL): controls max_iters, tol_abs, tol_rel.
 * @param out_final_resid If non-NULL, receives the final residual norm ||b - A x||; may be set to 0 when ||b|| == 0.
 *
 * @return Number of iterations actually performed (<= opts->max_iters). Returns 0 on invalid input or degenerate ||b|| == 0;
 *         a numerical failure (non-positive or non-finite p^T A p) also terminates the iteration early.
 */
ZX_API uint32_t ZX_CALL zx_pcg_solve(size_t n,
                                     const float* b,
                                     float* x,
                                     zx_apply_A_fn apply_A,
                                     zx_apply_prec_fn apply_prec,
                                     void* user,
                                     const zx_pcg_opts* opts,
                                     float* out_final_resid)
{
    if (!b || !x || !apply_A || !opts || n == 0) return 0;
    const uint32_t max_iters = std::max<uint32_t>(1, opts->max_iters);
    const float tol_abs = (opts->tol_abs > 0.0f) ? opts->tol_abs : 0.0f;
    const float tol_rel = (opts->tol_rel > 0.0f) ? opts->tol_rel : 0.0f;

    std::vector<float> r(n, 0.0f), z(n, 0.0f), p(n, 0.0f), Ap(n, 0.0f);

    // r = b - A x
    apply_A(x, Ap.data(), user);
    for (size_t i=0;i<n;++i) r[i] = b[i] - Ap[i];

    float b_norm = std::sqrt(std::max(0.0f, dot(n, b, b)));
    if (b_norm == 0.0f) {
        if (out_final_resid) *out_final_resid = 0.0f;
        return 0;
    }

    if (apply_prec) apply_prec(r.data(), z.data(), user); else copy(n, r.data(), z.data());
    copy(n, z.data(), p.data());
    float rz_old = dot(n, r.data(), z.data());

    uint32_t k = 0;
    for (; k < max_iters; ++k) {
        apply_A(p.data(), Ap.data(), user);
        float pAp = dot(n, p.data(), Ap.data());
        if (pAp <= 0.0f || std::isnan(pAp) || std::isinf(pAp)) {
            break; // not SPD or numerical failure
        }
        float alpha = rz_old / pAp;
        axpy(n, +alpha, p.data(), x);  // x = x + alpha p
        axpy(n, -alpha, Ap.data(), r.data()); // r = r - alpha A p

        float r_norm = std::sqrt(std::max(0.0f, dot(n, r.data(), r.data())));
        float thresh = std::max(tol_abs, tol_rel * b_norm);
        if (r_norm <= thresh) { if (out_final_resid) *out_final_resid = r_norm; ++k; break; }

        if (apply_prec) apply_prec(r.data(), z.data(), user); else copy(n, r.data(), z.data());
        float rz_new = dot(n, r.data(), z.data());
        float beta = rz_new / std::max(1e-30f, rz_old);
        for (size_t i=0;i<n;++i) p[i] = z[i] + beta * p[i];
        rz_old = rz_new;
    }

    return k;
}

} // extern "C"


