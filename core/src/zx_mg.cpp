/*!
 * \file zx_mg.cpp
 * \brief Geometric Multigrid preconditioner (1D Poisson V-cycle) implementation.
 * \author Colin Macritchie (Ripple Group, LLC)
 */

#include "zx/zx_mg.h"
#include <vector>
#include <algorithm>
#include <cmath>

struct zx_mg_level {
    size_t n;
    float h2inv;
    std::vector<float> x, r, z, tmp;
};

struct zx_mg_context {
    std::vector<zx_mg_level> levels;
    zx_mg_opts opts;
};

static inline void apply_A_1d(const float* x, float* y, size_t n, float h2inv){
    y[0] = h2inv * (2.0f * x[0] - x[1]);
    for (size_t i=1;i<n-1;++i) y[i] = h2inv * (-x[i-1] + 2.0f*x[i] - x[i+1]);
    y[n-1] = h2inv * (2.0f * x[n-1] - x[n-2]);
}

static void jacobi_smooth(size_t iters, float omega, zx_mg_level& L){
    const size_t n = L.n;
    std::vector<float> Ax(n, 0.0f);
    for (size_t k=0;k<iters;++k){
        apply_A_1d(L.x.data(), Ax.data(), n, L.h2inv);
        for (size_t i=0;i<n;++i){
            float diag = 2.0f * L.h2inv;
            float res = L.r[i] - Ax[i];
            L.x[i] += omega * res / diag;
        }
    }
}

static void restrict_full_weighting(const std::vector<float>& fine, std::vector<float>& coarse){
    const size_t nf = fine.size();
    const size_t nc = coarse.size();
    // Assume nf = 2*nc - 1 (standard 1D restriction)
    coarse[0] = 0.5f * (fine[0] + fine[1]);
    for (size_t i=1;i<nc-1;++i){
        size_t j = 2*i;
        coarse[i] = 0.25f * (fine[j-1] + 2.0f*fine[j] + fine[j+1]);
    }
    coarse[nc-1] = 0.5f * (fine[nf-2] + fine[nf-1]);
}

static void prolong_linear(const std::vector<float>& coarse, std::vector<float>& fine){
    const size_t nc = coarse.size();
    const size_t nf = fine.size();
    // nf = 2*nc - 1
    for (size_t i=0;i<nc;++i) fine[2*i] = coarse[i];
    for (size_t i=0;i<nc-1;++i) fine[2*i+1] = 0.5f * (coarse[i] + coarse[i+1]);
}

static void vcycle(zx_mg_context* ctx, size_t level_idx){
    zx_mg_level& L = ctx->levels[level_idx];
    const zx_mg_opts& o = ctx->opts;
    jacobi_smooth(o.pre_smooth, o.omega, L);

    if (level_idx + 1 < ctx->levels.size()){
        // compute residual r - A x
        std::vector<float> Ax(L.n, 0.0f);
        apply_A_1d(L.x.data(), Ax.data(), L.n, L.h2inv);
        for (size_t i=0;i<L.n;++i) L.tmp[i] = L.r[i] - Ax[i];
        // restrict to coarse
        zx_mg_level& C = ctx->levels[level_idx+1];
        restrict_full_weighting(L.tmp, C.r);
        std::fill(C.x.begin(), C.x.end(), 0.0f);
        vcycle(ctx, level_idx+1);
        // prolong and correct
        prolong_linear(C.x, L.tmp);
        for (size_t i=0;i<L.n;++i) L.x[i] += L.tmp[i];
        jacobi_smooth(o.post_smooth, o.omega, L);
    } else {
        jacobi_smooth(o.coarse_iters, o.omega, L);
    }
}

extern "C" {

ZX_API zx_mg_context* ZX_CALL zx_mg_create_poisson1d(size_t n, const zx_mg_opts* opts){
    if (n < 3 || !opts) return nullptr;
    zx_mg_context* ctx = new zx_mg_context();
    ctx->opts = *opts;
    size_t cur = n;
    float h2inv = 1.0f;
    for (uint32_t lvl=0; lvl<opts->max_levels && cur >= 3; ++lvl){
        zx_mg_level L; L.n = cur; L.h2inv = h2inv; L.x.assign(cur,0.0f); L.r.assign(cur,0.0f); L.z.assign(cur,0.0f); L.tmp.assign(cur,0.0f);
        ctx->levels.push_back(std::move(L));
        if (cur == 3) break;
        size_t next = (cur + 1) / 2; // ensure nf ≈ 2*nc - 1
        if (next < 3) break;
        cur = next;
        h2inv *= 4.0f; // grid spacing doubles each level -> 1/h^2 scales by 1/4; here we keep h^2 consistent per stencil scaling usage
    }
    if (ctx->levels.empty()) { delete ctx; return nullptr; }
    return ctx;
}

ZX_API void ZX_CALL zx_mg_destroy(zx_mg_context* ctx){ delete ctx; }

ZX_API void ZX_CALL zx_mg_prec_apply(const float* r, float* z, void* user_ctx){
    zx_mg_context* ctx = (zx_mg_context*)user_ctx;
    if (!ctx || !r || !z) return;
    // load RHS into finest level
    zx_mg_level& F = ctx->levels.front();
    std::copy(r, r + F.n, F.r.begin());
    std::fill(F.x.begin(), F.x.end(), 0.0f);
    vcycle(ctx, 0);
    std::copy(F.x.begin(), F.x.end(), z);
}

} // extern "C"


