/*!
 * \file pcg_smoke_test.cpp
 * \brief Smoke test for PCG on SPD diagonal system.
 */

#include "zx/zx_pcg.h"
#include <vector>
#include <cassert>
#include <cmath>

struct DiagCtx { std::vector<float> d; };

static void apply_A(const float* x, float* y, void* user){
    DiagCtx* c = (DiagCtx*)user;
    size_t n = c->d.size();
    for (size_t i=0;i<n;++i) y[i] = c->d[i] * x[i];
}

static void apply_prec(const float* r, float* z, void* user){
    DiagCtx* c = (DiagCtx*)user;
    size_t n = c->d.size();
    for (size_t i=0;i<n;++i) z[i] = r[i] / c->d[i];
}

int main(){
    const size_t n = 16;
    DiagCtx ctx; ctx.d.resize(n);
    for (size_t i=0;i<n;++i) ctx.d[i] = 1.0f + (float)i;
    std::vector<float> b(n, 1.0f), x(n, 0.0f);
    zx_pcg_opts opts{100, 1e-6f, 1e-6f};
    float final_res = 0.0f;
    uint32_t iters = zx_pcg_solve(n, b.data(), x.data(), apply_A, apply_prec, &ctx, &opts, &final_res);
    assert(iters > 0);
    // Check A x ≈ b
    std::vector<float> Ax(n, 0.0f);
    apply_A(x.data(), Ax.data(), &ctx);
    float err = 0.0f; for (size_t i=0;i<n;++i) err += (Ax[i] - b[i])*(Ax[i] - b[i]);
    assert(std::sqrt(err) < 1e-4f);
    (void)final_res;
    return 0;
}


