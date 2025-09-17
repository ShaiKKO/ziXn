/*!
 * \file mg_prec_test.cpp
 * \brief Test PCG with multigrid preconditioner on 1D Poisson.
 */

#include "zx/zx_pcg.h"
#include "zx/zx_mg.h"
#include <vector>
#include <cassert>
#include <cmath>

struct Poisson1D { size_t n; float h2inv; };
struct CombinedCtx { Poisson1D A; zx_mg_context* mg; };

static void apply_A(const float* x, float* y, void* user){
    CombinedCtx* c = (CombinedCtx*)user;
    size_t n = c->A.n; float h2i = c->A.h2inv;
    y[0] = h2i * (2.0f * x[0] - x[1]);
    for (size_t i=1;i<n-1;++i) y[i] = h2i * (-x[i-1] + 2.0f*x[i] - x[i+1]);
    y[n-1] = h2i * (2.0f * x[n-1] - x[n-2]);
}

static void apply_mg_prec(const float* r, float* z, void* user){
    CombinedCtx* c = (CombinedCtx*)user;
    zx_mg_prec_apply(r, z, c->mg);
}

int main(){
    const size_t n = 257; // power-of-two plus one for proper coarsening
    zx_mg_opts mgo{6, 2, 2, 8, 0.8f};
    CombinedCtx ctx; ctx.A = Poisson1D{n, 1.0f}; ctx.mg = zx_mg_create_poisson1d(n, &mgo);
    assert(ctx.mg);
    std::vector<float> b(n, 1.0f), x(n, 0.0f);
    zx_pcg_opts opts{200, 1e-6f, 1e-6f};
    float final_res = 0.0f;
    uint32_t iters = zx_pcg_solve(n, b.data(), x.data(), apply_A, apply_mg_prec, &ctx, &opts, &final_res);
    assert(iters > 0);
    std::vector<float> Ax(n, 0.0f); apply_A(x.data(), Ax.data(), &ctx);
    float err2 = 0.0f; for (size_t i=0;i<n;++i){ float e = Ax[i]-b[i]; err2 += e*e; }
    assert(std::sqrt(err2) < 1e-3f);
    zx_mg_destroy(ctx.mg);
    return 0;
}


