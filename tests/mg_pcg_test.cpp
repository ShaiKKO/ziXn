/*!
 * \file mg_pcg_test.cpp
 * \brief Simple geometric 1D Poisson system solved by PCG with Jacobi smoother as preconditioner.
 */

#include "zx/zx_pcg.h"
#include <vector>
#include <cassert>
#include <cmath>

struct Poisson1D {
    size_t n;
    float h2inv;
};

static void apply_A(const float* x, float* y, void* user){
    Poisson1D* p = (Poisson1D*)user;
    size_t n = p->n;
    float h2i = p->h2inv;
    y[0] = h2i * (2.0f * x[0] - x[1]);
    for (size_t i=1;i<n-1;++i) y[i] = h2i * (-x[i-1] + 2.0f*x[i] - x[i+1]);
    y[n-1] = h2i * (2.0f * x[n-1] - x[n-2]);
}

static void apply_jacobi_prec(const float* r, float* z, void* user){
    Poisson1D* p = (Poisson1D*)user;
    float diag = 2.0f * p->h2inv;
    size_t n = p->n;
    for (size_t i=0;i<n;++i) z[i] = r[i] / diag;
}

int main(){
    const size_t n = 64; // fine grid
    Poisson1D ctx{n, 1.0f};
    std::vector<float> b(n, 1.0f), x(n, 0.0f);
    zx_pcg_opts opts{200, 1e-6f, 1e-6f};
    float final_res = 0.0f;
    uint32_t iters = zx_pcg_solve(n, b.data(), x.data(), apply_A, apply_jacobi_prec, &ctx, &opts, &final_res);
    assert(iters > 0);
    // Residual check
    std::vector<float> Ax(n, 0.0f); apply_A(x.data(), Ax.data(), &ctx);
    float err2=0.0f; for (size_t i=0;i<n;++i){ float e = Ax[i]-b[i]; err2 += e*e; }
    assert(std::sqrt(err2) < 1e-3f);
    return 0;
}


