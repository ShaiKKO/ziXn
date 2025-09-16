/*!
 * \file determinism_tests.cpp
 * \brief Determinism: repeated runs must be identical.
 */

#include "zx/zx_integration.h"
#include "zx/zx_determinism.h"
#include <cassert>
#include <cstring>

static void run_dambreak_bytes(uint8_t* out, size_t out_size){
    std::memset(out, 0, out_size);
    zx_dambreak_params p{}; p.tiles=1; p.h=1.0f; p.dt=0.05f; p.steps=20; p.bed_k=2; p.init_head_u=1.0f; p.gamma_dot=1.0f;
    p.hbp = zx_hbp_params{1.0f,0.0f,1.0f,0.0f,1.0f}; p.permeability=1e-10f; p.k_min=1e-12f; p.mu_min=0.01f; p.mu_max=100.0f;
    p.beta_min=1e2f; p.beta_max=1e12f; p.policy=ZX_MU_CLAMP_HARD; p.softness_k=1.0f;
    zx_dambreak_metrics m = zx_integration_dambreak_run(&p);
    std::memcpy(out, &m, sizeof(m) < out_size ? sizeof(m) : out_size);
}

int main(){
    zx_set_determinism(1);
    zx_seed_rng(42);
    uint8_t a[64], b[64];
    run_dambreak_bytes(a, sizeof(a));
    run_dambreak_bytes(b, sizeof(b));
    assert(std::memcmp(a,b,sizeof(a)) == 0);
    return 0;
}


