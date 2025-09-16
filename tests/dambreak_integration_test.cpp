/*!
 * \file dambreak_integration_test.cpp
 * \brief Dam-break over porous bed integration proxy test.
 */

#include "zx/zx_integration.h"
#include <cassert>

int main(){
    zx_dambreak_params p{};
    p.tiles = 1; p.h = 1.0f; p.dt = 0.05f; p.steps = 20; p.bed_k = 2;
    p.init_head_u = 1.0f; p.gamma_dot = 1.0f;
    p.hbp = zx_hbp_params{1.0f, 0.0f, 1.0f, 0.0f, 1.0f};
    p.permeability = 1e-10f; p.k_min = 1e-12f; p.mu_min = 0.01f; p.mu_max = 100.0f; p.beta_min = 1e2f; p.beta_max = 1e12f;
    p.policy = ZX_MU_CLAMP_HARD; p.softness_k = 1.0f;

    zx_dambreak_metrics m = zx_integration_dambreak_run(&p);
    assert(m.front_x > 0.0f);
    assert(m.kinetic_j > 0.0f);
    return 0;
}


