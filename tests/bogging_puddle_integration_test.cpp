/*!
 * \file bogging_puddle_integration_test.cpp
 * \brief Wheel bogging and puddle creep integration proxy tests.
 */

#include "zx/zx_integration.h"
#include <cassert>

int main(){
    // Wheel bogging: expect some drag and non-negative sink depth
    zx_bogging_params bp{}; bp.h=1.0f; bp.dt=0.05f; bp.steps=20; bp.wheel_radius_nodes=4; bp.wheel_pull_u=1.0f; bp.wheel_push_w=0.1f;
    bp.hbp = zx_hbp_params{1.0f,0.0f,1.0f,0.0f,1.0f}; bp.mu_min=0.01f; bp.mu_max=100.0f; bp.permeability=1e-10f; bp.k_min=1e-12f;
    bp.beta_min=1e2f; bp.beta_max=1e12f; bp.gamma_dot=1.0f; bp.policy=ZX_MU_CLAMP_HARD; bp.softness_k=1.0f;
    zx_bogging_metrics bm = zx_integration_wheel_bogging_run(&bp);
    assert(bm.drag_N > 0.0f);
    assert(bm.sink_depth_m >= 0.0f);

    // Puddle creep: expect positive creep distance
    zx_puddle_params pp{}; pp.h=1.0f; pp.dt=0.05f; pp.steps=20; pp.init_head_u=1.0f;
    pp.hbp = zx_hbp_params{1.0f,0.0f,1.0f,0.0f,1.0f}; pp.mu_min=0.01f; pp.mu_max=100.0f; pp.permeability=1e-10f; pp.k_min=1e-12f;
    pp.beta_min=1e2f; pp.beta_max=1e12f; pp.gamma_dot=1.0f; pp.policy=ZX_MU_CLAMP_HARD; pp.softness_k=1.0f;
    zx_puddle_metrics pm = zx_integration_puddle_creep_run(&pp);
    assert(pm.creep_dist_x > 0.0f);
    return 0;
}


