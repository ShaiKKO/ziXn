/*!
 * \file hbp_clamp_policy_test.cpp
 * \brief Tests for μ_eff clamp policies.
 */

#include "zx/zx_hbp.h"
#include <cassert>
#include <cmath>

static bool le(float a, float b, float eps=1e-6f){ return a <= b + eps; }
static bool ge(float a, float b, float eps=1e-6f){ return a + eps >= b; }

int main(){
    zx_hbp_params p{0.5f, 10.0f, 0.5f, 5.0f, 10.0f};
    // Force very large mu via parameters and gamma_dot small
    float mu_raw = zx_hbp_mu_eff_policy(1e-9f, &p, -1e9f, 1e9f, ZX_MU_CLAMP_NONE, 0.0f);
    // Hard clamp should enforce bounds exactly
    float mu_hard = zx_hbp_mu_eff_policy(1e-9f, &p, 0.1f, 3.0f, ZX_MU_CLAMP_HARD, 0.0f);
    assert(ge(mu_hard, 0.1f) && le(mu_hard, 3.0f));

    // Smooth clamp should also respect bounds but allow soft approach
    float mu_smooth = zx_hbp_mu_eff_policy(1e-9f, &p, 0.1f, 3.0f, ZX_MU_CLAMP_SMOOTH_TANH, 2.0f);
    assert(ge(mu_smooth, 0.1f) && le(mu_smooth, 3.0f));
    (void)mu_raw;
    return 0;
}


