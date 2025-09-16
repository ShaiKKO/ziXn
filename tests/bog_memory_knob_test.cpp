/*!
 * \file bog_memory_knob_test.cpp
 * \brief Validate bog memory beta trend with shear rate and clamping.
 */

#include "zx/zx_mixture.h"
#include <cassert>
#include <cmath>

static bool le(float a, float b, float eps=1e-6f){ return a <= b + eps; }
static bool ge(float a, float b, float eps=1e-6f){ return a + eps >= b; }

int main(){
    zx_hbp_params hbp{0.5f, 2.0f, 0.6f, 5.0f, 8.0f};
    float k = 1e-10f; // m^2
    float k_min = 1e-12f;
    float mu_min = 0.01f, mu_max = 1000.0f;
    float beta_min = 1e3f, beta_max = 1e12f; // Pa·s/m^2

    float beta_lo = zx_bog_beta_hbp(1e-4f, &hbp, k, k_min, mu_min, mu_max, beta_min, beta_max, ZX_MU_CLAMP_HARD, 0.0f);
    float beta_hi = zx_bog_beta_hbp(10.0f, &hbp, k, k_min, mu_min, mu_max, beta_min, beta_max, ZX_MU_CLAMP_HARD, 0.0f);
    // For n<1 and tau_y>0, effective viscosity decreases with shear, so beta should not increase
    assert(le(beta_hi, beta_lo));

    // Clamping behavior
    float beta_smallk = zx_bog_beta_hbp(0.0f, &hbp, 0.0f, k_min, mu_min, mu_max, beta_min, beta_max, ZX_MU_CLAMP_SMOOTH_TANH, 2.0f);
    assert(ge(beta_smallk, beta_min) && le(beta_smallk, beta_max));
    return 0;
}


