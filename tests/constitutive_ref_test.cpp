/*!
 * \file constitutive_ref_test.cpp
 * \brief Basic tests for MC->DP mapping, elastic Lamé, and DP return mapping.
 */

#include "zx/zx_constitutive_ref.h"
#include <cassert>
#include <cmath>

int main(){
    zx_mc_params mc{ 30.0f, 2.0f };
    zx_dp_params dp{}; zx_mc_to_dp(&mc, &dp);
    assert(dp.alpha > 0 && dp.k > 0);

    zx_elastic_params ep{ 10.0e6f, 0.3f };
    float lam=0, mu=0; zx_elastic_lame(&ep, &lam, &mu);
    assert(lam > 0 && mu > 0);

    float sigma_trial[9] = { 5e3f,0,0, 0,5e3f,0, 0,0,5e3f };
    float sigma_out[9]; float dgam=0;
    zx_cap_params cap{1, -2.0e4f};
    zx_dp_return_map(&ep, &dp, &cap, sigma_trial, sigma_out, &dgam);
    // Expect output to lie on or within yield: alpha*I1 + sqrt(J2) - k <= 0
    float I1, J2; zx_stress_invariants(sigma_out, &I1, &J2);
    float f = dp.alpha * I1 + std::sqrt(J2) - dp.k;
    assert(f <= 1e-3f);
    return 0;
}


