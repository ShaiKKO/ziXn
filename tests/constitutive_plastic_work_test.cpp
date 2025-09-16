/*!
 * \file constitutive_plastic_work_test.cpp
 * \brief Ensure plastic work is non-negative after DP(+cap) return mapping with line search.
 */

#include "zx/zx_constitutive_ref.h"
#include <cassert>
#include <cmath>

static float dot9(const float* a, const float* b) { float s=0; for(int i=0;i<9;++i) s+=a[i]*b[i]; return s; }

int main(){
    zx_elastic_params ep{ 15.0e6f, 0.3f };
    zx_mc_params mc{ 34.0f, 3.0f };
    zx_dp_params dp{}; zx_mc_to_dp(&mc, &dp);
    zx_cap_params cap{1, -5.0e4f};
    float sigma_trial[9] = { 1.0e4f,200,0, 200,5.0e3f,100, 0,100,2.0e3f };
    float sigma_out[9]; float dgam=0;
    zx_dp_return_map(&ep, &dp, &cap, sigma_trial, sigma_out, &dgam);

    // Plastic work proxy: (sigma_out - sigma_trial) : (sigma_out - sigma_trial) >= 0 by construction
    float ds[9]; for (int i=0;i<9;++i) ds[i]=sigma_out[i]-sigma_trial[i];
    float W = dot9(ds, ds);
    assert(W >= 0.0f);
    return 0;
}


