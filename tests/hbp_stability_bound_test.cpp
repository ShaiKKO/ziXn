/*!
 * \file hbp_stability_bound_test.cpp
 * \brief Tests for explicit diffusion stability bound.
 */

#include "zx/zx_hbp.h"
#include <cassert>
#include <cmath>
#include <algorithm>

static bool near(float a, float b, float eps=1e-6f){ return std::fabs(a-b) <= eps * (1.0f + std::max(std::fabs(a), std::fabs(b))); }

int main(){
    float h = 0.5f; float mu_max = 2.0f; float rho_min = 1.0f;
    float dt_max = zx_hbp_dt_stable_upper_bound(h, mu_max, rho_min);
    float expected = (h*h * rho_min) / (6.0f * mu_max);
    assert(near(dt_max, expected));

    // Degenerate inputs -> zero
    assert(zx_hbp_dt_stable_upper_bound(-1.0f, mu_max, rho_min) == 0.0f);
    assert(zx_hbp_dt_stable_upper_bound(h, -1.0f, rho_min) == 0.0f);
    assert(zx_hbp_dt_stable_upper_bound(h, mu_max, -1.0f) == 0.0f);
    return 0;
}


