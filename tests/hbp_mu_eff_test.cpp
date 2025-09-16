/*!
 * \file hbp_mu_eff_test.cpp
 * \brief Unit tests for HBP effective viscosity function.
 */

#include "zx/zx_hbp.h"
#include <cassert>
#include <cmath>
#include <algorithm>

static bool near(float a, float b, float eps=1e-4f){ return std::fabs(a-b) <= eps * (1.0f + std::max(std::fabs(a), std::fabs(b))); }

int main(){
    zx_hbp_params p{0.5f, 0.0f, 0.7f, 3.0f, 10.0f};
    {
        float mu = zx_hbp_mu_eff(0.0f, &p, 0.0f, 100.0f);
        float expected = p.mu0 + p.tau_y * p.m; // K=0
        assert(near(mu, expected));
    }
    {
        float g = 1.0f;
        float mu = zx_hbp_mu_eff(g, &p, 0.0f, 100.0f);
        float expected = p.mu0 + p.tau_y * (1.0f - std::exp(-p.m * g)) / g;
        assert(near(mu, expected));
    }
    {
        float mu = zx_hbp_mu_eff(0.0f, &p, 0.0f, 1.0f);
        assert(near(mu, 1.0f)); // clamped to mu_max
    }
    return 0;
}


