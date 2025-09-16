/*!
 * \file mc_return_map_test.cpp
 * \brief Sanity test for MC principal-space return with dilatancy.
 */

#include "zx/zx_constitutive_ref.h"
#include <cassert>
#include <cmath>

int main(){
    zx_elastic_params ep{ 20.0e6f, 0.3f };
    zx_mc_params mc{ 32.0f, 2.0f };
    float sigma_trial[9] = { 8.0e3f,0,0, 0,3.0e3f,0, 0,0,1.0e3f };
    float sigma_out[9]; float dgam=0; float epsv=0;
    zx_mc_return_map(&ep, &mc, 5.0f, sigma_trial, sigma_out, &dgam, &epsv);
    // Yield check proxy
    float I1 = sigma_out[0]+sigma_out[4]+sigma_out[8];
    float p = I1/3.0f; float tau_max = 0.5f*(sigma_out[0]-sigma_out[8]);
    float alpha = std::sin(32.0f * 3.14159265358979323846f / 180.0f) / (1.0f - std::sin(32.0f * 3.14159265358979323846f / 180.0f));
    float c_pa = 2.0f * 1000.0f;
    float f = tau_max + alpha * (-p) - c_pa;
    assert(f <= 5.0e2f); // within small tolerance for simple proxy
    return 0;
}


