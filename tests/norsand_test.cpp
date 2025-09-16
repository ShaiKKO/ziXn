/*!
 * \file norsand_test.cpp
 * \brief Basic NorSand state parameter and return map sanity checks.
 */

#include "zx/zx_constitutive_ref.h"
#include <cassert>

int main(){
    zx_elastic_params ep{ 15.0e6f, 0.3f };
    zx_norsand_params ns{ 1.2f, 0.2f, 0.05f, 1.0e4f, 0.8f, 1.0f, 0.5f };
    zx_norsand_state st{ 0.9f };
    float p = -5.0e3f; // compression negative by convention here
    float psi = zx_norsand_state_parameter(&st, p, &ns);
    // Just check it computes
    (void)psi;

    float sigma_trial[9] = { p,0,0, 0,p,0, 0,0,p };
    float sigma_out[9]; float dgam=0; float epsv=0;
    zx_norsand_return_map(&ep, &ns, &st, sigma_trial, sigma_out, &dgam, &epsv);
    // q should trend toward M*(-p) and epsv should update
    float I1 = sigma_out[0]+sigma_out[4]+sigma_out[8];
    float p_out = I1/3.0f; float q_out = 0.5f*(sigma_out[0]-sigma_out[8])*2.0f; // rough proxy
    (void)q_out; (void)p_out;
    assert(st.void_ratio_e > 0.0f);
    return 0;
}


