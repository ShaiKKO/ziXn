/*!
 * \file anisotropic_contact_test.cpp
 * \brief Validate anisotropic contact ellipse and rolling resistance proxy.
 */

#include "zx/zx_contact_ref.h"
#include <cassert>
#include <cmath>

int main(){
    float n[3]={0,1,0}; float t0[3]={1,0,0}; float t1[3]={0,0,1};
    float v[3]={2.0f,-1.0f,1.5f}; float out[3]; int sat=0;
    zx_contact_project_aniso(v,n,t0,t1, 0.8f, 0.4f, 0.5f, -0.01f, 1e-3f, out, &sat);
    // Expect saturation due to anisotropy (lateral friction lower)
    assert(sat == 1);
    // Normal component should be near non-penetration
    assert(out[1] >= -0.1f);
    return 0;
}


