/*!
 * \file anisotropy_scene_test.cpp
 * \brief Scene-like validation: ellipse ratio affects longitudinal vs lateral clamp.
 */

#include "zx/zx_contact_ref.h"
#include <cassert>

static void proj(float mulong, float mulat, float v0_in, float v1_in, float& v0_out, float& v1_out)
{
    float n[3]={0,1,0}; float t0[3]={1,0,0}; float t1[3]={0,0,1};
    float vin[3]={v0_in, 1.0f, v1_in}; float vout[3]; int sat=0;
    zx_contact_project_aniso(vin, n, t0, t1, mulong, mulat, 0.0f, -0.01f, 1e-3f, vout, &sat);
    v0_out = vout[0]; v1_out = vout[2];
}

int main(){
    float v0a, v1a, v0b, v1b;
    proj(0.9f, 0.9f, 1.0f, 1.0f, v0a, v1a); // isotropic
    proj(0.9f, 0.3f, 1.0f, 1.0f, v0b, v1b); // anisotropic: lateral lower
    // Expect lateral tangential reduced more than longitudinal when mulat < mulong
    assert((v1b - v1a) < (v0b - v0a));
    return 0;
}


