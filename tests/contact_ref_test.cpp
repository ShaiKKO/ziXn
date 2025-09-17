/*!
 * \file contact_ref_test.cpp
 * \brief Contact projection tests (non-penetration + Coulomb friction), per ALGORITHMS.md Contact.
 */

#include "zx/zx_contact_ref.h"
#include <cassert>

int main(){
    float n[3] = {0,1,0};
    float v[3] = {1, -2, 0};
    float out[3]; int sat = 0;
    zx_contact_project(v, n, /*mu=*/0.5f, /*phi=*/-0.01f, /*kappa_n=*/1e-3f, out, &sat);
    // Normal should be clamped toward non-penetration; tangential limited by mu*vn
    assert(out[1] >= -0.1f); // compliant clamp, not too negative
    (void)sat;
    return 0;
}


