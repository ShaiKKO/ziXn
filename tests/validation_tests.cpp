/*!
 * \file validation_tests.cpp
 * \brief Tests for inclined plane and column collapse proxies.
 */

#include "zx/zx_validation.h"
#include <cassert>
#include <cmath>

int main(){
    zx_mc_params mc{ 30.0f, 0.0f };
    float theta_c = zx_validation_inclined_plane_theta_c(&mc);
    assert(std::fabs(theta_c - 30.0f * 3.14159265358979323846f/180.0f) < 1e-5f);

    float r1 = zx_validation_column_collapse_runout_ratio(30.0f, 1.0f);
    float r2 = zx_validation_column_collapse_runout_ratio(35.0f, 1.0f);
    assert(r1 > r2); // higher φ -> shorter runout

    float r3 = zx_validation_column_collapse_runout_ratio(30.0f, 2.0f);
    assert(r3 > r1); // higher aspect ratio -> longer runout
    return 0;
}


