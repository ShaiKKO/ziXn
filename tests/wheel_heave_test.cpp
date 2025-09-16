/*!
 * \file wheel_heave_test.cpp
 * \brief Wheel-like heave: narrower footprint yields higher peak heave for same load.
 */

#include "zx/zx_heave_validation.h"
#include <vector>
#include <algorithm>
#include <cassert>

static float max_value(const std::vector<float>& y) {
    return *std::max_element(y.begin(), y.end());
}

int main(){
    const float extent = 1.0f;
    zx_heave_params hp{ 1500.0f, 0.40f, 0.01f, 1.0f };
    uint32_t samples = (uint32_t)(2*extent/hp.grid_dx)+1;
    std::vector<float> y_wide(samples), y_narrow(samples);
    zx_heave_profile_compute(&hp, extent, y_wide.data(), samples);

    hp.footprint_width = 0.20f; // narrower wheel contact
    zx_heave_profile_compute(&hp, extent, y_narrow.data(), samples);

    float peak_wide = max_value(y_wide);
    float peak_narrow = max_value(y_narrow);
    assert(peak_narrow > peak_wide);
    return 0;
}


