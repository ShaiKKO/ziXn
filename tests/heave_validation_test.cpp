/*!
 * \file heave_validation_test.cpp
 * \brief Validate heave profile monotonicity with respect to load and dilatancy.
 */

#include "zx/zx_heave_validation.h"
#include <vector>
#include <cassert>

int main(){
    zx_heave_params hp{ 1000.0f, 0.3f, 0.01f, 1.0f };
    float extent = 1.0f; uint32_t samples = (uint32_t)(2*extent/hp.grid_dx)+1;
    std::vector<float> y1(samples), y2(samples);
    zx_heave_profile_compute(&hp, extent, y1.data(), samples);
    float V1 = zx_heave_berm_volume(y1.data(), samples, hp.grid_dx);
    hp.load_newton *= 2.0f; hp.dilatancy *= 1.5f;
    zx_heave_profile_compute(&hp, extent, y2.data(), samples);
    float V2 = zx_heave_berm_volume(y2.data(), samples, hp.grid_dx);
    assert(V2 > V1);
    return 0;
}


