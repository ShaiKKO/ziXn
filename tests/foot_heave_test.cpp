/*!
 * \file foot_heave_test.cpp
 * \brief Foot-like heave: increasing dilatancy increases berm volume under fixed load.
 */

#include "zx/zx_heave_validation.h"
#include <vector>
#include <cassert>

int main(){
    const float extent = 0.5f;
    zx_heave_params hp{ 800.0f, 0.10f, 0.005f, 0.8f };
    uint32_t samples = (uint32_t)(2*extent/hp.grid_dx)+1;
    std::vector<float> y0(samples), y1(samples);
    zx_heave_profile_compute(&hp, extent, y0.data(), samples);
    float V0 = zx_heave_berm_volume(y0.data(), samples, hp.grid_dx);
    hp.dilatancy = 1.6f;
    zx_heave_profile_compute(&hp, extent, y1.data(), samples);
    float V1 = zx_heave_berm_volume(y1.data(), samples, hp.grid_dx);
    assert(V1 > V0);
    return 0;
}


