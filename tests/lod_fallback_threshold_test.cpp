/*!
 * \file lod_fallback_threshold_test.cpp
 * \brief Threshold + hysteresis behavior for LOD fallback policy.
 */

#include "zx/zx_lod.h"
#include <cassert>

int main(){
    zx_lod_fallback_policy p{ /*active_tiles_max*/10, /*step_ms_max*/5.0f, /*enter*/2, /*exit*/2, /*blend*/3 };
    zx_lod_fallback_state s; zx_lod_fallback_init(&s);
    int use=0;
    // below thresholds -> remain off
    use = zx_lod_fallback_update(&p, 5, 1.0f, &s); assert(use==0);
    use = zx_lod_fallback_update(&p, 5, 1.0f, &s); assert(use==0);
    // cross tiles threshold -> need enter_frames to engage
    use = zx_lod_fallback_update(&p, 20, 1.0f, &s); assert(use==0);
    use = zx_lod_fallback_update(&p, 20, 1.0f, &s); assert(use==1);
    // cool down -> need exit_frames to disengage
    use = zx_lod_fallback_update(&p, 1, 1.0f, &s); // 1
    use = zx_lod_fallback_update(&p, 1, 1.0f, &s); assert(use==0);
    return 0;
}


