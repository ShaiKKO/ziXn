/*!
 * \file lod_fallback_blend_test.cpp
 * \brief Blend window behavior for LOD fallback policy.
 */

#include "zx/zx_lod.h"
#include <cassert>

int main(){
    zx_lod_fallback_policy p{ 0, 0.0f, 1, 1, 3 };
    zx_lod_fallback_state s; zx_lod_fallback_init(&s);
    int use=0;
    use = zx_lod_fallback_update(&p, 100, 10.0f, &s); assert(use==1); // activated
    assert(s.blend_remaining == 2 || s.blend_remaining == 3-1);
    zx_lod_fallback_update(&p, 100, 10.0f, &s);
    zx_lod_fallback_update(&p, 100, 10.0f, &s);
    // after three frames, blend window expires
    assert(s.blend_remaining == 0);
    return 0;
}


