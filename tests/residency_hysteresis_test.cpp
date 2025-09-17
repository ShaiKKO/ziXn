/*!
 * \file residency_hysteresis_test.cpp
 * \brief Hysteresis correctness for tile residency manager.
 */

#include "zx/zx_residency.h"
#include <cassert>

int main(){
    zx_residency_opts o{2,2,1};
    zx_residency* r = zx_residency_create(&o);
    uint32_t en=0, ex=0, pf=0;
    // First tick: one enter after 2 hot frames -> not yet
    zx_residency_tick(r, 0,0,0, 0, &en,&ex,&pf); assert(en==0 && ex==0);
    zx_residency_tick(r, 0,0,0, 0, &en,&ex,&pf); assert(en>=1);
    // Move away: need 2 cold frames to exit
    zx_residency_tick(r, 10,0,0, 0, &en,&ex,&pf);
    zx_residency_tick(r, 10,0,0, 0, &en,&ex,&pf); assert(ex>=1);
    zx_residency_destroy(r);
    return 0;
}


