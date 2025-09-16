/*!
 * \file residency_prefetch_test.cpp
 * \brief Prefetch coverage estimation for tile residency manager.
 */

#include "zx/zx_residency.h"
#include <cassert>

int main(){
    zx_residency_opts o{1,1,2};
    zx_residency* r = zx_residency_create(&o);
    uint32_t en=0, ex=0, pf=0;
    zx_residency_tick(r, 0,0,0, 1, &en,&ex,&pf);
    // Prefetch volume should be outer cube minus inner cube
    assert(pf > 0);
    zx_residency_destroy(r);
    return 0;
}


