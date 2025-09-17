/*!
 * \file residency_pinning_test.cpp
 * \brief Ensure pinned regions remain active regardless of hysteresis.
 */

#include "zx/zx_residency.h"
#include <cassert>

int main(){
    zx_residency_opts o{1,2,1};
    zx_residency* r = zx_residency_create(&o);
    zx_residency_pin_box(r, -1,-1,-1, 1,1,1);
    for (int t=0;t<5;++t){ zx_residency_tick(r, 100,100,100, 0, nullptr,nullptr,nullptr); }
    // All pinned tiles must be active (3*3*3 = 27)
    assert(zx_residency_get_active_count(r) == 27u);
    zx_residency_unpin_all(r);
    for (int t=0;t<5;++t){ zx_residency_tick(r, 100,100,100, 0, nullptr,nullptr,nullptr); }
    // After unpinning and with zero radius, all should decay to zero
    assert(zx_residency_get_active_count(r) == 0u);
    zx_residency_destroy(r);
    return 0;
}


