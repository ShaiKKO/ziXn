/*!
 * \file counters_smoke.cpp
 * \brief Validate zx_counters size/version and retrieval.
 */

#include "zx/zx_abi.h"
#include <cassert>

int main(){
    zx_procs P{}; auto st = zxGetProcTable(ZX_ABI_VERSION, &P); assert(st == ZX_OK);
    zx_counters ctr{}; ctr.size = sizeof(ctr);
    uint32_t count = 0; st = P.get_counters(0, &ctr, &count); assert(st == ZX_OK);
    assert(count == 1);
    assert(ctr.version != 0);
    return 0;
}


