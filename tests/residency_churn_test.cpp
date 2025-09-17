/*!
 * \file residency_churn_test.cpp
 * \brief Verify churn stats reflect enters/exits per tick.
 */

#include "zx/zx_residency.h"
#include <cassert>

int main(){
    zx_residency_opts o{1,1,0};
    zx_residency* r = zx_residency_create(&o);
    uint32_t en=0, ex=0, pf=0;
    zx_residency_tick(r, 0,0,0, 1, &en,&ex,&pf); // activate a cube 3^3 = 27
    uint32_t cen=0, cex=0, churn=0; zx_residency_get_last_churn(r, &cen,&cex,&churn);
    assert(cen == en && cex == ex && churn == (en+ex));
    // Move center far away to cause exits
    zx_residency_tick(r, 100,100,100, 0, &en,&ex,&pf);
    zx_residency_get_last_churn(r, &cen,&cex,&churn);
    assert(cen == en && cex == ex && churn == (en+ex));
    zx_residency_destroy(r);
    return 0;
}


