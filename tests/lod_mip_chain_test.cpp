/*!
 * \file lod_mip_chain_test.cpp
 * \brief Validate multi-level mip chain consistency (downsample then upsample).
 */

#include "zx/zx_lod.h"
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>

int main(){
    const uint32_t W=16,H=16; std::vector<float> base(W*H);
    for (uint32_t y=0;y<H;++y) for (uint32_t x=0;x<W;++x) base[y*W+x] = std::sin(0.1f*(x+y));

    std::vector<float> m1((W/2)*(H/2)), m2((W/4)*(H/4));
    zx_lod_downsample_2x(base.data(), W,H,W, m1.data(), W/2,H/2,W/2);
    zx_lod_downsample_2x(m1.data(), W/2,H/2,W/2, m2.data(), W/4,H/4,W/4);

    // Reconstruct to full res and compare
    std::vector<float> up1(W*H), up2(W*H);
    std::vector<float> tmp((W/2)*(H/2));
    zx_lod_upsample_2x(m2.data(), W/4,H/4,W/4, tmp.data(), W/2,H/2,W/2);
    zx_lod_upsample_2x(tmp.data(), W/2,H/2,W/2, up2.data(), W,H,W);
    zx_lod_upsample_2x(m1.data(), W/2,H/2,W/2, up1.data(), W,H,W);

    float maxd1=0, maxd2=0;
    for (uint32_t i=0;i<W*H;++i) { maxd1 = std::max(maxd1, std::fabs(base[i]-up1[i])); }
    for (uint32_t i=0;i<W*H;++i) { maxd2 = std::max(maxd2, std::fabs(base[i]-up2[i])); }
    // The error should not grow unbounded; coarse level error should be larger than single level
    assert(maxd1 < 3.5f);
    assert(maxd2 < 7.0f);
    assert(maxd2 >= maxd1);
    return 0;
}


