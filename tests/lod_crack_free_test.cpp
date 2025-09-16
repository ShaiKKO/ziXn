/*!
 * \file lod_crack_free_test.cpp
 * \brief Validate no cracks on borders after 2x downsample for smoothly varying tiles.
 */

#include "zx/zx_lod.h"
#include <vector>
#include <cassert>
#include <cmath>

static void make_smooth_tile(std::vector<float>& t, uint32_t W, uint32_t H, float xoffset)
{
    t.resize(W*H);
    for (uint32_t y=0;y<H;++y) for (uint32_t x=0;x<W;++x) {
        float xf = (x + xoffset) * 0.1f;
        float yf = y * 0.1f;
        t[y*W + x] = std::sin(xf) + std::cos(yf);
    }
}

int main(){
    const uint32_t W=8,H=8;
    std::vector<float> A, B; make_smooth_tile(A, W,H, 0.0f); make_smooth_tile(B, W,H, (float)W);
    // Ensure shared border continuity: right edge of A equals left edge of B
    for (uint32_t y=0;y<H;++y) B[y*W + 0] = A[y*W + (W-1)];

    std::vector<float> Ad(W/2 * H/2), Bd(W/2 * H/2);
    zx_lod_downsample_2x(A.data(), W,H,W, Ad.data(), W/2,H/2,W/2);
    zx_lod_downsample_2x(B.data(), W,H,W, Bd.data(), W/2,H/2,W/2);
    float border = zx_lod_border_consistency_check(Ad.data(), W/2,H/2,W/2, Bd.data(), W/2,H/2,W/2, 0);
    assert(border < 1e-3f);
    return 0;
}


