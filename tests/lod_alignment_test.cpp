/*!\
 * \file lod_alignment_test.cpp\
 * \brief Validate downsample/upsample correctness with non-16B and non-tight pitches.\
 */\

#include "zx/zx_lod.h"
#include <vector>
#include <cassert>
#include <cmath>

int main(){
    // Constant field downsample with extra pitch
    {
        const uint32_t W=8,H=8; const uint32_t SP=W+4; // non-tight
        std::vector<float> src(SP*H, 0.0f);
        for (uint32_t y=0;y<H;++y) for (uint32_t x=0;x<W;++x) src[y*SP+x] = 3.14f;
        std::vector<float> dst((W/2)*(H/2), -1.0f);
        zx_lod_downsample_2x(src.data(), W,H,SP, dst.data(), W/2,H/2,W/2);
        for (float v: dst) assert(std::fabs(v - 3.14f) < 1e-6f);
    }

    // Gradient field downsample with odd extra pitch
    {
        const uint32_t W=8,H=8; const uint32_t SP=W+5;
        std::vector<float> src(SP*H, 0.0f);
        for (uint32_t y=0;y<H;++y) for (uint32_t x=0;x<W;++x) src[y*SP+x] = (float)(x + 2*y);
        std::vector<float> dst((W/2)*(H/2), 0.0f);
        zx_lod_downsample_2x(src.data(), W,H,SP, dst.data(), W/2,H/2,W/2);
        // spot-check center
        assert(dst[(H/4)*(W/2) + (W/4)] > 0.0f);
    }

    // Upsample with non-tight destination pitch
    {
        const uint32_t SW=4,SH=4; const uint32_t DP=8+3; // non-tight dest pitch
        std::vector<float> src(SW*SH, 0.0f);
        for (uint32_t y=0;y<SH;++y) for (uint32_t x=0;x<SW;++x) src[y*SW+x] = (float)(x==y);
        std::vector<float> dst(DP*(SH*2), 0.0f);
        zx_lod_upsample_2x(src.data(), SW,SH,SW, dst.data(), SW*2,SH*2, DP);
        // sum should be positive
        double sum=0; for (uint32_t y=0;y<SH*2;++y) for (uint32_t x=0;x<SW*2;++x) sum += dst[y*DP+x];
        assert(sum > 0.0);
    }

    return 0;
}


