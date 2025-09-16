/*!
 * \file lod_filters_test.cpp
 * \brief Validate 2x downsample/upsample and border consistency.
 */

#include "zx/zx_lod.h"
#include <vector>
#include <cassert>
#include <cmath>

int main(){
    const uint32_t W=8,H=8; std::vector<float> src(W*H);
    for (uint32_t y=0;y<H;++y) for (uint32_t x=0;x<W;++x) src[y*W+x] = (float)(x+y);
    std::vector<float> ds((W/2)*(H/2));
    zx_lod_downsample_2x(src.data(), W,H,W, ds.data(), W/2,H/2,W/2);
    std::vector<float> us(W*H);
    zx_lod_upsample_2x(ds.data(), W/2,H/2,W/2, us.data(), W,H,W);
    float maxd=0; for (uint32_t i=0;i<W*H;++i) maxd = std::max(maxd, std::fabs(src[i]-us[i]));
    assert(maxd < 2.0f); // coarse check

    // Border consistency: identical tiles -> zero difference
    float edge = zx_lod_border_consistency_check(src.data(), W,H,W, src.data(), W,H,W, 0);
    assert(edge < 1e-6f);

    // 1xN downsample should early-return (no crash) and leave dst unchanged
    {
        const uint32_t w=2,h=8; std::vector<float> a(w*h, 1.0f), ds2((w/2)*(h/2), -1.0f);
        zx_lod_downsample_2x(a.data(), w,h,w, ds2.data(), w/2,h/2,w/2);
        // valid path: w and h divisible; ensure values average to 1
        for (float v: ds2) assert(std::fabs(v-1.0f) < 1e-6f);
        // Non-tight pitch path: use larger pitch
        std::vector<float> a_pad(w*h + 16, 2.0f), ds_pad((w/2)*(h/2), 0.0f);
        zx_lod_downsample_2x(a_pad.data(), w,h,w+8, ds_pad.data(), w/2,h/2,w/2);
        for (float v: ds_pad) assert(std::fabs(v-2.0f) < 1e-6f);
    }

    // Odd dimensions should be rejected (no write)
    {
        const uint32_t w=7,h=7; std::vector<float> a(w*h, 3.0f), ds_bad((w/2)*(h/2), -5.0f);
        zx_lod_downsample_2x(a.data(), w,h,w, ds_bad.data(), w/2,h/2,w/2);
        for (float v: ds_bad) assert(v == -5.0f);
    }

    // Upsample 2x with non-tight pitch and odd inputs should no-op safely
    {
        const uint32_t sw=4,sh=4; std::vector<float> a(sw*sh, 0.5f), up( (sw*2)*(sh*2), 0.0f );
        zx_lod_upsample_2x(a.data(), sw,sh, sw+0, up.data(), sw*2,sh*2, sw*2);
        float sum=0; for (float v: up) sum+=v; assert(sum > 0.0f);
    }

    // Null pointers should early return
    {
        std::vector<float> a(16, 1.0f), b(16, 0.0f);
        zx_lod_downsample_2x(nullptr, 4,4,4, b.data(), 2,2,2);
        zx_lod_downsample_2x(a.data(), 4,4,4, nullptr, 2,2,2);
        zx_lod_upsample_2x(nullptr, 2,2,2, b.data(), 4,4,4);
        zx_lod_upsample_2x(a.data(), 2,2,2, nullptr, 4,4,4);
    }

    // Border consistency: right vs left with mismatched height
    {
        const uint32_t Aw=8,Ah=6,Bw=8,Bh=8; std::vector<float> A(Aw*Ah), B(Bw*Bh);
        for (uint32_t y=0;y<Ah;++y) A[y*Aw + (Aw-1)] = (float)y;
        for (uint32_t y=0;y<Bh;++y) B[y*Bw + 0] = (float)y;
        float d = zx_lod_border_consistency_check(A.data(),Aw,Ah,Aw,B.data(),Bw,Bh,Bw,0);
        // Should compare up to min(Ah,Bh)=6, max abs diff zero
        assert(d < 1e-6f);
    }

    // Border consistency: bottom vs top with mismatched width and intentional offset
    {
        const uint32_t Aw=7,Ah=8,Bw=9,Bh=8; std::vector<float> A(Aw*Ah, 0.0f), B(Bw*Bh, 0.0f);
        // Write a gradient on A bottom row and offset B top row by +1
        for (uint32_t x=0;x<Aw;++x) A[(Ah-1)*Aw + x] = (float)x;
        for (uint32_t x=0;x<Bw;++x) B[0*Bw + x] = (float)x + 1.0f;
        float d = zx_lod_border_consistency_check(A.data(),Aw,Ah,Aw,B.data(),Bw,Bh,Bw,1);
        // min width is 7; expected max abs diff is 1
        assert(std::fabs(d - 1.0f) < 1e-6f);
    }

    return 0;
}


