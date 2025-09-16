/*!
 * \file zx_writeback.cpp
 * \brief CPU-side writeback helper implementing clip-safe displacement copy.
 */

#include "zx/zx_writeback.h"
#include <string.h>

void zx_writeback_copy_displacement(const float* src, uint32_t sw, uint32_t sh, zx_write_target* t, zx_rect r)
{
    if (!src || !t || !t->disp || sw==0 || sh==0 || t->width==0 || t->height==0) return;
    int32_t x0 = r.x; int32_t y0 = r.y; int32_t w = (int32_t)r.w; int32_t h = (int32_t)r.h;
    if (w <= 0 || h <= 0) return;
    int32_t x1 = x0 + w; int32_t y1 = y0 + h;
    int32_t cx0 = x0 < 0 ? 0 : x0;
    int32_t cy0 = y0 < 0 ? 0 : y0;
    int32_t cx1 = x1 > (int32_t)t->width ? (int32_t)t->width : x1;
    int32_t cy1 = y1 > (int32_t)t->height ? (int32_t)t->height : y1;
    if (cx1 <= cx0 || cy1 <= cy0) return;
    const uint32_t dst_pitch = t->row_pitch_bytes / sizeof(float);
    for (int32_t y = cy0; y < cy1; ++y) {
        int32_t sy = y - y0; if ((unsigned)sy >= sh) continue;
        float* dst_row = t->disp + y*dst_pitch;
        const float* src_row = src + sy*sw;
        int32_t sx0 = cx0 - x0; if (sx0 < 0) sx0 = 0; if ((unsigned)sx0 >= sw) continue;
        int32_t sx1 = sx0 + (cx1 - cx0); if ((unsigned)sx1 > sw) sx1 = (int32_t)sw;
        const uint32_t n = (uint32_t)(sx1 - sx0);
        memcpy(dst_row + cx0, src_row + sx0, n * sizeof(float));
    }
}


