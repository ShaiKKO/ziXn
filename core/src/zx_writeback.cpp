/*!
 * \file zx_writeback.cpp
 * \brief CPU-side writeback helper implementing clip-safe displacement copy.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_writeback.h"
#include <algorithm>
#include <cstring>

/**
 * @brief Copy a source displacement patch into a target displacement surface with clipping.
 *
 * Copies a sw-by-sh float patch pointed to by src into the target surface t at pixel
 * rectangle r, clipping both source and destination to t's bounds so all writes stay
 * in-range. Rows are copied with contiguous memcpy operations; rows of the source that
 * fall outside r or outside the source extent are skipped.
 *
 * @param src Pointer to the source patch (sw * sh floats). Must not be NULL.
 * @param sw Width of the source patch in floats.
 * @param sh Height of the source patch in rows.
 * @param t Target write surface; its disp pointer and dimensions are used and must be valid.
 * @param r Destination rectangle in target pixel coordinates (x,y,w,h).
 *
 * @note No return value. The function returns immediately without writes if inputs are
 * invalid (NULL pointers, zero dimensions, or if the clipped rectangle has non-positive size).
 */
void zx_writeback_copy_displacement(const float* src, uint32_t sw, uint32_t sh, zx_write_target* t,
                                    zx_rect r)
{
  if (src == nullptr || t == nullptr || t->disp == nullptr || sw == 0 || sh == 0 || t->width == 0 ||
      t->height == 0)
  {
    return;
  }
  auto x0 = r.x;
  auto y0 = r.y;
  auto w  = r.w;
  auto h  = r.h;
  if (w <= 0 || h <= 0)
  {
    return;
  }
  int32_t x1  = x0 + w;
  int32_t y1  = y0 + h;
  int32_t cx0 = x0 < 0 ? 0 : x0;
  int32_t cy0 = y0 < 0 ? 0 : y0;
  int32_t cx1 = x1 > (int32_t) t->width ? (int32_t) t->width : x1;
  int32_t cy1 = y1 > (int32_t) t->height ? (int32_t) t->height : y1;
  if (cx1 <= cx0 || cy1 <= cy0)
  {
    return;
  }
  const uint32_t dst_pitch = t->row_pitch_bytes / sizeof(float);
  for (int32_t y = cy0; y < cy1; ++y)
  {
    int32_t sy = y - y0;
    if (static_cast<unsigned>(sy) >= sh)
    {
      continue;
    }
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    float* dst_row = t->disp + (static_cast<size_t>(y) * static_cast<size_t>(dst_pitch));
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const float* src_row = src + (static_cast<size_t>(sy) * static_cast<size_t>(sw));
    int32_t sx0          = cx0 - x0;
    sx0                  = std::max(sx0, 0);
    if (static_cast<unsigned>(sx0) >= sw)
    {
      continue;
    }
    int32_t sx1 = sx0 + (cx1 - cx0);
    if (static_cast<unsigned>(sx1) > sw)
    {
      sx1 = static_cast<int32_t>(sw);
    }
    auto n = static_cast<uint32_t>(sx1 - sx0);
    std::memcpy(&dst_row[cx0], &src_row[sx0], n * sizeof(float));
  }
}
