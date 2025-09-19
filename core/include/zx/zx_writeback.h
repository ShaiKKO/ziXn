/*!
\file zx_writeback.h
\brief Shared terrain writeback surfaces and CPU encode helper.
\author Colin Macritchie (Ripple Group, LLC)
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.

Notes:
- Displacement surface is R32F; masks are optional and may be NULL.
- Row pitch is in bytes; width/height in elements.
*/

#ifndef ZX_WRITEBACK_H
#define ZX_WRITEBACK_H

#include "zx_abi.h"
#include <stdint.h>

typedef struct zx_write_target
{
  /* Displacement: R32F, width x height, row_pitch in bytes */
  float* disp;
  uint32_t width;
  uint32_t height;
  uint32_t row_pitch_bytes;
  /* Optional masks (can be NULL) */
  uint16_t* comp_moist;  /* R16F packed via FP16 bits or fixed10; test path uses uint16 */
  uint8_t* contact_mask; /* R8U */
  uint8_t* debug_mask;   /* R8U */
} zx_write_target;

typedef struct zx_rect
{
  int32_t x;
  int32_t y;
  int32_t w;
  int32_t h;
} zx_rect;

/** \brief Clip-safe displacement write into target at rect (x,y,w,h).
 * @param src_patch Source patch (size src_w*src_h)
 * @param src_w Source width
 * @param src_h Source height
 * @param target Target surface (disp must not be NULL)
 * @param rect Destination rectangle (pixels)
 */
ZX_API void ZX_CALL zx_writeback_copy_displacement(const float* src_patch, uint32_t src_w,
                                                   uint32_t src_h, zx_write_target* target,
                                                   zx_rect rect);

#endif /* ZX_WRITEBACK_H */
