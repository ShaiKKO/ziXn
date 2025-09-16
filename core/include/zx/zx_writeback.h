/*!
\file zx_writeback.h
\brief Shared terrain writeback surfaces and CPU encode helper.
\author Colin Macritchie (Ripple Group, LLC)
*/

#ifndef ZX_WRITEBACK_H
#define ZX_WRITEBACK_H

#include <stdint.h>
#include "zx_abi.h"

typedef struct zx_write_target {
    /* Displacement: R32F, width x height, row_pitch in bytes */
    float*   disp;
    uint32_t width;
    uint32_t height;
    uint32_t row_pitch_bytes;
    /* Optional masks (can be NULL) */
    uint16_t* comp_moist;   /* R16F packed via FP16 bits or fixed10; test path uses uint16 */
    uint8_t*  contact_mask; /* R8U */
    uint8_t*  debug_mask;   /* R8U */
} zx_write_target;

typedef struct zx_rect {
    int32_t x;
    int32_t y;
    int32_t w;
    int32_t h;
} zx_rect;

/* Clip-safe displacement write: copies src_patch[w*h] into target at (x,y). */
ZX_API void ZX_CALL zx_writeback_copy_displacement(
    const float* src_patch, uint32_t src_w, uint32_t src_h,
    zx_write_target* target, zx_rect rect);

#endif /* ZX_WRITEBACK_H */


