/*!
\file zx_lod.h
\brief Presentation LOD filters for displacement/masks with mip-consistent rules.
*/

#ifndef ZX_LOD_H
#define ZX_LOD_H

#include <stdint.h>
#include "zx_abi.h"

/* Downsample by 2× using average filter for displacement. Src/dst are row-major arrays.
 * src_pitch/dst_pitch are in elements (floats), not bytes.
 */
ZX_API void ZX_CALL zx_lod_downsample_2x(
    const float* src, uint32_t src_w, uint32_t src_h, uint32_t src_pitch,
    float* dst, uint32_t dst_w, uint32_t dst_h, uint32_t dst_pitch);

/* Upsample by 2× using bilinear filter for validation and crack checking. */
ZX_API void ZX_CALL zx_lod_upsample_2x(
    const float* src, uint32_t src_w, uint32_t src_h, uint32_t src_pitch,
    float* dst, uint32_t dst_w, uint32_t dst_h, uint32_t dst_pitch);

/* Check border consistency between two adjacent tiles after downsample.
 * side = 0: compare right edge of A with left edge of B.
 * side = 1: compare bottom edge of A with top edge of B.
 * Returns max absolute difference across compared edge.
 */
ZX_API float ZX_CALL zx_lod_border_consistency_check(
    const float* A, uint32_t Aw, uint32_t Ah, uint32_t Apitch,
    const float* B, uint32_t Bw, uint32_t Bh, uint32_t Bpitch,
    int side);

/* Fallback policy/state for present-LOD under pressure signals */
typedef struct zx_lod_fallback_policy {
    uint32_t active_tiles_max;
    float    step_ms_max;
    uint32_t enter_frames;
    uint32_t exit_frames;
    uint32_t blend_frames;
} zx_lod_fallback_policy;

typedef struct zx_lod_fallback_state {
    uint32_t active_frames;
    uint32_t inactive_frames;
    uint32_t blend_remaining;
    uint32_t activations;
} zx_lod_fallback_state;

ZX_API void ZX_CALL zx_lod_fallback_init(zx_lod_fallback_state* s);
ZX_API int  ZX_CALL zx_lod_fallback_update(const zx_lod_fallback_policy* p,
                                           uint32_t residency_active_tiles,
                                           float last_step_ms,
                                           zx_lod_fallback_state* s);
// returns 1 if coarse LOD should be used this frame

/* Global fallback configuration for integration/CLI wiring */
ZX_API void ZX_CALL zx_lod_set_enabled(int on);
ZX_API int  ZX_CALL zx_lod_is_enabled(void);
ZX_API void ZX_CALL zx_lod_set_default_policy(const zx_lod_fallback_policy* p);
ZX_API void ZX_CALL zx_lod_get_default_policy(zx_lod_fallback_policy* out);

/* SIMD override for LOD kernels: -1 auto, 0 scalar, 2 AVX2 (if supported) */
ZX_API void ZX_CALL zx_lod_set_simd_override(int mode);

#endif /* ZX_LOD_H */


