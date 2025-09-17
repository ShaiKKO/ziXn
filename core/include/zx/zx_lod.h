/*!
\file zx_lod.h
\brief Presentation LOD filters for displacement/masks with mip-consistent rules.
\author Colin Macritchie (Ripple Group, LLC)
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.

Level-of-Detail (LOD) utilities used at presentation time:
- Downsample by 2× using average filter (mip chain construction)
- Upsample by 2× using bilinear filter (validation, crack checks)
- Edge consistency checks across tile borders
- Frame-to-frame fallback policy with hysteresis and blend window

Units and conventions:
- Arrays are row-major; pitches are in elements (floats), not bytes
- Times are milliseconds (ms); frame counts are integral (frames)
- Tolerances are expressed as absolute differences in float units

SIMD and platform notes (Windows-first):
- Runtime dispatch selects AVX2 or scalar implementations via zx_lod_set_simd_override()
- Memory alignment is not required; intrinsics use unaligned loads
*/

#ifndef ZX_LOD_H
#define ZX_LOD_H

#include <stdint.h>
#include "zx_abi.h"

/** \brief Downsample a displacement field by 2× using a 2×2 average filter.
 *
 * Src/dst are row-major arrays. Pitches are in elements (floats), not bytes.
 * @param src Source pointer (must not be NULL)
 * @param src_w Source width (must equal 2×dst_w)
 * @param src_h Source height (must equal 2×dst_h)
 * @param src_pitch Source row pitch (elements)
 * @param dst Destination pointer (must not be NULL)
 * @param dst_w Destination width
 * @param dst_h Destination height
 * @param dst_pitch Destination row pitch (elements)
 * @note No-op on invalid preconditions.
 */
ZX_API void ZX_CALL zx_lod_downsample_2x(
    const float* src, uint32_t src_w, uint32_t src_h, uint32_t src_pitch,
    float* dst, uint32_t dst_w, uint32_t dst_h, uint32_t dst_pitch);

/** \brief Upsample a displacement field by 2× using separable bilinear filter.
 * @param src Source pointer (must not be NULL)
 * @param src_w Source width
 * @param src_h Source height
 * @param src_pitch Source row pitch (elements)
 * @param dst Destination pointer (must not be NULL)
 * @param dst_w Destination width (must equal 2×src_w)
 * @param dst_h Destination height (must equal 2×src_h)
 * @param dst_pitch Destination row pitch (elements)
 * @note No-op on invalid preconditions.
 */
ZX_API void ZX_CALL zx_lod_upsample_2x(
    const float* src, uint32_t src_w, uint32_t src_h, uint32_t src_pitch,
    float* dst, uint32_t dst_w, uint32_t dst_h, uint32_t dst_pitch);

/** \brief Check border consistency between two adjacent tiles.
 * Compares edges after downsample to detect cracks.
 * @param A Tile A pointer (must not be NULL)
 * @param Aw Width of A; @param Ah Height of A; @param Apitch Row pitch of A (elements)
 * @param B Tile B pointer (must not be NULL)
 * @param Bw Width of B; @param Bh Height of B; @param Bpitch Row pitch of B (elements)
 * @param side 0=compare A right to B left; 1=compare A bottom to B top
 * @return Max absolute difference across the compared edge
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

/** \brief Initialize fallback state to defaults. */
ZX_API void ZX_CALL zx_lod_fallback_init(zx_lod_fallback_state* s);
/** \brief Update fallback state based on current pressure signals.
 * Pressure if (residency_active_tiles > active_tiles_max) || (last_step_ms > step_ms_max).
 * @param p Policy thresholds (must not be NULL)
 * @param residency_active_tiles Active tile count
 * @param last_step_ms Last step duration in ms
 * @param s In/out state (must not be NULL)
 * @return 1 if coarse LOD should be used this frame, else 0
 */
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


