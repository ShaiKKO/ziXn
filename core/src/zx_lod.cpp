/*!
 * \file zx_lod.cpp
 * \brief Presentation LOD filters implementation (scalar + AVX2 dispatch).
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 *
 * Fallback state machine (hysteresis + blend):
 *
 *  State vars: active_frames, inactive_frames, blend_remaining, activations
 *
 *  Diagram (frames):
 *    [pressure=1]  active_frames++  -----------------------------\
 *                    if first activation: blend_remaining=blend  |
 *                    if active_frames >= enter -> use_coarse=1   |
 *    [pressure=0]  inactive_frames++                            |
 *                    if inactive_frames >= exit -> use_coarse=0  |
 *                    else retain use_coarse                      |
 *    Each frame: if blend_remaining>0 then blend_remaining--     |
 *  --------------------------------------------------------------/
 *  Pressure condition: (active_tiles > active_tiles_max) || (last_step_ms > step_ms_max)
 */

#include "zx/zx_lod.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <mutex>
#include <vector>
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

#if defined(__GNUC__) || defined(__clang__)
#define ZX_TARGET_AVX2 __attribute__((target("avx2")))
#else
#define ZX_TARGET_AVX2
#endif

// Internal SIMD constants
namespace
{
  constexpr int k_avx_lanes            = 8;
  constexpr int k_avx_align            = 32;
  constexpr uint32_t k_prefetch_offset = 64U;
  constexpr uint32_t k_block8          = 8U;
  constexpr uint32_t k_block16         = 16U;
#ifdef ZX_LOD_ENABLE_AVX2_BLOCK32
  constexpr uint32_t k_block32 = 32U;
#endif
  constexpr float k_half                           = 0.5F;
  constexpr float k_quarter                        = 0.25F;
  constexpr float k_neg_quarter                    = -0.25F;
  constexpr int k_gather_scale                     = static_cast<int>(sizeof(float));
  constexpr int k_cpuid_leaf_avx2                  = 7;       // CPUID leaf 7 (structured extended)
  constexpr int k_cpuid_ebx_avx2_bit               = 5;       // EBX bit 5 => AVX2
  constexpr int k_cpuid_ecx_osxsave_bit            = 27;      // ECX bit 27 => OSXSAVE
  constexpr unsigned long long k_xcr0_xmm_ymm_mask = 0x6ULL;  // XMM|YMM enabled in XCR0
  constexpr uint32_t k_default_active_tiles_max    = 1000000U;
  constexpr float k_default_step_ms_max            = 1.0e9F;
  constexpr uint32_t k_default_enter_frames        = 2U;
  constexpr uint32_t k_default_exit_frames         = 2U;
  constexpr uint32_t k_default_blend_frames        = 3U;
  constexpr uint32_t k_nt_min_row_width            = 128U;  // elements; enable NT on wide rows
}  // namespace

/* SIMD override: -1 auto, 0 scalar, 2 AVX2 */
static int g_simd_override = -1;  // -1 auto, 0 scalar, 2 AVX2
/* Store policy: -1 auto, 0 disabled, 1 enabled (when safe) */
static int g_store_policy = -1;

#if defined(__x86_64__) || defined(_M_X64)
static inline int cpu_supports_avx2(); /* forward */
#endif
/* SIMD-dispatched 2x downsample and upsample kernels */
typedef void (*pp_down2_fn)(const float*, uint32_t, uint32_t, uint32_t, float*, uint32_t, uint32_t,
                            uint32_t);
typedef void (*pp_up2_fn)(const float*, uint32_t, uint32_t, uint32_t, float*, uint32_t, uint32_t,
                          uint32_t);

static void down2_scalar(const float* src, uint32_t sw, uint32_t sh, uint32_t sp, float* dst,
                         uint32_t dw, uint32_t dh, uint32_t dp)
{
  if ((src == nullptr) || (dst == nullptr))
  {
    return;
  }
  if ((sw < 2U) || (sh < 2U) || (dw == 0U) || (dh == 0U) || (sp < sw) || (dp < dw))
  {
    return;
  }
  if (((dw * 2U) != sw) || ((dh * 2U) != sh))
  {
    return;
  }
  constexpr float k_quarter = 0.25F;
  for (uint32_t y = 0; y < dh; ++y)
  {
    // Scalar row accumulation; optional streaming store if policy enables and row wide
    const bool use_nt =
        (g_store_policy == 1) || ((g_store_policy == -1) && (dw >= k_nt_min_row_width));
    for (uint32_t x = 0; x < dw; ++x)
    {
      const float a = src[((2U * y + 0U) * sp) + (2U * x + 0U)];
      const float b = src[((2U * y + 0U) * sp) + (2U * x + 1U)];
      const float c = src[((2U * y + 1U) * sp) + (2U * x + 0U)];
      const float d = src[((2U * y + 1U) * sp) + (2U * x + 1U)];
      const float v = k_quarter * (a + b + c + d);
#if defined(__x86_64__) || defined(_M_X64)
      if (use_nt)
      {
        _mm_stream_ss(dst + (static_cast<size_t>(y) * dp) + x, _mm_set_ss(v));
      }
      else
#endif
      {
        dst[(static_cast<size_t>(y) * dp) + x] = v;
      }
    }
  }
}

#if defined(__x86_64__) || defined(_M_X64)
static ZX_TARGET_AVX2 void down2_avx2(const float* src, uint32_t sw, uint32_t sh, uint32_t sp,
                                      float* dst, uint32_t dw, uint32_t dh, uint32_t dp)
{
  if ((src == nullptr) || (dst == nullptr) || (sw < 2U) || (sh < 2U) || ((dw * 2U) != sw) ||
      ((dh * 2U) != sh))
  {
    return;
  }
  const uint32_t vec_w = (dw / k_block8) * k_block8;
  const __m256 quarter = _mm256_set1_ps(k_quarter);
  for (uint32_t y = 0; y < dh; ++y)
  {
    const float* row0 = src + static_cast<size_t>((2U * y + 0U) * sp);
    const float* row1 = src + static_cast<size_t>((2U * y + 1U) * sp);
    uint32_t x        = 0;
    const bool use_nt =
        (g_store_policy == 1) || ((g_store_policy == -1) && (dw >= k_nt_min_row_width));
    for (; x < vec_w; x += k_block8)
    {
      // Load 16 contiguous samples from each row
      const float* r0a_ptr = row0 + static_cast<size_t>((2U * x) + 0U);
      const float* r0b_ptr = row0 + static_cast<size_t>((2U * x) + 8U);
      const float* r1a_ptr = row1 + static_cast<size_t>((2U * x) + 0U);
      const float* r1b_ptr = row1 + static_cast<size_t>((2U * x) + 8U);
      __m256 r0a           = _mm256_loadu_ps(r0a_ptr);
      __m256 r0b           = _mm256_loadu_ps(r0b_ptr);
      __m256 r1a           = _mm256_loadu_ps(r1a_ptr);
      __m256 r1b           = _mm256_loadu_ps(r1b_ptr);
      // Horizontal add adjacent pairs within lanes
      __m256 s0            = _mm256_hadd_ps(r0a, r0b);
      __m256 s1            = _mm256_hadd_ps(r1a, r1b);
      __m256 sum           = _mm256_add_ps(s0, s1);
      __m256 out           = _mm256_mul_ps(sum, quarter);
      const size_t dst_row = static_cast<size_t>(y) * static_cast<size_t>(dp);
      if (use_nt)
      {
        _mm256_stream_ps(dst + dst_row + x, out);
      }
      else
      {
        _mm256_storeu_ps(dst + dst_row + x, out);
      }
    }
    for (; x < dw; ++x)
    {
      const float a = row0[(2U * x) + 0U];
      const float b = row0[(2U * x) + 1U];
      const float c = row1[(2U * x) + 0U];
      const float d = row1[(2U * x) + 1U];
      const float v = k_quarter * (a + b + c + d);
      if (use_nt)
      {
        _mm_stream_ss(dst + (static_cast<size_t>(y) * dp) + x, _mm_set_ss(v));
      }
      else
      {
        dst[(static_cast<size_t>(y) * dp) + x] = v;
      }
    }
  }
}
#endif

static pp_down2_fn g_down2_impl = down2_scalar;
static std::once_flag g_down2_once;
static void init_down2_impl()
{
  if (g_simd_override == 0)
  {
    g_down2_impl = down2_scalar;
    return;
  }
#if defined(__x86_64__) || defined(_M_X64)
  if (g_simd_override == 2)
  {
    g_down2_impl = (cpu_supports_avx2() != 0) ? down2_avx2 : down2_scalar;  // explicit int->bool
    return;
  }
  if (cpu_supports_avx2() != 0)
  {
    g_down2_impl = down2_avx2;
    return;
  }
#endif
  g_down2_impl = down2_scalar;
}

static void up2_scalar(const float* src, uint32_t sw, uint32_t sh, uint32_t sp, float* dst,
                       uint32_t dw, uint32_t dh, uint32_t dp)
{
  if ((src == nullptr) || (dst == nullptr))
  {
    return;
  }
  if ((sw == 0U) || (sh == 0U) || (dw == 0U) || (dh == 0U) || (sp < sw) || (dp < dw))
  {
    return;
  }
  if ((dw != (sw * 2U)) || (dh != (sh * 2U)))
  {
    return;
  }
  for (uint32_t y = 0; y < dh; ++y)
  {
    const float fy = ((static_cast<float>(y) + k_half) * k_half) - k_half;
    int y0         = 0;
    int y1         = 0;
    float ty       = 0.0F;
    {
      const int y_floor = static_cast<int>(std::floor(static_cast<double>(fy)));
      y0                = std::clamp(y_floor, 0, static_cast<int>(sh) - 1);
      y1                = std::min(y0 + 1, static_cast<int>(sh) - 1);
      ty                = fy - static_cast<float>(y0);
    }
    for (uint32_t x = 0; x < dw; ++x)
    {
      const float fx = ((static_cast<float>(x) + k_half) * k_half) - k_half;
      int x0         = static_cast<int>(std::floor(static_cast<double>(fx)));
      x0             = std::clamp(x0, 0, static_cast<int>(sw) - 1);
      int x1         = std::min(static_cast<int>(sw) - 1, x0 + 1);
      const float tx = fx - static_cast<float>(x0);
      const float a  = src[(y0 * sp) + x0];
      const float b  = src[(y0 * sp) + x1];
      const float c  = src[(y1 * sp) + x0];
      const float d  = src[(y1 * sp) + x1];
      const float ab = a + (tx * (b - a));
      const float cd = c + (tx * (d - c));
      dst[(static_cast<size_t>(y) * dp) + x] = ab + (ty * (cd - ab));
    }
  }
}

#if defined(__x86_64__) || defined(_M_X64)
/**
 * @brief AVX2-accelerated 2x bilinear upsampling of a single-channel float image.
 *
 * Performs a 2x bilinear upsample from src (sw x sh, stride sp) into dst (dw x dh, stride dp)
 * using AVX2 vector gathers and per-column precomputed indices/fractions. If inputs are
 * invalid (null pointers or dw != sw*2 or dh != sh*2) the function falls back to the scalar
 * implementation (up2_scalar).
 *
 * The implementation:
 * - Precomputes x index pairs and fractional x weights for all vector lanes in a row.
 * - Iterates rows, computes y index/weight per row, then evaluates bilinear interpolation
 *   across 8/16-column AVX2 blocks using _mm256_i32gather_ps and fused arithmetic.
 * - Handles remaining vector-width tails with 8-wide AVX2 loops and a final scalar tail.
 * - Optionally supports a 32-column block repetition when ZX_LOD_ENABLE_AVX2_BLOCK32 is defined.
 *
 * Notes:
 * - Targeted for AVX2; compiled with ZX_TARGET_AVX2 attribute.
 * - Writes dw*dh samples into dst; expects single-channel float data.
 * - No return value; side effect is writing the upsampled image into dst.
 */
static ZX_TARGET_AVX2 void up2_avx2(const float* src, uint32_t sw, uint32_t sh, uint32_t sp,
                                    float* dst, uint32_t dw, uint32_t dh, uint32_t dp)
{
  if ((src == nullptr) || (dst == nullptr) || (dw != sw * 2U) || (dh != sh * 2U))
  {
    up2_scalar(src, sw, sh, sp, dst, dw, dh, dp);
    return;
  }
  const uint32_t vec_w = (dw / k_block8) * k_block8;
  const uint32_t nv    = (vec_w / k_block8);
  // Precompute x indices and fractional tx for all vector lanes across the row
  std::vector<int> pre_ix0;
  std::vector<int> pre_ix1;
  std::vector<float> pre_tx;
  pre_ix0.resize(static_cast<size_t>(nv) * 8U);
  pre_ix1.resize(static_cast<size_t>(nv) * 8U);
  pre_tx.resize(static_cast<size_t>(nv) * 8U);
  __m256 half        = _mm256_set1_ps(k_half);
  __m256 neg_quarter = _mm256_set1_ps(k_neg_quarter);
  __m256i zero       = _mm256_set1_epi32(0);
  __m256i xmax       = _mm256_set1_epi32(static_cast<int>(sw) - 1);
  for (uint32_t x = 0, i = 0; x < vec_w; x += k_block8, ++i)
  {
    __m256i vx =
        _mm256_setr_epi32(static_cast<int>(x) + 0, static_cast<int>(x) + 1, static_cast<int>(x) + 2,
                          static_cast<int>(x) + 3, static_cast<int>(x) + 4, static_cast<int>(x) + 5,
                          static_cast<int>(x) + 6, static_cast<int>(x) + 7);
    __m256 vxf    = _mm256_cvtepi32_ps(vx);
    __m256 vfx    = _mm256_add_ps(_mm256_mul_ps(vxf, half), neg_quarter);
    __m256 vfloor = _mm256_floor_ps(vfx);
    __m256i vx0   = _mm256_cvttps_epi32(vfloor);
    vx0           = _mm256_max_epi32(zero, _mm256_min_epi32(vx0, xmax));
    __m256 vtx    = _mm256_sub_ps(vfx, vfloor);
    __m256i vx1   = _mm256_min_epi32(_mm256_add_epi32(vx0, _mm256_set1_epi32(1)), xmax);
    // store into scalar arrays to avoid template attribute warnings on __m256* vector elements
    size_t base = static_cast<size_t>(i) * static_cast<size_t>(k_avx_lanes);
    alignas(k_avx_align) std::array<int, k_avx_lanes> x0buf{};
    alignas(k_avx_align) std::array<int, k_avx_lanes> x1buf{};
    alignas(k_avx_align) std::array<float, k_avx_lanes> txbuf{};
    _mm256_store_si256(reinterpret_cast<__m256i*>(x0buf.data()), vx0);
    _mm256_store_si256(reinterpret_cast<__m256i*>(x1buf.data()), vx1);
    _mm256_store_ps(txbuf.data(), vtx);
    for (int k = 0; k < k_avx_lanes; ++k)
    {
      pre_ix0[base + k] = x0buf[k];
      pre_ix1[base + k] = x1buf[k];
      pre_tx[base + k]  = txbuf[k];
    }
  }
  for (uint32_t y = 0; y < dh; ++y)
  {
    auto fy = static_cast<float>(
        ((static_cast<double>(y) + static_cast<double>(k_half)) * static_cast<double>(k_half)) -
        static_cast<double>(k_half));
    float ty_scalar = 0.0F;
    __m256 vty;
    int y0i = 0;
    int y1i = 0;
    {
      const int y_floor = static_cast<int>(std::floor(static_cast<double>(fy)));
      y0i               = std::clamp(y_floor, 0, static_cast<int>(sh) - 1);
      y1i               = std::min(y0i + 1, static_cast<int>(sh) - 1);
      ty_scalar         = fy - static_cast<float>(y0i);
      vty               = _mm256_set1_ps(ty_scalar);
    }
    const float* row0 = src + (static_cast<size_t>(y0i) * sp);
    const float* row1 = src + (static_cast<size_t>(y1i) * sp);
    // process 16 columns per iteration where possible
    uint32_t x = 0;
    uint32_t i = 0;
    for (; (x + k_block16) <= vec_w; x += k_block16, i += 2U)
    {
      const int* px0a   = &pre_ix0[static_cast<size_t>(i) * static_cast<size_t>(k_avx_lanes)];
      const int* px1a   = &pre_ix1[static_cast<size_t>(i) * static_cast<size_t>(k_avx_lanes)];
      const float* ptxa = &pre_tx[static_cast<size_t>(i) * static_cast<size_t>(k_avx_lanes)];
      const int* px0b   = &pre_ix0[static_cast<size_t>(i + 1) * static_cast<size_t>(k_avx_lanes)];
      const int* px1b   = &pre_ix1[static_cast<size_t>(i + 1) * static_cast<size_t>(k_avx_lanes)];
      const float* ptxb = &pre_tx[static_cast<size_t>(i + 1) * static_cast<size_t>(k_avx_lanes)];
      __m256i vx0a      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px0a));
      __m256i vx1a      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px1a));
      __m256 vtxa       = _mm256_loadu_ps(ptxa);
      __m256i vx0b      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px0b));
      __m256i vx1b      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px1b));
      __m256 vtxb       = _mm256_loadu_ps(ptxb);
      // prefetch further ahead
      _mm_prefetch(reinterpret_cast<const char*>(row0 + x + k_prefetch_offset), _MM_HINT_T0);
      _mm_prefetch(reinterpret_cast<const char*>(row1 + x + k_prefetch_offset), _MM_HINT_T0);
      // Portability: AVX2 intrinsics
      __m256 a0   = _mm256_i32gather_ps(row0, vx0a, k_gather_scale);
      __m256 b0   = _mm256_i32gather_ps(row0, vx1a, k_gather_scale);
      __m256 c0   = _mm256_i32gather_ps(row1, vx0a, k_gather_scale);
      __m256 d0   = _mm256_i32gather_ps(row1, vx1a, k_gather_scale);
      __m256 ab0  = _mm256_add_ps(_mm256_mul_ps(vtxa, _mm256_sub_ps(b0, a0)), a0);
      __m256 cd0  = _mm256_add_ps(_mm256_mul_ps(vtxa, _mm256_sub_ps(d0, c0)), c0);
      __m256 out0 = _mm256_add_ps(_mm256_mul_ps(vty, _mm256_sub_ps(cd0, ab0)), ab0);
      _mm256_storeu_ps(dst + (static_cast<size_t>(y) * dp) + x, out0);

      __m256 a1   = _mm256_i32gather_ps(row0, vx0b, k_gather_scale);
      __m256 b1   = _mm256_i32gather_ps(row0, vx1b, k_gather_scale);
      __m256 c1   = _mm256_i32gather_ps(row1, vx0b, k_gather_scale);
      __m256 d1   = _mm256_i32gather_ps(row1, vx1b, k_gather_scale);
      __m256 ab1  = _mm256_add_ps(_mm256_mul_ps(vtxb, _mm256_sub_ps(b1, a1)), a1);
      __m256 cd1  = _mm256_add_ps(_mm256_mul_ps(vtxb, _mm256_sub_ps(d1, c1)), c1);
      __m256 out1 = _mm256_add_ps(_mm256_mul_ps(vty, _mm256_sub_ps(cd1, ab1)), ab1);
      _mm256_storeu_ps(dst + (static_cast<size_t>(y) * dp) + x + k_block8, out1);
    }
#ifdef ZX_LOD_ENABLE_AVX2_BLOCK32
    // Optional: process 32 columns per iteration when enough width remains
    while ((x + k_block32) <= vec_w)
    {
      // two iterations of the 16-col body
      for (int rep = 0; rep < 2; ++rep)
      {
        const int* px0a   = &pre_ix0[static_cast<size_t>(i) * 8U];
        const int* px1a   = &pre_ix1[static_cast<size_t>(i) * 8U];
        const float* ptxa = &pre_tx[static_cast<size_t>(i) * 8U];
        const int* px0b   = &pre_ix0[static_cast<size_t>(i + 1) * 8U];
        const int* px1b   = &pre_ix1[static_cast<size_t>(i + 1) * 8U];
        const float* ptxb = &pre_tx[static_cast<size_t>(i + 1) * 8U];
        __m256i vx0a      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px0a));
        __m256i vx1a      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px1a));
        __m256 vtxa       = _mm256_loadu_ps(ptxa);
        __m256i vx0b      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px0b));
        __m256i vx1b      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px1b));
        __m256 vtxb       = _mm256_loadu_ps(ptxb);
        _mm_prefetch(reinterpret_cast<const char*>(row0 + x + k_prefetch_offset), _MM_HINT_T0);
        _mm_prefetch(reinterpret_cast<const char*>(row1 + x + k_prefetch_offset), _MM_HINT_T0);
        __m256 a0   = _mm256_i32gather_ps(row0, vx0a, k_gather_scale);
        __m256 b0   = _mm256_i32gather_ps(row0, vx1a, k_gather_scale);
        __m256 c0   = _mm256_i32gather_ps(row1, vx0a, k_gather_scale);
        __m256 d0   = _mm256_i32gather_ps(row1, vx1a, k_gather_scale);
        __m256 ab0  = _mm256_add_ps(_mm256_mul_ps(vtxa, _mm256_sub_ps(b0, a0)), a0);
        __m256 cd0  = _mm256_add_ps(_mm256_mul_ps(vtxa, _mm256_sub_ps(d0, c0)), c0);
        __m256 out0 = _mm256_add_ps(_mm256_mul_ps(vty, _mm256_sub_ps(cd0, ab0)), ab0);
        _mm256_storeu_ps(dst + (size_t) y * dp + x, out0);

        __m256 a1   = _mm256_i32gather_ps(row0, vx0b, 4);
        __m256 b1   = _mm256_i32gather_ps(row0, vx1b, 4);
        __m256 c1   = _mm256_i32gather_ps(row1, vx0b, 4);
        __m256 d1   = _mm256_i32gather_ps(row1, vx1b, 4);
        __m256 ab1  = _mm256_add_ps(_mm256_mul_ps(vtxb, _mm256_sub_ps(b1, a1)), a1);
        __m256 cd1  = _mm256_add_ps(_mm256_mul_ps(vtxb, _mm256_sub_ps(d1, c1)), c1);
        __m256 out1 = _mm256_add_ps(_mm256_mul_ps(vty, _mm256_sub_ps(cd1, ab1)), ab1);
        _mm256_storeu_ps(dst + (size_t) y * dp + x + k_block8, out1);

        x += k_block16;
        i += 2;
      }
    }
#endif
    for (; x < vec_w; x += k_block8, ++i)
    {
      const int* px0   = &pre_ix0[static_cast<size_t>(i) * static_cast<size_t>(k_avx_lanes)];
      const int* px1   = &pre_ix1[static_cast<size_t>(i) * static_cast<size_t>(k_avx_lanes)];
      const float* ptx = &pre_tx[static_cast<size_t>(i) * static_cast<size_t>(k_avx_lanes)];
      __m256i vx0      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px0));
      __m256i vx1      = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(px1));
      __m256 vtx       = _mm256_loadu_ps(ptx);
      // gathers
      __m256 a   = _mm256_i32gather_ps(row0, vx0, k_gather_scale);
      __m256 b   = _mm256_i32gather_ps(row0, vx1, k_gather_scale);
      __m256 c   = _mm256_i32gather_ps(row1, vx0, k_gather_scale);
      __m256 d   = _mm256_i32gather_ps(row1, vx1, k_gather_scale);
      __m256 ab  = _mm256_add_ps(_mm256_mul_ps(vtx, _mm256_sub_ps(b, a)), a);
      __m256 cd  = _mm256_add_ps(_mm256_mul_ps(vtx, _mm256_sub_ps(d, c)), c);
      __m256 out = _mm256_add_ps(_mm256_mul_ps(vty, _mm256_sub_ps(cd, ab)), ab);
      _mm256_storeu_ps(dst + (static_cast<size_t>(y) * dp) + x, out);
    }
    // tail
    for (; x < dw; ++x)
    {
      auto fx = static_cast<float>(
          ((static_cast<double>(x) + static_cast<double>(k_half)) * static_cast<double>(k_half)) -
          static_cast<double>(k_half));
      int x0   = static_cast<int>(std::floor(static_cast<double>(fx)));
      x0       = std::clamp(x0, 0, static_cast<int>(sw) - 1);
      int x1   = std::min(x0 + 1, static_cast<int>(sw) - 1);
      float tx = fx - static_cast<float>(x0);
      float a  = row0[x0];
      float b  = row0[x1];
      float c  = row1[x0];
      float d  = row1[x1];
      float ab = a + (tx * (b - a));
      float cd = c + (tx * (d - c));
      dst[(static_cast<size_t>(y) * dp) + x] = ab + ty_scalar * (cd - ab);
    }
  }
}
#endif

static pp_up2_fn g_up2_impl = up2_scalar;
static std::once_flag g_up2_once;
static void init_up2_impl()
{
  if (g_simd_override == 0)
  {
    g_up2_impl = up2_scalar;
    return;
  }
#if defined(__x86_64__) || defined(_M_X64)
  if (g_simd_override == 2)
  {
    g_up2_impl = (cpu_supports_avx2() != 0) ? up2_avx2 : up2_scalar;
    return;
  }
  if (cpu_supports_avx2() != 0)
  {
    g_up2_impl = up2_avx2;
    return;
  }
#endif
  g_up2_impl = up2_scalar;
}

/**
 * @brief Downsample a single-channel float image by 2x using SIMD-dispatched implementation.
 *
 * Selects an AVX2-accelerated kernel when available (and not overridden), otherwise falls back
 * to a scalar implementation. Inputs must satisfy dw == sw/2 and dh == sh/2. No-op on invalid
 * arguments.
 *
 * @param src Source pointer (size sw x sh, stride sp in elements).
 * @param sw Source width in elements.
 * @param sh Source height in elements.
 * @param sp Source row stride in elements.
 * @param dst Destination pointer (size dw x dh, stride dp in elements).
 * @param dw Destination width in elements (must be sw/2).
 * @param dh Destination height in elements (must be sh/2).
 * @param dp Destination row stride in elements.
 */
void zx_lod_downsample_2x(const float* src, uint32_t sw, uint32_t sh, uint32_t sp, float* dst,
                          uint32_t dw, uint32_t dh, uint32_t dp)
{
  std::call_once(g_down2_once, init_down2_impl);
  g_down2_impl(src, sw, sh, sp, dst, dw, dh, dp);
}

/**
 * @brief Upsample a single-channel float image by 2x using SIMD-dispatched implementation.
 *
 * Selects an AVX2-accelerated kernel when available (and not overridden), otherwise falls back
 * to a scalar implementation. Inputs must satisfy dw == sw*2 and dh == sh*2. No-op on invalid
 * arguments.
 *
 * @param src Source pointer (size sw x sh, stride sp in elements).
 * @param sw Source width in elements.
 * @param sh Source height in elements.
 * @param sp Source row stride in elements.
 * @param dst Destination pointer (size dw x dh, stride dp in elements).
 * @param dw Destination width in elements (must be sw*2).
 * @param dh Destination height in elements (must be sh*2).
 * @param dp Destination row stride in elements.
 */
void zx_lod_upsample_2x(const float* src, uint32_t sw, uint32_t sh, uint32_t sp, float* dst,
                        uint32_t dw, uint32_t dh, uint32_t dp)
{
  std::call_once(g_up2_once, init_up2_impl);
  g_up2_impl(src, sw, sh, sp, dst, dw, dh, dp);
}

/*
 * SIMD dispatch for border consistency check (Windows-first, AVX2 primary)
 */

typedef float (*pp_border_fn)(const float*, uint32_t, uint32_t, uint32_t, const float*, uint32_t,
                              uint32_t, uint32_t, int);

static float border_scalar(const float* a, uint32_t a_w, uint32_t a_h, uint32_t a_p, const float* b,
                           uint32_t b_w, uint32_t b_h, uint32_t b_p, int side)
{
  if ((a == nullptr) || (b == nullptr) || (a_w == 0U) || (a_h == 0U) || (b_w == 0U) || (b_h == 0U))
  {
    return 0.0F;
  }
  float max_abs = 0.0F;
  if (side == 0)
  {  // A right edge vs B left edge
    uint32_t x_a = a_w - 1U;
    uint32_t x_b = 0U;
    uint32_t h   = std::min(a_h, b_h);
    for (uint32_t y = 0; y < h; ++y)
    {
      float da = a[(y * a_p) + x_a];
      float db = b[(y * b_p) + x_b];
      auto d   = static_cast<float>(std::fabs(static_cast<double>(da) - static_cast<double>(db)));
      max_abs  = std::max(max_abs, d);
    }
  }
  else
  {  // bottom vs top
    uint32_t y_a = a_h - 1U;
    uint32_t y_b = 0U;
    uint32_t w   = std::min(a_w, b_w);
    for (uint32_t x = 0; x < w; ++x)
    {
      float da = a[(y_a * a_p) + x];
      float db = b[(y_b * b_p) + x];
      auto d   = static_cast<float>(std::fabs(static_cast<double>(da) - static_cast<double>(db)));
      max_abs  = std::max(max_abs, d);
    }
  }
  return max_abs;
}

#if defined(__x86_64__) || defined(_M_X64)
static inline int cpu_supports_avx2()
{
#ifdef _MSC_VER
  int cpu_info[4] = {0};
  __cpuid(cpu_info, 0);
  int n_ids = cpu_info[0];
  int avx2  = 0;
  if (n_ids >= k_cpuid_leaf_avx2)
  {
    __cpuidex(cpu_info, k_cpuid_leaf_avx2, 0);
    avx2 = ((cpu_info[1] & (1 << k_cpuid_ebx_avx2_bit)) != 0) ? 1 : 0;
  }
  // Check OSXSAVE/AVX state
  __cpuid(cpu_info, 1);
  const bool osxsave = ((cpu_info[2] & (1 << k_cpuid_ecx_osxsave_bit)) != 0);
  if (osxsave)
  {
    unsigned long long x = _xgetbv(0);
    if ((x & k_xcr0_xmm_ymm_mask) != k_xcr0_xmm_ymm_mask)
    {
      avx2 = 0;
    }
  }
  return avx2;
#elif defined(__GNUC__) || defined(__clang__)
  return __builtin_cpu_supports("avx2");
#else
  return 0;
#endif
}

static ZX_TARGET_AVX2 float border_avx2(const float* a, uint32_t a_w, uint32_t a_h, uint32_t a_p,
                                        const float* b, uint32_t b_w, uint32_t b_h, uint32_t b_p,
                                        int side)
{
  if ((a == nullptr) || (b == nullptr) || (a_w == 0U) || (a_h == 0U) || (b_w == 0U) || (b_h == 0U))
  {
    return 0.0F;
  }
  float max_abs_scalar  = 0.0F;
  __m256 vmax           = _mm256_set1_ps(0.0F);
  const __m256 signmask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
  if (side == 0)
  {
    uint32_t x_a = a_w - 1U;
    uint32_t x_b = 0U;
    uint32_t h   = std::min(a_h, b_h);
    uint32_t y   = 0;
    while ((y + k_block8) <= h)
    {
      std::array<int, k_avx_lanes> idx_a{};
      std::array<int, k_avx_lanes> idx_b{};
      for (int i = 0; i < k_avx_lanes; ++i)
      {
        idx_a[static_cast<size_t>(i)] =
            static_cast<int>((static_cast<uint32_t>(y + i) * a_p) + x_a);
        idx_b[static_cast<size_t>(i)] =
            static_cast<int>((static_cast<uint32_t>(y + i) * b_p) + x_b);
      }
      __m256 va = _mm256_i32gather_ps(
          a, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(idx_a.data())), k_gather_scale);
      __m256 vb = _mm256_i32gather_ps(
          b, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(idx_b.data())), k_gather_scale);
      __m256 d = _mm256_sub_ps(va, vb);
      d        = _mm256_and_ps(d, signmask);
      vmax     = _mm256_max_ps(vmax, d);
      y += k_block8;
    }
    alignas(k_avx_align) std::array<float, k_avx_lanes> tmp{};
    _mm256_store_ps(tmp.data(), vmax);
    for (float v : tmp)
    {
      max_abs_scalar = std::max(max_abs_scalar, v);
    }
    for (; y < h; ++y)
    {
      float da = a[(y * a_p) + x_a];
      float db = b[(y * b_p) + x_b];
      auto d   = static_cast<float>(std::fabs(static_cast<double>(da) - static_cast<double>(db)));
      max_abs_scalar = std::max(max_abs_scalar, d);
    }
  }
  else
  {
    uint32_t y_a = a_h - 1U;
    uint32_t y_b = 0U;
    uint32_t w   = std::min(a_w, b_w);
    uint32_t x   = 0;
    while ((x + k_block8) <= w)
    {
      __m256 va = _mm256_loadu_ps(&a[(y_a * a_p) + x]);
      __m256 vb = _mm256_loadu_ps(&b[(y_b * b_p) + x]);
      __m256 d  = _mm256_sub_ps(va, vb);
      d         = _mm256_and_ps(d, signmask);
      vmax      = _mm256_max_ps(vmax, d);
      x += k_block8;
    }
    alignas(k_avx_align) std::array<float, k_avx_lanes> tmp{};
    _mm256_store_ps(tmp.data(), vmax);
    for (float v : tmp)
    {
      max_abs_scalar = std::max(max_abs_scalar, v);
    }
    for (; x < w; ++x)
    {
      float da = a[(y_a * a_p) + x];
      float db = b[(y_b * b_p) + x];
      auto d   = static_cast<float>(std::fabs(static_cast<double>(da) - static_cast<double>(db)));
      max_abs_scalar = std::max(max_abs_scalar, d);
    }
  }
  return max_abs_scalar;
}
#endif

static pp_border_fn g_border_impl = border_scalar;
static std::once_flag g_border_once;

static void init_border_impl()
{
  if (g_simd_override == 0)
  {
    g_border_impl = border_scalar;
    return;
  }
#if defined(__x86_64__) || defined(_M_X64)
  if (g_simd_override == 2)
  {
    g_border_impl = (cpu_supports_avx2() != 0) ? border_avx2 : border_scalar;
    return;
  }
  if (cpu_supports_avx2() != 0)
  {
    g_border_impl = border_avx2;
    return;
  }
#endif
  g_border_impl = border_scalar;
}

float zx_lod_border_consistency_check(const float* a, uint32_t a_w, uint32_t a_h, uint32_t a_p,
                                      const float* b, uint32_t b_w, uint32_t b_h, uint32_t b_p,
                                      int side)
{
  std::call_once(g_border_once, init_border_impl);
  return g_border_impl(a, a_w, a_h, a_p, b, b_w, b_h, b_p, side);
}

static void recompute_dispatch()
{
  // border
#if defined(__x86_64__) || defined(_M_X64)
  const bool has_avx2 = (cpu_supports_avx2() != 0);
  switch (g_simd_override)
  {
  case 0:
    g_border_impl = border_scalar;
    break;
  case 2:
    g_border_impl = has_avx2 ? border_avx2 : border_scalar;
    break;
  default:
    g_border_impl = border_scalar;
    break;
  }
#else
  g_border_impl = border_scalar;
#endif
  // downsample
#if defined(__x86_64__) || defined(_M_X64)
  switch (g_simd_override)
  {
  case 0:
    g_down2_impl = down2_scalar;
    break;
  case 2:
    g_down2_impl = has_avx2 ? down2_avx2 : down2_scalar;
    break;
  default:
    g_down2_impl = down2_scalar;
    break;
  }
#else
  g_down2_impl = down2_scalar;
#endif
  // upsample
#if defined(__x86_64__) || defined(_M_X64)
  switch (g_simd_override)
  {
  case 0:
    g_up2_impl = up2_scalar;
    break;
  case 2:
    g_up2_impl = has_avx2 ? up2_avx2 : up2_scalar;
    break;
  default:
    g_up2_impl = up2_scalar;
    break;
  }
#else
  g_up2_impl = up2_scalar;
#endif
}

void ZX_CALL zx_lod_set_simd_override(int mode)
{
  g_simd_override = mode;
  recompute_dispatch();
}

extern "C" void ZX_CALL zx_lod_set_store_policy(int mode)
{
  if ((mode < -1) || (mode > 1))
  {
    mode = -1;
  }
  g_store_policy = mode;
}

extern "C" int ZX_CALL zx_lod_get_store_policy(void)
{
  return g_store_policy;
}

void ZX_CALL zx_lod_fallback_init(zx_lod_fallback_state* s)
{
  if (s != nullptr)
  {
    s->active_frames   = 0;
    s->inactive_frames = 0;
    s->blend_remaining = 0;
    s->activations     = 0;
  }
}

int ZX_CALL zx_lod_fallback_update(const zx_lod_fallback_policy* p, uint32_t residency_active_tiles,
                                   float last_step_ms, zx_lod_fallback_state* s)
{
  if ((p == nullptr) || (s == nullptr))
  {
    return 0;
  }
  const int pressure =
      ((residency_active_tiles > p->active_tiles_max) || (last_step_ms > p->step_ms_max)) ? 1 : 0;
  int use_coarse = (s->active_frames > 0) ? 1 : 0;
  if (pressure != 0)
  {
    s->active_frames += 1;
    s->inactive_frames = 0;
    if (s->active_frames == 1)
    {
      s->blend_remaining = p->blend_frames;
      s->activations += 1;
    }
  }
  else
  {
    s->inactive_frames += 1;
    s->active_frames = (s->inactive_frames >= p->exit_frames) ? 0 : s->active_frames;
  }
  if (s->active_frames >= p->enter_frames)
  {
    use_coarse = 1;
  }
  if (s->inactive_frames >= p->exit_frames)
  {
    use_coarse = 0;
  }
  if (s->blend_remaining > 0)
  {
    s->blend_remaining -= 1;
  }
  return use_coarse;
}

/* Global fallback defaults */
static int g_lod_enabled                   = 0;
static zx_lod_fallback_policy g_lod_policy = {k_default_active_tiles_max, k_default_step_ms_max,
                                              k_default_enter_frames, k_default_exit_frames,
                                              k_default_blend_frames};

void ZX_CALL zx_lod_set_enabled(int on)
{
  g_lod_enabled = (on != 0) ? 1 : 0;
}
int ZX_CALL zx_lod_is_enabled()
{
  return g_lod_enabled;
}
void ZX_CALL zx_lod_set_default_policy(const zx_lod_fallback_policy* p)
{
  if (p != nullptr)
  {
    g_lod_policy = *p;
  }
}
void ZX_CALL zx_lod_get_default_policy(zx_lod_fallback_policy* out)
{
  if (out != nullptr)
  {
    *out = g_lod_policy;
  }
}
