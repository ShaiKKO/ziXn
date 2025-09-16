/*!
 * \file zx_lod.cpp
 * \brief Presentation LOD filters implementation.
 */

#include "zx/zx_lod.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <mutex>
#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

#if defined(__GNUC__) || defined(__clang__)
#define ZX_TARGET_AVX2 __attribute__((target("avx2")))
#else
#define ZX_TARGET_AVX2
#endif

/* SIMD override: -1 auto, 0 scalar, 2 AVX2 */
static int g_simd_override = -1;

#if defined(__x86_64__) || defined(_M_X64)
static inline int cpu_supports_avx2(); /* forward */
#endif
/* SIMD-dispatched 2x downsample and upsample kernels */
typedef void (*pp_down2_fn)(const float*,uint32_t,uint32_t,uint32_t,float*,uint32_t,uint32_t,uint32_t);
typedef void (*pp_up2_fn)(const float*,uint32_t,uint32_t,uint32_t,float*,uint32_t,uint32_t,uint32_t);

static void down2_scalar(const float* src, uint32_t sw, uint32_t sh, uint32_t sp,
                         float* dst, uint32_t dw, uint32_t dh, uint32_t dp)
{
    if (!src || !dst || sw < 2 || sh < 2 || dw*2 != sw || dh*2 != sh) return;
    for (uint32_t y = 0; y < dh; ++y) {
        for (uint32_t x = 0; x < dw; ++x) {
            const float a = src[(2*y+0)*sp + (2*x+0)];
            const float b = src[(2*y+0)*sp + (2*x+1)];
            const float c = src[(2*y+1)*sp + (2*x+0)];
            const float d = src[(2*y+1)*sp + (2*x+1)];
            dst[y*dp + x] = 0.25f * (a + b + c + d);
        }
    }
}

#if defined(__x86_64__) || defined(_M_X64)
static ZX_TARGET_AVX2 void down2_avx2(const float* src, uint32_t sw, uint32_t sh, uint32_t sp,
                                      float* dst, uint32_t dw, uint32_t dh, uint32_t dp)
{
    if (!src || !dst || sw < 2 || sh < 2 || dw*2 != sw || dh*2 != sh) return;
    const uint32_t vecW = (dw/8)*8;
    const __m256 quarter = _mm256_set1_ps(0.25f);
    for (uint32_t y=0;y<dh;++y){
        const float* row0 = src + (2*y+0)*sp;
        const float* row1 = src + (2*y+1)*sp;
        uint32_t x=0;
        for (; x<vecW; x+=8){
            // Load 16 contiguous samples from each row
            __m256 r0a = _mm256_loadu_ps(row0 + 2*x + 0);
            __m256 r0b = _mm256_loadu_ps(row0 + 2*x + 8);
            __m256 r1a = _mm256_loadu_ps(row1 + 2*x + 0);
            __m256 r1b = _mm256_loadu_ps(row1 + 2*x + 8);
            // Horizontal add adjacent pairs within lanes
            __m256 s0 = _mm256_hadd_ps(r0a, r0b);
            __m256 s1 = _mm256_hadd_ps(r1a, r1b);
            __m256 sum = _mm256_add_ps(s0, s1);
            __m256 out = _mm256_mul_ps(sum, quarter);
            _mm256_storeu_ps(dst + y*dp + x, out);
        }
        for (; x<dw; ++x){
            const float a = row0[2*x+0]; const float b = row0[2*x+1];
            const float c = row1[2*x+0]; const float d = row1[2*x+1];
            dst[y*dp + x] = 0.25f * (a+b+c+d);
        }
    }
}
#endif

static pp_down2_fn g_down2_impl = down2_scalar;
static std::once_flag g_down2_once;
static void init_down2_impl(){
    if (g_simd_override == 0){ g_down2_impl = down2_scalar; return; }
#if defined(__x86_64__) || defined(_M_X64)
    if (g_simd_override == 2){ g_down2_impl = cpu_supports_avx2()? down2_avx2 : down2_scalar; return; }
    if (cpu_supports_avx2()){ g_down2_impl = down2_avx2; return; }
#endif
    g_down2_impl = down2_scalar;
}

static void up2_scalar(const float* src, uint32_t sw, uint32_t sh, uint32_t sp,
                       float* dst, uint32_t dw, uint32_t dh, uint32_t dp)
{
    if (!src || !dst || dw != sw*2 || dh != sh*2) return;
    for (uint32_t y = 0; y < dh; ++y) {
        float fy = (y + 0.5f) * 0.5f - 0.5f;
        int y0 = (int)std::floor((double)fy); int y1 = std::min((int)sh-1, y0+1);
        float ty = fy - (float)y0; y0 = std::max(0, y0);
        for (uint32_t x = 0; x < dw; ++x) {
            float fx = (x + 0.5f) * 0.5f - 0.5f;
            int x0 = (int)std::floor((double)fx); int x1 = std::min((int)sw-1, x0+1);
            float tx = fx - (float)x0; x0 = std::max(0, x0);
            float a = src[y0*sp + x0]; float b = src[y0*sp + x1];
            float c = src[y1*sp + x0]; float d = src[y1*sp + x1];
            float ab = a + tx*(b - a);
            float cd = c + tx*(d - c);
            dst[y*dp + x] = ab + ty*(cd - ab);
        }
    }
}

#if defined(__x86_64__) || defined(_M_X64)
static ZX_TARGET_AVX2 void up2_avx2(const float* src, uint32_t sw, uint32_t sh, uint32_t sp,
                                    float* dst, uint32_t dw, uint32_t dh, uint32_t dp)
{
    // Placeholder: fallback to scalar until optimized AVX2 bilinear is implemented
    up2_scalar(src,sw,sh,sp,dst,dw,dh,dp);
}
#endif

static pp_up2_fn g_up2_impl = up2_scalar;
static std::once_flag g_up2_once;
static void init_up2_impl(){
    if (g_simd_override == 0){ g_up2_impl = up2_scalar; return; }
#if defined(__x86_64__) || defined(_M_X64)
    if (g_simd_override == 2){ g_up2_impl = cpu_supports_avx2()? up2_avx2 : up2_scalar; return; }
    if (cpu_supports_avx2()){ g_up2_impl = up2_avx2; return; }
#endif
    g_up2_impl = up2_scalar;
}

void zx_lod_downsample_2x(const float* src, uint32_t sw, uint32_t sh, uint32_t sp,
                          float* dst, uint32_t dw, uint32_t dh, uint32_t dp)
{
    std::call_once(g_down2_once, init_down2_impl);
    g_down2_impl(src,sw,sh,sp,dst,dw,dh,dp);
}

void zx_lod_upsample_2x(const float* src, uint32_t sw, uint32_t sh, uint32_t sp,
                        float* dst, uint32_t dw, uint32_t dh, uint32_t dp)
{
    std::call_once(g_up2_once, init_up2_impl);
    g_up2_impl(src,sw,sh,sp,dst,dw,dh,dp);
}

/*
 * SIMD dispatch for border consistency check (Windows-first, AVX2 primary)
 */

typedef float (*pp_border_fn)(const float*, uint32_t, uint32_t, uint32_t,
                              const float*, uint32_t, uint32_t, uint32_t,
                              int);

static float border_scalar(const float* A, uint32_t Aw, uint32_t Ah, uint32_t Ap,
                           const float* B, uint32_t Bw, uint32_t Bh, uint32_t Bp,
                           int side)
{
    if (!A || !B || Aw==0 || Ah==0 || Bw==0 || Bh==0) return 0.0f;
    float max_abs = 0.0f;
    if (side == 0) { // A right edge vs B left edge
        uint32_t xA = Aw - 1, xB = 0;
        uint32_t h = std::min(Ah, Bh);
        for (uint32_t y=0; y<h; ++y) {
            float da = A[y*Ap + xA]; float db = B[y*Bp + xB];
            float d = (float)std::fabs((double)da - (double)db);
            if (d > max_abs) max_abs = d;
        }
    } else { // bottom vs top
        uint32_t yA = Ah - 1, yB = 0;
        uint32_t w = std::min(Aw, Bw);
        for (uint32_t x=0; x<w; ++x) {
            float da = A[yA*Ap + x]; float db = B[yB*Bp + x];
            float d = (float)std::fabs((double)da - (double)db);
            if (d > max_abs) max_abs = d;
        }
    }
    return max_abs;
}

#if defined(__x86_64__) || defined(_M_X64)
static inline int cpu_supports_avx2(){
#if defined(_MSC_VER)
    int cpuInfo[4] = {0};
    __cpuid(cpuInfo, 0);
    int nIds = cpuInfo[0];
    int avx2 = 0;
    if (nIds >= 7){ __cpuidex(cpuInfo, 7, 0); avx2 = (cpuInfo[1] & (1<<5)) != 0; }
    // Check OSXSAVE/AVX state
    __cpuid(cpuInfo, 1);
    int osxsave = (cpuInfo[2] & (1<<27)) != 0;
    if (osxsave){ unsigned long long x = _xgetbv(0); if ((x & 0x6) != 0x6) avx2 = 0; }
    return avx2;
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_cpu_supports("avx2");
#else
    return 0;
#endif
}

static ZX_TARGET_AVX2 float border_avx2(const float* A, uint32_t Aw, uint32_t Ah, uint32_t Ap,
                         const float* B, uint32_t Bw, uint32_t Bh, uint32_t Bp,
                         int side)
{
    if (!A || !B || Aw==0 || Ah==0 || Bw==0 || Bh==0) return 0.0f;
    float max_abs_scalar = 0.0f;
    __m256 vmax = _mm256_set1_ps(0.0f);
    const __m256 signmask = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    if (side == 0){
        uint32_t xA = Aw - 1, xB = 0;
        uint32_t h = std::min(Ah, Bh);
        uint32_t y=0;
        for (; y+8<=h; y+=8){
            int idxA[8]; int idxB[8];
            for (int i=0;i<8;++i){ idxA[i] = (int)((y+i)*Ap + xA); idxB[i] = (int)((y+i)*Bp + xB); }
            __m256 a = _mm256_i32gather_ps(A, _mm256_loadu_si256((const __m256i*)idxA), 4);
            __m256 b = _mm256_i32gather_ps(B, _mm256_loadu_si256((const __m256i*)idxB), 4);
            __m256 d = _mm256_sub_ps(a,b);
            d = _mm256_and_ps(d, signmask);
            vmax = _mm256_max_ps(vmax, d);
        }
        alignas(32) float tmp[8]; _mm256_store_ps(tmp, vmax);
        for (int i=0;i<8;++i) if (tmp[i]>max_abs_scalar) max_abs_scalar = tmp[i];
        for (; y<h; ++y){
            float da = A[y*Ap + xA]; float db = B[y*Bp + xB];
            float d = (float)std::fabs((double)da - (double)db);
            if (d > max_abs_scalar) max_abs_scalar = d;
        }
    } else {
        uint32_t yA = Ah - 1, yB = 0;
        uint32_t w = std::min(Aw, Bw);
        uint32_t x=0;
        for (; x+8<=w; x+=8){
            __m256 a = _mm256_loadu_ps(&A[yA*Ap + x]);
            __m256 b = _mm256_loadu_ps(&B[yB*Bp + x]);
            __m256 d = _mm256_sub_ps(a,b);
            d = _mm256_and_ps(d, signmask);
            vmax = _mm256_max_ps(vmax, d);
        }
        alignas(32) float tmp[8]; _mm256_store_ps(tmp, vmax);
        for (int i=0;i<8;++i) if (tmp[i]>max_abs_scalar) max_abs_scalar = tmp[i];
        for (; x<w; ++x){
            float da = A[yA*Ap + x]; float db = B[yB*Bp + x];
            float d = (float)std::fabs((double)da - (double)db);
            if (d > max_abs_scalar) max_abs_scalar = d;
        }
    }
    return max_abs_scalar;
}
#endif

static pp_border_fn g_border_impl = border_scalar;
static std::once_flag g_border_once;

static void init_border_impl(){
    if (g_simd_override == 0) { g_border_impl = border_scalar; return; }
#if defined(__x86_64__) || defined(_M_X64)
    if (g_simd_override == 2){ g_border_impl = cpu_supports_avx2()? border_avx2 : border_scalar; return; }
    if (cpu_supports_avx2()) { g_border_impl = border_avx2; return; }
#endif
    g_border_impl = border_scalar;
}

float zx_lod_border_consistency_check(const float* A, uint32_t Aw, uint32_t Ah, uint32_t Ap,
                                      const float* B, uint32_t Bw, uint32_t Bh, uint32_t Bp,
                                      int side)
{
    std::call_once(g_border_once, init_border_impl);
    return g_border_impl(A,Aw,Ah,Ap,B,Bw,Bh,Bp,side);
}

static void recompute_dispatch(){
    // border
#if defined(__x86_64__) || defined(_M_X64)
    if (g_simd_override == 0) g_border_impl = border_scalar;
    else if (g_simd_override == 2) g_border_impl = cpu_supports_avx2()? border_avx2 : border_scalar;
    else g_border_impl = cpu_supports_avx2()? border_avx2 : border_scalar;
#else
    g_border_impl = border_scalar;
#endif
    // downsample
    extern pp_down2_fn g_down2_impl; extern void init_down2_impl();
#if defined(__x86_64__) || defined(_M_X64)
    if (g_simd_override == 0) g_down2_impl = down2_scalar;
    else if (g_simd_override == 2) g_down2_impl = cpu_supports_avx2()? down2_avx2 : down2_scalar;
    else g_down2_impl = cpu_supports_avx2()? down2_avx2 : down2_scalar;
#else
    g_down2_impl = down2_scalar;
#endif
    // upsample
    extern pp_up2_fn g_up2_impl; extern void init_up2_impl();
#if defined(__x86_64__) || defined(_M_X64)
    if (g_simd_override == 0) g_up2_impl = up2_scalar;
    else if (g_simd_override == 2) g_up2_impl = cpu_supports_avx2()? up2_avx2 : up2_scalar;
    else g_up2_impl = cpu_supports_avx2()? up2_avx2 : up2_scalar;
#else
    g_up2_impl = up2_scalar;
#endif
}

ZX_API void ZX_CALL zx_lod_set_simd_override(int mode){ g_simd_override = mode; recompute_dispatch(); }

ZX_API void ZX_CALL zx_lod_fallback_init(zx_lod_fallback_state* s){ if (s){ s->active_frames=0; s->inactive_frames=0; s->blend_remaining=0; s->activations=0; } }

ZX_API int ZX_CALL zx_lod_fallback_update(const zx_lod_fallback_policy* p,
                                          uint32_t residency_active_tiles,
                                          float last_step_ms,
                                          zx_lod_fallback_state* s)
{
    if (!p || !s) return 0;
    const int pressure = (residency_active_tiles > p->active_tiles_max) || (last_step_ms > p->step_ms_max);
    int use_coarse = (s->active_frames > 0);
    if (pressure) {
        s->active_frames += 1; s->inactive_frames = 0;
        if (s->active_frames == 1) { s->blend_remaining = p->blend_frames; s->activations += 1; }
    } else {
        s->inactive_frames += 1; s->active_frames = (s->inactive_frames >= p->exit_frames) ? 0 : s->active_frames;
    }
    if (s->active_frames >= p->enter_frames) use_coarse = 1;
    if (s->inactive_frames >= p->exit_frames) use_coarse = 0;
    if (s->blend_remaining > 0) s->blend_remaining -= 1;
    return use_coarse;
}

/* Global fallback defaults */
static int g_lod_enabled = 0;
static zx_lod_fallback_policy g_lod_policy = { 1000000u, 1e9f, 2u, 2u, 3u };

ZX_API void ZX_CALL zx_lod_set_enabled(int on){ g_lod_enabled = (on!=0); }
ZX_API int  ZX_CALL zx_lod_is_enabled(void){ return g_lod_enabled; }
ZX_API void ZX_CALL zx_lod_set_default_policy(const zx_lod_fallback_policy* p){ if (p) g_lod_policy = *p; }
ZX_API void ZX_CALL zx_lod_get_default_policy(zx_lod_fallback_policy* out){ if (out) *out = g_lod_policy; }


