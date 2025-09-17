/*!
\file zx_abi.h
\brief Public C ABI for ziXn: handles, status codes, proc table.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
*/

#ifndef ZX_ABI_H
#define ZX_ABI_H

#ifdef __cplusplus
extern "C"
{
#endif

#if defined(_WIN32)
#if defined(ZX_BUILD)
#define ZX_API __declspec(dllexport)
#else
#define ZX_API __declspec(dllimport)
#endif
#define ZX_CALL __cdecl
#else
#define ZX_API __attribute__((visibility("default")))
#define ZX_CALL
#endif

/* Pointer restrict qualifier for performance hints */
#if defined(_MSC_VER)
#define ZX_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
#define ZX_RESTRICT __restrict__
#else
#define ZX_RESTRICT
#endif

#include <stddef.h>
#include <stdint.h>

  /* ABI version 0x00010000 = 1.0.0. Bump on breaking struct or proc-table changes. */
  static constexpr uint32_t zx_abi_version = 0x00010000U;  // NOLINT(cppcoreguidelines-macro-usage)

  /* Common fixed sizes */
  static constexpr size_t zx_mat3_size = 9U;

  typedef uint64_t zx_handle;
  typedef zx_handle zx_context;
  typedef zx_handle zx_device;
  typedef zx_handle zx_scene;
  typedef zx_handle zx_material;
  typedef zx_handle zx_rig;

  typedef enum zx_status
  {
    ZX_OK            = 0,
    ZX_E_INVALID     = -1,
    ZX_E_UNSUPPORTED = -2,
    ZX_E_OOM         = -3,
    ZX_E_DEVICE      = -4,
    ZX_E_TIMEOUT     = -5,
    ZX_W_SOFT        = 1
  } zx_status;

  typedef struct zx_vec3
  {
    float x, y, z;
  } zx_vec3;
  typedef struct zx_mat3
  {
    float m[zx_mat3_size];
  } zx_mat3;

  typedef struct zx_context_desc
  {
    uint32_t size; /* must be set by caller to sizeof(zx_context_desc) */
    void* (*host_alloc)(size_t, void*);
    void (*host_free)(void*, void*);
    void* user;
    uint32_t feature_bits; /* e.g., ZX_FEATURE_F16, ZX_FEATURE_DET */
  } zx_context_desc;

  /* Forward-declared descriptor structs (opaque to clients until detailed headers exist) */
  typedef struct zx_device_desc
  {
    uint32_t size;
    uint32_t reserved;
  } zx_device_desc;
  typedef struct zx_scene_desc
  {
    uint32_t size;
    uint32_t reserved;
  } zx_scene_desc;
  typedef struct zx_material_desc
  {
    uint32_t size;
    uint32_t reserved;
  } zx_material_desc;
  typedef struct zx_rig_desc
  {
    uint32_t size;
    uint32_t reserved;
  } zx_rig_desc;

  typedef struct zx_frame_begin
  {
    uint32_t size;
    double time_seconds;
  } zx_frame_begin;
  typedef struct zx_sim_params
  {
    uint32_t size;
    float dt;
    uint32_t substeps;
  } zx_sim_params;
  typedef struct zx_writeback_desc
  {
    uint32_t size;
    uint32_t flags; /* future use */
    /* CPU writeback target (optional for GPU interop) */
    void* displacement; /* float* */
    uint32_t width;
    uint32_t height;
    uint32_t row_pitch_bytes;
  } zx_writeback_desc;
  typedef struct zx_frame_end
  {
    uint32_t size;
    uint64_t fence;
  } zx_frame_end;

  typedef struct zx_counters
  {
    uint32_t size;    /* sizeof(zx_counters) that caller allocated */
    uint32_t version; /* structure version for forward/back compat */
    /* Scene/runtime counters */
    uint32_t particle_count;
    uint32_t active_tiles;
    uint32_t contact_saturated; /* nodes clamped by friction cone */
    /* Timings (ms) */
    float last_step_ms;
    float p2g_ms;
    float grid_ms; /* constitutive + contact */
    float g2p_ms;
    float writeback_ms;
    /* Invariants / diagnostics */
    float mass_drift_abs;  /* absolute mass drift ratio */
    float max_penetration; /* grid-space */
    uint32_t substeps;
  } zx_counters;

  typedef struct zx_procs
  {
    uint32_t size;    /* must be set by caller to sizeof(zx_procs) */
    uint32_t version; /* semantic ABI version, e.g., 0x00010000 */
    /* creation */
    zx_status(ZX_CALL* create_context)(const zx_context_desc*, zx_context*);
    void(ZX_CALL* destroy_context)(zx_context);
    /* device/backends */
    zx_status(ZX_CALL* bind_device)(zx_context, const zx_device_desc*, zx_device*);
    zx_status(ZX_CALL* unbind_device)(zx_device);
    /* scene/material/rig lifecycle */
    zx_status(ZX_CALL* create_scene)(zx_context, const zx_scene_desc*, zx_scene*);
    void(ZX_CALL* destroy_scene)(zx_scene);
    zx_status(ZX_CALL* create_material)(zx_context, const zx_material_desc*, zx_material*);
    void(ZX_CALL* destroy_material)(zx_material);
    zx_status(ZX_CALL* create_rig)(zx_context, const zx_rig_desc*, zx_rig*);
    void(ZX_CALL* destroy_rig)(zx_rig);
    /* frame */
    zx_status(ZX_CALL* begin_frame)(zx_scene, const zx_frame_begin*);
    zx_status(ZX_CALL* simulate)(zx_scene, const zx_sim_params*);
    zx_status(ZX_CALL* writeback)(zx_scene, const zx_writeback_desc*);
    zx_status(ZX_CALL* end_frame)(zx_scene, const zx_frame_end*);
    /* telemetry */
    zx_status(ZX_CALL* get_counters)(zx_context, zx_counters*, uint32_t* count);
    /* errors/logging */
    const char*(ZX_CALL* error_string)(zx_status);
  } zx_procs;

  /** \brief Retrieve the procedure table for the requested ABI version.
   * @param abi_version Expected ABI version (e.g., ZX_ABI_VERSION)
   * @param out_procs Output pointer to proc table (size must be set)
   * @return ZX_OK on success; ZX_E_UNSUPPORTED if version mismatch; ZX_E_INVALID on bad args
   */
  ZX_API zx_status ZX_CALL zx_get_proc_table(uint32_t abi_version, zx_procs* out_procs);
/* Back-compat alias */
#define zxGetProcTable zx_get_proc_table

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* ZX_ABI_H */
