/*!
 * \file zx_api.cpp
 * \brief Minimal proc table export for ziXn.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \version 1.0.0
 * \date 2025-09-16
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_abi.h"
#include <string.h>

/* ABI version defined in public header */

/**
 * @brief Minimal stub for creating a context token.
 *
 * Validates the provided descriptor and output pointer, then writes a non-zero
 * token (1) to *out on success.
 *
 * @param d Descriptor for the context; must be non-null and have d->size >=
 * sizeof(zx_context_desc).
 * @param out Output pointer to receive the opaque context token; must be non-null.
 * @return ZX_OK on success, ZX_E_INVALID if inputs are null or the descriptor size is too small.
 */

static zx_status ZX_CALL create_context_stub(const zx_context_desc* d, zx_context* out)
{
  if (!d || !out || d->size < sizeof(zx_context_desc))
    return ZX_E_INVALID;
  *out = (zx_context) 1; /* non-zero token stub */
  return ZX_OK;
}

/**
 * @brief No-op stub that destroys a context token.
 *
 * This stub accepts a zx_context token and performs no action. It exists to
 * satisfy the API surface for implementations that do not require teardown.
 */
static void ZX_CALL destroy_context_stub(zx_context) {}
/**
 * @brief Stub implementation for binding a device to a context.
 *
 * Validates the device descriptor and output pointer, and on success returns
 * a non-zero opaque device token through `out`.
 *
 * @param dd Device descriptor; must be non-null and have `dd->size >= sizeof(zx_device_desc)`.
 * @param out Pointer to receive the created device token; must be non-null.
 * @return ZX_OK on success with `*out` set to a non-zero token (1).
 * @return ZX_E_INVALID if `dd` or `out` is null or `dd->size` is too small.
 */
static zx_status ZX_CALL bind_device_stub(zx_context, const zx_device_desc* dd, zx_device* out)
{
  if (!dd || !out || dd->size < sizeof(zx_device_desc))
    return ZX_E_INVALID;
  *out = (zx_device) 1;
  return ZX_OK;
}
/**
 * @brief No-op stub for unbinding a device.
 *
 * This stub is a placeholder for device unbinding logic; it performs no action
 * and always reports success.
 *
 * @return zx_status ZX_OK
 */
static zx_status ZX_CALL unbind_device_stub(zx_device)
{
  return ZX_OK;
}
/**
 * @brief Create a scene handle from a scene description.
 *
 * Validates the provided scene description and, on success, writes a non-zero scene
 * token to @p out and returns ZX_OK.
 *
 * The function returns ZX_E_INVALID if @p sd or @p out is null, or if @p sd->size
 * is smaller than sizeof(zx_scene_desc).
 *
 * @param sd Pointer to the scene description; must be non-null and have a valid
 *           size field (>= sizeof(zx_scene_desc)).
 * @param out Pointer to receive the created scene token; must be non-null.
 * @return ZX_OK on success, ZX_E_INVALID on invalid input.
 */
static zx_status ZX_CALL create_scene_stub(zx_context, const zx_scene_desc* sd, zx_scene* out)
{
  if (!sd || !out || sd->size < sizeof(zx_scene_desc))
    return ZX_E_INVALID;
  *out = (zx_scene) 1;
  return ZX_OK;
}
/**
 * @brief No-op stub for destroying a scene.
 *
 * Accepts a zx_scene token for API compatibility but performs no action. Present to populate the
 * proc table where a scene-destruction callback is required.
 */
static void ZX_CALL destroy_scene_stub(zx_scene) {}
/**
 * @brief Creates a minimal material handle stub.
 *
 * Validates the provided material description and writes a non-zero token to the output
 * material handle on success.
 *
 * @param md Pointer to a material description; must be non-null and have `md->size`
 *           at least `sizeof(zx_material_desc)`.
 * @param out Pointer to receive the created material token; set to a non-zero value on success.
 * @return ZX_OK on success, or ZX_E_INVALID if inputs are null or `md->size` is too small.
 */
static zx_status ZX_CALL create_material_stub(zx_context, const zx_material_desc* md,
                                              zx_material* out)
{
  if (!md || !out || md->size < sizeof(zx_material_desc))
    return ZX_E_INVALID;
  *out = (zx_material) 1;
  return ZX_OK;
}
/**
 * @brief Destroy a material token.
 *
 * Accepts a material token and performs no action in this stub implementation.
 *
 * @param material Opaque material token (ignored).
 */
static void ZX_CALL destroy_material_stub(zx_material) {}
/**
 * @brief Create a minimal rig handle from a descriptor.
 *
 * Validates the descriptor and output pointer, writes a non-zero opaque rig
 * token on success, and returns a ZX status code.
 *
 * @param rd Pointer to the rig descriptor; must be non-null and have `rd->size`
 *           >= sizeof(zx_rig_desc).
 * @param out Output pointer that receives the created rig token; must be non-null.
 * @return ZX_OK on success (and `*out` is set to a non-zero token).
 * @return ZX_E_INVALID if `rd` or `out` is null or `rd->size` is too small.
 */
static zx_status ZX_CALL create_rig_stub(zx_context, const zx_rig_desc* rd, zx_rig* out)
{
  if (!rd || !out || rd->size < sizeof(zx_rig_desc))
    return ZX_E_INVALID;
  *out = (zx_rig) 1;
  return ZX_OK;
}
/**
 * @brief Stub implementation that destroys a rig token.
 *
 * This is a no-op placeholder: the provided rig token is accepted but ignored,
 * and no resources are freed or modified.
 *
 * @param rig The rig token to destroy (ignored).
 */
static void ZX_CALL destroy_rig_stub(zx_rig) {}
/**
 * @brief Validate a frame-begin descriptor for the given scene.
 *
 * Checks that the provided `fb` pointer is non-null and that its `size`
 * field is at least `sizeof(zx_frame_begin)`.
 *
 * @param fb Pointer to a frame-begin descriptor to validate.
 * @return ZX_OK if `fb` is valid; ZX_E_INVALID if `fb` is null or too small.
 */
static zx_status ZX_CALL begin_frame_stub(zx_scene, const zx_frame_begin* fb)
{
  if (!fb || fb->size < sizeof(zx_frame_begin))
    return ZX_E_INVALID;
  return ZX_OK;
}
/**
 * @brief Validate simulation parameters for a scene and report status.
 *
 * Checks that the provided zx_sim_params pointer is non-null and its
 * size field is at least sizeof(zx_sim_params). The zx_scene parameter
 * is not inspected by this stub.
 *
 * @param sp Pointer to simulation parameters; must be non-null and have
 *           sp->size >= sizeof(zx_sim_params).
 * @return zx_status ZX_OK if `sp` is valid; ZX_E_INVALID if `sp` is null
 *         or its size is too small.
 */
static zx_status ZX_CALL simulate_stub(zx_scene, const zx_sim_params* sp)
{
  if (!sp || sp->size < sizeof(zx_sim_params))
    return ZX_E_INVALID;
  return ZX_OK;
}
/**
 * @brief Perform writeback operations for a scene; clears a CPU displacement buffer if provided.
 *
 * Validates the writeback descriptor and, if a CPU-side displacement buffer and non-zero
 * width/height/row_pitch are provided, interprets the buffer as an array of floats and
 * writes 0.0f into every element inside the width×height region using row_pitch_bytes
 * as the byte stride between rows.
 *
 * @param wb Pointer to a zx_writeback_desc structure describing the writeback; must be non-null
 *           and have wb->size >= sizeof(zx_writeback_desc).
 * @return ZX_OK on success. Returns ZX_E_INVALID if `wb` is null or `wb->size` is too small.
 */
static zx_status ZX_CALL writeback_stub(zx_scene, const zx_writeback_desc* wb)
{
  if (!wb || wb->size < sizeof(zx_writeback_desc))
    return ZX_E_INVALID;
  /* If a CPU displacement buffer is provided, perform a trivial clear for now */
  if (wb->displacement && wb->width && wb->height && wb->row_pitch_bytes)
  {
    float* disp          = (float*) wb->displacement;
    const uint32_t pitch = wb->row_pitch_bytes / sizeof(float);
    for (uint32_t y = 0; y < wb->height; ++y)
    {
      for (uint32_t x = 0; x < wb->width; ++x)
        disp[y * pitch + x] = 0.0f;
    }
  }
  return ZX_OK;
}
/**
 * @brief Validate and finalize a frame.
 *
 * Checks that the provided frame-end descriptor is non-null and sized at least
 * sizeof(zx_frame_end); if valid, indicates success.
 *
 * @param fe Pointer to the frame-end descriptor to validate.
 * @return zx_status ZX_OK if the descriptor is valid; ZX_E_INVALID if `fe` is null
 * or its `size` is smaller than sizeof(zx_frame_end).
 */
static zx_status ZX_CALL end_frame_stub(zx_scene, const zx_frame_end* fe)
{
  if (!fe || fe->size < sizeof(zx_frame_end))
    return ZX_E_INVALID;
  return ZX_OK;
}
/**
 * @brief Populate a counters structure with default/zeroed values.
 *
 * Validates the provided counters pointer and its size, writes a single
 * counters entry and fills it with default values (version 0x00010000,
 * zeroed numeric fields, and substeps = 1).
 *
 * @param c Output counters structure; must be non-null and have c->size >= sizeof(zx_counters).
 * @param count Out parameter; set to the number of entries written (always 1 on success).
 * @return ZX_OK on success; ZX_E_INVALID if inputs are null or c->size is too small.
 */
static zx_status ZX_CALL get_counters_stub(zx_context, zx_counters* c, uint32_t* count)
{
  if (!c || !count || c->size < sizeof(zx_counters))
    return ZX_E_INVALID;
  *count               = 1;
  c->version           = 0x00010000u;
  c->particle_count    = 0;
  c->active_tiles      = 0;
  c->contact_saturated = 0;
  c->last_step_ms      = 0.0f;
  c->p2g_ms = c->grid_ms = c->g2p_ms = c->writeback_ms = 0.0f;
  c->mass_drift_abs                                    = 0.0f;
  c->max_penetration                                   = 0.0f;
  c->substeps                                          = 1;
  return ZX_OK;
}
/**
 * @brief Convert a zx_status value to a human-readable string literal.
 *
 * Returns a static, null-terminated string describing the provided ZX status code
 * (e.g. "ZX_OK", "ZX_E_INVALID"). Unrecognized codes map to "ZX_UNKNOWN".
 *
 * @param s ZX status code to describe.
 * @return const char* Static string literal; caller must not free or modify it.
 */
static const char* ZX_CALL error_string_impl(zx_status s)
{
  switch (s)
  {
  case ZX_OK:
    return "ZX_OK";
  case ZX_E_INVALID:
    return "ZX_E_INVALID";
  case ZX_E_UNSUPPORTED:
    return "ZX_E_UNSUPPORTED";
  case ZX_E_OOM:
    return "ZX_E_OOM";
  case ZX_E_DEVICE:
    return "ZX_E_DEVICE";
  case ZX_E_TIMEOUT:
    return "ZX_E_TIMEOUT";
  case ZX_W_SOFT:
    return "ZX_W_SOFT";
  default:
    return "ZX_UNKNOWN";
  }
}

/**
 * @brief Populate a zx_procs table with the implementation pointers for the supported ABI.
 *
 * Validates inputs and ABI version, clears the output structure to avoid exposing
 * stale pointers on failure, and fills the table with stub implementations and
 * metadata when the requested ABI matches ZX_ABI_VERSION.
 *
 * @param abi_version ABI version requested by the caller. Must equal ZX_ABI_VERSION.
 * @param out_procs Pointer to a zx_procs structure to be filled; must be non-null.
 * @return zx_status ZX_OK on success.
 *         ZX_E_INVALID if out_procs is null.
 *         ZX_E_UNSUPPORTED if abi_version does not match ZX_ABI_VERSION.
 */
zx_status ZX_CALL zx_get_proc_table(uint32_t abi_version, zx_procs* out_procs)
{
  if (!out_procs)
    return ZX_E_INVALID;
  /* Zero on entry so callers never observe stale pointers on failure */
  memset(out_procs, 0, sizeof(*out_procs));
  if (abi_version != zx_abi_version)
    return ZX_E_UNSUPPORTED;
  out_procs->size             = sizeof(zx_procs);
  out_procs->version          = zx_abi_version;
  out_procs->create_context   = &create_context_stub;
  out_procs->destroy_context  = &destroy_context_stub;
  out_procs->bind_device      = &bind_device_stub;
  out_procs->unbind_device    = &unbind_device_stub;
  out_procs->create_scene     = &create_scene_stub;
  out_procs->destroy_scene    = &destroy_scene_stub;
  out_procs->create_material  = &create_material_stub;
  out_procs->destroy_material = &destroy_material_stub;
  out_procs->create_rig       = &create_rig_stub;
  out_procs->destroy_rig      = &destroy_rig_stub;
  out_procs->begin_frame      = &begin_frame_stub;
  out_procs->simulate         = &simulate_stub;
  out_procs->writeback        = &writeback_stub;
  out_procs->end_frame        = &end_frame_stub;
  out_procs->get_counters     = &get_counters_stub;
  out_procs->error_string     = &error_string_impl;
  return ZX_OK;
}
