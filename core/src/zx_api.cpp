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

/* Remove unused stub to silence -Wunused-function; error_string_impl is the public API */

static zx_status ZX_CALL create_context_stub(const zx_context_desc* d, zx_context* out)
{
    if (!d || !out || d->size < sizeof(zx_context_desc)) return ZX_E_INVALID;
    *out = (zx_context)1; /* non-zero token stub */
    return ZX_OK;
}

static void ZX_CALL destroy_context_stub(zx_context) {}
static zx_status ZX_CALL bind_device_stub(zx_context, const zx_device_desc* dd, zx_device* out)
{ if (!dd || !out || dd->size < sizeof(zx_device_desc)) return ZX_E_INVALID; *out=(zx_device)1; return ZX_OK; }
static zx_status ZX_CALL unbind_device_stub(zx_device){ return ZX_OK; }
static zx_status ZX_CALL create_scene_stub(zx_context, const zx_scene_desc* sd, zx_scene* out)
{ if (!sd || !out || sd->size < sizeof(zx_scene_desc)) return ZX_E_INVALID; *out=(zx_scene)1; return ZX_OK; }
static void ZX_CALL destroy_scene_stub(zx_scene) {}
static zx_status ZX_CALL create_material_stub(zx_context, const zx_material_desc* md, zx_material* out)
{ if (!md || !out || md->size < sizeof(zx_material_desc)) return ZX_E_INVALID; *out=(zx_material)1; return ZX_OK; }
static void ZX_CALL destroy_material_stub(zx_material) {}
static zx_status ZX_CALL create_rig_stub(zx_context, const zx_rig_desc* rd, zx_rig* out)
{ if (!rd || !out || rd->size < sizeof(zx_rig_desc)) return ZX_E_INVALID; *out=(zx_rig)1; return ZX_OK; }
static void ZX_CALL destroy_rig_stub(zx_rig) {}
static zx_status ZX_CALL begin_frame_stub(zx_scene, const zx_frame_begin* fb)
{ if (!fb || fb->size < sizeof(zx_frame_begin)) return ZX_E_INVALID; return ZX_OK; }
static zx_status ZX_CALL simulate_stub(zx_scene, const zx_sim_params* sp)
{ if (!sp || sp->size < sizeof(zx_sim_params)) return ZX_E_INVALID; return ZX_OK; }
static zx_status ZX_CALL writeback_stub(zx_scene, const zx_writeback_desc* wb)
{
    if (!wb || wb->size < sizeof(zx_writeback_desc)) return ZX_E_INVALID;
    /* If a CPU displacement buffer is provided, perform a trivial clear for now */
    if (wb->displacement && wb->width && wb->height && wb->row_pitch_bytes) {
        float* disp = (float*)wb->displacement;
        const uint32_t pitch = wb->row_pitch_bytes / sizeof(float);
        for (uint32_t y=0; y<wb->height; ++y) {
            for (uint32_t x=0; x<wb->width; ++x) disp[y*pitch + x] = 0.0f;
        }
    }
    return ZX_OK;
}
static zx_status ZX_CALL end_frame_stub(zx_scene, const zx_frame_end* fe)
{ if (!fe || fe->size < sizeof(zx_frame_end)) return ZX_E_INVALID; return ZX_OK; }
static zx_status ZX_CALL get_counters_stub(zx_context, zx_counters* c, uint32_t* count)
{
    if (!c || !count || c->size < sizeof(zx_counters)) return ZX_E_INVALID;
    *count = 1;
    c->version = 0x00010000u;
    c->particle_count = 0;
    c->active_tiles = 0;
    c->contact_saturated = 0;
    c->last_step_ms = 0.0f;
    c->p2g_ms = c->grid_ms = c->g2p_ms = c->writeback_ms = 0.0f;
    c->mass_drift_abs = 0.0f;
    c->max_penetration = 0.0f;
    c->substeps = 1;
    return ZX_OK;
}
static const char* ZX_CALL error_string_impl(zx_status s)
{
    switch (s) {
        case ZX_OK: return "ZX_OK";
        case ZX_E_INVALID: return "ZX_E_INVALID";
        case ZX_E_UNSUPPORTED: return "ZX_E_UNSUPPORTED";
        case ZX_E_OOM: return "ZX_E_OOM";
        case ZX_E_DEVICE: return "ZX_E_DEVICE";
        case ZX_E_TIMEOUT: return "ZX_E_TIMEOUT";
        case ZX_W_SOFT: return "ZX_W_SOFT";
        default: return "ZX_UNKNOWN";
    }
}

ZX_API zx_status ZX_CALL zxGetProcTable(uint32_t abi_version, zx_procs* out_procs)
{
    if (!out_procs) return ZX_E_INVALID;
    /* Zero on entry so callers never observe stale pointers on failure */
    memset(out_procs, 0, sizeof(*out_procs));
    if (abi_version != ZX_ABI_VERSION) return ZX_E_UNSUPPORTED;
    out_procs->size = sizeof(zx_procs);
    out_procs->version = ZX_ABI_VERSION;
    out_procs->create_context = &create_context_stub;
    out_procs->destroy_context = &destroy_context_stub;
    out_procs->bind_device = &bind_device_stub;
    out_procs->unbind_device = &unbind_device_stub;
    out_procs->create_scene = &create_scene_stub;
    out_procs->destroy_scene = &destroy_scene_stub;
    out_procs->create_material = &create_material_stub;
    out_procs->destroy_material = &destroy_material_stub;
    out_procs->create_rig = &create_rig_stub;
    out_procs->destroy_rig = &destroy_rig_stub;
    out_procs->begin_frame = &begin_frame_stub;
    out_procs->simulate = &simulate_stub;
    out_procs->writeback = &writeback_stub;
    out_procs->end_frame = &end_frame_stub;
    out_procs->get_counters = &get_counters_stub;
    out_procs->error_string = &error_string_impl;
    return ZX_OK;
}


