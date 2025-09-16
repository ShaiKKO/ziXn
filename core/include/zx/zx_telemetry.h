/*!
\file zx_telemetry.h
\brief Thread-safe telemetry buffer and CSV/JSON exporters.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
*/

#ifndef ZX_TELEMETRY_H
#define ZX_TELEMETRY_H

#include <stdint.h>

#include "zx_abi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct zx_telemetry zx_telemetry; /* opaque */

/* Create telemetry with ring capacity (samples). */
ZX_API zx_telemetry* ZX_CALL zx_telemetry_create(uint32_t capacity);
ZX_API void ZX_CALL zx_telemetry_destroy(zx_telemetry* ctx);

/* Begin and end a step sample. scene is a short identifier (e.g., "dambreak"). */
ZX_API void ZX_CALL zx_telemetry_begin_step(zx_telemetry* ctx, const char* scene, uint32_t step_index);
ZX_API void ZX_CALL zx_telemetry_set_counter(zx_telemetry* ctx, const char* name, float value);
ZX_API void ZX_CALL zx_telemetry_end_step(zx_telemetry* ctx);

/* Exporters return 0 on success, non-zero on error. */
ZX_API int ZX_CALL zx_telemetry_export_csv(zx_telemetry* ctx, const char* path);
ZX_API int ZX_CALL zx_telemetry_export_json(zx_telemetry* ctx, const char* path);

#ifdef __cplusplus
}
#endif

#endif /* ZX_TELEMETRY_H */


