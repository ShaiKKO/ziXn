/*!
\file zx_telemetry.h
\brief Thread-safe telemetry buffer and CSV/JSON exporters.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0

\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_TELEMETRY_H
#define ZX_TELEMETRY_H

#include <cstdint>

#include "zx_abi.h"

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct zx_telemetry zx_telemetry; /* opaque */

  /** \brief Create telemetry buffer; capacity==0 selects a reasonable default. */
  ZX_API zx_telemetry* ZX_CALL zx_telemetry_create(uint32_t capacity);
  ZX_API void ZX_CALL zx_telemetry_destroy(zx_telemetry* ctx);

  /** \brief Begin a step (no-op if ctx==NULL). scene may be NULL. */
  ZX_API void ZX_CALL zx_telemetry_begin_step(zx_telemetry* ctx, const char* scene,
                                              uint32_t step_index);
  /** \brief Set a counter during an active step (no-op if ctx==NULL or not in a step). name must be
   * non-NULL. */
  ZX_API void ZX_CALL zx_telemetry_set_counter(zx_telemetry* ctx, const char* name, float value);
  /** \brief End current step (no-op if ctx==NULL or not in a step). */
  ZX_API void ZX_CALL zx_telemetry_end_step(zx_telemetry* ctx);

  /** \brief Exporters return 0 on success, negative on error; no-op if ctx==NULL. */
  /* Exporters return 0 on success, non-zero on error. */
  ZX_API int ZX_CALL zx_telemetry_export_csv(zx_telemetry* ctx, const char* path);
  ZX_API int ZX_CALL zx_telemetry_export_json(zx_telemetry* ctx, const char* path);

  /** \brief Add an error entry with frame context for triage (no-op if ctx==NULL).
   * @param ctx Telemetry handle
   * @param scene Scene identifier (may be NULL)
   * @param step_index Frame/step index
   * @param code Short machine code (non-NULL)
   * @param message Human-readable message (non-NULL)
   */
  ZX_API void ZX_CALL zx_telemetry_add_error(zx_telemetry* ctx, const char* scene,
                                             uint32_t step_index, const char* code,
                                             const char* message);

#ifdef __cplusplus
}
#endif

#endif /* ZX_TELEMETRY_H */
