/*!
\file zx_residency.h
\brief Tile residency manager with hysteresis and prefetch hints.
\author Colin Macritchie (Ripple Group, LLC)
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_RESIDENCY_H
#define ZX_RESIDENCY_H

#include "zx_abi.h"
#include <cstdint>

#ifdef __cplusplus
extern "C"
{
#endif

  typedef struct zx_residency zx_residency; /* opaque */

  typedef struct
  {
    /* require this many consecutive hot frames to mark active */
    uint32_t enter_frames;
    /* require this many consecutive cold frames to mark inactive */
    uint32_t exit_frames;
    /* number of rings beyond active radius to prefetch */
    uint32_t prefetch_rings;
  } zx_residency_opts;

  /** \brief Create residency context.
   * @param opts Options (may have zero values; must not be NULL)
   * @return Residency handle; caller owns and must destroy
   */
  ZX_API zx_residency* ZX_CALL zx_residency_create(const zx_residency_opts* opts);
  ZX_API void ZX_CALL zx_residency_destroy(zx_residency* ctx);

  /** \brief Advance residency one frame centered at (cx,cy,cz) with active radius (tiles).
   * @param ctx Residency handle (no-op if NULL)
   * @param cx Center X; @param cy Center Y; @param cz Center Z (tiles)
   * @param active_radius Active radius in tiles
   * @param enters Out: newly entered tiles this tick (may be NULL)
   * @param exits Out: newly exited tiles this tick (may be NULL)
   * @param prefetch_count Out: current prefetch ring size (may be NULL)
   */
  /* Advance residency; no-op if ctx is NULL. Outputs may be NULL. */
  ZX_API void ZX_CALL zx_residency_tick(zx_residency* ctx, int cx, int cy, int cz,
                                        uint32_t active_radius, uint32_t* enters, uint32_t* exits,
                                        uint32_t* prefetch_count);

  /** \brief Get current active tile count (0 if ctx==NULL). */
  ZX_API uint32_t ZX_CALL zx_residency_get_active_count(const zx_residency* ctx);

  /** \brief Pin an axis-aligned box of tiles so they remain active regardless of hysteresis.
   * Coordinates are inclusive; caller may pass x0>x1 etc.; the implementation will normalize.
   */
  ZX_API void ZX_CALL zx_residency_pin_box(zx_residency* ctx, int x0, int y0, int z0, int x1,
                                           int y1, int z1);

  /** \brief Clear all pinned regions. */
  ZX_API void ZX_CALL zx_residency_unpin_all(zx_residency* ctx);

  /** \brief Update prefetch rings at runtime. */
  ZX_API void ZX_CALL zx_residency_set_prefetch_rings(zx_residency* ctx, uint32_t rings);

  /** \brief Retrieve last-tick churn stats (enters/exits and their sum). */
  ZX_API void ZX_CALL zx_residency_get_last_churn(const zx_residency* ctx, uint32_t* enters,
                                                  uint32_t* exits, uint32_t* churn);

#ifdef __cplusplus
}
#endif

#endif /* ZX_RESIDENCY_H */
