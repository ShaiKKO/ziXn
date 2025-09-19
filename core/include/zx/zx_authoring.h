/**
 * @file zx_authoring.h
 * @brief Authoring UI mappings and safe clamps for material presets.
 * @details Provides normalized UI models and C-ABI mapping functions that clamp and translate
 *          UI slider inputs into solver parameter structures. Thread-safe and side-effect free.
 * @copyright
 *   (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
 *   Licensed for use within the ziXn project under project terms.
 */

#ifndef ZX_AUTHORING_H
#define ZX_AUTHORING_H

#include "zx_abi.h"
#include "zx_constitutive_ref.h"
#include <cstdint>

/**
 * @brief UI models are normalized [0,1] sliders; mapping functions clamp to safe ranges.
 */

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * @brief Normalized authoring model for sand.
   */
  typedef struct ZxUiSand
  {
    float dryness; /**< 0 = very wet, 1 = very dry. */
    float grain;   /**< 0 = fine, 1 = coarse. */
  } zx_ui_sand;

  /**
   * @brief Normalized authoring model for snow.
   */
  typedef struct ZxUiSnow
  {
    float temperature; /**< proxy η: 0=cold/powder, 1=warm/near slush. */
    float packing;     /**< 0 = loose, 1 = packed. */
  } zx_ui_snow;

  /**
   * @brief Map sand UI to elastic, MC, and NorSand parameters with safe clamps.
   * @param[in] ui Sand UI input (must not be NULL).
   * @param[out] out_elastic Optional elastic parameters.
   * @param[out] out_mc Optional Mohr–Coulomb parameters.
   * @param[out] out_ns Optional NorSand parameters.
   * @param[out] out_ns_state Optional NorSand state.
   */
  ZX_API void ZX_CALL zx_authoring_map_sand(const zx_ui_sand* ui, zx_elastic_params* out_elastic,
                                            zx_mc_params* out_mc, zx_norsand_params* out_ns,
                                            zx_norsand_state* out_ns_state);

  /**
   * @brief Map snow UI to elastic and MC(+cap via cohesion) with safe clamps.
   * @param[in] ui Snow UI input (must not be NULL).
   * @param[out] out_elastic Optional elastic parameters.
   * @param[out] out_mc Optional Mohr–Coulomb parameters.
   */
  ZX_API void ZX_CALL zx_authoring_map_snow(const zx_ui_snow* ui, zx_elastic_params* out_elastic,
                                            zx_mc_params* out_mc);

#ifdef __cplusplus
}
#endif

#endif /* ZX_AUTHORING_H */
