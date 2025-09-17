/*!
\file zx_presets.h
\brief Authoring presets for materials: sands, gravel, soils, clays, snow.
\author Colin Macritchie (Ripple Group, LLC)
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_PRESETS_H
#define ZX_PRESETS_H

#include <stdint.h>
#include "zx_abi.h"
#include "zx_constitutive_ref.h"

/** \brief Returns number of available presets. */
ZX_API int ZX_CALL zx_preset_count(void);

/** \brief Returns preset name by index (stable). */
ZX_API const char* ZX_CALL zx_preset_name(int index);

/** \brief Populate physics parameters for a preset by name.
 * @param name Preset identifier
 * @param out_elastic Out elastic params (may be NULL)
 * @param out_mc Out Mohr–Coulomb params (may be NULL)
 * @param out_ns_params Out NorSand params (may be NULL)
 * @param out_ns_state Out NorSand state (may be NULL)
 * @return 1 on success, 0 otherwise
 */
ZX_API int ZX_CALL zx_preset_get(
    const char* name,
    zx_elastic_params* out_elastic,
    zx_mc_params* out_mc,
    zx_norsand_params* out_ns_params,
    zx_norsand_state* out_ns_state);

#endif /* ZX_PRESETS_H */


