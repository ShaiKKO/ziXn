/*!
\file zx_presets.h
\brief Authoring presets for materials: sands, gravel, soils, clays, snow.
*/

#ifndef ZX_PRESETS_H
#define ZX_PRESETS_H

#include <stdint.h>
#include "zx_abi.h"
#include "zx_constitutive_ref.h"

/* Returns number of available presets. */
ZX_API int ZX_CALL zx_preset_count(void);

/* Returns preset name by index (stable). */
ZX_API const char* ZX_CALL zx_preset_name(int index);

/* Populate physics parameters for a preset by name. Returns 1 on success, 0 otherwise.
 * Any of the out pointers may be NULL if not needed. NorSand is provided for sands; others may be NULL.
 */
ZX_API int ZX_CALL zx_preset_get(
    const char* name,
    zx_elastic_params* out_elastic,
    zx_mc_params* out_mc,
    zx_norsand_params* out_ns_params,
    zx_norsand_state* out_ns_state);

#endif /* ZX_PRESETS_H */


