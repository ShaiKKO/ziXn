/*!
\file zx_authoring.h
\brief Authoring UI mappings and safe clamps for material presets.
*/

#ifndef ZX_AUTHORING_H
#define ZX_AUTHORING_H

#include <stdint.h>
#include "zx_abi.h"
#include "zx_constitutive_ref.h"

/* UI models are normalized [0,1] sliders; functions clamp and map to solver parameters. */

typedef struct zx_ui_sand {
    float dryness;   /* 0=very wet, 1=very dry */
    float grain;     /* 0=fine, 1=coarse */
} zx_ui_sand;

typedef struct zx_ui_snow {
    float temperature; /* proxy η: 0=cold/powder, 1=warm/near slush */
    float packing;     /* 0=loose, 1=packed */
} zx_ui_snow;

/* Map sand UI to elastic, MC, and NorSand parameters with safe clamps. */
ZX_API void ZX_CALL zx_authoring_map_sand(
    const zx_ui_sand* ui,
    zx_elastic_params* out_elastic,
    zx_mc_params* out_mc,
    zx_norsand_params* out_ns,
    zx_norsand_state* out_ns_state);

/* Map snow UI to elastic and MC(+cap via cohesion) with safe clamps. */
ZX_API void ZX_CALL zx_authoring_map_snow(
    const zx_ui_snow* ui,
    zx_elastic_params* out_elastic,
    zx_mc_params* out_mc);

#endif /* ZX_AUTHORING_H */


