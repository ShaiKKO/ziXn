/*!
 * \file zx_presets.cpp
 * \brief Authoring presets mapped to solver parameters.
 */

#include "zx/zx_presets.h"
#include <cstring>

struct PresetRow {
    const char* name;
    zx_elastic_params elastic;
    zx_mc_params mc;
    zx_norsand_params ns;
    zx_norsand_state ns_state;
};

static const PresetRow kPresets[] = {
    { "Fine Sand", { 1.5e7f, 0.3f }, { 32.0f, 0.5f }, { 1.2f, 0.18f, 0.05f, 1.0e4f, 0.85f, 1.0f, 0.6f }, { 0.90f } },
    { "Coarse Sand", { 2.0e7f, 0.28f }, { 35.0f, 0.8f }, { 1.3f, 0.16f, 0.05f, 1.0e4f, 0.80f, 1.0f, 0.6f }, { 0.88f } },
    { "Gravel", { 5.0e7f, 0.25f }, { 42.0f, 0.5f }, { 1.5f, 0.12f, 0.04f, 1.0e4f, 0.75f, 1.0f, 0.5f }, { 0.75f } },
    { "Loam", { 1.0e7f, 0.30f }, { 28.0f, 2.0f }, { 1.1f, 0.22f, 0.06f, 1.0e4f, 0.95f, 1.0f, 0.7f }, { 1.0f } },
    { "Clay", { 2.0e7f, 0.32f }, { 20.0f, 50.0f }, { 0.9f, 0.25f, 0.08f, 1.0e4f, 1.10f, 1.0f, 0.4f }, { 1.10f } },
    { "Powder Snow", { 0.8e6f, 0.25f }, { 20.0f, 1.0f }, { 0.6f, 0.30f, 0.09f, 1.0e4f, 1.20f, 1.0f, 0.3f }, { 1.20f } },
    { "Packed Snow", { 5.0e6f, 0.30f }, { 25.0f, 3.0f }, { 0.8f, 0.25f, 0.07f, 1.0e4f, 1.00f, 1.0f, 0.5f }, { 1.00f } }
};

int zx_preset_count(void) { return (int)(sizeof(kPresets)/sizeof(kPresets[0])); }

const char* zx_preset_name(int index)
{
    if (index < 0 || index >= zx_preset_count()) return "";
    return kPresets[index].name;
}

int zx_preset_get(const char* name, zx_elastic_params* out_elastic, zx_mc_params* out_mc, zx_norsand_params* out_ns_params, zx_norsand_state* out_ns_state)
{
    for (int i=0;i<zx_preset_count();++i) {
        if (std::strcmp(name, kPresets[i].name) == 0) {
            if (out_elastic) *out_elastic = kPresets[i].elastic;
            if (out_mc) *out_mc = kPresets[i].mc;
            if (out_ns_params) *out_ns_params = kPresets[i].ns;
            if (out_ns_state) *out_ns_state = kPresets[i].ns_state;
            return 1;
        }
    }
    return 0;
}


