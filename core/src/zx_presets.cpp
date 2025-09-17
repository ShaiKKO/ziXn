/*!
 * \file zx_presets.cpp
 * \brief Authoring presets mapped to solver parameters.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_presets.h"
#include <array>
#include <cstring>

struct PresetRow
{
  const char* name;
  zx_elastic_params elastic;
  zx_mc_params mc;
  zx_norsand_params ns;
  zx_norsand_state ns_state;
};

static const std::array<PresetRow, 7> kPresets = {{{"Fine Sand",
                                                    {1.5e7f, 0.3f},
                                                    {32.0f, 0.5f},
                                                    {1.2f, 0.18f, 0.05f, 1.0e4f, 0.85f, 1.0f, 0.6f},
                                                    {0.90f}},
                                                   {"Coarse Sand",
                                                    {2.0e7f, 0.28f},
                                                    {35.0f, 0.8f},
                                                    {1.3f, 0.16f, 0.05f, 1.0e4f, 0.80f, 1.0f, 0.6f},
                                                    {0.88f}},
                                                   {"Gravel",
                                                    {5.0e7f, 0.25f},
                                                    {42.0f, 0.5f},
                                                    {1.5f, 0.12f, 0.04f, 1.0e4f, 0.75f, 1.0f, 0.5f},
                                                    {0.75f}},
                                                   {"Loam",
                                                    {1.0e7f, 0.30f},
                                                    {28.0f, 2.0f},
                                                    {1.1f, 0.22f, 0.06f, 1.0e4f, 0.95f, 1.0f, 0.7f},
                                                    {1.0f}},
                                                   {"Clay",
                                                    {2.0e7f, 0.32f},
                                                    {20.0f, 50.0f},
                                                    {0.9f, 0.25f, 0.08f, 1.0e4f, 1.10f, 1.0f, 0.4f},
                                                    {1.10f}},
                                                   {"Powder Snow",
                                                    {0.8e6f, 0.25f},
                                                    {20.0f, 1.0f},
                                                    {0.6f, 0.30f, 0.09f, 1.0e4f, 1.20f, 1.0f, 0.3f},
                                                    {1.20f}},
                                                   {"Packed Snow",
                                                    {5.0e6f, 0.30f},
                                                    {25.0f, 3.0f},
                                                    {0.8f, 0.25f, 0.07f, 1.0e4f, 1.00f, 1.0f, 0.5f},
                                                    {1.00f}}}};

/** \brief Number of available presets. */
int zx_preset_count(void)
{
  return static_cast<int>(kPresets.size());
}

/**
 * @brief Return the stable name of a preset by index.
 *
 * Returns a pointer to the internal, null-terminated preset name for the
 * given index. The returned string is owned by the library and remains valid
 * for the lifetime of the program; do not attempt to free it.
 *
 * @param index Preset index (0 .. zx_preset_count()-1). If out of range, an empty
 *              string ("") is returned.
 * @return const char* Pointer to the preset name or an empty string on invalid index.
 */
const char* zx_preset_name(int index)
{
  if (index < 0 || index >= zx_preset_count())
    return "";
  return kPresets.at(static_cast<size_t>(index)).name;
}

/**
 * @brief Populate solver parameter structures from a named authoring preset.
 *
 * Looks up a preset by exact, case-sensitive name and copies its parameter
 * values into any non-null output pointers.
 *
 * @param name Preset identifier (must not be NULL). Matching is exact and case-sensitive.
 * @param out_elastic If non-null, receives the preset's elastic parameters.
 * @param out_mc If non-null, receives the preset's Mohr–Coulomb parameters.
 * @param out_ns_params If non-null, receives the preset's NorSand parameters.
 * @param out_ns_state If non-null, receives the preset's NorSand state.
 * @return int 1 if a preset with the given name was found and values copied, 0 otherwise.
 */
int zx_preset_get(const char* name, zx_elastic_params* out_elastic, zx_mc_params* out_mc,
                  zx_norsand_params* out_ns_params, zx_norsand_state* out_ns_state)
{
  for (int i = 0; i < zx_preset_count(); ++i)
  {
    if (std::strcmp(name, kPresets.at(static_cast<size_t>(i)).name) == 0)
    {
      if (out_elastic)
        *out_elastic = kPresets.at(static_cast<size_t>(i)).elastic;
      if (out_mc)
        *out_mc = kPresets.at(static_cast<size_t>(i)).mc;
      if (out_ns_params)
        *out_ns_params = kPresets.at(static_cast<size_t>(i)).ns;
      if (out_ns_state)
        *out_ns_state = kPresets.at(static_cast<size_t>(i)).ns_state;
      return 1;
    }
  }
  return 0;
}
