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

static const std::array<PresetRow, 7> k_presets = {
    {{"Fine Sand",
      {1.5e7F, 0.30F},
      {32.0F, 0.50F},
      {1.2F, 0.18F, 0.05F, 1.0e4F, 0.85F, 1.0F, 0.6F},
      {0.90F}},
     {"Coarse Sand",
      {2.0e7F, 0.28F},
      {35.0F, 0.80F},
      {1.3F, 0.16F, 0.05F, 1.0e4F, 0.80F, 1.0F, 0.6F},
      {0.88F}},
     {"Gravel",
      {5.0e7F, 0.25F},
      {42.0F, 0.50F},
      {1.5F, 0.12F, 0.04F, 1.0e4F, 0.75F, 1.0F, 0.5F},
      {0.75F}},
     {"Loam",
      {1.0e7F, 0.30F},
      {28.0F, 2.0F},
      {1.1F, 0.22F, 0.06F, 1.0e4F, 0.95F, 1.0F, 0.7F},
      {1.0F}},
     {"Clay",
      {2.0e7F, 0.32F},
      {20.0F, 50.0F},
      {0.9F, 0.25F, 0.08F, 1.0e4F, 1.10F, 1.0F, 0.4F},
      {1.10F}},
     {"Powder Snow",
      {0.8e6F, 0.25F},
      {20.0F, 1.0F},
      {0.6F, 0.30F, 0.09F, 1.0e4F, 1.20F, 1.0F, 0.3F},
      {1.20F}},
     {"Packed Snow",
      {5.0e6F, 0.30F},
      {25.0F, 3.0F},
      {0.8F, 0.25F, 0.07F, 1.0e4F, 1.00F, 1.0F, 0.5F},
      {1.00F}}}};

/** \brief Number of available presets. */
int zx_preset_count()
{
  return static_cast<int>(k_presets.size());
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
  if ((index < 0) || (index >= zx_preset_count()))
  {
    return "";
  }
  return k_presets.at(static_cast<size_t>(index)).name;
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
    if (std::strcmp(name, k_presets.at(static_cast<size_t>(i)).name) == 0)
    {
      if (out_elastic != nullptr)
      {
        *out_elastic = k_presets.at(static_cast<size_t>(i)).elastic;
      }
      if (out_mc != nullptr)
      {
        *out_mc = k_presets.at(static_cast<size_t>(i)).mc;
      }
      if (out_ns_params != nullptr)
      {
        *out_ns_params = k_presets.at(static_cast<size_t>(i)).ns;
      }
      if (out_ns_state != nullptr)
      {
        *out_ns_state = k_presets.at(static_cast<size_t>(i)).ns_state;
      }
      return 1;
    }
  }
  return 0;
}
