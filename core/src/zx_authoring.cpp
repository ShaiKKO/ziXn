/*!
 * \file zx_authoring.cpp
 * \brief UI mapping and clamps for presets.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_authoring.h"
#include <algorithm>

// Local named constants to avoid magic numbers and improve readability
namespace
{
  constexpr float k_zero = 0.0F;
  constexpr float k_one  = 1.0F;

  // Sand mapping coefficients
  constexpr float k_sand_young_e_base          = 1.5e7F;
  constexpr float k_sand_young_e_gain_grain    = 3.5e7F;
  constexpr float k_sand_poisson_base          = 0.28F;
  constexpr float k_sand_poisson_slope_dryness = -0.06F;
  constexpr float k_sand_friction_deg_base     = 28.0F;
  constexpr float k_sand_friction_gain_grain   = 14.0F;
  constexpr float k_sand_friction_gain_dryness = 6.0F;
  constexpr float k_sand_cohesion_kpa_base     = 3.0F;
  constexpr float k_sand_ns_m_base             = 1.1F;
  constexpr float k_sand_ns_m_gain_grain       = 0.4F;
  constexpr float k_sand_lambda_cs_base        = 0.18F;
  constexpr float k_sand_lambda_cs_slope_grain = -0.06F;
  constexpr float k_sand_kappa                 = 0.05F;
  constexpr float k_sand_p_ref                 = 1.0e4F;
  constexpr float k_sand_e_ref_base            = 0.85F;
  constexpr float k_sand_e_ref_slope_grain     = -0.1F;
  constexpr float k_sand_dilatancy_base        = 0.5F;
  constexpr float k_sand_dilatancy_gain_dry    = 0.3F;
  constexpr float k_sand_void_ratio_base       = 0.9F;
  constexpr float k_sand_void_ratio_slope_g    = -0.1F;

  // Snow mapping coefficients
  constexpr float k_snow_young_e_base          = 0.8e6F;
  constexpr float k_snow_young_e_gain_packing  = 5.0e6F;
  constexpr float k_snow_poisson_base          = 0.25F;
  constexpr float k_snow_poisson_gain_packing  = 0.05F;
  constexpr float k_snow_friction_deg_base     = 20.0F;
  constexpr float k_snow_friction_gain_packing = 10.0F;
  constexpr float k_snow_friction_slope_temp   = -5.0F;
  constexpr float k_snow_cohesion_kpa_base     = 1.0F;
  constexpr float k_snow_cohesion_gain_packing = 4.0F;
  constexpr float k_snow_cohesion_slope_temp   = -2.0F;
}  // namespace

static float clamp01(float x)
{
  if (x < k_zero)
  {
    return k_zero;
  }
  if (x > k_one)
  {
    return k_one;
  }
  return x;
}

void zx_authoring_map_sand(const zx_ui_sand* ui, zx_elastic_params* ep, zx_mc_params* mc,
                           zx_norsand_params* ns, zx_norsand_state* st)
{
  const float dryness = clamp01(ui->dryness);
  const float grain   = clamp01(ui->grain);
  // Elastic
  ep->young_E    = k_sand_young_e_base + k_sand_young_e_gain_grain * grain;  // stiffer for coarse
  ep->poisson_nu = k_sand_poisson_base + k_sand_poisson_slope_dryness * dryness;
  // MC
  mc->friction_deg = k_sand_friction_deg_base + k_sand_friction_gain_grain * grain +
                     k_sand_friction_gain_dryness * dryness;
  mc->cohesion_kpa = (k_one - dryness) * k_sand_cohesion_kpa_base;
  // NorSand
  ns->M               = k_sand_ns_m_base + k_sand_ns_m_gain_grain * grain;
  ns->lambda_cs       = k_sand_lambda_cs_base + k_sand_lambda_cs_slope_grain * grain;
  ns->kappa           = k_sand_kappa;
  ns->p_ref           = k_sand_p_ref;
  ns->e_ref           = k_sand_e_ref_base + k_sand_e_ref_slope_grain * grain;
  ns->n_exp           = k_one;
  ns->dilatancy_scale = k_sand_dilatancy_base + k_sand_dilatancy_gain_dry * dryness;
  st->void_ratio_e    = k_sand_void_ratio_base + k_sand_void_ratio_slope_g * grain;
}

void zx_authoring_map_snow(const zx_ui_snow* ui, zx_elastic_params* ep, zx_mc_params* mc)
{
  const float t  = clamp01(ui->temperature);
  const float p  = clamp01(ui->packing);
  ep->young_E    = k_snow_young_e_base + k_snow_young_e_gain_packing * p;
  ep->poisson_nu = k_snow_poisson_base + k_snow_poisson_gain_packing * p;
  mc->friction_deg =
      k_snow_friction_deg_base + k_snow_friction_gain_packing * p + k_snow_friction_slope_temp * t;
  mc->cohesion_kpa =
      k_snow_cohesion_kpa_base + k_snow_cohesion_gain_packing * p + k_snow_cohesion_slope_temp * t;
}
