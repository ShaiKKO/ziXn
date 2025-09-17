/*!
 * \file zx_authoring.cpp
 * \brief UI mapping and clamps for presets.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_authoring.h"
#include <algorithm>

static float clamp01(float x){ return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x); }

void zx_authoring_map_sand(const zx_ui_sand* ui, zx_elastic_params* ep, zx_mc_params* mc, zx_norsand_params* ns, zx_norsand_state* st)
{
    const float dryness = clamp01(ui->dryness);
    const float grain   = clamp01(ui->grain);
    // Elastic
    ep->young_E = 1.5e7f + 3.5e7f * grain; // stiffer for coarse
    ep->poisson_nu = 0.28f - 0.06f * dryness;
    // MC
    mc->friction_deg = 28.0f + 14.0f * grain + 6.0f * dryness;
    mc->cohesion_kpa = (1.0f - dryness) * 3.0f;
    // NorSand
    ns->M = 1.1f + 0.4f * grain;
    ns->lambda_cs = 0.18f - 0.06f * grain;
    ns->kappa = 0.05f;
    ns->p_ref = 1.0e4f; ns->e_ref = 0.85f - 0.1f * grain; ns->n_exp = 1.0f;
    ns->dilatancy_scale = 0.5f + 0.3f * dryness;
    st->void_ratio_e = 0.9f - 0.1f * grain;
}

void zx_authoring_map_snow(const zx_ui_snow* ui, zx_elastic_params* ep, zx_mc_params* mc)
{
    const float t = clamp01(ui->temperature);
    const float p = clamp01(ui->packing);
    ep->young_E = 0.8e6f + 5.0e6f * p;
    ep->poisson_nu = 0.25f + 0.05f * p;
    mc->friction_deg = 20.0f + 10.0f * p - 5.0f * t;
    mc->cohesion_kpa = 1.0f + 4.0f * p - 2.0f * t;
}


