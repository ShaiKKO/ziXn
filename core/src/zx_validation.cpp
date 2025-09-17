/*!
 * \file zx_validation.cpp
 * \brief Simple analytical proxies for validation harnesses.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_validation.h"
#include <cmath>
#include <algorithm>

float zx_validation_inclined_plane_theta_c(const zx_mc_params* mc)
{
    // θ_c ≈ φ for cohesionless; if cohesion present, this is a lower bound.
    return mc->friction_deg * 3.14159265358979323846f / 180.0f;
}

float zx_validation_column_collapse_runout_ratio(float friction_deg, float ar)
{
    // Empirical monotone: larger φ -> shorter runout; larger aspect ratio -> longer.
    // Use a simple proxy L/H ≈ a * ar^b / tan(φ), with clamps for stability.
    const float a = 1.2f, b = 0.8f;
    const float phi = friction_deg * 3.14159265358979323846f / 180.0f;
    const float t = std::max(0.2f, std::tan(std::max(5.0f*3.14159265358979323846f/180.0f, phi)));
    const float ratio = a * std::pow(std::max(0.2f, ar), b) / t;
    return ratio;
}


