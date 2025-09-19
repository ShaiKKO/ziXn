/*!
 * \file zx_validation.cpp
 * \brief Simple analytical proxies for validation harnesses.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_validation.h"
#include <algorithm>
#include <cmath>

namespace
{
  constexpr float zx_pi        = 3.14159265358979323846F;
  constexpr float zx_deg2rad   = zx_pi / 180.0F;
  constexpr float zx_tan_floor = 0.2F;
  constexpr float zx_min_deg   = 5.0F;
}  // namespace

/**
 * @brief Compute the critical slope angle theta_c (radians) from material friction.
 * @param mc Input Mohr–Coulomb parameters; uses friction_deg (degrees).
 * @return float theta_c in radians (lower bound when cohesion present).
 */
float zx_validation_inclined_plane_theta_c(const zx_mc_params* mc)
{
  // θ_c ≈ φ for cohesionless; if cohesion present, this is a lower bound.
  return mc->friction_deg * zx_deg2rad;
}

/**
 * @brief Estimate collapse runout ratio (L/H) for a granular column.
 * @param friction_deg Friction angle in degrees.
 * @param ar Initial aspect ratio (height/width).
 * @return float Estimated runout ratio L/H using a clamped proxy.
 */
float zx_validation_column_collapse_runout_ratio(float friction_deg, float ar)
{
  // Empirical monotone: larger φ -> shorter runout; larger aspect ratio -> longer.
  // Use a simple proxy L/H ≈ a * ar^b / tan(φ), with clamps for stability.
  const float a     = 1.2F;
  const float b     = 0.8F;
  const float phi   = friction_deg * zx_deg2rad;
  const float t     = std::max(zx_tan_floor, std::tan(std::max(zx_min_deg * zx_deg2rad, phi)));
  const float ratio = a * std::pow(std::max(zx_tan_floor, ar), b) / t;
  return ratio;
}
