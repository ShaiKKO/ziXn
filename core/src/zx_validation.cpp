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
  constexpr float ZX_PI        = 3.14159265358979323846F;  // NOLINT(readability-magic-numbers)
  constexpr float ZX_DEG2RAD   = ZX_PI / 180.0F;           // NOLINT(readability-magic-numbers)
  constexpr float ZX_TAN_FLOOR = 0.2F;                     // NOLINT(readability-magic-numbers)
  constexpr float ZX_MIN_DEG   = 5.0F;                     // NOLINT(readability-magic-numbers)
}  // namespace

/**
 * @brief Compute the critical slope angle theta_c (radians) from material friction.
 *
 * Converts the material friction angle stored in mc->friction_deg (degrees) to radians.
 * For cohesionless materials θ_c ≈ φ; when cohesion is present this value should be
 * treated as a lower bound.
 *
 * @param mc Pointer to zx_mc_params containing the field `friction_deg` (degrees).
 *           Caller must ensure `mc` is non-null.
 * @return float theta_c in radians.
 */
float zx_validation_inclined_plane_theta_c(const zx_mc_params* mc)
{
  // θ_c ≈ φ for cohesionless; if cohesion present, this is a lower bound.
  return mc->friction_deg * ZX_DEG2RAD;
}

/**
 * @brief Estimate collapse runout ratio (L/H) for a granular column.
 *
 * Computes an empirical proxy for the runout ratio L/H as a monotone function
 * of internal friction angle and column aspect ratio:
 * L/H ≈ a * ar^b / tan(phi), with clamping for numerical stability.
 *
 * friction_deg is interpreted in degrees. ar is the initial column aspect
 * ratio (height/width). To avoid extreme values the function clamps:
 * - phi is floored at 5° when computing tan(phi),
 * - tan(...) is floored at 0.2,
 * - ar is floored at 0.2 before exponentiation.
 *
 * @param friction_deg Friction angle in degrees (φ).
 * @param ar Initial aspect ratio (height/width).
 * @return float Estimated runout ratio L/H.
 */
float zx_validation_column_collapse_runout_ratio(
    float friction_deg, float ar)  // NOLINT(bugprone-easily-swappable-parameters)
{
  // Empirical monotone: larger φ -> shorter runout; larger aspect ratio -> longer.
  // Use a simple proxy L/H ≈ a * ar^b / tan(φ), with clamps for stability.
  const float a     = 1.2F;
  const float b     = 0.8F;
  const float phi   = friction_deg * ZX_DEG2RAD;
  const float t     = std::max(ZX_TAN_FLOOR, std::tan(std::max(ZX_MIN_DEG * ZX_DEG2RAD, phi)));
  const float ratio = a * std::pow(std::max(ZX_TAN_FLOOR, ar), b) / t;
  return ratio;
}
