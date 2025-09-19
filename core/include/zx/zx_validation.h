/*!
\file zx_validation.h
\brief Analytical/empirical validation utilities: inclined plane and column collapse proxies.
\author Colin Macritchie (Ripple Group, LLC)
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_VALIDATION_H
#define ZX_VALIDATION_H

#include "zx_abi.h"
#include "zx_constitutive_ref.h"
#include <stdint.h>

/** \brief Critical angle for onset of sliding on an inclined plane for dry granular MC.
 * For MC with friction φ and negligible cohesion, θ_c ≈ φ. Returns radians.
 */
ZX_API float ZX_CALL zx_validation_inclined_plane_theta_c(const zx_mc_params* mc);

/** \brief Column collapse runout proxy L/H for aspect ratio ar and friction.
 * Provides a monotonic estimate to be refined by empirical fit during calibration.
 */
ZX_API float ZX_CALL zx_validation_column_collapse_runout_ratio(float friction_deg,
                                                                float aspect_ratio);

#endif /* ZX_VALIDATION_H */
