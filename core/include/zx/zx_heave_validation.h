/*!
\file zx_heave_validation.h
\brief Heave/berm profile validation utilities under wheel/foot loading.
*/

#ifndef ZX_HEAVE_VALIDATION_H
#define ZX_HEAVE_VALIDATION_H

#include "zx_abi.h"
#include <cstdint>

typedef struct ZxHeaveParams
{
  float load_newton;     /* applied vertical load */
  float footprint_width; /* contact width (m) */
  float grid_dx;         /* sampling spacing (m) */
  float dilatancy;       /* effective dilatancy scalar (>0 increases heave) */
} zx_heave_params;

/* Compute a 1D heave profile y(x) sampled on [-extent, extent] with step grid_dx.
 * The model distributes volumetric plastic strain around the footprint and integrates
 * it to a surface displacement using a symmetric kernel.
 * out_y must have length samples = floor(2*extent/grid_dx)+1
 */
ZX_API void ZX_CALL zx_heave_profile_compute(const zx_heave_params* hp, float extent, float* out_y,
                                             uint32_t samples);

/* Integrate berm (heave) volume by trapezoidal rule over the positive part of y. */
ZX_API float ZX_CALL zx_heave_berm_volume(const float* y, uint32_t samples, float dx);

#endif /* ZX_HEAVE_VALIDATION_H */
