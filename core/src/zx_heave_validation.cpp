/**
 * @file zx_heave_validation.cpp
 * @brief Heave profile proxy under localized loading.
 * @details Provides simple analytical/proxy routines for heave computation used by tests.
 *          Not part of performance-critical paths. Thread-safe and pure.
 */

#include "zx/zx_heave_validation.h"
#include <algorithm>
#include <cmath>

void zx_heave_profile_compute(const zx_heave_params* hp, float extent, float* out_y,
                              uint32_t samples)
{
  if ((hp == nullptr) || (out_y == nullptr) || (samples == 0U))
  {
    return;
  }
  const float dx    = hp->grid_dx;
  const float half  = extent;
  const float sigma = 0.5F * hp->footprint_width;  // spread
  const float a     = (hp->load_newton * hp->dilatancy) / (std::max(1.0e-3F, hp->footprint_width));
  for (uint32_t i = 0; i < samples; ++i)
  {
    const float x = (-half) + (static_cast<float>(i) * dx);
    const float w = std::exp(-(x * x) / (2.0F * sigma * sigma));
    out_y[i]      = (a * w) * dx;  // proxy displacement
  }
}

float zx_heave_berm_volume(const float* y, uint32_t samples, float dx)
{
  if ((y == nullptr) || (samples < 2U))
  {
    return 0.0F;
  }
  float v = 0.0F;
  for (uint32_t i = 1; i < samples; ++i)
  {
    const float a = std::max(0.0F, y[i - 1]);
    const float b = std::max(0.0F, y[i]);
    v += 0.5F * (a + b) * dx;
  }
  return v;
}
