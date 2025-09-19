/*!
 * \file zx_heave_validation.cpp
 * \brief Heave profile proxy under localized loading.
 */

#include "zx/zx_heave_validation.h"
#include <algorithm>
#include <cmath>

void zx_heave_profile_compute(const zx_heave_params* hp, float extent, float* out_y,
                              uint32_t samples)
{
  if (!hp || !out_y || samples == 0)
    return;
  float dx    = hp->grid_dx;
  float half  = extent;
  float sigma = 0.5f * hp->footprint_width;  // spread
  float A     = (hp->load_newton * hp->dilatancy) / (std::max(1e-3f, hp->footprint_width));
  for (uint32_t i = 0; i < samples; ++i)
  {
    float x  = -half + i * dx;
    float w  = std::exp(-(x * x) / (2.0f * sigma * sigma));
    out_y[i] = (A * w) * dx;  // proxy displacement
  }
}

float zx_heave_berm_volume(const float* y, uint32_t samples, float dx)
{
  if (!y || samples < 2)
    return 0.0f;
  float V = 0.0f;
  for (uint32_t i = 1; i < samples; ++i)
  {
    float a = std::max(0.0f, y[i - 1]);
    float b = std::max(0.0f, y[i]);
    V += 0.5f * (a + b) * dx;
  }
  return V;
}
