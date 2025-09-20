/*!\
 * \file lod_mip_random_test.cpp\
 * \brief Random-texture mip chain test to check border crack-free property over cycles.\
 */

#include "zx/zx_lod.h"
#include <vector>
#include <cassert>
#include <cstdint>

static uint32_t lcg(uint32_t& s)
{
  s = s * 1664525U + 1013904223U;
  return s;
}

int main()
{
  const uint32_t W = 16, H = 16;  // larger tile to stress borders
  for (uint32_t seed = 1; seed <= 5; ++seed)
  {
    uint32_t rng = seed;
    std::vector<float> A(W * H), B(W * H);
    for (uint32_t y = 0; y < H; ++y)
    {
      for (uint32_t x = 0; x < W; ++x)
      {
        float r      = (float) (lcg(rng) & 0xFFFF) / 65535.0f;
        A[y * W + x] = r;
      }
    }
    // Copy right edge of A to left edge of B to enforce continuity
    for (uint32_t y = 0; y < H; ++y)
    {
      B[y * W + 0] = A[y * W + (W - 1)];
    }
    // Fill remainder of B with random values
    for (uint32_t y = 0; y < H; ++y)
    {
      for (uint32_t x = 1; x < W; ++x)
      {
        float r      = (float) (lcg(rng) & 0xFFFF) / 65535.0f;
        B[y * W + x] = r;
      }
    }

    // Build two levels of mip and test borders at each level
    std::vector<float> A1((W / 2) * (H / 2)), B1((W / 2) * (H / 2));
    zx_lod_downsample_2x(A.data(), W, H, W, A1.data(), W / 2, H / 2, W / 2);
    zx_lod_downsample_2x(B.data(), W, H, W, B1.data(), W / 2, H / 2, W / 2);
    float d0 = zx_lod_border_consistency_check(A1.data(), W / 2, H / 2, W / 2, B1.data(), W / 2,
                                               H / 2, W / 2, 0);
    assert(d0 < 1e-3f);

    std::vector<float> A2((W / 4) * (H / 4)), B2((W / 4) * (H / 4));
    zx_lod_downsample_2x(A1.data(), W / 2, H / 2, W / 2, A2.data(), W / 4, H / 4, W / 4);
    zx_lod_downsample_2x(B1.data(), W / 2, H / 2, W / 2, B2.data(), W / 4, H / 4, W / 4);
    float d1 = zx_lod_border_consistency_check(A2.data(), W / 4, H / 4, W / 4, B2.data(), W / 4,
                                               H / 4, W / 4, 0);
    assert(d1 < 1e-3f);
  }
  return 0;
}
