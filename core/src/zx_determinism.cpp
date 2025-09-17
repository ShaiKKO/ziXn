/*!
 * \file zx_determinism.cpp
 * \brief Global determinism controls (CPU reference path).
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_determinism.h"
#include <atomic>

static std::atomic<int> g_det{0};
static std::atomic<uint64_t> g_seed{0x9e3779b97f4a7c15ULL};

extern "C"
{

  void ZX_CALL zx_set_determinism(int enable)
  {
    g_det.store(enable ? 1 : 0, std::memory_order_relaxed);
  }
  int ZX_CALL zx_get_determinism(void)
  {
    return g_det.load(std::memory_order_relaxed);
  }
  void ZX_CALL zx_seed_rng(uint64_t seed)
  {
    if (seed == 0)
      seed = 1;
    g_seed.store(seed, std::memory_order_relaxed);
  }
  uint64_t ZX_CALL zx_get_rng_seed(void)
  {
    return g_seed.load(std::memory_order_relaxed);
  }
}
