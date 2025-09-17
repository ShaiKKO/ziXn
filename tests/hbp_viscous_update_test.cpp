/*!
 * \file hbp_viscous_update_test.cpp
 * \brief Unit test for coarse-grid viscous diffusion step.
 */

#include "zx/zx_hbp.h"
#include "zx/zx_tiles.h"
#include "zx/zx_tiles_api.h"
#include <cassert>
#include <cstdlib>

static void* zx_alloc(size_t n, void*)
{
  return std::malloc(n);
}
static void zx_free(void* p, void*)
{
  std::free(p);
}

static inline uint32_t nidx(uint32_t ix, uint32_t iy, uint32_t iz)
{
  return ix + (uint32_t) ZX_TILE_B * (iy + (uint32_t) ZX_TILE_B * iz);
}

int main()
{
  zx_tile* pool = zx_create_tile_pool(1, zx_alloc, nullptr);
  assert(pool);
  const uint32_t B = (uint32_t) ZX_TILE_B;
  const uint32_t N = B * B * B;
  for (uint32_t i = 0; i < N; ++i)
  {
    pool[0].nodes[i].mass  = 1.0f;
    pool[0].nodes[i].mom_x = 0.0f;
    pool[0].nodes[i].mom_y = 0.0f;
    pool[0].nodes[i].mom_z = 0.0f;
  }

  const uint32_t cx = B / 2, cy = B / 2, cz = B / 2;
  pool[0].nodes[nidx(cx, cy, cz)].mom_x = 1.0f;  // unit velocity at center in x

  zx_hbp_params p;
  p.mu0   = 1.0f;
  p.K     = 0.0f;
  p.n     = 1.0f;
  p.tau_y = 0.0f;
  p.m     = 1.0f;
  float h = 1.0f, dt = 0.05f;  // stable for nu=1 -> dt<~0.166
  zx_hbp_update_coarse_grid(pool, 1, h, dt, &p, 1e-5f, 1e3f);

  float vc_after = pool[0].nodes[nidx(cx, cy, cz)].mom_x / pool[0].nodes[nidx(cx, cy, cz)].mass;
  assert(vc_after < 1.0f);
  float vxp = pool[0].nodes[nidx(cx + 1, cy, cz)].mom_x;
  assert(vxp > 0.0f);
  float vxm = pool[0].nodes[nidx(cx - 1, cy, cz)].mom_x;
  assert(vxm > 0.0f);
  float vyp = pool[0].nodes[nidx(cx, cy + 1, cz)].mom_x;
  assert(vyp > 0.0f);
  float vym = pool[0].nodes[nidx(cx, cy - 1, cz)].mom_x;
  assert(vym > 0.0f);
  float vzp = pool[0].nodes[nidx(cx, cy, cz + 1)].mom_x;
  assert(vzp > 0.0f);
  float vzm = pool[0].nodes[nidx(cx, cy, cz - 1)].mom_x;
  assert(vzm > 0.0f);

  zx_destroy_tile_pool(pool, zx_free, nullptr);
  return 0;
}
