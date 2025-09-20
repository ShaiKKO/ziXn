/*!
 * \file determinism_strict_test.cpp
 * \brief Strict determinism: fixed seed, identical telemetry JSON across runs.
 * \author Colin Macritchie (Ripple Group, LLC)
 */

#include "zx/zx_determinism.h"
#include "zx/zx_integration.h"
#include "zx/zx_telemetry.h"
#include <cassert>
#include <cstdio>

static void run_scene_to_json(const char* path)
{
  zx_set_determinism(1);
  zx_seed_rng(0x0123456789ABCDEFULL);
  zx_telemetry* telem = zx_telemetry_create(1024);
  // mirror CLI scene 'dambreak'
  zx_dambreak_params p{};
  p.tiles        = 1;
  p.h            = 1.0f;
  p.dt           = 0.05f;
  p.steps        = 20;
  p.bed_k        = 2;
  p.init_head_u  = 1.0f;
  p.gamma_dot    = 1.0f;
  p.hbp          = zx_hbp_params{1.0f, 0.0f, 1.0f, 0.0f, 1.0f};
  p.permeability = 1e-10f;
  p.k_min        = 1e-12f;
  p.mu_min       = 0.01f;
  p.mu_max       = 100.0f;
  p.beta_min     = 1e2f;
  p.beta_max     = 1e12f;
  p.policy       = ZX_MU_CLAMP_HARD;
  p.softness_k   = 1.0f;
  zx_telemetry_begin_step(telem, "dambreak", 0);
  zx_dambreak_metrics m = zx_integration_dambreak_run(&p);
  zx_telemetry_set_counter(telem, "front_x", m.front_x);
  zx_telemetry_set_counter(telem, "kinetic_j", m.kinetic_j);
  zx_telemetry_end_step(telem);
  int rc = zx_telemetry_export_json(telem, path);
  assert(rc == 0);
  zx_telemetry_destroy(telem);
}

static bool files_equal(const char* a, const char* b)
{
  FILE* fa = std::fopen(a, "rb");
  if (!fa)
    return false;
  FILE* fb = std::fopen(b, "rb");
  if (!fb)
  {
    std::fclose(fa);
    return false;
  }
  int ca = 0, cb = 0;
  bool ok = true;
  while (ok)
  {
    ca = std::fgetc(fa);
    cb = std::fgetc(fb);
    if (ca != cb)
    {
      ok = false;
      break;
    }
    if (ca == EOF || cb == EOF)
      break;
  }
  std::fclose(fa);
  std::fclose(fb);
  return ok;
}

int main()
{
  const char* A = "det_a.json";
  const char* B = "det_b.json";
  run_scene_to_json(A);
  run_scene_to_json(B);
  assert(files_equal(A, B));
  // Validate seed remained as set
  assert(zx_get_rng_seed() == 0x0123456789ABCDEFULL);
  return 0;
}
