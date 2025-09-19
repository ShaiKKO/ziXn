/*!
 * \file zx_cli.cpp
 * \brief Minimal CLI harness to query validation proxies (inclined plane, column collapse).
 */

#include "zx/zx_contact_ref.h"
#include "zx/zx_determinism.h"
#include "zx/zx_integration.h"
#include "zx/zx_lod.h"
#include "zx/zx_presets.h"
#include "zx/zx_residency.h"
#include "zx/zx_telemetry.h"
#include "zx/zx_validation.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

/**
 * @brief Print command-line usage and available modes for the ziXn CLI.
 *
 * Prints usage patterns and option summaries for all supported modes:
 * "inclined", "collapse", "anisotropy", "presets", "scene", and "bench"
 * (including both bench-lod and bench-up2 variants), along with their
 * expected command-line arguments and optional flags.
 */
static void print_usage()
{
  const char* lines[] = {
      "ziXn CLI\n",
      "Usage:\n",
      "  zx_cli --mode inclined --phi <deg>\n",
      "  zx_cli --mode collapse --phi <deg> --ar <H_over_L0>\n",
      "  zx_cli --mode anisotropy --mulong <mu0> --mulat <mu1> --vn <vn> --v0 <v_t0> --v1 <v_t1>\n",
      "  zx_cli --mode presets [--out <path>]\n",
      "  zx_cli --mode scene --scene {dambreak|bogging|puddle} [--telemetry-out <path>] "
      "[--deterministic 0|1] [--seed <u64>] [--fallback on|off|auto] [--fallback-thresholds "
      "<tiles> <ms> <enter> <exit> <blend>]\n",
      "  zx_cli --mode bench --bench-lod [--size N] [--iters K] [--simd auto|scalar|avx2]\n",
      "  zx_cli --mode bench --bench-up2  [--size N] [--iters K]     (runs scalar and AVX2)\n",
  };
  for (const char* s : lines)
  {
    std::fputs(s, stdout);
  }
}

/**
 * @brief Command-line entry point for the ZX validation and benchmarking CLI.
 *
 * Parses command-line options and dispatches one of several modes:
 * - "inclined": compute critical angle for an inclined-plane validation.
 * - "collapse": compute column-collapse runout ratio.
 * - "anisotropy": project a velocity vector using anisotropic contact model.
 * - "presets": emit JSON describing built-in presets (optionally to a file).
 * - "bench": run microbenchmarks for LOD downsample/upsample and border checks.
 * - "scene": run a small predefined integration scene (dambreak, bogging, puddle)
 *   and optionally export telemetry JSON.
 *
 * Recognized flags include (--mode, --phi, --ar, --mulong, --mulat, --vn, --v0, --v1,
 * --out, --telemetry-out, --scene, --deterministic, --seed, --fallback,
 * --fallback-thresholds, --bench-lod, --bench-up2, --size, --iters, --simd, --help).
 *
 * Behavior notes:
 * - When running "presets" with --out, the program attempts to open the named file;
 *   failure to open the file results in an error (exit code 2).
 * - "bench" requires either --bench-lod or --bench-up2; a too-small --size (<16)
 *   prints an error and exits with code 2.
 * - "scene" requires --scene and supports optional determinism/seed, LOD fallback
 *   policy, and telemetry export via --telemetry-out.
 *
 * Return:
 * - 0 on successful completion of the selected mode.
 * - 1 for missing/invalid usage or unrecognized mode.
 * - 2 for filesystem or benchmark size errors (e.g., fopen failure, size too small).
 */
int main(int argc, char** argv)
{
  const char* mode          = nullptr;
  float phi_deg             = 30.0F;
  float ar                  = 1.0F;
  float mu_long             = 0.8F;
  float mu_lat              = 0.4F;
  float vn                  = 1.0F;
  float v0                  = 1.0F;
  float v1                  = 0.0F;
  const char* out_path      = nullptr;
  const char* telemetry_out = nullptr;
  const char* scene_name    = nullptr;
  int deterministic         = 0;
  unsigned long long seed   = 0ULL;
  const char* fallback_mode = "off";  // off|on|auto
  zx_lod_fallback_policy fbp{1000000U, 1.0e9F, 2U, 2U, 3U};
  int prefetch_rings_cli = -1;
  int pin_args[6];
  int have_pin          = 0;
  int bench_lod         = 0;
  int bench_up2         = 0;
  int bench_size        = 1024;
  int bench_iters       = 100;
  const char* simd_mode = "auto";
  for (int i = 1; i < argc; ++i)
  {
    if (std::strcmp(argv[i], "--mode") == 0 && i + 1 < argc)
    {
      mode = argv[++i];
    }
    else if (std::strcmp(argv[i], "--phi") == 0 && i + 1 < argc)
    {
      phi_deg = static_cast<float>(std::atof(argv[++i]));
    }
    else if (std::strcmp(argv[i], "--ar") == 0 && i + 1 < argc)
    {
      ar = static_cast<float>(std::atof(argv[++i]));
    }
    else if (std::strcmp(argv[i], "--mulong") == 0 && i + 1 < argc)
    {
      mu_long = static_cast<float>(std::atof(argv[++i]));
    }
    else if (std::strcmp(argv[i], "--mulat") == 0 && i + 1 < argc)
    {
      mu_lat = static_cast<float>(std::atof(argv[++i]));
    }
    else if (std::strcmp(argv[i], "--vn") == 0 && i + 1 < argc)
    {
      vn = static_cast<float>(std::atof(argv[++i]));
    }
    else if (std::strcmp(argv[i], "--v0") == 0 && i + 1 < argc)
    {
      v0 = static_cast<float>(std::atof(argv[++i]));
    }
    else if (std::strcmp(argv[i], "--v1") == 0 && i + 1 < argc)
    {
      v1 = static_cast<float>(std::atof(argv[++i]));
    }
    else if (std::strcmp(argv[i], "--out") == 0 && i + 1 < argc)
    {
      out_path = argv[++i];
    }
    else if (std::strcmp(argv[i], "--help") == 0)
    {
      print_usage();
      return 0;
    }
    else if (std::strcmp(argv[i], "--telemetry-out") == 0 && i + 1 < argc)
    {
      telemetry_out = argv[++i];
    }
    else if (std::strcmp(argv[i], "--scene") == 0 && i + 1 < argc)
    {
      scene_name = argv[++i];
    }
    else if (std::strcmp(argv[i], "--deterministic") == 0 && i + 1 < argc)
    {
      deterministic = std::atoi(argv[++i]);
    }
    else if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc)
    {
      seed = strtoull(argv[++i], nullptr, 10);
    }
    else if (std::strcmp(argv[i], "--fallback") == 0 && i + 1 < argc)
    {
      fallback_mode = argv[++i];
    }
    else if (std::strcmp(argv[i], "--fallback-thresholds") == 0 && i + 5 < argc)
    {
      fbp.active_tiles_max = (uint32_t) std::strtoul(argv[++i], nullptr, 10);
      fbp.step_ms_max      = (float) std::atof(argv[++i]);
      fbp.enter_frames     = (uint32_t) std::strtoul(argv[++i], nullptr, 10);
      fbp.exit_frames      = (uint32_t) std::strtoul(argv[++i], nullptr, 10);
      fbp.blend_frames     = (uint32_t) std::strtoul(argv[++i], nullptr, 10);
    }
    else if (std::strcmp(argv[i], "--bench-lod") == 0)
    {
      bench_lod = 1;
    }
    else if (std::strcmp(argv[i], "--bench-up2") == 0)
    {
      bench_up2 = 1;
    }
    else if (std::strcmp(argv[i], "--size") == 0 && i + 1 < argc)
    {
      bench_size = std::atoi(argv[++i]);
    }
    else if (std::strcmp(argv[i], "--iters") == 0 && i + 1 < argc)
    {
      bench_iters = std::atoi(argv[++i]);
    }
    else if (std::strcmp(argv[i], "--simd") == 0 && i + 1 < argc)
    {
      simd_mode = argv[++i];
    }
    else if (std::strcmp(argv[i], "--prefetch-rings") == 0 && i + 1 < argc)
    {
      prefetch_rings_cli = std::atoi(argv[++i]);
    }
    else if (std::strcmp(argv[i], "--pin-box") == 0 && i + 6 < argc)
    {
      for (int k = 0; k < 6; ++k)
      {
        pin_args[k] = std::atoi(argv[++i]);
      }
      have_pin = 1;
    }
  }
  if (mode == nullptr)
  {
    print_usage();
    return 1;
  }

  if (std::strcmp(mode, "inclined") == 0)
  {
    zx_mc_params mc{phi_deg, 0.0F};
    float theta_c = zx_validation_inclined_plane_theta_c(&mc);
    {
      char buf[256];
      std::snprintf(buf, sizeof(buf), "inclined_plane: phi=%.3f deg -> theta_c=%.3f deg\n", phi_deg,
                    theta_c * 180.0 / 3.14159265358979323846);
      std::fputs(buf, stdout);
    }
    return 0;
  }
  else if (std::strcmp(mode, "collapse") == 0)
  {
    float ratio = zx_validation_column_collapse_runout_ratio(phi_deg, ar);
    {
      char buf[256];
      std::snprintf(buf, sizeof(buf), "column_collapse: phi=%.3f deg, ar=%.3f -> L/H=%.3f\n",
                    phi_deg, ar, ratio);
      std::fputs(buf, stdout);
    }
    return 0;
  }
  else if (std::strcmp(mode, "anisotropy") == 0)
  {
    float n[3]   = {0, 1, 0};
    float t0[3]  = {1, 0, 0};
    float t1[3]  = {0, 0, 1};
    float vin[3] = {v0, vn, v1}, vout[3];
    int sat      = 0;
    zx_contact_project_aniso(vin, n, t0, t1, mu_long, mu_lat, /*rr*/ 0.5F, /*phi*/ -0.01F,
                             /*kappa*/ 1e-3F, vout, &sat);
    {
      char buf[256];
      std::snprintf(buf, sizeof(buf),
                    "anisotropy: mulong=%.3f mulat=%.3f vn=%.3f v0=%.3f v1=%.3f -> "
                    "out=(%.3f,%.3f,%.3f) sat=%d\n",
                    mu_long, mu_lat, vn, v0, v1, vout[0], vout[1], vout[2], sat);
      std::fputs(buf, stdout);
    }
    return 0;
  }
  else if (std::strcmp(mode, "presets") == 0)
  {
    FILE* f = stdout;
    if (out_path != nullptr)
    {
      f = std::fopen(out_path, "wb");
      if (f == nullptr)
      {
        std::perror("fopen");
        return 2;
      }
    }
    std::fputs("{\n  \"presets\": [\n", f);
    for (int i = 0; i < zx_preset_count(); ++i)
    {
      const char* nm = zx_preset_name(i);
      zx_elastic_params ep{};
      zx_mc_params mc{};
      zx_norsand_params ns{};
      zx_norsand_state st{};
      zx_preset_get(nm, &ep, &mc, &ns, &st);
      char buf[512];
      std::snprintf(
          buf, sizeof(buf),
          "    {\n      \"name\": \"%s\",\n      \"elastic\": {\"E\": %.6g, \"nu\": %.6g},\n      "
          "\"mc\": {\"phi_deg\": %.3f, \"cohesion_kpa\": %.3f},\n      \"norsand\": {\"M\": %.3f, "
          "\"lambda_cs\": %.3f, \"kappa\": %.3f, \"p_ref\": %.3g, \"e_ref\": %.3f, \"n_exp\": "
          "%.3f, \"dilatancy_scale\": %.3f, \"void_ratio_e\": %.3f}\n    }%s\n",
          nm, ep.young_E, ep.poisson_nu, mc.friction_deg, mc.cohesion_kpa, ns.M, ns.lambda_cs,
          ns.kappa, ns.p_ref, ns.e_ref, ns.n_exp, ns.dilatancy_scale, st.void_ratio_e,
          (i + 1 < zx_preset_count() ? "," : ""));
      std::fputs(buf, f);
    }
    std::fputs("  ]\n}\n", f);
    if (out_path != nullptr)
    {
      std::fclose(f);
    }
    return 0;
  }
  else if (std::strcmp(mode, "bench") == 0)
  {
    if ((bench_lod == 0) && (bench_up2 == 0))
    {
      print_usage();
      return 1;
    }
    if (std::strcmp(simd_mode, "scalar") == 0)
    {
      zx_lod_set_simd_override(0);
    }
    else if (std::strcmp(simd_mode, "avx2") == 0)
    {
      zx_lod_set_simd_override(2);
    }
    else
    {
      zx_lod_set_simd_override(-1);
    }
    const int w = bench_size;
    const int h = bench_size;
    if ((w < 16) || (h < 16))
    {
      std::fputs("size too small\n", stderr);
      return 2;
    }
    std::vector<float> a((size_t) w * h);
    std::vector<float> b((size_t) w * h);
    for (int y = 0; y < h; ++y)
    {
      for (int x = 0; x < w; ++x)
      {
        a[(size_t) y * w + x] = static_cast<float>(((x * 13 + y * 7) & 255)) * 0.01F;
      }
    }
    for (int y = 0; y < h; ++y)
    {
      b[(size_t) y * w + 0] = a[(size_t) y * w + (w - 1)];
      for (int x = 1; x < w; ++x)
      {
        b[(size_t) y * w + x] = static_cast<float>(((x * 3 + y * 5) & 255)) * 0.02F;
      }
    }
    int dw = w / 2;
    int dh = h / 2;
    std::vector<float> d((size_t) dw * dh);
    if (bench_up2)
    {
      // Prepare a synthetic downsampled image as input for upsample
      for (int y = 0; y < dh; ++y)
      {
        for (int x = 0; x < dw; ++x)
        {
          d[(size_t) y * dw + x] = static_cast<float>((((x * 17 + y * 11) ^ 0x55) & 255)) * 0.007F;
        }
      }
      std::vector<float> u((size_t) w * h);
      // Warm-up
      zx_lod_set_simd_override(0);
      zx_lod_upsample_2x(d.data(), dw, dh, dw, u.data(), w, h, w);
      zx_lod_set_simd_override(2);
      zx_lod_upsample_2x(d.data(), dw, dh, dw, u.data(), w, h, w);
      // Scalar timing
      zx_lod_set_simd_override(0);
      auto ts0 = std::chrono::steady_clock::now();
      for (int it = 0; it < bench_iters; ++it)
      {
        zx_lod_upsample_2x(d.data(), dw, dh, dw, u.data(), w, h, w);
      }
      auto ts1         = std::chrono::steady_clock::now();
      double ms_scalar = std::chrono::duration<double, std::milli>(ts1 - ts0).count();
      // AVX2 timing (may fall back to scalar if unsupported)
      zx_lod_set_simd_override(2);
      auto tv0 = std::chrono::steady_clock::now();
      for (int it = 0; it < bench_iters; ++it)
      {
        zx_lod_upsample_2x(d.data(), dw, dh, dw, u.data(), w, h, w);
      }
      auto tv1       = std::chrono::steady_clock::now();
      double ms_avx2 = std::chrono::duration<double, std::milli>(tv1 - tv0).count();
      double sumU    = 0.0;
      for (size_t i = 0; i < u.size(); ++i)
      {
        sumU += u[i];
      }
      {
        char buf[256];
        std::snprintf(buf, sizeof(buf), "bench_up2 size=%dx%d iters=%d\n", w, h, bench_iters);
        std::fputs(buf, stdout);
        std::snprintf(buf, sizeof(buf), "upsample_2x scalar: %.3f ms\n", ms_scalar);
        std::fputs(buf, stdout);
        std::snprintf(buf, sizeof(buf), "upsample_2x  avx2*: %.3f ms\n", ms_avx2);
        std::fputs(buf, stdout);
      }
      if (ms_avx2 > 0.0)
      {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "speedup (scalar/avx2): %.2fx\n", ms_scalar / ms_avx2);
        std::fputs(buf, stdout);
      }
      {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "checksum=%.6f\n", sumU);
        std::fputs(buf, stdout);
      }
      return 0;
    }
    auto t0 = std::chrono::steady_clock::now();
    for (int it = 0; it < bench_iters; ++it)
    {
      zx_lod_downsample_2x(a.data(), w, h, w, d.data(), dw, dh, dw);
    }
    auto t1      = std::chrono::steady_clock::now();
    double ms_ds = std::chrono::duration<double, std::milli>(t1 - t0).count();
    float accum  = 0.0F;
    auto t2      = std::chrono::steady_clock::now();
    for (int it = 0; it < bench_iters; ++it)
    {
      accum += zx_lod_border_consistency_check(a.data(), w, h, w, b.data(), w, h, w, 0);
    }
    auto t3      = std::chrono::steady_clock::now();
    double ms_bc = std::chrono::duration<double, std::milli>(t3 - t2).count();
    // Upsample microbench: from dw x dh back to W x H
    std::vector<float> u2((size_t) w * h);
    auto t4 = std::chrono::steady_clock::now();
    for (int it = 0; it < bench_iters; ++it)
    {
      zx_lod_upsample_2x(d.data(), dw, dh, dw, u2.data(), w, h, w);
    }
    auto t5       = std::chrono::steady_clock::now();
    double ms_us  = std::chrono::duration<double, std::milli>(t5 - t4).count();
    double us_sum = 0.0;
    for (size_t i = 0; i < u2.size(); ++i)
    {
      us_sum += u2[i];
    }
    {
      char buf[256];
      std::snprintf(buf, sizeof(buf), "bench_lod size=%dx%d iters=%d simd=%s\n", w, h, bench_iters,
                    simd_mode);
      std::fputs(buf, stdout);
      std::snprintf(buf, sizeof(buf), "downsample_2x: %.3f ms (dst0=%.3f)\n", ms_ds,
                    d.empty() ? 0.0F : d[0]);
      std::fputs(buf, stdout);
      std::snprintf(buf, sizeof(buf), "border_check:  %.3f ms (accum=%.6f)\n", ms_bc, accum);
      std::fputs(buf, stdout);
      std::snprintf(buf, sizeof(buf), "upsample_2x:   %.3f ms (sum=%.3f)\n", ms_us, us_sum);
      std::fputs(buf, stdout);
    }
    return 0;
  }
  else if (std::strcmp(mode, "scene") == 0)
  {
    if (scene_name == nullptr)
    {
      print_usage();
      return 1;
    }
    if (deterministic != 0)
    {
      zx_set_determinism(1);
    }
    if (seed != 0ULL)
    {
      zx_seed_rng((uint64_t) seed);
    }
    // Apply LOD global policy and enablement
    zx_lod_set_default_policy(&fbp);
    zx_lod_set_enabled((std::strcmp(fallback_mode, "on") == 0) ||
                       (std::strcmp(fallback_mode, "auto") == 0));
    zx_telemetry* telem = zx_telemetry_create(1024);
    auto t_scene_begin  = std::chrono::steady_clock::now();
    if (std::strcmp(scene_name, "dambreak") == 0)
    {
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
      auto t_fb_begin       = std::chrono::steady_clock::now();
      zx_dambreak_metrics m = zx_integration_dambreak_run(&p);
      auto t_fb_end         = std::chrono::steady_clock::now();
      zx_telemetry_set_counter(telem, "front_x", m.front_x);
      zx_telemetry_set_counter(telem, "kinetic_j", m.kinetic_j);
      double fb_ms = std::chrono::duration<double, std::milli>(t_fb_end - t_fb_begin).count();
      if (!deterministic)
      {
        zx_telemetry_set_counter(telem, "fallback_eval_ms", (float) fb_ms);
      }
      // Residency metrics and optional fallback counters
      {
        zx_residency_opts ro{1U, 2U, 1U};
        if (prefetch_rings_cli >= 0)
        {
          ro.prefetch_rings = (uint32_t) prefetch_rings_cli;
        }
        zx_residency* rez = zx_residency_create(&ro);
        if (have_pin != 0)
        {
          zx_residency_pin_box(rez, pin_args[0], pin_args[1], pin_args[2], pin_args[3], pin_args[4],
                               pin_args[5]);
        }
        zx_lod_fallback_state fbs;
        zx_lod_fallback_init(&fbs);
        uint32_t active_frames = 0, churn_en = 0, churn_ex = 0, pf_sum = 0, active_max = 0,
                 active_sum = 0;
        for (uint32_t s = 0; s < p.steps; ++s)
        {
          uint32_t en = 0U;
          uint32_t ex = 0U;
          uint32_t pf = 0U;
          zx_residency_tick(rez, (int) s, 0, 0, 0, &en, &ex, &pf);
          uint32_t active = zx_residency_get_active_count(rez);
          churn_en += en;
          churn_ex += ex;
          pf_sum += pf;
          active_max = std::max(active_max, active);
          active_sum += active;
          if (zx_lod_is_enabled() != 0)
          {
            if (zx_lod_fallback_update(&fbp, active, 1.0F, &fbs) != 0)
            {
              active_frames++;
            }
          }
        }
        float active_mean = (p.steps > 0) ? ((float) active_sum / (float) p.steps) : 0.0F;
        zx_telemetry_set_counter(telem, "residency_active_max", (float) active_max);
        zx_telemetry_set_counter(telem, "residency_active_mean", active_mean);
        zx_telemetry_set_counter(telem, "residency_churn_enter", (float) churn_en);
        zx_telemetry_set_counter(telem, "residency_churn_exit", (float) churn_ex);
        zx_telemetry_set_counter(telem, "residency_prefetch_sum", (float) pf_sum);
        // Uniform counters across scenes
        zx_telemetry_set_counter(telem, "particle_count", 0.0F);
        zx_telemetry_set_counter(telem, "active_tiles", (float) active_max);
        if (zx_lod_is_enabled() != 0)
        {
          zx_telemetry_set_counter(telem, "fallback_activations", (float) fbs.activations);
          zx_telemetry_set_counter(telem, "fallback_active_frames", (float) active_frames);
          zx_telemetry_set_counter(telem, "fallback_blend_frames", (float) fbp.blend_frames);
        }
        zx_residency_destroy(rez);
      }
      zx_telemetry_end_step(telem);
    }
    else if (std::strcmp(scene_name, "bogging") == 0)
    {
      zx_bogging_params p{};
      p.h                  = 1.0f;
      p.dt                 = 0.05f;
      p.steps              = 20;
      p.wheel_radius_nodes = 4;
      p.wheel_pull_u       = 1.0f;
      p.wheel_push_w       = 0.1f;
      p.hbp                = zx_hbp_params{1.0f, 0.0f, 1.0f, 0.0f, 1.0f};
      p.mu_min             = 0.01f;
      p.mu_max             = 100.0f;
      p.permeability       = 1e-10f;
      p.k_min              = 1e-12f;
      p.beta_min           = 1e2f;
      p.beta_max           = 1e12f;
      p.gamma_dot          = 1.0f;
      p.policy             = ZX_MU_CLAMP_HARD;
      p.softness_k         = 1.0f;
      zx_telemetry_begin_step(telem, "bogging", 0);
      auto t_fb_begin      = std::chrono::steady_clock::now();
      zx_bogging_metrics m = zx_integration_wheel_bogging_run(&p);
      auto t_fb_end        = std::chrono::steady_clock::now();
      zx_telemetry_set_counter(telem, "drag_N", m.drag_N);
      zx_telemetry_set_counter(telem, "sink_depth_m", m.sink_depth_m);
      double fb_ms = std::chrono::duration<double, std::milli>(t_fb_end - t_fb_begin).count();
      if (!deterministic)
      {
        zx_telemetry_set_counter(telem, "fallback_eval_ms", (float) fb_ms);
      }
      {
        zx_residency_opts ro{1U, 2U, 1U};
        if (prefetch_rings_cli >= 0)
        {
          ro.prefetch_rings = (uint32_t) prefetch_rings_cli;
        }
        zx_residency* rez = zx_residency_create(&ro);
        if (have_pin != 0)
        {
          zx_residency_pin_box(rez, pin_args[0], pin_args[1], pin_args[2], pin_args[3], pin_args[4],
                               pin_args[5]);
        }
        zx_lod_fallback_state fbs;
        zx_lod_fallback_init(&fbs);
        uint32_t active_frames = 0, churn_en = 0, churn_ex = 0, pf_sum = 0, active_max = 0,
                 active_sum = 0;
        for (uint32_t s = 0; s < p.steps; ++s)
        {
          uint32_t en = 0U;
          uint32_t ex = 0U;
          uint32_t pf = 0U;
          zx_residency_tick(rez, (int) s, 0, 0, 0, &en, &ex, &pf);
          uint32_t active = zx_residency_get_active_count(rez);
          churn_en += en;
          churn_ex += ex;
          pf_sum += pf;
          active_max = std::max(active_max, active);
          active_sum += active;
          if (zx_lod_is_enabled() != 0)
          {
            if (zx_lod_fallback_update(&fbp, active, 1.0F, &fbs) != 0)
            {
              active_frames++;
            }
          }
        }
        float active_mean = (p.steps > 0) ? ((float) active_sum / (float) p.steps) : 0.0F;
        zx_telemetry_set_counter(telem, "residency_active_max", (float) active_max);
        zx_telemetry_set_counter(telem, "residency_active_mean", active_mean);
        zx_telemetry_set_counter(telem, "residency_churn_enter", (float) churn_en);
        zx_telemetry_set_counter(telem, "residency_churn_exit", (float) churn_ex);
        zx_telemetry_set_counter(telem, "residency_prefetch_sum", (float) pf_sum);
        zx_telemetry_set_counter(telem, "particle_count", 0.0F);
        zx_telemetry_set_counter(telem, "active_tiles", (float) active_max);
        if (zx_lod_is_enabled() != 0)
        {
          zx_telemetry_set_counter(telem, "fallback_activations", (float) fbs.activations);
          zx_telemetry_set_counter(telem, "fallback_active_frames", (float) active_frames);
          zx_telemetry_set_counter(telem, "fallback_blend_frames", (float) fbp.blend_frames);
        }
        zx_residency_destroy(rez);
      }
      zx_telemetry_end_step(telem);
    }
    else if (std::strcmp(scene_name, "puddle") == 0)
    {
      zx_puddle_params p{};
      p.h            = 1.0f;
      p.dt           = 0.05f;
      p.steps        = 20;
      p.init_head_u  = 1.0f;
      p.hbp          = zx_hbp_params{1.0f, 0.0f, 1.0f, 0.0f, 1.0f};
      p.mu_min       = 0.01f;
      p.mu_max       = 100.0f;
      p.permeability = 1e-10f;
      p.k_min        = 1e-12f;
      p.beta_min     = 1e2f;
      p.beta_max     = 1e12f;
      p.gamma_dot    = 1.0f;
      p.policy       = ZX_MU_CLAMP_HARD;
      p.softness_k   = 1.0f;
      zx_telemetry_begin_step(telem, "puddle", 0);
      auto t_fb_begin     = std::chrono::steady_clock::now();
      zx_puddle_metrics m = zx_integration_puddle_creep_run(&p);
      auto t_fb_end       = std::chrono::steady_clock::now();
      zx_telemetry_set_counter(telem, "creep_dist_x", m.creep_dist_x);
      double fb_ms = std::chrono::duration<double, std::milli>(t_fb_end - t_fb_begin).count();
      if (!deterministic)
      {
        zx_telemetry_set_counter(telem, "fallback_eval_ms", (float) fb_ms);
      }
      {
        zx_residency_opts ro{1U, 2U, 1U};
        if (prefetch_rings_cli >= 0)
        {
          ro.prefetch_rings = (uint32_t) prefetch_rings_cli;
        }
        zx_residency* rez = zx_residency_create(&ro);
        if (have_pin != 0)
        {
          zx_residency_pin_box(rez, pin_args[0], pin_args[1], pin_args[2], pin_args[3], pin_args[4],
                               pin_args[5]);
        }
        zx_lod_fallback_state fbs;
        zx_lod_fallback_init(&fbs);
        uint32_t active_frames = 0, churn_en = 0, churn_ex = 0, pf_sum = 0, active_max = 0,
                 active_sum = 0;
        for (uint32_t s = 0; s < p.steps; ++s)
        {
          uint32_t en = 0U;
          uint32_t ex = 0U;
          uint32_t pf = 0U;
          zx_residency_tick(rez, (int) s, 0, 0, 0, &en, &ex, &pf);
          uint32_t active = zx_residency_get_active_count(rez);
          churn_en += en;
          churn_ex += ex;
          pf_sum += pf;
          active_max = std::max(active_max, active);
          active_sum += active;
          if (zx_lod_is_enabled() != 0)
          {
            if (zx_lod_fallback_update(&fbp, active, 1.0F, &fbs) != 0)
            {
              active_frames++;
            }
          }
        }
        float active_mean = (p.steps > 0) ? ((float) active_sum / (float) p.steps) : 0.0F;
        zx_telemetry_set_counter(telem, "residency_active_max", (float) active_max);
        zx_telemetry_set_counter(telem, "residency_active_mean", active_mean);
        zx_telemetry_set_counter(telem, "residency_churn_enter", (float) churn_en);
        zx_telemetry_set_counter(telem, "residency_churn_exit", (float) churn_ex);
        zx_telemetry_set_counter(telem, "residency_prefetch_sum", (float) pf_sum);
        zx_telemetry_set_counter(telem, "particle_count", 0.0F);
        zx_telemetry_set_counter(telem, "active_tiles", (float) active_max);
        if (zx_lod_is_enabled())
        {
          zx_telemetry_set_counter(telem, "fallback_activations", (float) fbs.activations);
          zx_telemetry_set_counter(telem, "fallback_active_frames", (float) active_frames);
          zx_telemetry_set_counter(telem, "fallback_blend_frames", (float) fbp.blend_frames);
        }
        zx_residency_destroy(rez);
      }
      zx_telemetry_end_step(telem);
    }
    else
    {
      print_usage();
      return 1;
    }
    auto t_scene_end = std::chrono::steady_clock::now();
    double scene_ms =
        std::chrono::duration<double, std::milli>(t_scene_end - t_scene_begin).count();
    if (deterministic == 0)
    {
      zx_telemetry_set_counter(telem, "scene_ms", (float) scene_ms);
    }
    if (telemetry_out != nullptr)
    {
      zx_telemetry_export_json(telem, telemetry_out);
    }
    zx_telemetry_destroy(telem);
    return 0;
  }
  else
  {
    print_usage();
    return 1;
  }
}
