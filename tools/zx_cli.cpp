/*!
 * \file zx_cli.cpp
 * \brief Minimal CLI harness to query validation proxies (inclined plane, column collapse).
 */

#include "zx/zx_validation.h"
#include "zx/zx_presets.h"
#include "zx/zx_contact_ref.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include "zx/zx_integration.h"
#include "zx/zx_telemetry.h"
#include "zx/zx_residency.h"
#include "zx/zx_lod.h"
#include "zx/zx_determinism.h"

static void print_usage() {
    std::printf("ziXn CLI\n");
    std::printf("Usage:\n");
    std::printf("  zx_cli --mode inclined --phi <deg>\n");
    std::printf("  zx_cli --mode collapse --phi <deg> --ar <H_over_L0>\n");
    std::printf("  zx_cli --mode anisotropy --mulong <mu0> --mulat <mu1> --vn <vn> --v0 <v_t0> --v1 <v_t1>\n");
    std::printf("  zx_cli --mode presets [--out <path>]\n");
    std::printf("  zx_cli --mode scene --scene {dambreak|bogging|puddle} [--telemetry-out <path>] [--deterministic 0|1] [--seed <u64>] [--fallback on|off|auto] [--fallback-thresholds <tiles> <ms> <enter> <exit> <blend>]\n");
}

int main(int argc, char** argv) {
    const char* mode = nullptr;
    float phi_deg = 30.0f;
    float ar = 1.0f;
    float mu_long = 0.8f, mu_lat = 0.4f, vn = 1.0f, v0 = 1.0f, v1 = 0.0f;
    const char* out_path = nullptr;
    const char* telemetry_out = nullptr;
    const char* scene_name = nullptr;
    int deterministic = 0; unsigned long long seed = 0ULL;
    const char* fallback_mode = "off"; // off|on|auto
    zx_lod_fallback_policy fbp{1000000u, 1e9f, 2u, 2u, 3u};
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--mode") == 0 && i+1 < argc) { mode = argv[++i]; }
        else if (std::strcmp(argv[i], "--phi") == 0 && i+1 < argc) { phi_deg = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--ar") == 0 && i+1 < argc) { ar = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--mulong") == 0 && i+1 < argc) { mu_long = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--mulat") == 0 && i+1 < argc) { mu_lat = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--vn") == 0 && i+1 < argc) { vn = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--v0") == 0 && i+1 < argc) { v0 = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--v1") == 0 && i+1 < argc) { v1 = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--out") == 0 && i+1 < argc) { out_path = argv[++i]; }
        else if (std::strcmp(argv[i], "--help") == 0) { print_usage(); return 0; }
        else if (std::strcmp(argv[i], "--telemetry-out") == 0 && i+1 < argc) { telemetry_out = argv[++i]; }
        else if (std::strcmp(argv[i], "--scene") == 0 && i+1 < argc) { scene_name = argv[++i]; }
        else if (std::strcmp(argv[i], "--deterministic") == 0 && i+1 < argc) { deterministic = std::atoi(argv[++i]); }
        else if (std::strcmp(argv[i], "--seed") == 0 && i+1 < argc) { seed = strtoull(argv[++i], nullptr, 10); }
        else if (std::strcmp(argv[i], "--fallback") == 0 && i+1 < argc) { fallback_mode = argv[++i]; }
        else if (std::strcmp(argv[i], "--fallback-thresholds") == 0 && i+5 < argc) {
            fbp.active_tiles_max = (uint32_t)std::strtoul(argv[++i], nullptr, 10);
            fbp.step_ms_max = (float)std::atof(argv[++i]);
            fbp.enter_frames = (uint32_t)std::strtoul(argv[++i], nullptr, 10);
            fbp.exit_frames  = (uint32_t)std::strtoul(argv[++i], nullptr, 10);
            fbp.blend_frames = (uint32_t)std::strtoul(argv[++i], nullptr, 10);
        }
    }
    if (!mode) { print_usage(); return 1; }

    if (std::strcmp(mode, "inclined") == 0) {
        zx_mc_params mc{ phi_deg, 0.0f };
        float theta_c = zx_validation_inclined_plane_theta_c(&mc);
        std::printf("inclined_plane: phi=%.3f deg -> theta_c=%.3f deg\n", phi_deg, theta_c*180.0/3.14159265358979323846);
        return 0;
    } else if (std::strcmp(mode, "collapse") == 0) {
        float ratio = zx_validation_column_collapse_runout_ratio(phi_deg, ar);
        std::printf("column_collapse: phi=%.3f deg, ar=%.3f -> L/H=%.3f\n", phi_deg, ar, ratio);
        return 0;
    } else if (std::strcmp(mode, "anisotropy") == 0) {
        float n[3]={0,1,0}; float t0[3]={1,0,0}; float t1[3]={0,0,1};
        float vin[3]={v0, vn, v1}, vout[3]; int sat=0;
        zx_contact_project_aniso(vin, n, t0, t1, mu_long, mu_lat, /*rr*/0.5f, /*phi*/-0.01f, /*kappa*/1e-3f, vout, &sat);
        std::printf("anisotropy: mulong=%.3f mulat=%.3f vn=%.3f v0=%.3f v1=%.3f -> out=(%.3f,%.3f,%.3f) sat=%d\n",
                    mu_long, mu_lat, vn, v0, v1, vout[0], vout[1], vout[2], sat);
        return 0;
    } else if (std::strcmp(mode, "presets") == 0) {
        FILE* f = stdout;
        if (out_path) {
            f = std::fopen(out_path, "wb");
            if (!f) { std::perror("fopen"); return 2; }
        }
        std::fprintf(f, "{\n  \"presets\": [\n");
        for (int i=0;i<zx_preset_count();++i) {
            const char* nm = zx_preset_name(i);
            zx_elastic_params ep{}; zx_mc_params mc{}; zx_norsand_params ns{}; zx_norsand_state st{};
            zx_preset_get(nm, &ep, &mc, &ns, &st);
            std::fprintf(f,
                "    {\n      \"name\": \"%s\",\n      \"elastic\": {\"E\": %.6g, \"nu\": %.6g},\n      \"mc\": {\"phi_deg\": %.3f, \"cohesion_kpa\": %.3f},\n      \"norsand\": {\"M\": %.3f, \"lambda_cs\": %.3f, \"kappa\": %.3f, \"p_ref\": %.3g, \"e_ref\": %.3f, \"n_exp\": %.3f, \"dilatancy_scale\": %.3f, \"void_ratio_e\": %.3f}\n    }%s\n",
                nm, ep.young_E, ep.poisson_nu, mc.friction_deg, mc.cohesion_kpa,
                ns.M, ns.lambda_cs, ns.kappa, ns.p_ref, ns.e_ref, ns.n_exp, ns.dilatancy_scale, st.void_ratio_e,
                (i+1<zx_preset_count()?",":""));
        }
        std::fprintf(f, "  ]\n}\n");
        if (out_path) std::fclose(f);
        return 0;
    } else if (std::strcmp(mode, "scene") == 0) {
        if (!scene_name) { print_usage(); return 1; }
        if (deterministic) { zx_set_determinism(1); }
        if (seed) { zx_seed_rng((uint64_t)seed); }
        // Apply LOD global policy and enablement
        zx_lod_set_default_policy(&fbp);
        zx_lod_set_enabled((std::strcmp(fallback_mode, "on")==0) || (std::strcmp(fallback_mode, "auto")==0));
        zx_telemetry* telem = zx_telemetry_create(1024);
        auto t_scene_begin = std::chrono::steady_clock::now();
        if (std::strcmp(scene_name, "dambreak") == 0) {
            zx_dambreak_params p{}; p.tiles=1; p.h=1.0f; p.dt=0.05f; p.steps=20; p.bed_k=2; p.init_head_u=1.0f; p.gamma_dot=1.0f;
            p.hbp = zx_hbp_params{1.0f,0.0f,1.0f,0.0f,1.0f}; p.permeability=1e-10f; p.k_min=1e-12f; p.mu_min=0.01f; p.mu_max=100.0f;
            p.beta_min=1e2f; p.beta_max=1e12f; p.policy=ZX_MU_CLAMP_HARD; p.softness_k=1.0f;
            zx_telemetry_begin_step(telem, "dambreak", 0);
            auto t_fb_begin = std::chrono::steady_clock::now();
            zx_dambreak_metrics m = zx_integration_dambreak_run(&p);
            auto t_fb_end = std::chrono::steady_clock::now();
            zx_telemetry_set_counter(telem, "front_x", m.front_x);
            zx_telemetry_set_counter(telem, "kinetic_j", m.kinetic_j);
            double fb_ms = std::chrono::duration<double, std::milli>(t_fb_end - t_fb_begin).count();
            zx_telemetry_set_counter(telem, "fallback_eval_ms", (float)fb_ms);
            // Residency metrics and optional fallback counters
            {
                zx_residency_opts ro{1,2,1}; zx_residency* rez = zx_residency_create(&ro);
                zx_lod_fallback_state fbs; zx_lod_fallback_init(&fbs);
                uint32_t active_frames=0, churn_en=0, churn_ex=0, pf_sum=0, active_max=0, active_sum=0;
                for (uint32_t s=0; s<p.steps; ++s){
                    uint32_t en=0, ex=0, pf=0; zx_residency_tick(rez,(int)s,0,0,0,&en,&ex,&pf);
                    uint32_t active=zx_residency_get_active_count(rez);
                    churn_en+=en; churn_ex+=ex; pf_sum+=pf; if (active>active_max) active_max=active; active_sum+=active;
                    if (zx_lod_is_enabled()) { if (zx_lod_fallback_update(&fbp, active, 1.0f, &fbs)) active_frames++; }
                }
                float active_mean = (p.steps>0)? ((float)active_sum/(float)p.steps) : 0.0f;
                zx_telemetry_set_counter(telem, "residency_active_max", (float)active_max);
                zx_telemetry_set_counter(telem, "residency_active_mean", active_mean);
                zx_telemetry_set_counter(telem, "residency_churn_enter", (float)churn_en);
                zx_telemetry_set_counter(telem, "residency_churn_exit", (float)churn_ex);
                zx_telemetry_set_counter(telem, "residency_prefetch_sum", (float)pf_sum);
                if (zx_lod_is_enabled()){
                    zx_telemetry_set_counter(telem, "fallback_activations", (float)fbs.activations);
                    zx_telemetry_set_counter(telem, "fallback_active_frames", (float)active_frames);
                    zx_telemetry_set_counter(telem, "fallback_blend_frames", (float)fbp.blend_frames);
                }
                zx_residency_destroy(rez);
            }
            zx_telemetry_end_step(telem);
        } else if (std::strcmp(scene_name, "bogging") == 0) {
            zx_bogging_params p{}; p.h=1.0f; p.dt=0.05f; p.steps=20; p.wheel_radius_nodes=4; p.wheel_pull_u=1.0f; p.wheel_push_w=0.1f;
            p.hbp = zx_hbp_params{1.0f,0.0f,1.0f,0.0f,1.0f}; p.mu_min=0.01f; p.mu_max=100.0f; p.permeability=1e-10f; p.k_min=1e-12f;
            p.beta_min=1e2f; p.beta_max=1e12f; p.gamma_dot=1.0f; p.policy=ZX_MU_CLAMP_HARD; p.softness_k=1.0f;
            zx_telemetry_begin_step(telem, "bogging", 0);
            auto t_fb_begin = std::chrono::steady_clock::now();
            zx_bogging_metrics m = zx_integration_wheel_bogging_run(&p);
            auto t_fb_end = std::chrono::steady_clock::now();
            zx_telemetry_set_counter(telem, "drag_N", m.drag_N);
            zx_telemetry_set_counter(telem, "sink_depth_m", m.sink_depth_m);
            double fb_ms = std::chrono::duration<double, std::milli>(t_fb_end - t_fb_begin).count();
            zx_telemetry_set_counter(telem, "fallback_eval_ms", (float)fb_ms);
            {
                zx_residency_opts ro{1,2,1}; zx_residency* rez = zx_residency_create(&ro);
                zx_lod_fallback_state fbs; zx_lod_fallback_init(&fbs);
                uint32_t active_frames=0, churn_en=0, churn_ex=0, pf_sum=0, active_max=0, active_sum=0;
                for (uint32_t s=0; s<p.steps; ++s){ uint32_t en=0, ex=0, pf=0; zx_residency_tick(rez,(int)s,0,0,0,&en,&ex,&pf);
                    uint32_t active=zx_residency_get_active_count(rez); churn_en+=en; churn_ex+=ex; pf_sum+=pf; if (active>active_max) active_max=active; active_sum+=active;
                    if (zx_lod_is_enabled()) { if (zx_lod_fallback_update(&fbp, active, 1.0f, &fbs)) active_frames++; } }
                float active_mean = (p.steps>0)? ((float)active_sum/(float)p.steps) : 0.0f;
                zx_telemetry_set_counter(telem, "residency_active_max", (float)active_max);
                zx_telemetry_set_counter(telem, "residency_active_mean", active_mean);
                zx_telemetry_set_counter(telem, "residency_churn_enter", (float)churn_en);
                zx_telemetry_set_counter(telem, "residency_churn_exit", (float)churn_ex);
                zx_telemetry_set_counter(telem, "residency_prefetch_sum", (float)pf_sum);
                if (zx_lod_is_enabled()){
                    zx_telemetry_set_counter(telem, "fallback_activations", (float)fbs.activations);
                    zx_telemetry_set_counter(telem, "fallback_active_frames", (float)active_frames);
                    zx_telemetry_set_counter(telem, "fallback_blend_frames", (float)fbp.blend_frames);
                }
                zx_residency_destroy(rez);
            }
            zx_telemetry_end_step(telem);
        } else if (std::strcmp(scene_name, "puddle") == 0) {
            zx_puddle_params p{}; p.h=1.0f; p.dt=0.05f; p.steps=20; p.init_head_u=1.0f;
            p.hbp = zx_hbp_params{1.0f,0.0f,1.0f,0.0f,1.0f}; p.mu_min=0.01f; p.mu_max=100.0f; p.permeability=1e-10f; p.k_min=1e-12f;
            p.beta_min=1e2f; p.beta_max=1e12f; p.gamma_dot=1.0f; p.policy=ZX_MU_CLAMP_HARD; p.softness_k=1.0f;
            zx_telemetry_begin_step(telem, "puddle", 0);
            auto t_fb_begin = std::chrono::steady_clock::now();
            zx_puddle_metrics m = zx_integration_puddle_creep_run(&p);
            auto t_fb_end = std::chrono::steady_clock::now();
            zx_telemetry_set_counter(telem, "creep_dist_x", m.creep_dist_x);
            double fb_ms = std::chrono::duration<double, std::milli>(t_fb_end - t_fb_begin).count();
            zx_telemetry_set_counter(telem, "fallback_eval_ms", (float)fb_ms);
            {
                zx_residency_opts ro{1,2,1}; zx_residency* rez = zx_residency_create(&ro);
                zx_lod_fallback_state fbs; zx_lod_fallback_init(&fbs);
                uint32_t active_frames=0, churn_en=0, churn_ex=0, pf_sum=0, active_max=0, active_sum=0;
                for (uint32_t s=0; s<p.steps; ++s){ uint32_t en=0, ex=0, pf=0; zx_residency_tick(rez,(int)s,0,0,0,&en,&ex,&pf);
                    uint32_t active=zx_residency_get_active_count(rez); churn_en+=en; churn_ex+=ex; pf_sum+=pf; if (active>active_max) active_max=active; active_sum+=active;
                    if (zx_lod_is_enabled()) { if (zx_lod_fallback_update(&fbp, active, 1.0f, &fbs)) active_frames++; } }
                float active_mean = (p.steps>0)? ((float)active_sum/(float)p.steps) : 0.0f;
                zx_telemetry_set_counter(telem, "residency_active_max", (float)active_max);
                zx_telemetry_set_counter(telem, "residency_active_mean", active_mean);
                zx_telemetry_set_counter(telem, "residency_churn_enter", (float)churn_en);
                zx_telemetry_set_counter(telem, "residency_churn_exit", (float)churn_ex);
                zx_telemetry_set_counter(telem, "residency_prefetch_sum", (float)pf_sum);
                if (zx_lod_is_enabled()){
                    zx_telemetry_set_counter(telem, "fallback_activations", (float)fbs.activations);
                    zx_telemetry_set_counter(telem, "fallback_active_frames", (float)active_frames);
                    zx_telemetry_set_counter(telem, "fallback_blend_frames", (float)fbp.blend_frames);
                }
                zx_residency_destroy(rez);
            }
            zx_telemetry_end_step(telem);
        } else { print_usage(); return 1; }
        auto t_scene_end = std::chrono::steady_clock::now();
        double scene_ms = std::chrono::duration<double, std::milli>(t_scene_end - t_scene_begin).count();
        zx_telemetry_set_counter(telem, "scene_ms", (float)scene_ms);
        if (telemetry_out) { zx_telemetry_export_json(telem, telemetry_out); }
        zx_telemetry_destroy(telem);
        return 0;
    } else {
        print_usage();
        return 1;
    }
}


