/*!
 * \file telemetry_depth_test.cpp
 * \brief Integration: run scenes via CLI to JSON and assert required keys exist.
 */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

static bool file_contains(const char* path, const char* needle){
    FILE* f = std::fopen(path, "rb"); if (!f) return false;
    std::fseek(f, 0, SEEK_END); long n = std::ftell(f); std::fseek(f, 0, SEEK_SET);
    if (n <= 0) { std::fclose(f); return false; }
    char* buf = (char*)std::malloc((size_t)n+1); if (!buf) { std::fclose(f); return false; }
    size_t rd = std::fread(buf, 1, (size_t)n, f); buf[rd] = '\0'; std::fclose(f);
    bool ok = std::strstr(buf, needle) != nullptr; std::free(buf); return ok;
}

static void run_scene_and_assert(const char* scene){
    const char* out = "wsl-build/telem_depth.json";
    // Use fallback 'on' to ensure fallback_* counters are emitted
    char cmd[512];
    std::snprintf(cmd, sizeof(cmd), "./zx_cli --mode scene --scene %s --deterministic 1 --seed 42 --fallback on --telemetry-out %s > /dev/null 2>&1", scene, out);
    int rc = std::system(cmd);
    assert(rc == 0);
    // Required keys
    assert(file_contains(out, "\"particle_count\""));
    assert(file_contains(out, "\"active_tiles\""));
    assert(file_contains(out, "\"residency_active_max\""));
    assert(file_contains(out, "\"residency_churn_enter\""));
    assert(file_contains(out, "\"residency_prefetch_sum\""));
    // Fallback keys (since fallback on)
    assert(file_contains(out, "\"fallback_activations\""));
    assert(file_contains(out, "\"fallback_active_frames\""));
    assert(file_contains(out, "\"fallback_blend_frames\""));
}

int main(){
    run_scene_and_assert("dambreak");
    run_scene_and_assert("bogging");
    run_scene_and_assert("puddle");
    return 0;
}


