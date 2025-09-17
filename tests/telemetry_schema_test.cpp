/*!
 * \file telemetry_schema_test.cpp
 * \brief Telemetry schema and counter correctness tests.
 */

#include "zx/zx_telemetry.h"
#include <cassert>
#include <cstdio>

int main(){
    zx_telemetry* t = zx_telemetry_create(8);
    zx_telemetry_begin_step(t, "dambreak", 0);
    zx_telemetry_set_counter(t, "front_x", 1.0f);
    zx_telemetry_set_counter(t, "kinetic_j", 2.5f);
    zx_telemetry_set_counter(t, "particle_count", 0.0f);
    zx_telemetry_set_counter(t, "active_tiles", 1.0f);
    zx_telemetry_set_counter(t, "fallback_activations", 0.0f);
    zx_telemetry_set_counter(t, "fallback_active_frames", 0.0f);
    zx_telemetry_set_counter(t, "fallback_blend_frames", 0.0f);
    zx_telemetry_end_step(t);
    zx_telemetry_add_error(t, "dambreak", 0, "E_SAMPLE", "sample error triage");
    int rc1 = zx_telemetry_export_csv(t, "telemetry_test.csv");
    int rc2 = zx_telemetry_export_json(t, "telemetry_test.json");
    assert(rc1==0 && rc2==0);
    zx_telemetry_destroy(t);
    return 0;
}


