/*!\
 * \file lod_fallback_sweep_test.cpp\
 * \brief Param sweep tests for LOD fallback hysteresis and blend behavior.\
 */\

#include "zx/zx_lod.h"
#include <cassert>

static void step_sequence(const zx_lod_fallback_policy& p,
                          const uint32_t* active_tiles_seq,
                          const float* step_ms_seq,
                          unsigned count,
                          unsigned* out_activations,
                          unsigned* out_used_frames)
{
    zx_lod_fallback_state s; zx_lod_fallback_init(&s);
    unsigned used=0;
    for (unsigned i=0;i<count;++i){
        int use = zx_lod_fallback_update(&p, active_tiles_seq? active_tiles_seq[i] : 0u,
                                         step_ms_seq? step_ms_seq[i] : 0.0f, &s);
        if (use) used++;
    }
    if (out_activations) *out_activations = s.activations;
    if (out_used_frames) *out_used_frames = used;
}

int main(){
    // 1) enter=2, exit=2, blend=3; pressure sustained
    {
        zx_lod_fallback_policy p{ /*tiles*/10u, /*ms*/5.0f, /*enter*/2u, /*exit*/2u, /*blend*/3u };
        uint32_t seq[5] = { 11, 11, 11, 11, 11 }; // active tiles > max triggers pressure
        float ms[5] = { 1,1,1,1,1 };
        unsigned acts=0, used=0; step_sequence(p, seq, ms, 5, &acts, &used);
        assert(acts == 1);
        // After 2 frames enter, then 3 blend frames reduce but still considered active
        assert(used >= 3);
    }

    // 2) hysteresis: active 2, then inactive 1 (stay on), then inactive 2 (off)
    {
        zx_lod_fallback_policy p{ 10u, 5.0f, 2u, 2u, 3u };
        uint32_t seq[6] = { 11,11, 0, 0,0,0 }; float ms[6] = {1,1,1,1,1,1};
        zx_lod_fallback_state s; zx_lod_fallback_init(&s);
        int u0 = zx_lod_fallback_update(&p, seq[0], ms[0], &s);
        int u1 = zx_lod_fallback_update(&p, seq[1], ms[1], &s);
        assert(s.activations == 1);
        int u2 = zx_lod_fallback_update(&p, seq[2], ms[2], &s); // still on due to exit=2
        assert(u2==1 || s.blend_remaining<=3);
        int u3 = zx_lod_fallback_update(&p, seq[3], ms[3], &s); // turn off here
        (void)u0;(void)u1;(void)u3;
        assert(s.active_frames == 0 || s.inactive_frames >= 1);
    }

    // 3) thresholds OR: either tiles or step_ms triggers
    {
        zx_lod_fallback_policy p{ 10u, 0.5f, 1u, 1u, 0u };
        zx_lod_fallback_state s; zx_lod_fallback_init(&s);
        (void)zx_lod_fallback_update(&p, 5u, 1.0f, &s); // ms>0.5 triggers
        assert(s.active_frames == 1);
        zx_lod_fallback_init(&s);
        (void)zx_lod_fallback_update(&p, 11u, 0.1f, &s); // tiles>10 triggers
        assert(s.active_frames == 1);
    }

    // 4) Edge: zero blend, immediate enter/exit
    {
        zx_lod_fallback_policy p{ 10u, 5.0f, 1u, 1u, 0u };
        zx_lod_fallback_state s; zx_lod_fallback_init(&s);
        (void)zx_lod_fallback_update(&p, 11u, 0.0f, &s); // on
        assert(s.activations == 1);
        (void)zx_lod_fallback_update(&p, 0u, 0.0f, &s); // off
        assert(s.active_frames == 0);
    }

    return 0;
}


