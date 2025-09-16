/*!
 * \file zx_residency.cpp
 * \brief Tile residency manager with hysteresis and prefetch hints.
 */

#include "zx/zx_residency.h"
#include <unordered_map>
#include <vector>
#include <algorithm>

struct TileState { int hot_frames=0; int cold_frames=0; bool active=false; };

struct zx_residency {
    zx_residency_opts opts;
    std::unordered_map<long long, TileState> states;
    uint32_t active_count=0;
};

static inline long long key(int x,int y,int z){ return ( ( (long long)(unsigned int)x) << 42) ^ ( ((long long)(unsigned int)y) << 21) ^ (long long)(unsigned int)z; }

extern "C" {

ZX_API zx_residency* ZX_CALL zx_residency_create(const zx_residency_opts* opts){
    zx_residency* r = new zx_residency(); r->opts = *opts; return r;
}

ZX_API void ZX_CALL zx_residency_destroy(zx_residency* ctx){ delete ctx; }

ZX_API void ZX_CALL zx_residency_tick(zx_residency* ctx,
                                      int cx, int cy, int cz,
                                      uint32_t active_radius,
                                      uint32_t* enters,
                                      uint32_t* exits,
                                      uint32_t* prefetch_count)
{
    if (!ctx) return; if (enters) *enters=0; if (exits) *exits=0; if (prefetch_count) *prefetch_count=0;
    const int R = (int)active_radius;
    // Mark hot tiles within radius
    for (int z=-R; z<=R; ++z) for (int y=-R; y<=R; ++y) for (int x=-R; x<=R; ++x){
        long long k = key(cx+x, cy+y, cz+z);
        TileState& s = ctx->states[k];
        s.hot_frames += 1; s.cold_frames = 0;
        if (!s.active && s.hot_frames >= (int)ctx->opts.enter_frames) { s.active=true; if (enters) (*enters)++; ctx->active_count++; }
    }
    // Cool everything else
    std::vector<long long> to_erase;
    for (auto& it : ctx->states){
        // if not touched this frame (not within radius)
        long long k = it.first; TileState& s = it.second;
        // We detect untouched by checking if hot_frames advanced this tick; we can't directly know, so we keep a decay model:
        if (s.hot_frames > 0) s.hot_frames -= 1; else s.cold_frames += 1;
        if (s.active && s.cold_frames >= (int)ctx->opts.exit_frames) { s.active=false; if (exits) (*exits)++; if (ctx->active_count>0) ctx->active_count--; }
        if (!s.active && s.cold_frames > (int)ctx->opts.exit_frames*2) to_erase.push_back(k);
    }
    for (auto k : to_erase) ctx->states.erase(k);

    if (prefetch_count) {
        const int PR = (int)ctx->opts.prefetch_rings; *prefetch_count = (uint32_t)((2*(R+PR)+1)*(2*(R+PR)+1)*(2*(R+PR)+1) - (2*R+1)*(2*R+1)*(2*R+1));
    }
}

ZX_API uint32_t ZX_CALL zx_residency_get_active_count(const zx_residency* ctx){ return ctx?ctx->active_count:0; }

}


