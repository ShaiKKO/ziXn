/*!
 * \file zx_telemetry.cpp
 * \brief Thread-safe telemetry buffer and CSV/JSON exporters (minimal core).
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_telemetry.h"
#include <vector>
#include <string>
#include <mutex>
#include <cstdio>

struct zx_sample {
    std::string scene;
    uint32_t step;
    std::vector<std::pair<std::string, float>> counters;
};

struct zx_telemetry {
    std::mutex mtx;
    std::vector<zx_sample> samples;
    uint32_t capacity;
    zx_sample current;
    bool in_step;
    struct err { std::string scene; uint32_t step; std::string code; std::string msg; };
    std::vector<err> errors;
};

extern "C" {

/** \brief Create a telemetry buffer; capacity==0 selects a reasonable default.
 * @param capacity Ring capacity (samples); 0 for default
 * @return New telemetry handle
 */
ZX_API zx_telemetry* ZX_CALL zx_telemetry_create(uint32_t capacity){
    zx_telemetry* t = new zx_telemetry();
    t->capacity = capacity ? capacity : 1024;
    t->in_step = false;
    return t;
}

/** \brief Destroy telemetry and release memory. */
ZX_API void ZX_CALL zx_telemetry_destroy(zx_telemetry* ctx){ delete ctx; }

/** \brief Begin a step (no-op if ctx==NULL). scene may be NULL.
 * @param ctx Telemetry handle
 * @param scene Scene identifier (may be NULL)
 * @param step_index Step index
 */
ZX_API void ZX_CALL zx_telemetry_begin_step(zx_telemetry* ctx, const char* scene, uint32_t step_index){
    if (!ctx) return;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    ctx->current = zx_sample{}; ctx->current.scene = scene ? scene : ""; ctx->current.step = step_index; ctx->in_step = true;
}

/** \brief Set a counter during an active step (no-op if ctx==NULL/not in step).
 * @param ctx Telemetry handle
 * @param name Counter name (must not be NULL)
 * @param value Counter value
 */
ZX_API void ZX_CALL zx_telemetry_set_counter(zx_telemetry* ctx, const char* name, float value){
    if (!ctx || !name) return;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    if (!ctx->in_step) return;
    ctx->current.counters.emplace_back(std::string(name), value);
}

/** \brief End current step (no-op if ctx==NULL or not in a step). */
ZX_API void ZX_CALL zx_telemetry_end_step(zx_telemetry* ctx){
    if (!ctx) return;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    if (!ctx->in_step) return;
    if (ctx->samples.size() >= ctx->capacity) ctx->samples.erase(ctx->samples.begin());
    ctx->samples.push_back(std::move(ctx->current));
    ctx->in_step = false;
}

/** \brief Add an error entry with scene/step context. */
ZX_API void ZX_CALL zx_telemetry_add_error(zx_telemetry* ctx, const char* scene, uint32_t step_index, const char* code, const char* message){
    if (!ctx || !code || !message) return;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    zx_telemetry::err e; e.scene = scene ? scene : ""; e.step = step_index; e.code = code; e.msg = message;
    ctx->errors.push_back(std::move(e));
}

static void export_header(FILE* f, const zx_sample& s){
    std::fprintf(f, "scene,step");
    for (const auto& kv : s.counters) std::fprintf(f, ",%s", kv.first.c_str());
    std::fprintf(f, "\n");
}

/** \brief Export recent samples to CSV.
 * @param ctx Telemetry handle
 * @param path Output file path
 * @return 0 on success, negative on error
 */
ZX_API int ZX_CALL zx_telemetry_export_csv(zx_telemetry* ctx, const char* path){
    if (!ctx || !path) return -1;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    FILE* f = std::fopen(path, "wb"); if (!f) return -2;
    if (ctx->samples.empty()) { std::fclose(f); return 0; }
    export_header(f, ctx->samples.front());
    for (const auto& s : ctx->samples){
        std::fprintf(f, "%s,%u", s.scene.c_str(), s.step);
        for (const auto& kv : s.counters) std::fprintf(f, ",%.9g", kv.second);
        std::fprintf(f, "\n");
    }
    std::fclose(f); return 0;
}

/** \brief Export recent samples to JSON.
 * @param ctx Telemetry handle
 * @param path Output file path
 * @return 0 on success, negative on error
 */
ZX_API int ZX_CALL zx_telemetry_export_json(zx_telemetry* ctx, const char* path){
    if (!ctx || !path) return -1;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    FILE* f = std::fopen(path, "wb"); if (!f) return -2;
    std::fprintf(f, "{\n  \"samples\": [\n");
    for (size_t i=0;i<ctx->samples.size();++i){
        const auto& s = ctx->samples[i];
        std::fprintf(f, "    { \"scene\": \"%s\", \"step\": %u, \"counters\": {", s.scene.c_str(), s.step);
        for (size_t k=0;k<s.counters.size();++k){
            const auto& kv = s.counters[k];
            std::fprintf(f, "\"%s\": %.9g%s", kv.first.c_str(), kv.second, (k+1<s.counters.size()?", ":""));
        }
        std::fprintf(f, "} }%s\n", (i+1<ctx->samples.size()?",":""));
    }
    std::fprintf(f, "  ],\n  \"errors\": [\n");
    for (size_t i=0;i<ctx->errors.size(); ++i){
        const auto& e = ctx->errors[i];
        std::fprintf(f, "    { \"scene\": \"%s\", \"step\": %u, \"code\": \"%s\", \"message\": \"%s\" }%s\n",
                    e.scene.c_str(), e.step, e.code.c_str(), e.msg.c_str(), (i+1<ctx->errors.size()?",":""));
    }
    std::fprintf(f, "  ]\n}\n");
    std::fclose(f); return 0;
}

} // extern "C"


