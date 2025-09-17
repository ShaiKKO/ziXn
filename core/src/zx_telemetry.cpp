/*!
 * \file zx_telemetry.cpp
 * \brief Thread-safe telemetry buffer and CSV/JSON exporters (minimal core).
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_telemetry.h"
#include <cstdio>
#include <mutex>
#include <string>
#include <vector>

struct zx_sample
{
  std::string scene{};
  uint32_t step{0};
  std::vector<std::pair<std::string, float>> counters{};
};

struct zx_telemetry
{
  std::mutex mtx{};
  std::vector<zx_sample> samples{};
  uint32_t capacity{0};
  zx_sample current{};
  bool in_step{false};
  struct err
  {
    std::string scene{};
    uint32_t step{0};
    std::string code{};
    std::string msg{};
  };
  std::vector<err> errors{};
};

extern "C"
{

  /** \brief Create a telemetry buffer; capacity==0 selects a reasonable default.
   * @param capacity Ring capacity (samples); 0 for default
   * @return New telemetry handle
   */
  /**
   * @brief Allocate and initialize a telemetry context.
   *
   * Creates a new zx_telemetry instance with the given sample capacity. If
   * `capacity` is zero, a default capacity of 1024 is used. The returned object
   * is heap-allocated and the caller is responsible for destroying it.
   *
   * @param capacity Maximum number of samples to retain; if zero, defaults to 1024.
   * @return zx_telemetry* Pointer to the newly allocated telemetry object (owner must
   * call zx_telemetry_destroy).
   */
  zx_telemetry* ZX_CALL zx_telemetry_create(uint32_t capacity)
  {
    auto* t                                    = new zx_telemetry();
    static constexpr uint32_t kDefaultCapacity = 1024u;
    t->capacity                                = (capacity != 0u) ? capacity : kDefaultCapacity;
    t->in_step                                 = false;
    return t;
  }

  /**
   * @brief Destroy a telemetry context and free its resources.
   *
   * Deletes the provided zx_telemetry object, releasing memory and any owned resources.
   * It is safe to pass nullptr (no action is performed).
   *
   * @param ctx Pointer to the telemetry context to destroy.
   */
  void ZX_CALL zx_telemetry_destroy(zx_telemetry* ctx)
  {
    delete ctx;
  }

  /**
   * @brief Start a new telemetry step and begin collecting counters.
   *
   * Begins a step by resetting the active sample, setting its scene (empty if NULL)
   * and step index, and marking the telemetry context as being "in step".
   * This function is a no-op if the provided telemetry handle is NULL.
   *
   * @param scene Optional scene identifier; if NULL the scene is set to an empty string.
   * @param step_index Index of the step to begin.
   */
  void ZX_CALL zx_telemetry_begin_step(zx_telemetry* ctx, const char* scene, uint32_t step_index)
  {
    if (!ctx)
      return;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    ctx->current       = zx_sample{};
    ctx->current.scene = (scene != nullptr) ? scene : "";
    ctx->current.step  = step_index;
    ctx->in_step       = true;
  }

  /**
   * @brief Record a named counter for the currently active telemetry step.
   *
   * If a step is active, appends the (name, value) pair to the in-progress sample in a
   * thread-safe manner. This function is a no-op if there is no active step or if the
   * telemetry context or `name` is null.
   *
   * @param name Counter name (must not be null).
   * @param value Counter value to record for the active step.
   */
  void ZX_CALL zx_telemetry_set_counter(zx_telemetry* ctx, const char* name, float value)
  {
    if ((ctx == nullptr) || (name == nullptr))
      return;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    if (!ctx->in_step)
      return;
    ctx->current.counters.emplace_back(std::string(name), value);
  }

  /**
   * @brief Finalize and commit the active telemetry step.
   *
   * If a step is active, moves the accumulated current sample into the internal
   * samples buffer and marks the telemetry as not in a step. If the buffer has
   * reached its configured capacity, the oldest sample is removed to make room.
   *
   * This function is a no-op when called with a null context or when no step is
   * currently active. The operation is performed under the telemetry mutex.
   */
  void ZX_CALL zx_telemetry_end_step(zx_telemetry* ctx)
  {
    if (ctx == nullptr)
      return;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    if (!ctx->in_step)
      return;
    if (ctx->samples.size() >= ctx->capacity)
      ctx->samples.erase(ctx->samples.begin());
    ctx->samples.push_back(std::move(ctx->current));
    ctx->in_step = false;
  }

  /**
   * @brief Record an error entry associated with a scene and step.
   *
   * Appends an error record into the telemetry's internal error log in a thread-safe manner.
   * If `scene` is null it is stored as an empty string. If `ctx`, `code`, or `message` is null
   * the function is a no-op.
   *
   * @param scene Optional scene context for the error (nullable).
   * @param step_index Step index associated with the error.
   * @param code Short error code string (must not be null).
   * @param message Human-readable error message (must not be null).
   */
  void ZX_CALL zx_telemetry_add_error(zx_telemetry* ctx, const char* scene, uint32_t step_index,
                                      const char* code, const char* message)
  {
    if ((ctx == nullptr) || (code == nullptr) || (message == nullptr))
      return;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    zx_telemetry::err e;
    e.scene = (scene != nullptr) ? scene : "";
    e.step  = step_index;
    e.code  = code;
    e.msg   = message;
    ctx->errors.push_back(std::move(e));
  }

  /**
   * @brief Writes a CSV header row for a sample to the given file.
   *
   * Writes "scene,step" followed by each counter name from the sample's
   * counters (in order) as CSV columns, then a newline.
   *
   * @param f Destination file stream (must be non-null and open for writing).
   * @param s Sample whose counter names are used to build the header row.
   */
  static void export_header(FILE* f, const zx_sample& s)
  {
    std::fprintf(f, "scene,step");
    for (const auto& kv : s.counters)
      std::fprintf(f, ",%s", kv.first.c_str());
    std::fprintf(f, "\n");
  }

  /**
   * @brief Export buffered telemetry samples to a CSV file.
   *
   * Writes the stored samples as CSV rows with a header derived from the first sample.
   * Each row contains: scene, step, then counter values in the same order as the header.
   * The function is thread-safe (acquires the telemetry mutex) and will create the output
   * file even when there are no samples (resulting file will be empty aside from being created).
   *
   * @param ctx Telemetry handle.
   * @param path Filesystem path to write the CSV output.
   * @return int 0 on success;
   *             -1 if either `ctx` or `path` is null;
   *             -2 if the output file cannot be opened for writing.
   */
  int ZX_CALL zx_telemetry_export_csv(zx_telemetry* ctx, const char* path)
  {
    if (!ctx || !path)
      return -1;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    FILE* f = std::fopen(path, "wb");
    if (!f)
      return -2;
    if (ctx->samples.empty())
    {
      std::fclose(f);
      return 0;
    }
    export_header(f, ctx->samples.front());
    for (const auto& s : ctx->samples)
    {
      std::fprintf(f, "%s,%u", s.scene.c_str(), s.step);
      for (const auto& kv : s.counters)
        std::fprintf(f, ",%.9g", kv.second);
      std::fprintf(f, "\n");
    }
    std::fclose(f);
    return 0;
  }

  /**
   * @brief Export the in-memory telemetry buffer to a JSON file.
   *
   * Writes a JSON object with two top-level arrays: "samples" and "errors".
   * - "samples" contains objects of the form:
   *   { "scene": "<string>", "step": <uint>, "counters": { "<name>": <value>, ... } }
   *   Numeric counter values are written with up to 9 significant digits.
   * - "errors" contains objects of the form:
   *   { "scene": "<string>", "step": <uint>, "code": "<string>", "message": "<string>" }
   *
   * @param ctx Telemetry handle whose buffered samples and errors will be exported.
   * @param path Filesystem path to write the JSON output to.
   * @return int 0 on success; -1 if an argument is null; -2 if the output file could not be opened.
   */
  int ZX_CALL zx_telemetry_export_json(zx_telemetry* ctx, const char* path)
  {
    if (!ctx || !path)
      return -1;
    std::lock_guard<std::mutex> lock(ctx->mtx);
    FILE* f = std::fopen(path, "wb");
    if (!f)
      return -2;
    std::fprintf(f, "{\n  \"samples\": [\n");
    for (size_t i = 0; i < ctx->samples.size(); ++i)
    {
      const auto& s = ctx->samples[i];
      std::fprintf(f, "    { \"scene\": \"%s\", \"step\": %u, \"counters\": {", s.scene.c_str(),
                   s.step);
      for (size_t k = 0; k < s.counters.size(); ++k)
      {
        const auto& kv = s.counters[k];
        std::fprintf(f, "\"%s\": %.9g%s", kv.first.c_str(), kv.second,
                     (k + 1 < s.counters.size() ? ", " : ""));
      }
      std::fprintf(f, "} }%s\n", (i + 1 < ctx->samples.size() ? "," : ""));
    }
    std::fprintf(f, "  ],\n  \"errors\": [\n");
    for (size_t i = 0; i < ctx->errors.size(); ++i)
    {
      const auto& e = ctx->errors[i];
      std::fprintf(
          f, "    { \"scene\": \"%s\", \"step\": %u, \"code\": \"%s\", \"message\": \"%s\" }%s\n",
          e.scene.c_str(), e.step, e.code.c_str(), e.msg.c_str(),
          (i + 1 < ctx->errors.size() ? "," : ""));
    }
    std::fprintf(f, "  ]\n}\n");
    std::fclose(f);
    return 0;
  }

}  // extern "C"
