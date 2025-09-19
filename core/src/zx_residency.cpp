/*!
 * \file zx_residency.cpp
 * \brief Tile residency manager with hysteresis and prefetch hints.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_residency.h"
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>
struct TileState
{
  int hot_frames  = 0;
  int cold_frames = 0;
  bool active     = false;
};

struct zx_residency
{
  zx_residency_opts opts{};
  std::unordered_map<long long, TileState> states;
  uint32_t active_count = 0;
  std::unordered_map<long long, bool> pinned;
  uint32_t last_enters = 0, last_exits = 0;
};

static inline long long key(int x, int y, int z)
{
  constexpr int shift_x = 42;
  constexpr int shift_y = 21;
  return ((static_cast<long long>(static_cast<unsigned int>(x)) << shift_x) ^
          (static_cast<long long>(static_cast<unsigned int>(y)) << shift_y) ^
          static_cast<long long>(static_cast<unsigned int>(z)));
}

extern "C"
{

  /**
   * @brief Create residency context.
   * @param opts Options for enter/exit frames and prefetch rings (must not be NULL).
   * @return Residency handle; caller owns and must destroy.
   */
  zx_residency*
      ZX_CALL /**
               * @brief Allocate and initialize a new residency context.
               *
               * Creates a new zx_residency instance on the heap and copies the supplied
               * options into the new context.
               *
               * @param opts Pointer to options used to initialize the context. Must be non-null.
               * @return zx_residency* Pointer to the newly allocated residency context.
               *         The caller takes ownership and must free it with zx_residency_destroy().
               */
      zx_residency_create(const zx_residency_opts* opts)
  {
    if (opts == nullptr)
    {
      return nullptr;
    }
    auto* r = new zx_residency();
    r->opts = *opts;
    return r;
  }

  /**
   * @brief Destroy and free a residency context.
   * @param ctx Residency context to destroy (nullptr is safe).
   */
  void ZX_CALL zx_residency_destroy(zx_residency* ctx)
  {
    delete ctx;
  }

  /**
   * @brief Advance residency state by one frame around a tile center.
   * @param ctx Residency handle (no-op if null).
   * @param center_x Center X coordinate in tiles.
   * @param center_y Center Y coordinate in tiles.
   * @param center_z Center Z coordinate in tiles.
   * @param active_radius Radius (in tiles) defining the cubic hot region.
   * @param enters Optional; incremented by tiles that became active this tick.
   * @param exits Optional; incremented by tiles that became inactive this tick.
   * @param prefetch_count Optional; receives ring size derived from active radius.
   */
  void ZX_CALL zx_residency_tick(zx_residency* ctx, int center_x, int center_y, int center_z,
                                 uint32_t active_radius, uint32_t* enters, uint32_t* exits,
                                 uint32_t* prefetch_count)
  {
    if (ctx == nullptr)
    {
      return;
    }
    ctx->last_enters = 0;
    ctx->last_exits  = 0;
    if (enters != nullptr)
    {
      *enters = 0;
    }
    if (exits != nullptr)
    {
      *exits = 0;
    }
    if (prefetch_count != nullptr)
    {
      *prefetch_count = 0;
    }
    const int r = static_cast<int>(active_radius);
    // Mark hot tiles within radius
    for (int z = -r; z <= r; ++z)
    {
      for (int y = -r; y <= r; ++y)
      {
        for (int x = -r; x <= r; ++x)
        {
          long long k  = key(center_x + x, center_y + y, center_z + z);
          TileState& s = ctx->states[k];
          s.hot_frames += 1;
          s.cold_frames = 0;
          if (!s.active && s.hot_frames >= static_cast<int>(ctx->opts.enter_frames))
          {
            s.active = true;
            ctx->last_enters++;
            if (enters != nullptr)
            {
              (*enters)++;
            }
            ctx->active_count++;
          }
        }
      }
    }
    // Cool everything else
    std::vector<long long> to_erase;
    for (auto& it : ctx->states)
    {
      // if not touched this frame (not within radius)
      long long k  = it.first;
      TileState& s = it.second;
      // We detect untouched by checking if hot_frames advanced this tick; we can't directly know,
      // so we keep a decay model:
      if (s.hot_frames > 0)
      {
        s.hot_frames -= 1;
      }
      else
      {
        s.cold_frames += 1;
      }
      if (s.active && s.cold_frames >= static_cast<int>(ctx->opts.exit_frames))
      {
        s.active = false;
        ctx->last_exits++;
        if (exits != nullptr)
        {
          (*exits)++;
        }
        if (ctx->active_count > 0)
        {
          ctx->active_count--;
        }
      }
      if (!s.active && s.cold_frames > static_cast<int>(ctx->opts.exit_frames) * 2)
      {
        to_erase.push_back(k);
      }
    }
    for (const auto k : to_erase)
    {
      ctx->states.erase(k);
    }

    // Force pinned tiles to remain active
    for (const auto& kv : ctx->pinned)
    {
      long long k  = kv.first;
      TileState& s = ctx->states[k];
      if (!s.active)
      {
        s.active = true;
        ctx->active_count++;
      }
      s.hot_frames  = std::max(s.hot_frames, static_cast<int>(ctx->opts.enter_frames));
      s.cold_frames = 0;
    }
    if (prefetch_count != nullptr)
    {
      const int pr             = static_cast<int>(ctx->opts.prefetch_rings);
      const int64_t a          = (2LL * (r + pr)) + 1;
      const int64_t b          = (2LL * r) + 1;
      const uint64_t vol_shell = static_cast<uint64_t>(a) * a * a;
      const uint64_t vol_core  = static_cast<uint64_t>(b) * b * b;
      const uint64_t diff      = (vol_shell > vol_core) ? (vol_shell - vol_core) : 0ULL;
      *prefetch_count          = (diff > std::numeric_limits<uint32_t>::max())
                                     ? std::numeric_limits<uint32_t>::max()
                                     : static_cast<uint32_t>(diff);
    }
  }

  /**
   * @brief Pin all tiles in an axis-aligned box so they remain active.
   * @param ctx Residency context (no-op if null).
   * @param x0 Min x, @param y0 Min y, @param z0 Min z
   * @param x1 Max x, @param y1 Max y, @param z1 Max z
   */
  void ZX_CALL zx_residency_pin_box(zx_residency* ctx, int x0, int y0, int z0, int x1, int y1,
                                    int z1)
  {
    if (ctx == nullptr)
    {
      return;
    }
    if (x0 > x1)
    {
      std::swap(x0, x1);
    }
    if (y0 > y1)
    {
      std::swap(y0, y1);
    }
    if (z0 > z1)
    {
      std::swap(z0, z1);
    }
    for (int z = z0; z <= z1; ++z)
    {
      for (int y = y0; y <= y1; ++y)
      {
        for (int x = x0; x <= x1; ++x)
        {
          ctx->pinned[key(x, y, z)] = true;
        }
      }
    }
  }

  /**
   * @brief Clear all pinned tiles.
   * @param ctx Residency context (no-op if null).
   */
  void ZX_CALL zx_residency_unpin_all(zx_residency* ctx)
  {
    if (ctx == nullptr)
    {
      return;
    }
    ctx->pinned.clear();
  }

  /**
   * @brief Set number of prefetch rings used to compute prefetch surface.
   * @param ctx Residency context (no-op if null).
   * @param rings Prefetch rings (non-negative integer)
   */
  void ZX_CALL zx_residency_set_prefetch_rings(zx_residency* ctx, uint32_t rings)
  {
    if (ctx == nullptr)
    {
      return;
    }
    ctx->opts.prefetch_rings = rings;
  }

  /**
   * @brief Get last tick enter/exit counts and their sum (churn).
   * @param ctx Residency context (may be null; outputs become zero).
   * @param enters Optional; receives last enters.
   * @param exits Optional; receives last exits.
   * @param churn Optional; receives enters+exits.
   */
  void ZX_CALL zx_residency_get_last_churn(const zx_residency* ctx, uint32_t* enters,
                                           uint32_t* exits, uint32_t* churn)
  {
    if (ctx == nullptr)
    {
      if (enters != nullptr)
      {
        *enters = 0;
      }
      if (exits != nullptr)
      {
        *exits = 0;
      }
      if (churn != nullptr)
      {
        *churn = 0;
      }
      return;
    }
    if (enters != nullptr)
    {
      *enters = ctx->last_enters;
    }
    if (exits != nullptr)
    {
      *exits = ctx->last_exits;
    }
    if (churn != nullptr)
    {
      *churn = ctx->last_enters + ctx->last_exits;
    }
  }

  /**
   * @brief Current active tile count (from the last tick).
   * @return Active tile count (0 if ctx is null).
   */
  uint32_t ZX_CALL zx_residency_get_active_count(const zx_residency* ctx)
  {
    return (ctx != nullptr) ? ctx->active_count : 0;
  }
#ifdef __cplusplus
}  // extern "C"
#endif
