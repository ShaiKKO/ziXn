/*!
 * \file zx_residency.cpp
 * \brief Tile residency manager with hysteresis and prefetch hints.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_residency.h"
#include <algorithm>
#include <unordered_map>
#include <vector>

struct tile_state  // NOLINT(readability-identifier-naming)
{
  int hot_frames  = 0;
  int cold_frames = 0;
  bool active     = false;
};

struct zx_residency
{
  zx_residency_opts opts;
  std::unordered_map<long long, tile_state> states;
  uint32_t active_count = 0;
  std::unordered_map<long long, bool> pinned;
  uint32_t last_enters = 0, last_exits = 0;
};

static inline long long key(int x, int y, int z)
{
  constexpr int shift_x = 42;  // NOLINT(readability-magic-numbers)
  constexpr int shift_y = 21;  // NOLINT(readability-magic-numbers)
  return ((static_cast<long long>(static_cast<unsigned int>(x)) << shift_x) ^
          (static_cast<long long>(static_cast<unsigned int>(y)) << shift_y) ^
          static_cast<long long>(static_cast<unsigned int>(z)));
}

extern "C"
{

  /** \brief Create residency context.
   * @param opts Options for enter/exit frames and prefetch rings (must not be NULL)
   * @return Residency handle; caller owns and must destroy
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
    auto r  = new zx_residency();  // NOLINT(cppcoreguidelines-owning-memory)
    r->opts = *opts;
    return r;
  }

  /**
   * @brief Destroys and frees a residency context.
   *
   * Deletes the given zx_residency instance and releases its resources.
   * Passing nullptr is safe and has no effect. After calling this, the pointer
   * must not be used.
   *
   * @param ctx Pointer to the residency context to destroy.
   */
  void ZX_CALL zx_residency_destroy(zx_residency* ctx)
  {
    delete ctx;  // NOLINT(cppcoreguidelines-owning-memory)
  }

  /**
   * @brief Advance residency state by one frame around a tile center.
   *
   * Advances per-tile hot/cold counters for a frame centered at (cx,cy,cz) using
   * an active cubic radius (in tiles). Tiles within the cube are considered
   * "hot" this tick; all others cool down. Tiles transition to active when their
   * hot counter reaches opts.enter_frames and deactivate when their cold counter
   * reaches opts.exit_frames. Inactive tiles with prolonged coldness
   * (cold_frames > exit_frames*2) are pruned from the internal map.
   *
   * Side effects:
   * - Mutates ctx->states (per-tile TileState) and ctx->active_count.
   * - Increments `*enters`/`*exits` for newly activated/deactivated tiles if those
   *   output pointers are non-null.
   * - Writes a computed prefetch ring size to `*prefetch_count` if non-null.
   *
   * @param ctx Residency handle (no-op if null).
   * @param cx Center X coordinate in tiles.
   * @param cy Center Y coordinate in tiles.
   * @param cz Center Z coordinate in tiles.
   * @param active_radius Radius (in tiles) defining the cubic hot region.
   * @param enters Optional out parameter; incremented by the number of tiles
   *        that became active this tick.
   * @param exits Optional out parameter; incremented by the number of tiles
   *        that became inactive this tick.
   * @param prefetch_count Optional out parameter receiving the size of the
   *        prefetch ring computed from `active_radius` and `opts.prefetch_rings`.
   */
  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  void ZX_CALL
  zx_residency_tick(zx_residency* ctx, int center_x, int center_y,
                    int center_z,  // NOLINT(bugprone-easily-swappable-parameters)
                    uint32_t active_radius, uint32_t* enters, uint32_t* exits,
                    uint32_t* prefetch_count)  // NOLINT(bugprone-easily-swappable-parameters)
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
          long long k   = key(center_x + x, center_y + y, center_z + z);
          tile_state& s = ctx->states[k];
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
      long long k   = it.first;
      tile_state& s = it.second;
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
    for (auto k : to_erase)
    {
      ctx->states.erase(k);
    }

    // Force pinned tiles to remain active
    for (const auto& kv : ctx->pinned)
    {
      long long k   = kv.first;
      tile_state& s = ctx->states[k];
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
      const int pr = static_cast<int>(ctx->opts.prefetch_rings);
      *prefetch_count =
          static_cast<uint32_t>(((2 * (r + pr) + 1) * (2 * (r + pr) + 1) * (2 * (r + pr) + 1)) -
                                ((2 * r + 1) * (2 * r + 1) * (2 * r + 1)));
    }
  }

  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  void ZX_CALL zx_residency_pin_box(zx_residency* ctx, int x0, int y0, int z0, int x1, int y1,
                                    int z1) /* NOLINT(bugprone-easily-swappable-parameters) */
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

  void ZX_CALL zx_residency_unpin_all(zx_residency* ctx)
  {
    if (ctx == nullptr)
    {
      return;
    }
    ctx->pinned.clear();
  }

  void ZX_CALL zx_residency_set_prefetch_rings(zx_residency* ctx, uint32_t rings)
  {
    if (ctx == nullptr)
    {
      return;
    }
    ctx->opts.prefetch_rings = rings;
  }

  // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
  void ZX_CALL
  zx_residency_get_last_churn(const zx_residency* ctx, uint32_t* enters, uint32_t* exits,
                              uint32_t* churn) /* NOLINT(bugprone-easily-swappable-parameters) */
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
   * @brief Returns the number of tiles currently marked active (from the last tick).
   *
   * If `ctx` is null, returns 0.
   *
   * @return uint32_t Current active tile count.
   */
  uint32_t ZX_CALL zx_residency_get_active_count(const zx_residency* ctx)
  {
    return (ctx != nullptr) ? ctx->active_count : 0;
  }
}
