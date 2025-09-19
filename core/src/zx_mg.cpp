/*!
 * \file zx_mg.cpp
 * \brief Geometric Multigrid preconditioner (1D Poisson V-cycle) implementation.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_mg.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace
{
  constexpr float k_two     = 2.0F;
  constexpr float k_half    = 0.5F;
  constexpr float k_quarter = 0.25F;
  constexpr float k_four    = 4.0F;
}  // namespace

struct ZxMgLevel
{
  size_t n    = 0;
  float h2inv = 0.0F;
  std::vector<float> x, r, z, tmp;
};

struct zx_mg_context
{
  std::vector<ZxMgLevel> levels;
  zx_mg_opts opts{};
};

static inline void apply_a_1d(const float* x, float* y, size_t n, float h2inv)
{
  y[0] = h2inv * (k_two * x[0] - x[1]);
  for (size_t i = 1; i < n - 1; ++i)
  {
    y[i] = h2inv * (-x[i - 1] + k_two * x[i] - x[i + 1]);
  }
  y[n - 1] = h2inv * (k_two * x[n - 1] - x[n - 2]);
}

static void jacobi_smooth(size_t iters, float omega, ZxMgLevel& level)
{
  const size_t n = level.n;
  auto& ax       = level.tmp;
  ax.resize(n);
  for (size_t k = 0; k < iters; ++k)
  {
    apply_a_1d(level.x.data(), ax.data(), n, level.h2inv);
    for (size_t i = 0; i < n; ++i)
    {
      float diag = k_two * level.h2inv;
      float res  = level.r[i] - ax[i];
      level.x[i] += omega * res / diag;
    }
  }
}

static void restrict_full_weighting(const std::vector<float>& fine, std::vector<float>& coarse)
{
  const size_t nf = fine.size();
  const size_t nc = coarse.size();
  coarse[0]       = k_half * (fine[0] + fine[1]);
  for (size_t i = 1; i < nc - 1; ++i)
  {
    size_t j  = 2 * i;
    coarse[i] = k_quarter * (fine[j - 1] + k_two * fine[j] + fine[j + 1]);
  }
  coarse[nc - 1] = k_half * (fine[nf - 2] + fine[nf - 1]);
}

/**
 * @brief Prolongates a coarse-grid vector to a fine grid using linear interpolation.
 *
 * Copies coarse values into the even indices of fine and fills odd indices with the
 * average of neighboring coarse values: fine[2*i] = coarse[i], fine[2*i+1] = 0.5*(coarse[i] +
 * coarse[i+1]).
 *
 * @param coarse Source vector on the coarse grid.
 * @param fine Destination vector on the fine grid; must be preallocated with size = 2*coarse.size()
 * - 1.
 */
static void prolong_linear(const std::vector<float>& coarse, std::vector<float>& fine)
{
  const size_t nc = coarse.size();
  // nf = 2*nc - 1
  for (size_t i = 0; i < nc; ++i)
  {
    fine[2 * i] = coarse[i];
  }
  for (size_t i = 0; i < nc - 1; ++i)
  {
    fine[(2 * i) + 1] = k_half * (coarse[i] + coarse[i + 1]);
  }
}

static void vcycle(zx_mg_context* ctx, size_t level_idx)
{
  ZxMgLevel& level    = ctx->levels[level_idx];
  const zx_mg_opts& o = ctx->opts;
  jacobi_smooth(o.pre_smooth, o.omega, level);

  if (level_idx + 1 < ctx->levels.size())
  {
    // compute residual r - A x
    auto& ax = level.tmp;
    ax.resize(level.n);
    apply_a_1d(level.x.data(), ax.data(), level.n, level.h2inv);
    for (size_t i = 0; i < level.n; ++i)
    {
      level.tmp[i] = level.r[i] - ax[i];
    }
    // restrict to coarse
    ZxMgLevel& coarse = ctx->levels[level_idx + 1];
    restrict_full_weighting(level.tmp, coarse.r);
    std::fill(coarse.x.begin(), coarse.x.end(), 0.0F);
    vcycle(ctx, level_idx + 1);
    // prolong and correct
    prolong_linear(coarse.x, level.tmp);
    for (size_t i = 0; i < level.n; ++i)
    {
      level.x[i] += level.tmp[i];
    }
    jacobi_smooth(o.post_smooth, o.omega, level);
  }
  else
  {
    jacobi_smooth(o.coarse_iters, o.omega, level);
  }
}

extern "C"
{

  /** \brief Build a 1D Poisson hierarchy for size n (n>=3).
   * @param n Problem size (>=3)
   * @param opts Options (must not be NULL)
   * @return Context pointer or NULL on failure
   */
  /**
   * @brief Create a 1D Poisson V-cycle hierarchy for PCG preconditioning.
   * @param n Finest grid size (>= 3).
   * @param opts Options (must not be NULL).
   * @return New multigrid context or NULL on failure.
   */
  zx_mg_context* ZX_CALL zx_mg_create_poisson1d(size_t n, const zx_mg_opts* opts)
  {
    if ((n < 3) || (opts == nullptr))
    {
      return nullptr;
    }
    // C-ABI factory with RAII during construction
    auto ctx    = std::make_unique<zx_mg_context>();
    ctx->opts   = *opts;
    size_t cur  = n;
    float h2inv = 1.0F;
    for (uint32_t lvl = 0; lvl < opts->max_levels && cur >= 3; ++lvl)
    {
      ZxMgLevel level;
      level.n     = cur;
      level.h2inv = h2inv;
      level.x.assign(cur, 0.0F);
      level.r.assign(cur, 0.0F);
      level.z.assign(cur, 0.0F);
      level.tmp.assign(cur, 0.0F);
      ctx->levels.push_back(std::move(level));
      if (cur == 3)
      {
        break;
      }
      size_t next = (cur + 1) / 2;  // ensure nf ≈ 2*nc - 1
      if (next < 3)
      {
        break;
      }
      cur = next;
      h2inv *= k_four;  // grid spacing doubles each level -> 1/h^2 scales by 1/4; here we keep h^2
                        // consistent per stencil scaling usage
    }
    if (ctx->levels.empty())
    {
      return nullptr;
    }
    return ctx.release();
  }

  /** \brief Destroy multigrid context. */
  /** @brief Destroy multigrid context. */
  void ZX_CALL zx_mg_destroy(zx_mg_context* ctx)
  {
    delete ctx;
  }

  /**
   * @brief Apply the multigrid preconditioner: compute z = M^{-1} r for a 1D Poisson hierarchy.
   *
   * Loads the input residual into the finest-level right-hand side, runs a V-cycle, and writes
   * the computed finest-level solution into z.
   *
   * The caller must provide arrays r and z of length equal to the finest level size configured
   * in the multigrid context. If any of r, z, or the context pointer is null, the function
   * returns immediately without modifying z.
   *
   * @param r Input residual vector on the finest grid (length = finest level n).
   * @param z Output vector that receives the preconditioned result (length = finest level n).
   */
  /**
   * @brief Apply multigrid preconditioner: compute z = M^{-1} r.
   * @param r Input residual (length finest n).
   * @param z Output vector (length finest n).
   * @param user_ctx Opaque zx_mg_context* passed from PCG caller.
   */
  void ZX_CALL zx_mg_prec_apply(const float* r, float* z, void* user_ctx)
  {
    auto* ctx = static_cast<zx_mg_context*>(user_ctx);
    if ((ctx == nullptr) || (r == nullptr) || (z == nullptr))
    {
      return;
    }
    // load RHS into finest level
    ZxMgLevel& f = ctx->levels.front();
    std::copy_n(r, f.n, f.r.begin());
    std::fill(f.x.begin(), f.x.end(), 0.0F);
    vcycle(ctx, 0);
    std::copy(f.x.begin(), f.x.end(), z);
  }

}  // extern "C"
