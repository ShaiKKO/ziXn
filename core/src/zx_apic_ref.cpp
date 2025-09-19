/*!
 * \file zx_apic_ref.cpp
 * \brief CPU reference APIC/MLS transfers used for validation.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_apic_ref.h"
#include "zx/zx_abi.h"
#include "zx/zx_determinism.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <vector>

// Minimal array view to centralize pointer indexing while preserving the C ABI
// in public headers.
template <typename T> class ArrayView
{
public:
  explicit ArrayView(T* ptr, size_t count) : ptr_(ptr), count_(count) {}

  T& operator[](size_t idx) const
  {
    // Bounds are validated by callers where applicable; this centralizes pointer
    // arithmetic to a single location.
    return ptr_[idx];
  }

  [[nodiscard]] T* data() const
  {
    return ptr_;
  }
  [[nodiscard]] size_t size() const
  {
    return count_;
  }

private:
  T* ptr_;
  size_t count_;
};

// File-scope constants and indices used across helpers
constexpr size_t k_vec3 = 3U;
constexpr size_t k_cdim = 9U;
constexpr float k_half  = 0.5F;
constexpr float k_four  = 4.0F;

// Indices into 3x3 affine C matrix (row-major)
constexpr size_t k_cxx = 0U;
constexpr size_t k_cxy = 1U;
constexpr size_t k_cxz = 2U;
constexpr size_t k_cyx = 3U;
constexpr size_t k_cyy = 4U;
constexpr size_t k_cyz = 5U;
constexpr size_t k_czx = 6U;
constexpr size_t k_czy = 7U;
constexpr size_t k_czz = 8U;

struct Contrib
{
  int gi;
  int pid;
  float m;
  float mx;
  float my;
  float mz;
};

// Deterministic finalize: sort and accumulate contributions to grid
static inline void finalize_p2g_contribs(std::vector<Contrib>& contribs,
                                         ArrayView<float>& m_grid_view,
                                         ArrayView<float>& p_grid_view)
{
  std::stable_sort(contribs.begin(), contribs.end(),
                   [](const Contrib& a, const Contrib& b)
                   {
                     if (a.gi != b.gi)
                     {
                       return a.gi < b.gi;
                     }
                     return a.pid < b.pid;
                   });
  int cur     = -1;
  float msum  = 0.0F;
  float mxsum = 0.0F;
  float mysum = 0.0F;
  float mzsum = 0.0F;
  auto flush  = [&]()
  {
    if (cur >= 0)
    {
      m_grid_view[cur] += msum;
      const size_t base3 = static_cast<size_t>(cur) * k_vec3;
      p_grid_view[base3 + 0U] += mxsum;
      p_grid_view[base3 + 1U] += mysum;
      p_grid_view[base3 + 2U] += mzsum;
    }
  };
  for (const auto& c : contribs)
  {
    if (c.gi != cur)
    {
      flush();
      cur   = c.gi;
      msum  = 0.0F;
      mxsum = 0.0F;
      mysum = 0.0F;
      mzsum = 0.0F;
    }
    msum += c.m;
    mxsum += c.mx;
    mysum += c.my;
    mzsum += c.mz;
  }
  flush();
}

static inline int idx3(int x, int y, int z, int nx, int ny)
{
  // (((z * ny) + y) * nx) + x
  return (((z * ny) + y) * nx) + x;
}

/**
 * @brief Quadratic B-spline weights (1D).
 * @param x Position relative to cell center in [0,1].
 * @param w Output 3 weights array (non-null).
 */
void zx_bspline_w3(float x, float w[3])
{
  // Public ABI requires pointer-based array. Localize the pointer indexing.
  constexpr float k_half_f           = 0.5F;
  constexpr float k_three_quarters_f = 0.75F;
  const float x0                     = k_half_f - x;
  const float x1                     = k_half_f + x;
  w[0]                               = k_half_f * x0 * x0;
  w[1]                               = k_three_quarters_f - x * x;
  w[2]                               = k_half_f * x1 * x1;
}

/**
 * @brief Derivative of quadratic B-spline (1D).
 * @param x Position relative to cell center in [0,1].
 * @param g Output 3 derivative weights array (non-null).
 */
void zx_bspline_dw3(float x, float g[3])
{
  // Public ABI requires pointer-based array. Localize the pointer indexing.
  constexpr float k_half_f = 0.5F;
  constexpr float k_two_f  = 2.0F;
  const float x0           = k_half_f - x;
  const float x1           = k_half_f + x;
  g[0]                     = -x0;
  g[1]                     = -k_two_f * x;
  g[2]                     = x1;
}

// Helper to compute grid base coordinates and B-spline weights for a particle
static inline void compute_support_and_weights(float gx, float gy, float gz, int& base_x,
                                               int& base_y, int& base_z, std::array<float, 3>& wx,
                                               std::array<float, 3>& wy, std::array<float, 3>& wz)
{
  base_x = static_cast<int>(std::floor(gx - k_half));
  base_y = static_cast<int>(std::floor(gy - k_half));
  base_z = static_cast<int>(std::floor(gz - k_half));
  zx_bspline_w3(gx - (static_cast<float>(base_x) + 1.0F), wx.data());
  zx_bspline_w3(gy - (static_cast<float>(base_y) + 1.0F), wy.data());
  zx_bspline_w3(gz - (static_cast<float>(base_z) + 1.0F), wz.data());
}

// Helper to accumulate a single particle's P2G contributions (deterministic or not)
static inline void accumulate_p2g_particle(size_t i, const ArrayView<const float>& pos_view,
                                           const ArrayView<const float>& vel_view, const float* C,
                                           size_t N, float mp,
                                           const std::array<float, 3>& origin_arr, float h, int nx,
                                           int ny, int nz, ArrayView<float>& m_grid_view,
                                           ArrayView<float>& p_grid_view,
                                           std::vector<Contrib>* contribs)
{
  const size_t pbase = static_cast<size_t>(k_vec3) * i;
  const float px     = pos_view[pbase + 0U];
  const float py     = pos_view[pbase + 1U];
  const float pz     = pos_view[pbase + 2U];
  const float vx0    = vel_view[pbase + 0U];
  const float vy0    = vel_view[pbase + 1U];
  const float vz0    = vel_view[pbase + 2U];

  const float invh = 1.0F / h;
  const float gx   = (px - origin_arr[0]) * invh;
  const float gy   = (py - origin_arr[1]) * invh;
  const float gz   = (pz - origin_arr[2]) * invh;

  int base_x = 0;
  int base_y = 0;
  int base_z = 0;
  std::array<float, 3> wx{};
  std::array<float, 3> wy{};
  std::array<float, 3> wz{};
  compute_support_and_weights(gx, gy, gz, base_x, base_y, base_z, wx, wy, wz);

  std::array<float, k_cdim> c_local{};
  bool has_c = false;
  if (C != nullptr)
  {
    ArrayView<const float> c_view(C, static_cast<size_t>(k_cdim) * N);
    const size_t cbase = static_cast<size_t>(k_cdim) * i;
    std::memcpy(c_local.data(), &c_view[cbase], sizeof(float) * k_cdim);
    has_c = true;
  }

  constexpr int k_support = 3;
  for (int dz = 0; dz < k_support; ++dz)
  {
    for (int dy = 0; dy < k_support; ++dy)
    {
      for (int dx = 0; dx < k_support; ++dx)
      {
        const int ix = base_x + dx;
        const int iy = base_y + dy;
        const int iz = base_z + dz;
        if ((static_cast<unsigned>(ix) >= static_cast<unsigned>(nx)) ||
            (static_cast<unsigned>(iy) >= static_cast<unsigned>(ny)) ||
            (static_cast<unsigned>(iz) >= static_cast<unsigned>(nz)))
        {
          continue;
        }
        const float w = wx.at(static_cast<size_t>(dx)) * wy.at(static_cast<size_t>(dy)) *
                        wz.at(static_cast<size_t>(dz));
        const int gi   = idx3(ix, iy, iz, nx, ny);
        const float mm = w * mp;
        if (contribs == nullptr)
        {
          m_grid_view[gi] += mm;
        }

        float vx_aff = vx0;
        float vy_aff = vy0;
        float vz_aff = vz0;
        if (has_c)
        {
          const float gx_rel = (static_cast<float>(ix) - gx);
          const float gy_rel = (static_cast<float>(iy) - gy);
          const float gz_rel = (static_cast<float>(iz) - gz);
          // APIC affine: v + C*(x_i - x_p) scaled by h
          vx_aff +=
              (c_local[k_cxx] * gx_rel + c_local[k_cxy] * gy_rel + c_local[k_cxz] * gz_rel) * h;
          vy_aff +=
              (c_local[k_cyx] * gx_rel + c_local[k_cyy] * gy_rel + c_local[k_cyz] * gz_rel) * h;
          vz_aff +=
              (c_local[k_czx] * gx_rel + c_local[k_czy] * gy_rel + c_local[k_czz] * gz_rel) * h;
        }

        if (contribs == nullptr)
        {
          const size_t base3 = static_cast<size_t>(gi) * k_vec3;
          p_grid_view[base3 + 0U] += mm * vx_aff;
          p_grid_view[base3 + 1U] += mm * vy_aff;
          p_grid_view[base3 + 2U] += mm * vz_aff;
        }
        else
        {
          contribs->push_back(
              Contrib{gi, static_cast<int>(i), mm, mm * vx_aff, mm * vy_aff, mm * vz_aff});
        }
      }
    }
  }
}

// Helper to accumulate a single particle's G2P gather
static inline void accumulate_g2p_particle(size_t i, const ArrayView<const float>& pos_view,
                                           ArrayView<float>& out_vel_view, ArrayView<float>* out_c,
                                           const std::array<float, 3>& origin_arr, float h, int nx,
                                           int ny, int nz,
                                           const ArrayView<const float>& m_grid_view,
                                           const ArrayView<const float>& v_grid_view)
{
  const size_t pbase = i * static_cast<size_t>(k_vec3);
  const float px     = pos_view[pbase + 0U];
  const float py     = pos_view[pbase + 1U];
  const float pz     = pos_view[pbase + 2U];

  const float invh = 1.0F / h;
  const float gx   = (px - origin_arr[0]) * invh;
  const float gy   = (py - origin_arr[1]) * invh;
  const float gz   = (pz - origin_arr[2]) * invh;

  int base_x = 0;
  int base_y = 0;
  int base_z = 0;
  std::array<float, 3> wx{};
  std::array<float, 3> wy{};
  std::array<float, 3> wz{};
  compute_support_and_weights(gx, gy, gz, base_x, base_y, base_z, wx, wy, wz);

  float vx = 0.0F;
  float vy = 0.0F;
  float vz = 0.0F;
  std::array<float, 3> cx{};
  std::array<float, 3> cy{};
  std::array<float, 3> cz{};

  constexpr int k_support = 3;
  for (int dz = 0; dz < k_support; ++dz)
  {
    for (int dy = 0; dy < k_support; ++dy)
    {
      for (int dx = 0; dx < k_support; ++dx)
      {
        const int ix = base_x + dx;
        const int iy = base_y + dy;
        const int iz = base_z + dz;
        if ((static_cast<unsigned>(ix) >= static_cast<unsigned>(nx)) ||
            (static_cast<unsigned>(iy) >= static_cast<unsigned>(ny)) ||
            (static_cast<unsigned>(iz) >= static_cast<unsigned>(nz)))
        {
          continue;
        }
        const int gi  = idx3(ix, iy, iz, nx, ny);
        const float m = m_grid_view[static_cast<size_t>(gi)];
        if (m <= 0.0F)
        {
          continue;
        }
        const float w = wx.at(static_cast<size_t>(dx)) * wy.at(static_cast<size_t>(dy)) *
                        wz.at(static_cast<size_t>(dz));
        const size_t base3 = static_cast<size_t>(gi) * k_vec3;
        const float vx_i   = v_grid_view[base3 + 0U];
        const float vy_i   = v_grid_view[base3 + 1U];
        const float vz_i   = v_grid_view[base3 + 2U];
        vx += w * vx_i;
        vy += w * vy_i;
        vz += w * vz_i;

        if (out_c != nullptr)
        {
          const float gx_rel = (static_cast<float>(ix) - gx);
          const float gy_rel = (static_cast<float>(iy) - gy);
          const float gz_rel = (static_cast<float>(iz) - gz);
          cx[0] += w * vx_i * gx_rel;
          cx[1] += w * vx_i * gy_rel;
          cx[2] += w * vx_i * gz_rel;
          cy[0] += w * vy_i * gx_rel;
          cy[1] += w * vy_i * gy_rel;
          cy[2] += w * vy_i * gz_rel;
          cz[0] += w * vz_i * gx_rel;
          cz[1] += w * vz_i * gy_rel;
          cz[2] += w * vz_i * gz_rel;
        }
      }
    }
  }

  const size_t out_base       = i * static_cast<size_t>(k_vec3);
  out_vel_view[out_base + 0U] = vx;
  out_vel_view[out_base + 1U] = vy;
  out_vel_view[out_base + 2U] = vz;
  if (out_c != nullptr)
  {
    const float s           = k_four * invh * invh;  // common APIC scaling
    const size_t cbase      = i * static_cast<size_t>(k_cdim);
    (*out_c)[cbase + k_cxx] = s * cx[0];
    (*out_c)[cbase + k_cxy] = s * cx[1];
    (*out_c)[cbase + k_cxz] = s * cx[2];
    (*out_c)[cbase + k_cyx] = s * cy[0];
    (*out_c)[cbase + k_cyy] = s * cy[1];
    (*out_c)[cbase + k_cyz] = s * cy[2];
    (*out_c)[cbase + k_czx] = s * cz[0];
    (*out_c)[cbase + k_czy] = s * cz[1];
    (*out_c)[cbase + k_czz] = s * cz[2];
  }
}

/**
 * @brief APIC/MLS particle-to-grid transfer (CPU reference).
 * @param N Number of particles.
 * @param pos Positions (xyz triplets, length 3N).
 * @param vel Velocities (xyz triplets, length 3N).
 * @param C Affine matrix per particle (row-major, length 9N) or nullptr.
 * @param mass Particle masses (length N).
 * @param origin Grid origin (min corner), length-3 array.
 * @param h Grid spacing.
 * @param nx Grid x size; @param ny Grid y size; @param nz Grid z size.
 * @param m_grid Out mass grid (length nx*ny*nz).
 * @param p_grid Out momentum grid (xyz per node, length 3*nx*ny*nz).
 */
void zx_apic_p2g_ref(size_t N, const float* pos, const float* vel, const float* C,
                     const float* mass, const float origin[3], float h, int nx, int ny, int nz,
                     float* ZX_RESTRICT m_grid, float* ZX_RESTRICT p_grid)
{
  std::memset(m_grid, 0, sizeof(float) * nx * ny * nz);
  std::memset(p_grid, 0, sizeof(float) * 3 * nx * ny * nz);
  constexpr size_t k_neighbors_3d  = 27U;
  const bool deterministic_enabled = (zx_get_determinism() != 0);
  std::vector<Contrib> contribs;
  if (deterministic_enabled)
  {
    contribs.reserve(N * static_cast<size_t>(k_neighbors_3d));
  }
  // Wrap raw arrays in views to avoid scattered pointer arithmetic and keep
  // indexing readable and tidy-friendly.
  // Use non-const view for unified API; underlying pointers remain treated as read-only.
  ArrayView<const float> pos_view(pos, static_cast<size_t>(k_vec3) * N);
  ArrayView<const float> vel_view(vel, static_cast<size_t>(k_vec3) * N);
  ArrayView<const float> mass_view(mass, N);
  ArrayView<float> m_grid_view(m_grid, static_cast<size_t>(nx) * static_cast<size_t>(ny) *
                                           static_cast<size_t>(nz));
  ArrayView<float> p_grid_view(p_grid, static_cast<size_t>(k_vec3) * static_cast<size_t>(nx) *
                                           static_cast<size_t>(ny) * static_cast<size_t>(nz));

  // Copy origin to a fixed-size array to avoid pointer indexing on a pointer parameter.
  std::array<float, 3> origin_arr{};
  std::memcpy(origin_arr.data(), origin, sizeof(float) * 3U);

  std::vector<Contrib>* contribs_ptr = deterministic_enabled ? &contribs : nullptr;
  for (size_t i = 0; i < N; ++i)
  {
    const float mp = mass_view[i];
    accumulate_p2g_particle(i, pos_view, vel_view, C, N, mp, origin_arr, h, nx, ny, nz, m_grid_view,
                            p_grid_view, contribs_ptr);
  }
  if (deterministic_enabled)
  {
    finalize_p2g_contribs(contribs, m_grid_view, p_grid_view);
  }
}

/**
 * @brief APIC/MLS grid-to-particle gather (CPU reference).
 * @param N Number of particles.
 * @param pos Positions (xyz triplets, length 3N).
 * @param out_vel Out velocities (xyz triplets, length 3N).
 * @param out_C Optional out affine matrices (row-major, length 9N) or nullptr.
 * @param origin Grid origin (min corner), length-3 array.
 * @param h Grid spacing.
 * @param nx Grid x size; @param ny Grid y size; @param nz Grid z size.
 * @param m_grid Mass grid (length nx*ny*nz).
 * @param v_grid Velocity grid (xyz per node, length 3*nx*ny*nz).
 */
void zx_apic_g2p_ref(size_t N, const float* pos, float* out_vel, float* out_C,
                     const float origin[3], float h, int nx, int ny, int nz,
                     const float* ZX_RESTRICT m_grid, const float* ZX_RESTRICT v_grid)
{
  // Wrap arrays in views and copy origin to avoid pointer indexing on parameters.
  ArrayView<const float> pos_view(pos, static_cast<size_t>(k_vec3) * N);
  ArrayView<float> out_vel_view(out_vel, static_cast<size_t>(k_vec3) * N);
  ArrayView<float> out_c_view(out_C, (out_C != nullptr) ? static_cast<size_t>(k_cdim) * N : 0U);
  ArrayView<const float> m_grid_view(m_grid, static_cast<size_t>(nx) * static_cast<size_t>(ny) *
                                                 static_cast<size_t>(nz));
  ArrayView<const float> v_grid_view(v_grid, static_cast<size_t>(k_vec3) * static_cast<size_t>(nx) *
                                                 static_cast<size_t>(ny) * static_cast<size_t>(nz));
  std::array<float, 3> origin_arr{};
  std::memcpy(origin_arr.data(), origin, sizeof(float) * 3U);

  ArrayView<float>* out_c_ptr = (out_C != nullptr) ? &out_c_view : nullptr;
  for (size_t i = 0; i < N; ++i)
  {
    accumulate_g2p_particle(i, pos_view, out_vel_view, out_c_ptr, origin_arr, h, nx, ny, nz,
                            m_grid_view, v_grid_view);
  }
}
