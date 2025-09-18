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

// Minimal array view to centralize pointer indexing and localize NOLINT usage.
// This avoids scattering pointer arithmetic across call sites while preserving
// the C ABI in public headers.
template <typename T>
class ArrayView
{
public:
  ArrayView(T* ptr, size_t count) : ptr_(ptr), count_(count) {}
  // NOLINTNEXTLINE(modernize-pass-by-value)
  ArrayView(const T* ptr, size_t count) : ptr_(const_cast<T*>(ptr)), count_(count) {}

  T& operator[](size_t idx) const
  {
    // Bounds are validated by callers where applicable; this centralizes pointer
    // arithmetic to a single location with a focused suppression.
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    return ptr_[idx];
  }

  T* data() const { return ptr_; }
  size_t size() const { return count_; }

private:
  T* ptr_;
  size_t count_;
};

static inline int idx3(int x, int y, int z, int nx, int ny)
{
  // ((z * ny) + y) * nx + x
  return ((z * ny) + y) * nx + x;
}

void zx_bspline_w3(float x, float w[3])
{
  // Public ABI requires pointer-based array. Localize the pointer indexing.
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const float x0 = 0.5F - x;
  const float x1 = 0.5F + x;
  w[0]           = 0.5F * x0 * x0;
  w[1]           = 0.75F - x * x;
  w[2]           = 0.5F * x1 * x1;
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

void zx_bspline_dw3(float x, float g[3])
{
  // Public ABI requires pointer-based array. Localize the pointer indexing.
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const float x0 = 0.5F - x;
  const float x1 = 0.5F + x;
  g[0]           = -x0;
  g[1]           = -2.0F * x;
  g[2]           = x1;
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

void zx_apic_p2g_ref(size_t N, const float* pos, const float* vel, const float* C,
                     const float* mass, const float origin[3], float h, int nx, int ny, int nz,
                     float* ZX_RESTRICT m_grid, float* ZX_RESTRICT p_grid)
{
  std::memset(m_grid, 0, sizeof(float) * nx * ny * nz);
  std::memset(p_grid, 0, sizeof(float) * 3 * nx * ny * nz);
  const float invh        = 1.0F / h;
  constexpr float k_half  = 0.5F;   // NOLINT(cppcoreguidelines-avoid-magic-numbers)
  const int deterministic = zx_get_determinism();
  struct Contrib
  {
    int gi;
    int pid;
    float m;
    float mx;
    float my;
    float mz;
  };
  std::vector<Contrib> contribs;
  if (deterministic)
  {
    contribs.reserve(N * 27);
  }
  // Wrap raw arrays in views to avoid scattered pointer arithmetic and keep
  // indexing readable and tidy-friendly.
  ArrayView<const float> pos_view(pos, static_cast<size_t>(3U) * N);
  ArrayView<const float> vel_view(vel, static_cast<size_t>(3U) * N);
  ArrayView<const float> mass_view(mass, N);
  ArrayView<float>       m_grid_view(m_grid, static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz));
  ArrayView<float>       p_grid_view(p_grid, static_cast<size_t>(3U) * static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz));

  // Copy origin to a fixed-size array to avoid pointer indexing on a pointer parameter.
  std::array<float, 3> origin_arr{};
  std::memcpy(origin_arr.data(), origin, sizeof(float) * 3U);

  for (size_t i = 0; i < N; ++i)
  {
    const size_t pbase = static_cast<size_t>(3U) * i;
    const float px     = pos_view[pbase + 0U];
    const float py     = pos_view[pbase + 1U];
    const float pz     = pos_view[pbase + 2U];
    const float vx0    = vel_view[pbase + 0U];
    const float vy0    = vel_view[pbase + 1U];
    const float vz0    = vel_view[pbase + 2U];

    std::array<float, 9> C_local{};
    bool has_C = false;
    if (C != nullptr)
    {
      ArrayView<const float> C_view(C, static_cast<size_t>(9U) * N);
      const size_t cbase = static_cast<size_t>(9U) * i;
      std::memcpy(C_local.data(), &C_view[cbase], sizeof(float) * 9U);
      has_C = true;
    }
    const float mp = mass_view[i];

    const float gx = (px - origin_arr[0]) * invh;
    const float gy = (py - origin_arr[1]) * invh;
    const float gz = (pz - origin_arr[2]) * invh;

    const int base_x = static_cast<int>(std::floor(gx - k_half));
    const int base_y = static_cast<int>(std::floor(gy - k_half));
    const int base_z = static_cast<int>(std::floor(gz - k_half));

    float wx[3], wy[3], wz[3];
    zx_bspline_w3(gx - (base_x + 1), wx);
    zx_bspline_w3(gy - (base_y + 1), wy);
    zx_bspline_w3(gz - (base_z + 1), wz);

    for (int dz = 0; dz < 3; ++dz)
      for (int dy = 0; dy < 3; ++dy)
        for (int dx = 0; dx < 3; ++dx)
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
          const float w = wx[dx] * wy[dy] * wz[dz];
          const int gi  = idx3(ix, iy, iz, nx, ny);
          if (!deterministic)
          {
            m_grid_view[gi] += w * mp;
          }

          float vx_aff = vx0;
          float vy_aff = vy0;
          float vz_aff = vz0;
          if (has_C)
          {
            const float gx_rel = (static_cast<float>(ix) - gx);
            const float gy_rel = (static_cast<float>(iy) - gy);
            const float gz_rel = (static_cast<float>(iz) - gz);
            // APIC affine: v + C*(x_i - x_p) scaled by h
            vx_aff += (C_local[0] * gx_rel + C_local[1] * gy_rel + C_local[2] * gz_rel) * h;
            vy_aff += (C_local[3] * gx_rel + C_local[4] * gy_rel + C_local[5] * gz_rel) * h;
            vz_aff += (C_local[6] * gx_rel + C_local[7] * gy_rel + C_local[8] * gz_rel) * h;
          }
          if (!deterministic)
          {
            const size_t base3 = static_cast<size_t>(gi) * 3U;
            p_grid_view[base3 + 0U] += w * mp * vx_aff;
            p_grid_view[base3 + 1U] += w * mp * vy_aff;
            p_grid_view[base3 + 2U] += w * mp * vz_aff;
          }
          else
          {
            const float mm = w * mp;
            contribs.push_back(Contrib{gi, (int) i, mm, mm * vx_aff, mm * vy_aff, mm * vz_aff});
          }
        }
  }
  if (deterministic)
  {
    std::stable_sort(contribs.begin(), contribs.end(),
                     [](const Contrib& a, const Contrib& b)
                     {
                       if (a.gi != b.gi)
                         return a.gi < b.gi;
                       return a.pid < b.pid;
                     });
    int cur    = -1;
    float msum = 0.0f, mxsum = 0.0f, mysum = 0.0f, mzsum = 0.0f;
    auto flush = [&]()
    {
      if (cur >= 0)
      {
        m_grid_view[cur] += msum;
        const size_t base3 = static_cast<size_t>(cur) * 3U;
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
        cur  = c.gi;
        msum = mxsum = mysum = mzsum = 0.0f;
      }
      msum += c.m;
      mxsum += c.mx;
      mysum += c.my;
      mzsum += c.mz;
    }
    flush();
  }
}

void zx_apic_g2p_ref(size_t N, const float* pos, float* out_vel, float* out_C,
                     const float origin[3], float h, int nx, int ny, int nz,
                     const float* ZX_RESTRICT m_grid, const float* ZX_RESTRICT v_grid)
{
  const float invh       = 1.0F / h;
  constexpr float k_half = 0.5F;   // NOLINT(cppcoreguidelines-avoid-magic-numbers)
  // Wrap arrays in views and copy origin to avoid pointer indexing on parameters.
  ArrayView<const float> pos_view(pos, static_cast<size_t>(3U) * N);
  ArrayView<float>       out_vel_view(out_vel, static_cast<size_t>(3U) * N);
  ArrayView<float>       out_C_view(out_C, out_C ? static_cast<size_t>(9U) * N : 0U);
  ArrayView<const float> m_grid_view(m_grid, static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz));
  ArrayView<const float> v_grid_view(v_grid, static_cast<size_t>(3U) * static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz));
  std::array<float, 3> origin_arr{};
  std::memcpy(origin_arr.data(), origin, sizeof(float) * 3U);

  for (size_t i = 0; i < N; ++i)
  {
    const size_t pbase = static_cast<size_t>(3U) * i;
    const float px     = pos_view[pbase + 0U];
    const float py     = pos_view[pbase + 1U];
    const float pz     = pos_view[pbase + 2U];

    const float gx = (px - origin_arr[0]) * invh;
    const float gy = (py - origin_arr[1]) * invh;
    const float gz = (pz - origin_arr[2]) * invh;

     const int base_x = static_cast<int>(std::floor(gx - k_half));
     const int base_y = static_cast<int>(std::floor(gy - k_half));
     const int base_z = static_cast<int>(std::floor(gz - k_half));

    float wx[3], wy[3], wz[3];
    zx_bspline_w3(gx - (base_x + 1), wx);
    zx_bspline_w3(gy - (base_y + 1), wy);
    zx_bspline_w3(gz - (base_z + 1), wz);

    float vx = 0, vy = 0, vz = 0;
    float Cx[3] = {0, 0, 0}, Cy[3] = {0, 0, 0}, Cz[3] = {0, 0, 0};

    for (int dz = 0; dz < 3; ++dz)
      for (int dy = 0; dy < 3; ++dy)
        for (int dx = 0; dx < 3; ++dx)
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
          const float w     = wx[dx] * wy[dy] * wz[dz];
          const size_t base3 = static_cast<size_t>(gi) * 3U;
          const float vx_i   = v_grid_view[base3 + 0U];
          const float vy_i   = v_grid_view[base3 + 1U];
          const float vz_i   = v_grid_view[base3 + 2U];
          vx += w * vx_i;
          vy += w * vy_i;
          vz += w * vz_i;

          if (out_C)
          {
             const float gx_rel = (static_cast<float>(ix) - gx);
             const float gy_rel = (static_cast<float>(iy) - gy);
             const float gz_rel = (static_cast<float>(iz) - gz);
            Cx[0] += w * vx_i * gx_rel;
            Cx[1] += w * vx_i * gy_rel;
            Cx[2] += w * vx_i * gz_rel;
            Cy[0] += w * vy_i * gx_rel;
            Cy[1] += w * vy_i * gy_rel;
            Cy[2] += w * vy_i * gz_rel;
            Cz[0] += w * vz_i * gx_rel;
            Cz[1] += w * vz_i * gy_rel;
            Cz[2] += w * vz_i * gz_rel;
          }
        }
    const size_t out_base = static_cast<size_t>(i) * 3U;
    out_vel_view[out_base + 0U] = vx;
    out_vel_view[out_base + 1U] = vy;
    out_vel_view[out_base + 2U] = vz;
    if (out_C)
    {
      const float s      = 4.0F * invh * invh;  // common APIC scaling
      const size_t cbase = static_cast<size_t>(i) * 9U;
      out_C_view[cbase + 0U] = s * Cx[0];
      out_C_view[cbase + 1U] = s * Cx[1];
      out_C_view[cbase + 2U] = s * Cx[2];
      out_C_view[cbase + 3U] = s * Cy[0];
      out_C_view[cbase + 4U] = s * Cy[1];
      out_C_view[cbase + 5U] = s * Cy[2];
      out_C_view[cbase + 6U] = s * Cz[0];
      out_C_view[cbase + 7U] = s * Cz[1];
      out_C_view[cbase + 8U] = s * Cz[2];
    }
  }
}
