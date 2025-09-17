/*!
 * \file zx_contact_ref.cpp
 * \brief CPU reference for grid-node contact projection.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_contact_ref.h"
#include <algorithm>
#include <cmath>

static inline float dot3(const float a[3], const float b[3])
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static inline void sub3(const float a[3], const float b[3], float out[3])
{
  out[0] = a[0] - b[0];
  out[1] = a[1] - b[1];
  out[2] = a[2] - b[2];
}
/**
 * @brief Compute out = a + s * b component-wise for 3D vectors.
 *
 * Performs a scaled addition of the 3-element vector b to a:
 * out[i] = a[i] + s * b[i] for i = 0,1,2.
 *
 * @param a Left-hand 3-element vector.
 * @param s Scalar multiplier applied to b.
 * @param b Right-hand 3-element vector to be scaled and added.
 * @param out Destination 3-element vector receiving the result.
 */
static inline void add3s(const float a[3], const float s, const float b[3], float out[3])
{
  out[0] = a[0] + s * b[0];
  out[1] = a[1] + s * b[1];
  out[2] = a[2] + s * b[2];
}
/**
 * @brief Returns the Euclidean length (L2 norm) of a 3-element vector.
 *
 * @param a 3-element array representing a 3D vector.
 * @return float Euclidean norm (sqrt of the dot product of the vector with itself).
 */
static inline float len3(const float a[3])
{
  return std::sqrt(dot3(a, a));
}

extern "C" /**
            * @brief Projects an incoming velocity onto a contact constraint (normal + Coulomb
            * friction).
            *
            * Projects v_in into a contact frame defined by surface normal n, clamps the normal
            * component to prevent penetration (with optional compliance), applies isotropic
            * Coulomb friction to the tangential component, and writes the resulting velocity
            * to v_out.
            *
            * If phi < 0, the normal component is clamped to at least -(-phi * kappa_n) to
            * allow a small compliant penetration-dependent allowance. Tangential velocity is
            * limited to mu * max(0, vn_post_clamp); if the tangential magnitude exceeds
            * this limit it is scaled down and saturation is indicated.
            *
            * @param v_in  Input velocity in world coordinates (length-3 array).
            * @param n     Contact normal (unit or non-unit; used to split normal/tangent).
            * @param mu    Coefficient of (isotropic) Coulomb friction.
            * @param phi   Signed penetration distance; when < 0 enables penetration-based clamp.
            * @param kappa_n Normal stiffness scale used with phi to compute allowed penetration
            * clamp.
            * @param v_out Output velocity after projection (length-3 array), overwritten by the
            * function.
            * @param out_sat If non-null, set to 1 when tangential friction saturates (velocity was
            * clamped), otherwise 0.
            */
    void ZX_CALL
    zx_contact_project(const float v_in[3], const float n[3], float mu, float phi, float kappa_n,
                       float v_out[3], int* out_sat)
{
  // Split into normal and tangential components
  const float vn = dot3(v_in, n);
  float vt[3]    = {v_in[0] - vn * n[0], v_in[1] - vn * n[1], v_in[2] - vn * n[2]};

  // Non-penetration: clamp normal velocity; allow small compliance proportional to penetration
  // depth
  float vn_proj = vn;
  if (phi < 0.0f)
  {
    const float allowance = -phi * kappa_n;  // more penetration => stronger clamp
    if (vn_proj < -allowance)
      vn_proj = -allowance;
  }

  // Coulomb friction on tangent
  const float vt_len = len3(vt);
  int saturated      = 0;
  if (vt_len > 0.0f)
  {
    const float max_t = mu * std::max(0.0f, vn_proj);  // use post-clamp normal magnitude
    if (vt_len > max_t)
    {
      const float s = max_t / vt_len;
      vt[0] *= s;
      vt[1] *= s;
      vt[2] *= s;
      saturated = 1;
    }
  }

  v_out[0] = vn_proj * n[0] + vt[0];
  v_out[1] = vn_proj * n[1] + vt[1];
  v_out[2] = vn_proj * n[2] + vt[2];
  if (out_sat)
    *out_sat = saturated;
}

static inline float dot(const float a[3], const float b[3])
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
static inline void madd3(float* o, const float a[3], float s)
{
  o[0] += s * a[0];
  o[1] += s * a[1];
  o[2] += s * a[2];
}

extern "C" void ZX_CALL zx_contact_project_aniso(const float v_in[3], const float n[3],
                                                 const float t0[3], const float t1[3],
                                                 float mu_long, float mu_lat, float rr_scale,
                                                 float phi, float kappa_n, float v_out[3],
                                                 int* out_sat)
{
  float vn      = dot(v_in, n);
  float v0      = dot(v_in, t0);
  float v1      = dot(v_in, t1);
  float vn_proj = vn;
  if (phi < 0.0f)
  {
    float allowance = -phi * kappa_n;
    if (vn_proj < -allowance)
      vn_proj = -allowance;
  }
  // Elliptical Coulomb: (v0/L0)^2 + (v1/L1)^2 <= 1 with L0=mu_long*vn_proj, L1=mu_lat*vn_proj
  float L0      = mu_long * std::max(0.0f, vn_proj);
  float L1      = mu_lat * std::max(0.0f, vn_proj);
  int saturated = 0;
  if (L0 > 0.0f && L1 > 0.0f)
  {
    float ell = (v0 * v0) / (L0 * L0) + (v1 * v1) / (L1 * L1);
    if (ell > 1.0f)
    {
      float s = 1.0f / std::sqrt(ell);
      v0 *= s;
      v1 *= s;
      saturated = 1;
    }
  }
  else
  {
    v0        = 0.0f;
    v1        = 0.0f;
    saturated = 1;
  }
  // Rolling resistance proxy: shrink tangential limit as |v_t| grows
  float vt_mag = std::sqrt(v0 * v0 + v1 * v1);
  float rr     = 1.0f / (1.0f + rr_scale * vt_mag);
  v0 *= rr;
  v1 *= rr;

  // Recompose
  v_out[0] = vn_proj * n[0];
  v_out[1] = vn_proj * n[1];
  v_out[2] = vn_proj * n[2];
  madd3(v_out, t0, v0);
  madd3(v_out, t1, v1);
  if (out_sat)
    *out_sat = saturated;
}
