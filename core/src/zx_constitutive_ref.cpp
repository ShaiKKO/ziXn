/*!
 * \file zx_constitutive_ref.cpp
 * \brief CPU reference constitutive: Elastic trial and DP(+cap) return mapping.
 */

#include "zx/zx_constitutive_ref.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>

static inline float deg2rad(float d)
{
  return d * 3.14159265358979323846f / 180.0f;
}

void zx_mc_to_dp(const zx_mc_params* mc, zx_dp_params* out_dp)
{
  const float phi   = deg2rad(mc->friction_deg);
  const float c_pa  = mc->cohesion_kpa * 1000.0f;
  const float sinp  = std::sin(phi);
  const float sqrt3 = 1.7320508075688772f;
  // Common MC->DP mapping (triaxial compression fit)
  out_dp->alpha = (2.0f * sinp) / (sqrt3 * (3.0f - sinp));
  out_dp->k     = (6.0f * c_pa * std::cos(phi)) / (sqrt3 * (3.0f - sinp));
}

void zx_elastic_lame(const zx_elastic_params* ep, float* lambda, float* mu)
{
  const float E  = ep->young_E;
  const float nu = ep->poisson_nu;
  *mu            = E / (2.0f * (1.0f + nu));
  *lambda        = (E * nu) / ((1.0f + nu) * (1.0f - 2.0f * nu));
}

void zx_stress_invariants(const float s[ZX_MAT3_SIZE], float* I1, float* J2)
{
  const float I1v = s[0] + s[4] + s[8];
  float p         = I1v / 3.0f;
  // deviatoric = s - pI
  float dev[ZX_MAT3_SIZE];
  std::memcpy(dev, s, sizeof(dev));
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  // J2 = 1/2 * dev:dev
  float dd = 0.0f;
  for (int i = 0; i < 9; ++i)
    dd += dev[i] * dev[i];
  *I1 = I1v;
  *J2 = 0.5f * dd;
}

void zx_dp_return_map(const zx_elastic_params* ep, const zx_dp_params* dp, const zx_cap_params* cap,
                      const float sigma_trial[ZX_MAT3_SIZE], float sigma_out[ZX_MAT3_SIZE],
                      float* delta_gamma_out)
{
  float I1, J2;
  zx_stress_invariants(sigma_trial, &I1, &J2);
  float sqrtJ2 = std::sqrt(std::max(0.0f, J2));
  float f      = dp->alpha * I1 + sqrtJ2 - dp->k;

  // Cap check (compression): if I1 < i1_min, project to cap line in I1
  float I1c = I1;
  if (cap && cap->enabled && I1 < cap->i1_min)
  {
    I1c = cap->i1_min;
  }
  if (f <= 0.0f && I1 == I1c)
  {
    std::memcpy(sigma_out, sigma_trial, sizeof(float) * 9);
    if (delta_gamma_out)
      *delta_gamma_out = 0.0f;
    return;
  }

  // Deviatoric decomposition
  float dev[9];
  float p = I1 / 3.0f;
  for (int i = 0; i < (int) ZX_MAT3_SIZE; ++i)
    dev[i] = sigma_trial[i];
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  float s2 = 0.0f;
  for (int i = 0; i < (int) ZX_MAT3_SIZE; ++i)
    s2 += dev[i] * dev[i];
  float s_norm = std::sqrt(std::max(1e-30f, s2));

  // Target invariants on the yield surface after cap clamp
  float I1_new     = I1c;
  float sqrtJ2_new = std::max(0.0f, dp->k - dp->alpha * I1_new);
  float J2_new     = sqrtJ2_new * sqrtJ2_new;

  // Scale deviatoric part to match J2_new
  float scale = 0.0f;
  if (s_norm > 0.0f)
  {
    float target_s_norm = std::sqrt(2.0f * J2_new);
    scale               = target_s_norm / s_norm;
  }

  float p_new = I1_new / 3.0f;
  float sigma_proj[ZX_MAT3_SIZE];
  for (int i = 0; i < (int) ZX_MAT3_SIZE; ++i)
    sigma_proj[i] = scale * dev[i];
  sigma_proj[0] += p_new;
  sigma_proj[4] += p_new;
  sigma_proj[8] += p_new;

  // Line search to ensure non-negative plastic work with associated flow
  // Flow direction N = ∂f/∂σ ≈ alpha*I + dev/(2*sqrt(J2)+eps)
  const float eps       = 1e-6f;
  float N[ZX_MAT3_SIZE] = {0};
  // alpha*I contribution
  N[0] += dp->alpha;
  N[4] += dp->alpha;
  N[8] += dp->alpha;
  // deviatoric contribution from trial
  float denom = 2.0f * std::max(sqrtJ2, eps);
  for (int i = 0; i < (int) ZX_MAT3_SIZE; ++i)
    N[i] += dev[i] / denom;

  // Plastic multiplier proxy from distance to yield
  float dgamma = std::max(0.0f, f);

  auto dot9 = [](const float* a, const float* b)
  {
    float s = 0;
    for (int i = 0; i < (int) ZX_MAT3_SIZE; ++i)
      s += a[i] * b[i];
    return s;
  };
  auto sub9 = [](const float* a, const float* b, float* o)
  {
    for (int i = 0; i < (int) ZX_MAT3_SIZE; ++i)
      o[i] = a[i] - b[i];
  };

  float step = 1.0f;
  float sigma_candidate[ZX_MAT3_SIZE];
  for (int iter = 0; iter < 10; ++iter)
  {
    for (int i = 0; i < (int) ZX_MAT3_SIZE; ++i)
      sigma_candidate[i] = sigma_trial[i] + step * (sigma_proj[i] - sigma_trial[i]);
    float d_eps_p[ZX_MAT3_SIZE];
    for (int i = 0; i < (int) ZX_MAT3_SIZE; ++i)
      d_eps_p[i] = (dgamma * step) * N[i];
    float d_sigma[ZX_MAT3_SIZE];
    sub9(sigma_candidate, sigma_trial, d_sigma);
    float W = dot9(d_sigma, d_eps_p);  // plastic work increment
    if (W >= -1e-6f)
    {
      std::memcpy(sigma_out, sigma_candidate, sizeof(float) * 9);
      if (delta_gamma_out)
        *delta_gamma_out = dgamma * step;
      return;
    }
    step *= 0.5f;  // backtrack
  }
  // Fallback if line search fails (should not): accept projection
  std::memcpy(sigma_out, sigma_proj, sizeof(float) * 9);
  if (delta_gamma_out)
    *delta_gamma_out = dgamma * step;
}

static void eig3_sym(const float s[9], float eval[3])
{
  // Very small symmetric 3x3 eigenvalue approximation: use invariants + Newton for cubic
  // For validation harnesses only; production path should use robust solver.
  const float I1 = s[0] + s[4] + s[8];
  const float I2 =
      s[0] * s[4] + s[4] * s[8] + s[8] * s[0] - (s[1] * s[3] + s[2] * s[6] + s[5] * s[7]);
  const float I3 = s[0] * (s[4] * s[8] - s[5] * s[7]) - s[1] * (s[3] * s[8] - s[5] * s[6]) +
                   s[2] * (s[3] * s[7] - s[4] * s[6]);
  // Solve cubic λ^3 - I1 λ^2 + I2 λ - I3 = 0 via Cardano (assuming distinct real roots typical for
  // stresses)
  const float a = -I1;
  const float b = I2;
  const float c = -I3;
  const float Q = (3 * b - a * a) / 9.0f;
  const float R = (9 * a * b - 27 * c - 2 * a * a * a) / 54.0f;
  const float D = Q * Q * Q + R * R;
  if (D >= 0.0f)
  {
    float sqrtD = std::sqrt(D);
    float S     = std::cbrt(R + sqrtD);
    float T     = std::cbrt(R - sqrtD);
    eval[0]     = -a / 3.0f + (S + T);
    eval[1]     = -a / 3.0f - (S + T) / 2.0f;
    eval[2]     = eval[1];
  }
  else
  {
    float theta = std::acos(R / std::sqrt(-Q * Q * Q));
    float r     = 2.0f * std::sqrt(-Q);
    eval[0]     = -a / 3.0f + r * std::cos(theta / 3.0f);
    eval[1]     = -a / 3.0f + r * std::cos((theta + 2.0f * 3.14159265358979323846f) / 3.0f);
    eval[2]     = -a / 3.0f + r * std::cos((theta + 4.0f * 3.14159265358979323846f) / 3.0f);
  }
}

void zx_mc_return_map(const zx_elastic_params* ep, const zx_mc_params* mc, float dilatancy_deg,
                      const float sigma_trial[ZX_MAT3_SIZE], float sigma_out[ZX_MAT3_SIZE],
                      float* delta_gamma_out, float* eps_v_pl)
{
  (void) ep;  // reserved for scaling in full implementation
  // Principal stresses (approximate) and invariants
  float lam[3];
  eig3_sym(sigma_trial, lam);
  // Sort descending
  std::sort(lam, lam + 3, std::greater<float>());

  const float phi  = deg2rad(mc->friction_deg);
  const float psi  = deg2rad(dilatancy_deg);
  const float c_pa = mc->cohesion_kpa * 1000.0f;

  // MC yield in principal space: tau_max + alpha p - c <= 0
  const float p       = (lam[0] + lam[1] + lam[2]) / 3.0f;
  const float tau_max = 0.5f * (lam[0] - lam[2]);
  const float alpha   = std::sin(phi) / (1.0f - std::sin(phi));
  float f =
      tau_max + alpha * (-p) - c_pa;  // note p is compression negative if using sign convention
  if (f <= 0.0f)
  {
    std::memcpy(sigma_out, sigma_trial, sizeof(float) * 9);
    if (delta_gamma_out)
      *delta_gamma_out = 0.0f;
    return;
  }

  // Non-associated flow: use dilatancy ψ for plastic potential
  const float alpha_psi = std::sin(psi) / (1.0f - std::sin(psi));

  // Simple radial reduction of deviatoric gap and volumetric shift to satisfy f=0 (approx)
  float scale = std::max(0.0f, (tau_max + alpha * (-p) - c_pa) / std::max(1e-6f, tau_max));
  float lam_new[3];
  lam_new[0] = lam[0] - scale * 0.5f * (lam[0] - p);
  lam_new[2] = lam[2] + scale * 0.5f * (p - lam[2]);
  lam_new[1] = lam[1];

  // Adjust volumetric part using dilatancy to reduce f; p_new = p - gamma * alpha_psi
  float gamma = f;  // proxy magnitude
  float p_new = p - gamma * alpha_psi;
  float shift = p_new - p;
  lam_new[0] += shift;
  lam_new[1] += shift;
  lam_new[2] += shift;

  // Recompose as diagonal Cauchy (no rotation change for reference path)
  std::memset(sigma_out, 0, sizeof(float) * 9);
  sigma_out[0] = lam_new[0];
  sigma_out[4] = lam_new[1];
  sigma_out[8] = lam_new[2];
  if (eps_v_pl)
    *eps_v_pl += -gamma * 3.0f * alpha_psi;
  if (delta_gamma_out)
    *delta_gamma_out = gamma;
}

float zx_norsand_e_cs(float p_mean, const zx_norsand_params* ns)
{
  // Ensure pressure positive for log; use absolute
  float p     = std::max(1.0f, std::fabs(p_mean) + ns->p_ref);
  float ratio = p / ns->p_ref;
  float term  = std::pow(std::log(ratio), ns->n_exp);
  return ns->e_ref - ns->lambda_cs * term;
}

float zx_norsand_state_parameter(const zx_norsand_state* st, float p_mean,
                                 const zx_norsand_params* ns)
{
  return st->void_ratio_e - zx_norsand_e_cs(p_mean, ns);
}

void zx_norsand_return_map(const zx_elastic_params* ep, const zx_norsand_params* ns,
                           zx_norsand_state* state, const float sigma_trial[ZX_MAT3_SIZE],
                           float sigma_out[ZX_MAT3_SIZE], float* delta_gamma_out, float* eps_v_pl)
{
  (void) ep;
  float I1, J2;
  zx_stress_invariants(sigma_trial, &I1, &J2);
  float p   = I1 / 3.0f;
  float q   = std::sqrt(std::max(0.0f, 3.0f * J2));
  float psi = zx_norsand_state_parameter(state, p, ns);

  // Target critical stress ratio q = M * (-p)
  float q_target = ns->M * (-p);
  float t        = 0.5f;  // relaxation factor
  float q_new    = (1.0f - t) * q + t * q_target;
  float J2_new   = (q_new * q_new) / 3.0f;

  // Deviatoric scaling
  float dev[9];
  for (int i = 0; i < 9; ++i)
    dev[i] = sigma_trial[i];
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  float s2 = 0.0f;
  for (int i = 0; i < 9; ++i)
    s2 += dev[i] * dev[i];
  float s_norm        = std::sqrt(std::max(1e-30f, s2));
  float target_s_norm = std::sqrt(2.0f * J2_new);
  float scale         = target_s_norm / s_norm;
  for (int i = 0; i < 9; ++i)
    sigma_out[i] = scale * dev[i];
  sigma_out[0] += p;
  sigma_out[4] += p;
  sigma_out[8] += p;

  // Volumetric update from dilatancy: dε_v^pl ≈ k * ψ * dγ
  float dgamma    = std::fabs(q - q_new);
  float deps_v_pl = ns->dilatancy_scale * psi * dgamma;
  if (eps_v_pl)
    *eps_v_pl += deps_v_pl;
  // Update void ratio with volumetric change proxy: e_new = e * exp(deps_v)
  state->void_ratio_e = state->void_ratio_e * std::exp(deps_v_pl);
  if (delta_gamma_out)
    *delta_gamma_out = dgamma;
}
