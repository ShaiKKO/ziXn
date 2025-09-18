/*!
 * \file zx_constitutive_ref.cpp
 * \brief CPU reference constitutive: Elastic trial and DP(+cap) return mapping.
 */

#include "zx/zx_constitutive_ref.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>

namespace
{
  constexpr float k_pi    = 3.14159265358979323846F;
  constexpr float k_deg   = 180.0F;
  constexpr float k_sqrt3 = 1.7320508075688772F;
  constexpr float k_eps30 = 1e-30F;
  constexpr float k_eps6  = 1e-6F;
  constexpr float k_half  = 0.5F;
  constexpr float k_two   = 2.0F;
  constexpr float k_three = 3.0F;
  constexpr float k_six   = 6.0F;
  constexpr int k_nine    = 9;
  constexpr int k_eight   = 8;
  constexpr int k_ten     = 10;
}  // namespace

static inline float deg2rad(float d)
{
  return d * k_pi / k_deg;
}

void zx_mc_to_dp(const zx_mc_params* mc, zx_dp_params* out_dp)
{
  const float phi   = deg2rad(mc->friction_deg);
  const float c_pa  = mc->cohesion_kpa * 1000.0F;
  const float sinp  = std::sin(phi);
  const float sqrt3 = k_sqrt3;
  // Common MC->DP mapping (triaxial compression fit)
  out_dp->alpha = (k_two * sinp) / (sqrt3 * (k_three - sinp));
  out_dp->k     = (k_six * c_pa * std::cos(phi)) / (sqrt3 * (k_three - sinp));
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void zx_elastic_lame(const zx_elastic_params* ep, float* lambda, float* mu)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
  const float young_e = ep->young_E;
  const float nu      = ep->poisson_nu;
  *mu                 = young_e / (k_two * (1.0F + nu));
  *lambda             = (young_e * nu) / ((1.0F + nu) * (1.0F - k_two * nu));
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void zx_stress_invariants(const float s[zx_mat3_size], float* I1, float* J2)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
  const float i1v = s[0] + s[4] + s[8];
  float p         = i1v / k_three;
  // deviatoric = s - pI
  float dev[zx_mat3_size];
  std::memcpy(dev, s, sizeof(dev));
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  // J2 = 1/2 * dev:dev
  float dd = 0.0F;
  for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
  {
    dd += dev[idx] * dev[idx];
  }
  *I1 = i1v;
  *J2 = k_half * dd;
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void zx_dp_return_map(const zx_elastic_params* ep, const zx_dp_params* dp, const zx_cap_params* cap,
                      const float sigma_trial[zx_mat3_size], float sigma_out[zx_mat3_size],
                      float* delta_gamma_out)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
  float I1, J2;
  zx_stress_invariants(sigma_trial, &I1, &J2);
  float sqrtJ2 = std::sqrt(std::max(0.0F, J2));
  float f      = dp->alpha * I1 + sqrtJ2 - dp->k;

  // Cap check (compression): if I1 < i1_min, project to cap line in I1
  float I1c = I1;
  if ((cap != nullptr) && (cap->enabled != 0) && (I1 < cap->i1_min))
  {
    I1c = cap->i1_min;
  }
  if (f <= 0.0F && I1 == I1c)
  {
    std::memcpy(sigma_out, sigma_trial, sizeof(float) * 9);
    if (delta_gamma_out != nullptr)
    {
      *delta_gamma_out = 0.0F;
    }
    return;
  }

  // Deviatoric decomposition
  float dev[9];
  float p = I1 / k_three;
  for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
  {
    dev[idx] = sigma_trial[idx];
  }
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  float s2 = 0.0F;
  for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
  {
    s2 += dev[idx] * dev[idx];
  }
  float s_norm = std::sqrt(std::max(k_eps30, s2));

  // Target invariants on the yield surface after cap clamp
  float I1_new     = I1c;
  float sqrtJ2_new = std::max(0.0F, dp->k - dp->alpha * I1_new);
  float J2_new     = sqrtJ2_new * sqrtJ2_new;

  // Scale deviatoric part to match J2_new
  float scale = 0.0F;
  if (s_norm > 0.0F)
  {
    float target_s_norm = std::sqrt(k_two * J2_new);
    scale               = target_s_norm / s_norm;
  }

  float p_new = I1_new / k_three;
  float sigma_proj[zx_mat3_size];
  for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
  {
    sigma_proj[idx] = scale * dev[idx];
  }
  sigma_proj[0] += p_new;
  sigma_proj[4] += p_new;
  sigma_proj[8] += p_new;

  // Line search to ensure non-negative plastic work with associated flow
  // Flow direction N = ∂f/∂σ ≈ alpha*I + dev/(2*sqrt(J2)+eps)
  const float eps       = k_eps6;
  float N[zx_mat3_size] = {0.0F};
  // alpha*I contribution
  N[0] += dp->alpha;
  N[4] += dp->alpha;
  N[8] += dp->alpha;
  // deviatoric contribution from trial
  float denom = k_two * std::max(sqrtJ2, eps);
  for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
  {
    N[idx] += dev[idx] / denom;
  }

  // Plastic multiplier proxy from distance to yield
  float dgamma = std::max(0.0F, f);

  auto dot9 = [](const float* a, const float* b)
  {
    float s = 0.0F;
    for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
    {
      s += a[idx] * b[idx];
    }
    return s;
  };
  auto sub9 = [](const float* a, const float* b, float* o)
  {
    for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
    {
      o[idx] = a[idx] - b[idx];
    }
  };

  float step = 1.0F;
  float sigma_candidate[zx_mat3_size];
  constexpr int k_max_iters = 10;
  for (int iter = 0; iter < k_max_iters; ++iter)
  {
    for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
    {
      sigma_candidate[idx] = sigma_trial[idx] + step * (sigma_proj[idx] - sigma_trial[idx]);
    }
    float d_eps_p[zx_mat3_size];
    for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
    {
      d_eps_p[idx] = (dgamma * step) * N[idx];
    }
    float d_sigma[zx_mat3_size];
    sub9(sigma_candidate, sigma_trial, d_sigma);
    float W = dot9(d_sigma, d_eps_p);  // plastic work increment
    if (W >= -k_eps6)
    {
      std::memcpy(sigma_out, sigma_candidate, sizeof(float) * 9);
      if (delta_gamma_out != nullptr)
      {
        *delta_gamma_out = dgamma * step;
      }
      return;
    }
    step *= 0.5F;  // backtrack
  }
  // Fallback if line search fails (should not): accept projection
  std::memcpy(sigma_out, sigma_proj, sizeof(float) * 9);
  if (delta_gamma_out != nullptr)
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
  const float Q = (3 * b - a * a) / 9.0F;
  const float R = (9 * a * b - 27 * c - 2 * a * a * a) / 54.0F;
  const float D = Q * Q * Q + R * R;
  if (D >= 0.0F)
  {
    float sqrtD = std::sqrt(D);
    float S     = std::cbrt(R + sqrtD);
    float T     = std::cbrt(R - sqrtD);
    eval[0]     = -a / 3.0F + (S + T);
    eval[1]     = -a / 3.0F - (S + T) / 2.0F;
    eval[2]     = eval[1];
  }
  else
  {
    float theta = std::acos(R / std::sqrt(-Q * Q * Q));
    float r     = 2.0F * std::sqrt(-Q);
    eval[0]     = -a / 3.0F + r * std::cos(theta / 3.0F);
    eval[1]     = -a / 3.0F + r * std::cos((theta + 2.0F * k_pi) / 3.0F);
    eval[2]     = -a / 3.0F + r * std::cos((theta + 4.0F * k_pi) / 3.0F);
  }
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void zx_mc_return_map(const zx_elastic_params* ep, const zx_mc_params* mc, float dilatancy_deg,
                      const float sigma_trial[zx_mat3_size], float sigma_out[zx_mat3_size],
                      float* delta_gamma_out, float* eps_v_pl)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
  (void) ep;  // reserved for scaling in full implementation
  // Principal stresses (approximate) and invariants
  float lam[3];
  eig3_sym(sigma_trial, lam);
  // Sort descending
  std::sort(lam, lam + 3, std::greater<float>());

  const float phi  = deg2rad(mc->friction_deg);
  const float psi  = deg2rad(dilatancy_deg);
  const float c_pa = mc->cohesion_kpa * 1000.0F;

  // MC yield in principal space: tau_max + alpha p - c <= 0
  const float p       = (lam[0] + lam[1] + lam[2]) / 3.0F;
  const float tau_max = 0.5F * (lam[0] - lam[2]);
  const float alpha   = std::sin(phi) / (1.0F - std::sin(phi));
  float f =
      tau_max + alpha * (-p) - c_pa;  // note p is compression negative if using sign convention
  if (f <= 0.0F)
  {
    std::memcpy(sigma_out, sigma_trial, sizeof(float) * 9);
    if (delta_gamma_out)
      *delta_gamma_out = 0.0F;
    return;
  }

  // Non-associated flow: use dilatancy ψ for plastic potential
  const float alpha_psi = std::sin(psi) / (1.0F - std::sin(psi));

  // Simple radial reduction of deviatoric gap and volumetric shift to satisfy f=0 (approx)
  float scale = std::max(0.0F, (tau_max + alpha * (-p) - c_pa) / std::max(k_eps6, tau_max));
  float lam_new[3];
  lam_new[0] = lam[0] - scale * 0.5F * (lam[0] - p);
  lam_new[2] = lam[2] + scale * 0.5F * (p - lam[2]);
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
  if (eps_v_pl != nullptr)
  {
    *eps_v_pl += -gamma * 3.0F * alpha_psi;
  }
  if (delta_gamma_out != nullptr)
  {
    *delta_gamma_out = gamma;
  }
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
float zx_norsand_e_cs(float p_mean, const zx_norsand_params* ns)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
  // Ensure pressure positive for log; use absolute
  float p     = std::max(1.0F, std::fabs(p_mean) + ns->p_ref);
  float ratio = p / ns->p_ref;
  float term  = std::pow(std::log(ratio), ns->n_exp);
  return ns->e_ref - ns->lambda_cs * term;
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
float zx_norsand_state_parameter(const zx_norsand_state* st, float p_mean,
                                 const zx_norsand_params* ns)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
  return st->void_ratio_e - zx_norsand_e_cs(p_mean, ns);
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
void zx_norsand_return_map(const zx_elastic_params* ep, const zx_norsand_params* ns,
                           zx_norsand_state* state, const float sigma_trial[zx_mat3_size],
                           float sigma_out[zx_mat3_size], float* delta_gamma_out, float* eps_v_pl)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
  (void) ep;
  float I1, J2;
  zx_stress_invariants(sigma_trial, &I1, &J2);
  float p   = I1 / 3.0F;
  float q   = std::sqrt(std::max(0.0F, 3.0F * J2));
  float psi = zx_norsand_state_parameter(state, p, ns);

  // Target critical stress ratio q = M * (-p)
  float q_target = ns->M * (-p);
  float t        = 0.5F;  // relaxation factor
  float q_new    = (1.0F - t) * q + t * q_target;
  float J2_new   = (q_new * q_new) / 3.0F;

  // Deviatoric scaling
  float dev[9];
  for (int i = 0; i < 9; ++i)
    dev[i] = sigma_trial[i];
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  float s2 = 0.0F;
  for (int i = 0; i < 9; ++i)
  {
    s2 += dev[i] * dev[i];
  }
  float s_norm        = std::sqrt(std::max(k_eps30, s2));
  float target_s_norm = std::sqrt(2.0F * J2_new);
  float scale         = target_s_norm / s_norm;
  for (int i = 0; i < 9; ++i)
    sigma_out[i] = scale * dev[i];
  sigma_out[0] += p;
  sigma_out[4] += p;
  sigma_out[8] += p;

  // Volumetric update from dilatancy: dε_v^pl ≈ k * ψ * dγ
  float dgamma    = std::fabs(q - q_new);
  float deps_v_pl = ns->dilatancy_scale * psi * dgamma;
  if (eps_v_pl != nullptr)
    *eps_v_pl += deps_v_pl;
  // Update void ratio with volumetric change proxy: e_new = e * exp(deps_v)
  state->void_ratio_e = state->void_ratio_e * std::exp(deps_v_pl);
  if (delta_gamma_out != nullptr)
    *delta_gamma_out = dgamma;
}
