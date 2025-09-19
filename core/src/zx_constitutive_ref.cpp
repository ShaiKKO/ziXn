/*!
 * \file zx_constitutive_ref.cpp
 * \brief CPU reference constitutive: Elastic trial and DP(+cap) return mapping.
 */

#include "zx/zx_constitutive_ref.h"
#include <algorithm>
#include <array>
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
}  // namespace

static inline float deg2rad(float d)
{
  return d * k_pi / k_deg;
}

/**
 * @brief Map Mohr–Coulomb parameters to Drucker–Prager (triaxial compression fit).
 * @param mc Input Mohr–Coulomb parameters.
 * @param out_dp Output Drucker–Prager parameters (alpha, k) populated.
 */
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

/**
 * @brief Compute Lamé parameters from Young's modulus and Poisson ratio.
 * @param ep Elastic parameters.
 * @param lambda Out: first Lamé parameter.
 * @param mu Out: shear modulus.
 */
void zx_elastic_lame(const zx_elastic_params* ep, float* lambda, float* mu)
{
  const float young_e = ep->young_E;
  const float nu      = ep->poisson_nu;
  *mu                 = young_e / (k_two * (1.0F + nu));
  *lambda             = (young_e * nu) / ((1.0F + nu) * (1.0F - k_two * nu));
}

/**
 * @brief Compute the first invariant I1 and the second deviatoric invariant J2.
 * @param s Symmetric Cauchy stress (9 elements, row-major).
 * @param i1_out Out: I1 = trace(s).
 * @param j2_out Out: J2 = 1/2 dev:dev.
 */
void zx_stress_invariants(const float s[zx_mat3_size], float* i1_out, float* j2_out)
{
  // Copy to std::array to avoid pointer arithmetic on C arrays
  std::array<float, zx_mat3_size> s_arr{};
  std::memcpy(s_arr.data(), s, sizeof(float) * zx_mat3_size);
  const float i1v = s_arr[0] + s_arr[4] + s_arr[8];
  const float p   = i1v / k_three;
  // deviatoric = s - pI
  std::array<float, zx_mat3_size> dev{};
  std::memcpy(dev.data(), s, sizeof(float) * zx_mat3_size);
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  // J2 = 1/2 * dev:dev
  float dd = 0.0F;
  for (const float v : dev)
  {
    dd += v * v;
  }
  *i1_out = i1v;
  *j2_out = k_half * dd;
}

/**
 * @brief Drucker–Prager (+ optional cap) return map (reference implementation).
 * @param ep Elastic parameters (reserved for future scaling).
 * @param dp Drucker–Prager parameters.
 * @param cap Optional cap parameters; may be null.
 * @param sigma_trial Trial Cauchy stress (9 elements).
 * @param sigma_out Out: returned stress on yield surface (9 elements).
 * @param delta_gamma_out Out: plastic multiplier magnitude (may be null).
 */
void zx_dp_return_map(const zx_elastic_params* ep, const zx_dp_params* dp, const zx_cap_params* cap,
                      const float sigma_trial[zx_mat3_size], float sigma_out[zx_mat3_size],
                      float* delta_gamma_out)
{
  (void) ep;
  float i1 = 0.0F;
  float j2 = 0.0F;
  zx_stress_invariants(sigma_trial, &i1, &j2);
  float sqrt_j2 = std::sqrt(std::max(0.0F, j2));
  float f       = (dp->alpha * i1) + sqrt_j2 - dp->k;

  // Cap check (compression): if I1 < i1_min, project to cap line in I1
  float i1c = i1;
  if ((cap != nullptr) && (cap->enabled != 0) && (i1 < cap->i1_min))
  {
    i1c = cap->i1_min;
  }
  if ((f <= 0.0F) && (i1 == i1c))
  {
    std::memcpy(sigma_out, sigma_trial, sizeof(float) * k_nine);
    if (delta_gamma_out != nullptr)
    {
      *delta_gamma_out = 0.0F;
    }
    return;
  }

  // Deviatoric decomposition
  std::array<float, zx_mat3_size> dev{};
  const float p = i1 / k_three;
  std::memcpy(dev.data(), sigma_trial, sizeof(float) * zx_mat3_size);
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  float s2 = 0.0F;
  for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
  {
    s2 += dev[idx] * dev[idx];
  }
  const float s_norm = std::sqrt(std::max(k_eps30, s2));

  // Target invariants on the yield surface after cap clamp
  float i1_new      = i1c;
  float sqrt_j2_new = std::max(0.0F, dp->k - (dp->alpha * i1_new));
  float j2_new      = sqrt_j2_new * sqrt_j2_new;

  // Scale deviatoric part to match J2_new
  float scale = 0.0F;
  if (s_norm > 0.0F)
  {
    const float target_s_norm = std::sqrt(k_two * j2_new);
    scale                     = target_s_norm / s_norm;
  }

  const float p_new = i1_new / k_three;
  std::array<float, zx_mat3_size> sigma_proj{};
  for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
  {
    sigma_proj[idx] = scale * dev[idx];
  }
  sigma_proj[0] += p_new;
  sigma_proj[4] += p_new;
  sigma_proj[8] += p_new;

  // Line search to ensure non-negative plastic work with associated flow
  // Flow direction N = ∂f/∂σ ≈ alpha*I + dev/(2*sqrt(J2)+eps)
  const float eps = k_eps6;
  std::array<float, zx_mat3_size> n_vec{};
  // alpha*I contribution
  n_vec[0] += dp->alpha;
  n_vec[4] += dp->alpha;
  n_vec[8] += dp->alpha;
  // deviatoric contribution from trial
  const float denom = k_two * std::max(sqrt_j2, eps);
  for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
  {
    n_vec[idx] += dev[idx] / denom;
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
  std::array<float, zx_mat3_size> sigma_candidate{};
  constexpr int k_max_iters = 10;
  for (int iter = 0; iter < k_max_iters; ++iter)
  {
    for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
    {
      sigma_candidate[idx] = sigma_trial[idx] + step * (sigma_proj[idx] - sigma_trial[idx]);
    }
    std::array<float, zx_mat3_size> d_eps_p{};
    for (int idx = 0; idx < (int) zx_mat3_size; ++idx)
    {
      d_eps_p[idx] = (dgamma * step) * n_vec[idx];
    }
    std::array<float, zx_mat3_size> d_sigma{};
    sub9(sigma_candidate.data(), sigma_trial, d_sigma.data());
    const float plastic_work = dot9(d_sigma.data(), d_eps_p.data());  // plastic work increment
    if (plastic_work >= -k_eps6)
    {
      std::memcpy(sigma_out, sigma_candidate.data(), sizeof(float) * k_nine);
      if (delta_gamma_out != nullptr)
      {
        *delta_gamma_out = dgamma * step;
      }
      return;
    }
    step *= 0.5F;  // backtrack
  }
  // Fallback if line search fails (should not): accept projection
  std::memcpy(sigma_out, sigma_proj.data(), sizeof(float) * k_nine);
  if (delta_gamma_out != nullptr)
  {
    *delta_gamma_out = dgamma * step;
  }
}

/**
 * @brief Approximate eigenvalues of a real symmetric 3x3 matrix using invariants.
 * @param s Input symmetric matrix (row-major 3x3, 9 elements).
 * @param eval Out eigenvalues (unordered).
 */
static void eig3_sym(const float s[9], float eval[3])
{
  // Very small symmetric 3x3 eigenvalue approximation: use invariants + Newton for cubic
  // For validation harnesses only; production path should use robust solver.
  const float i1 = s[0] + s[4] + s[8];
  const float i2 = ((s[0] * s[4]) + (s[4] * s[8]) + (s[8] * s[0])) -
                   ((s[1] * s[3]) + (s[2] * s[6]) + (s[5] * s[7]));
  const float i3 = (s[0] * ((s[4] * s[8]) - (s[5] * s[7]))) -
                   (s[1] * ((s[3] * s[8]) - (s[5] * s[6]))) +
                   (s[2] * ((s[3] * s[7]) - (s[4] * s[6])));
  // Solve cubic λ^3 - I1 λ^2 + I2 λ - I3 = 0 via Cardano (assuming distinct real roots typical for
  // stresses)
  const float a = -i1;
  const float b = i2;
  const float c = -i3;
  const float q = (3 * b - a * a) / 9.0F;
  const float r = (9 * a * b - 27 * c - 2 * a * a * a) / 54.0F;
  const float d = (q * q * q) + (r * r);
  if (d >= 0.0F)
  {
    float sqrt_d = std::sqrt(d);
    float s_val  = std::cbrt(r + sqrt_d);
    float t_val  = std::cbrt(r - sqrt_d);
    eval[0]      = -a / 3.0F + (s_val + t_val);
    eval[1]      = -a / 3.0F - (s_val + t_val) / 2.0F;
    eval[2]      = eval[1];
  }
  else
  {
    float theta = std::acos(r / std::sqrt(-(q * q * q)));
    float rad   = 2.0F * std::sqrt(-q);
    eval[0]     = -a / 3.0F + rad * std::cos(theta / 3.0F);
    eval[1]     = -a / 3.0F + rad * std::cos((theta + 2.0F * k_pi) / 3.0F);
    eval[2]     = -a / 3.0F + rad * std::cos((theta + 4.0F * k_pi) / 3.0F);
  }
}

void zx_mc_return_map(const zx_elastic_params* ep, const zx_mc_params* mc, float dilatancy_deg,
                      const float sigma_trial[zx_mat3_size], float sigma_out[zx_mat3_size],
                      float* delta_gamma_out, float* eps_v_pl)
{
  (void) ep;  // reserved for scaling in full implementation
  // Principal stresses (approximate) and invariants
  float lam[3];
  eig3_sym(sigma_trial, lam);
  // Sort descending
  std::sort(lam, lam + 3, std::greater<>());

  const float phi  = deg2rad(mc->friction_deg);
  const float psi  = deg2rad(dilatancy_deg);
  const float c_pa = mc->cohesion_kpa * 1000.0F;

  // MC yield in principal space: tau_max + alpha p - c <= 0
  const float p       = (lam[0] + lam[1] + lam[2]) / 3.0F;
  const float tau_max = 0.5F * (lam[0] - lam[2]);
  const float alpha   = std::sin(phi) / (1.0F - std::sin(phi));
  float f =
      (tau_max + (alpha * (-p))) - c_pa;  // note p is compression negative if using sign convention
  if (f <= 0.0F)
  {
    std::memcpy(sigma_out, sigma_trial, sizeof(float) * 9);
    if (delta_gamma_out != nullptr)
    {
      *delta_gamma_out = 0.0F;
    }
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
  float p_new = p - (gamma * alpha_psi);
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

float zx_norsand_e_cs(float p_mean, const zx_norsand_params* ns)
{
  // Ensure pressure positive for log; use absolute
  float p     = std::max(1.0F, std::fabs(p_mean) + ns->p_ref);
  float ratio = p / ns->p_ref;
  float term  = std::pow(std::log(ratio), ns->n_exp);
  return ns->e_ref - (ns->lambda_cs * term);
}

float zx_norsand_state_parameter(const zx_norsand_state* st, float p_mean,
                                 const zx_norsand_params* ns)
{
  return st->void_ratio_e - zx_norsand_e_cs(p_mean, ns);
}

void zx_norsand_return_map(const zx_elastic_params* ep, const zx_norsand_params* ns,
                           zx_norsand_state* state, const float sigma_trial[zx_mat3_size],
                           float sigma_out[zx_mat3_size], float* delta_gamma_out, float* eps_v_pl)
{
  (void) ep;
  float i1_val = 0.0F;
  float j2_val = 0.0F;
  zx_stress_invariants(sigma_trial, &i1_val, &j2_val);
  float p   = i1_val / 3.0F;
  float q   = std::sqrt(std::max(0.0F, 3.0F * j2_val));
  float psi = zx_norsand_state_parameter(state, p, ns);

  // Target critical stress ratio q = M * (-p)
  float q_target = ns->M * (-p);
  float t        = 0.5F;  // relaxation factor
  float q_new    = ((1.0F - t) * q) + (t * q_target);
  float j2_new   = (q_new * q_new) / 3.0F;

  // Deviatoric scaling
  float dev[9];
  for (int i = 0; i < 9; ++i)
  {
    dev[i] = sigma_trial[i];
  }
  dev[0] -= p;
  dev[4] -= p;
  dev[8] -= p;
  float s2 = 0.0F;
  for (float v : dev)
  {
    s2 += v * v;
  }
  float s_norm        = std::sqrt(std::max(k_eps30, s2));
  float target_s_norm = std::sqrt(2.0F * j2_new);
  float scale         = target_s_norm / s_norm;
  for (int i = 0; i < 9; ++i)
  {
    sigma_out[i] = scale * dev[i];
  }
  sigma_out[0] += p;
  sigma_out[4] += p;
  sigma_out[8] += p;

  // Volumetric update from dilatancy: dε_v^pl ≈ k * ψ * dγ
  float dgamma    = std::fabs(q - q_new);
  float deps_v_pl = ns->dilatancy_scale * psi * dgamma;
  if (eps_v_pl != nullptr)
  {
    *eps_v_pl += deps_v_pl;
  }
  // Update void ratio with volumetric change proxy: e_new = e * exp(deps_v)
  state->void_ratio_e = state->void_ratio_e * std::exp(deps_v_pl);
  if (delta_gamma_out != nullptr)
  {
    *delta_gamma_out = dgamma;
  }
}
