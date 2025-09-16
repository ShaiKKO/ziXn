/*!
 * \file zx_contact_ref.cpp
 * \brief CPU reference for grid-node contact projection.
 */

#include "zx/zx_contact_ref.h"
#include <cmath>
#include <algorithm>

static inline float dot3(const float a[3], const float b[3]) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static inline void sub3(const float a[3], const float b[3], float out[3]) { out[0]=a[0]-b[0]; out[1]=a[1]-b[1]; out[2]=a[2]-b[2]; }
static inline void add3s(const float a[3], const float s, const float b[3], float out[3]) { out[0]=a[0]+s*b[0]; out[1]=a[1]+s*b[1]; out[2]=a[2]+s*b[2]; }
static inline float len3(const float a[3]) { return std::sqrt(dot3(a,a)); }

void zx_contact_project(const float v_in[3], const float n[3], float mu, float phi, float kappa_n, float v_out[3], int* out_sat)
{
    // Split into normal and tangential components
    const float vn = dot3(v_in, n);
    float vt[3] = { v_in[0] - vn*n[0], v_in[1] - vn*n[1], v_in[2] - vn*n[2] };

    // Non-penetration: clamp normal velocity; allow small compliance proportional to penetration depth
    float vn_proj = vn;
    if (phi < 0.0f) {
        const float allowance = -phi * kappa_n; // more penetration => stronger clamp
        if (vn_proj < -allowance) vn_proj = -allowance;
    }

    // Coulomb friction on tangent
    const float vt_len = len3(vt);
    int saturated = 0;
    if (vt_len > 0.0f) {
        const float max_t = mu * std::max(0.0f, vn_proj); // use post-clamp normal magnitude
        if (vt_len > max_t) {
            const float s = max_t / vt_len;
            vt[0] *= s; vt[1] *= s; vt[2] *= s;
            saturated = 1;
        }
    }

    v_out[0] = vn_proj*n[0] + vt[0];
    v_out[1] = vn_proj*n[1] + vt[1];
    v_out[2] = vn_proj*n[2] + vt[2];
    if (out_sat) *out_sat = saturated;
}

static inline float dot(const float a[3], const float b[3]) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
static inline void madd3(float* o, const float a[3], float s) { o[0]+=s*a[0]; o[1]+=s*a[1]; o[2]+=s*a[2]; }

extern "C" ZX_API void ZX_CALL zx_contact_project_aniso(
    const float v_in[3], const float n[3], const float t0[3], const float t1[3],
    float mu_long, float mu_lat, float rr_scale, float phi, float kappa_n,
    float v_out[3], int* out_sat)
{
    float vn = dot(v_in, n);
    float v0 = dot(v_in, t0);
    float v1 = dot(v_in, t1);
    float vn_proj = vn;
    if (phi < 0.0f) {
        float allowance = -phi * kappa_n;
        if (vn_proj < -allowance) vn_proj = -allowance;
    }
    // Elliptical Coulomb: (v0/L0)^2 + (v1/L1)^2 <= 1 with L0=mu_long*vn_proj, L1=mu_lat*vn_proj
    float L0 = mu_long * std::max(0.0f, vn_proj);
    float L1 = mu_lat * std::max(0.0f, vn_proj);
    int saturated = 0;
    if (L0 > 0.0f && L1 > 0.0f) {
        float ell = (v0*v0)/(L0*L0) + (v1*v1)/(L1*L1);
        if (ell > 1.0f) {
            float s = 1.0f / std::sqrt(ell);
            v0 *= s; v1 *= s; saturated = 1;
        }
    } else {
        v0 = 0.0f; v1 = 0.0f; saturated = 1;
    }
    // Rolling resistance proxy: shrink tangential limit as |v_t| grows
    float vt_mag = std::sqrt(v0*v0 + v1*v1);
    float rr = 1.0f / (1.0f + rr_scale * vt_mag);
    v0 *= rr; v1 *= rr;

    // Recompose
    v_out[0] = vn_proj*n[0]; v_out[1] = vn_proj*n[1]; v_out[2] = vn_proj*n[2];
    madd3(v_out, t0, v0);
    madd3(v_out, t1, v1);
    if (out_sat) *out_sat = saturated;
}


