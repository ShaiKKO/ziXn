/*!
 * \file zx_contact_ref.cpp
 * \brief CPU reference for grid-node contact projection.
 */

#include "zx/zx_contact_ref.h"
#include <cmath>

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


