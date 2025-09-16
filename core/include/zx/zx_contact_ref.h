/*!
\file zx_contact_ref.h
\brief CPU reference contact projection: non-penetration + Coulomb friction.
\author Colin Macritchie (Ripple Group, LLC)
*/

#ifndef ZX_CONTACT_REF_H
#define ZX_CONTACT_REF_H

#include "zx_abi.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Project velocity v by contact normal n (unit), friction mu, and signed distance phi (phi<0 penetrates).
 * Compliance kappa_n >= 0 allows small negative normal velocity proportional to penetration.
 */
void zx_contact_project(
    const float v_in[3],     /* input velocity */
    const float n[3],        /* unit normal (out of surface) */
    float mu,                /* Coulomb coefficient */
    float phi,               /* signed distance; <0 penetrates */
    float kappa_n,           /* normal compliance (s^2/m) proxy; 0 => hard clamp */
    float v_out[3],          /* projected velocity */
    int* out_saturated       /* 1 if friction cone saturated */
);

/* Anisotropic tangential contact: ellipse limits along two tangent axes t0/t1 (orthonormal).
 * Rolling resistance proxy reduces admissible tangential magnitude by (1 / (1 + rr_scale*|vt|)).
 */
ZX_API void ZX_CALL zx_contact_project_aniso(
    const float v_in[3],     /* input velocity */
    const float n[3],        /* unit normal */
    const float t0[3],       /* tangent axis 0 (longitudinal) */
    const float t1[3],       /* tangent axis 1 (lateral) */
    float mu_long,           /* friction along t0 */
    float mu_lat,            /* friction along t1 */
    float rr_scale,          /* rolling resistance scale (>=0) */
    float phi,               /* signed distance; <0 penetrates */
    float kappa_n,           /* compliance */
    float v_out[3],          /* projected velocity */
    int* out_saturated
);
#ifdef __cplusplus
}
#endif

#endif /* ZX_CONTACT_REF_H */


