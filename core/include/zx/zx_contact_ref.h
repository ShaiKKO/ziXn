/*!
\file zx_contact_ref.h
\brief CPU reference contact projection: non-penetration + Coulomb friction.
\author Colin Macritchie (Ripple Group, LLC)
*/

#ifndef ZX_CONTACT_REF_H
#define ZX_CONTACT_REF_H

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

#endif /* ZX_CONTACT_REF_H */


