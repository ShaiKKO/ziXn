/*!
\file zx_contact_ref.h
\brief CPU reference contact projection: non-penetration + Coulomb friction.
\author Colin Macritchie (Ripple Group, LLC)
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

#ifndef ZX_CONTACT_REF_H
#define ZX_CONTACT_REF_H

#include "zx_abi.h"
#ifdef __cplusplus
extern "C"
{
#endif
  /** \brief Project velocity by contact normal and Coulomb friction with compliance.
   * @param v_in Input velocity
   * @param n Unit normal (out of surface)
   * @param mu Coulomb friction coefficient
   * @param phi Signed distance (<0 penetrates)
   * @param kappa_n Normal compliance (s^2/m); 0=hard clamp
   * @param v_out Out projected velocity
   * @param out_saturated Out flag: 1 if friction cone saturated
   */
  ZX_API void ZX_CALL zx_contact_project(
      const float v_in[3], /* input velocity */
      const float n[3],    /* unit normal (out of surface) */
      float mu,            /* Coulomb coefficient */
      float phi,           /* signed distance; <0 penetrates */
      float kappa_n,       /* normal compliance (s^2/m) proxy; 0 => hard clamp */
      float v_out[3],      /* projected velocity */
      int* out_saturated   /* 1 if friction cone saturated */
  );

  /** \brief Anisotropic tangential contact with ellipse limits along two tangent axes.
   * Rolling resistance reduces admissible tangential magnitude as |vt| grows.
   */
  ZX_API void ZX_CALL
  zx_contact_project_aniso(const float v_in[3], /* input velocity */
                           const float n[3],    /* unit normal */
                           const float t0[3],   /* tangent axis 0 (longitudinal) */
                           const float t1[3],   /* tangent axis 1 (lateral) */
                           float mu_long,       /* friction along t0 */
                           float mu_lat,        /* friction along t1 */
                           float rr_scale,      /* rolling resistance scale (>=0) */
                           float phi,           /* signed distance; <0 penetrates */
                           float kappa_n,       /* compliance */
                           float v_out[3],      /* projected velocity */
                           int* out_saturated);
#ifdef __cplusplus
}
#endif

#endif /* ZX_CONTACT_REF_H */
