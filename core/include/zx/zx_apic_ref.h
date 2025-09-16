/*!
\file zx_apic_ref.h
\brief CPU reference APIC/MLS P2G and G2P for validation.
\author Colin Macritchie (Ripple Group, LLC)
*/

#ifndef ZX_APIC_REF_H
#define ZX_APIC_REF_H

#include <stddef.h>
#include "zx_abi.h"

/* Quadratic B-spline basis (1D) */
ZX_API void ZX_CALL zx_bspline_w3(float x, float w[3]);
ZX_API void ZX_CALL zx_bspline_dw3(float x, float g[3]);

/* CPU reference transfers (3D grid) */
/* Grid is row-major [nz][ny][nx]; origin is min corner; h is spacing. */
ZX_API void ZX_CALL zx_apic_p2g_ref(
    size_t num_particles,
    const float* pos,         /* xyz triplets, length 3N */
    const float* vel,         /* xyz triplets, length 3N */
    const float* C,           /* 9 floats per particle (row-major), length 9N; can be null => zeros */
    const float* mass,        /* length N */
    const float origin[3],
    float h,
    int nx, int ny, int nz,
    float* ZX_RESTRICT m_grid, /* length nx*ny*nz */
    float* ZX_RESTRICT p_grid  /* xyz triplets per node, length 3*nx*ny*nz */
);

ZX_API void ZX_CALL zx_apic_g2p_ref(
    size_t num_particles,
    const float* pos,         /* 3N */
    float*       out_vel,     /* 3N */
    float*       out_C,       /* 9N, may be null to skip */
    const float origin[3],
    float h,
    int nx, int ny, int nz,
    const float* ZX_RESTRICT m_grid, /* nx*ny*nz */
    const float* ZX_RESTRICT v_grid  /* 3*nx*ny*nz (velocity per node) */
);

#endif /* ZX_APIC_REF_H */


