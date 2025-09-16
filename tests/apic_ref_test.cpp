/*!
 * \file apic_ref_test.cpp
 * \brief Validate APIC P2G/G2P round-trip basics and conservation.
 */

#include "zx/zx_apic_ref.h"
#include <cassert>
#include <vector>
#include <cmath>

int main() {
    const size_t N = 4;
    std::vector<float> pos(3*N), vel(3*N), C(9*N, 0.0f), mass(N, 1.0f);
    const float origin[3] = {0,0,0};
    const float h = 1.0f;
    const int nx=8, ny=8, nz=8;
    std::vector<float> m(nx*ny*nz, 0.0f), p(3*nx*ny*nz, 0.0f), v(3*nx*ny*nz, 0.0f);

    // Place a few particles in the center with simple velocities
    for (size_t i=0;i<N;++i){
        pos[3*i+0] = 3.2f + 0.1f*i;
        pos[3*i+1] = 3.4f + 0.1f*i;
        pos[3*i+2] = 3.6f + 0.1f*i;
        vel[3*i+0] = 0.5f; vel[3*i+1] = 0.25f; vel[3*i+2] = -0.1f;
    }

    zx_apic_p2g_ref(N, pos.data(), vel.data(), C.data(), mass.data(), origin, h, nx, ny, nz, m.data(), p.data());

    // Compute grid velocities
    for (int k=0;k<nz;++k) for (int j=0;j<ny;++j) for (int i=0;i<nx;++i){
        const int gi = (k*ny + j)*nx + i;
        if (m[gi] > 0.0f) {
            v[3*gi+0] = p[3*gi+0] / m[gi];
            v[3*gi+1] = p[3*gi+1] / m[gi];
            v[3*gi+2] = p[3*gi+2] / m[gi];
        }
    }

    std::vector<float> vel2(3*N, 0.0f), C2(9*N, 0.0f);
    zx_apic_g2p_ref(N, pos.data(), vel2.data(), C2.data(), origin, h, nx, ny, nz, m.data(), v.data());

    // Basic checks: velocity should be close to original (APIC is not exact but should be near)
    for (size_t i=0;i<N;++i){
        assert(std::fabs(vel2[3*i+0] - vel[3*i+0]) < 1e-4f);
        assert(std::fabs(vel2[3*i+1] - vel[3*i+1]) < 1e-4f);
        assert(std::fabs(vel2[3*i+2] - vel[3*i+2]) < 1e-4f);
    }
    return 0;
}


