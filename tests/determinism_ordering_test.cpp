/*!
 * \file determinism_ordering_test.cpp
 * \brief Show nondeterministic vs deterministic P2G accumulation behavior.
 */

#include "zx/zx_apic_ref.h"
#include "zx/zx_determinism.h"
#include <vector>
#include <cassert>
#include <cstring>

static void run_once(int det, std::vector<float>& m, std::vector<float>& p){
    zx_set_determinism(det);
    const int nx=8, ny=8, nz=8;
    std::vector<float> pos, vel, C, mass; pos.resize(3*2); vel.resize(3*2); mass.resize(2);
    // Two particles that map to overlapping 27-node stencils
    pos[0]=1.2f; pos[1]=1.2f; pos[2]=1.2f; vel[0]=1.0f; vel[1]=0.0f; vel[2]=0.0f; mass[0]=1.0f;
    pos[3]=1.25f; pos[4]=1.15f; pos[5]=1.22f; vel[3]=0.8f; vel[4]=0.1f; vel[5]=0.0f; mass[1]=1.0f;
    float origin[3]={0,0,0}; float h=1.0f;
    m.assign(nx*ny*nz, 0.0f); p.assign(3*nx*ny*nz, 0.0f);
    zx_apic_p2g_ref(2, pos.data(), vel.data(), nullptr, mass.data(), origin, h, nx,ny,nz, m.data(), p.data());
}

int main(){
    std::vector<float> m0, p0, m1, p1;
    run_once(0, m0, p0); // nondeterministic
    run_once(1, m1, p1); // deterministic
    // Deterministic should match itself when run twice
    std::vector<float> m2, p2; run_once(1, m2, p2);
    assert(m1 == m2);
    assert(p1 == p2);
    return 0;
}


