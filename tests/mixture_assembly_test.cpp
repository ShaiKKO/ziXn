/*!
 * \file mixture_assembly_test.cpp
 * \brief Tests for two-phase mixture helpers.
 */

#include "zx/zx_mixture.h"
#include <cassert>
#include <cmath>
#include <algorithm>

static bool near(float a, float b, float eps=1e-5f){ return std::fabs(a-b) <= eps * (1.0f + std::max(std::fabs(a), std::fabs(b))); }

int main(){
    // Darcy beta
    {
        float beta = zx_beta_darcy(1.0f, 1e-10f, 1e-12f);
        assert(beta > 0.0f);
        float beta2 = zx_beta_darcy(1.0f, 0.0f, 1e-9f);
        assert(near(beta2, 1.0e9f));
    }

    // Porosity update clamps into [phi_min,phi_max]
    {
        float phi = 0.4f;
        float phi_next = zx_update_porosity(phi, +5.0f, 0.1f, 0.05f, 0.6f);
        assert(phi_next <= 0.6f);
        phi_next = zx_update_porosity(phi, -5.0f, 0.1f, 0.05f, 0.6f);
        assert(phi_next >= 0.05f);
    }

    // Effective pressure
    {
        float p_eff = zx_effective_pressure(100.0f, 30.0f, 0.8f);
        assert(near(p_eff, 100.0f - 0.8f*30.0f));
    }

    // Drag symmetry
    {
        float Fs_x=0,Fs_y=0,Fs_z=0,Ff_x=0,Ff_y=0,Ff_z=0;
        zx_mixture_drag_force(5.0f, /*vf*/1,0,0, /*vs*/0,0,0, /*vol*/2.0f, /*dt*/0.5f,
                              &Fs_x,&Fs_y,&Fs_z,&Ff_x,&Ff_y,&Ff_z);
        assert(near(Fs_x, -Ff_x) && near(Fs_y, -Ff_y) && near(Fs_z, -Ff_z));
        assert(Fs_x > 0.0f); // fluid faster than solid -> positive force on solid
    }
    return 0;
}


