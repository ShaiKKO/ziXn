/*!
 * \file zx_apic_ref.cpp
 * \brief CPU reference APIC/MLS transfers used for validation.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
 */

#include "zx/zx_apic_ref.h"
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include "zx/zx_abi.h"
#include "zx/zx_determinism.h"

static inline int idx3(int x, int y, int z, int nx, int ny) { return (z*ny + y)*nx + x; }

void zx_bspline_w3(float x, float w[3]) {
    const float x0 = 0.5f - x;
    const float x1 = 0.5f + x;
    w[0] = 0.5f * x0 * x0;
    w[1] = 0.75f - x * x;
    w[2] = 0.5f * x1 * x1;
}

void zx_bspline_dw3(float x, float g[3]) {
    const float x0 = 0.5f - x;
    const float x1 = 0.5f + x;
    g[0] = -x0;
    g[1] = -2.0f * x;
    g[2] =  x1;
}

void zx_apic_p2g_ref(
    size_t N,
    const float* pos,
    const float* vel,
    const float* C,
    const float* mass,
    const float origin[3],
    float h,
    int nx, int ny, int nz,
    float* ZX_RESTRICT m_grid,
    float* ZX_RESTRICT p_grid)
{
    std::memset(m_grid, 0, sizeof(float)*nx*ny*nz);
    std::memset(p_grid, 0, sizeof(float)*3*nx*ny*nz);
    const float invh = 1.0f / h;
    const int deterministic = zx_get_determinism();
    struct Contrib { int gi; int pid; float m; float mx; float my; float mz; };
    std::vector<Contrib> contribs;
    if (deterministic) {
        contribs.reserve(N * 27);
    }
    for (size_t i = 0; i < N; ++i) {
        const float px = pos[3*i+0];
        const float py = pos[3*i+1];
        const float pz = pos[3*i+2];
        const float vx = vel[3*i+0];
        const float vy = vel[3*i+1];
        const float vz = vel[3*i+2];
        const float* Cp = C ? &C[9*i] : nullptr;
        const float mp = mass[i];

        const float gx = (px - origin[0]) * invh;
        const float gy = (py - origin[1]) * invh;
        const float gz = (pz - origin[2]) * invh;

        const int base_x = (int)std::floor(gx - 0.5f);
        const int base_y = (int)std::floor(gy - 0.5f);
        const int base_z = (int)std::floor(gz - 0.5f);

        float wx[3], wy[3], wz[3];
        zx_bspline_w3(gx - (base_x + 1), wx);
        zx_bspline_w3(gy - (base_y + 1), wy);
        zx_bspline_w3(gz - (base_z + 1), wz);

        for (int dz = 0; dz < 3; ++dz) for (int dy = 0; dy < 3; ++dy) for (int dx = 0; dx < 3; ++dx) {
            const int ix = base_x + dx;
            const int iy = base_y + dy;
            const int iz = base_z + dz;
            if ((unsigned)ix >= (unsigned)nx || (unsigned)iy >= (unsigned)ny || (unsigned)iz >= (unsigned)nz) continue;
            const float w = wx[dx] * wy[dy] * wz[dz];
            const int gi = idx3(ix,iy,iz,nx,ny);
            if (!deterministic) {
                m_grid[gi] += w * mp;
            }

            float vx_aff = vx;
            float vy_aff = vy;
            float vz_aff = vz;
            if (Cp) {
                const float gx_rel = (ix + 0.0f - gx);
                const float gy_rel = (iy + 0.0f - gy);
                const float gz_rel = (iz + 0.0f - gz);
                // APIC affine: v + C*(x_i - x_p) scaled by h
                vx_aff += (Cp[0]*gx_rel + Cp[1]*gy_rel + Cp[2]*gz_rel) * h;
                vy_aff += (Cp[3]*gx_rel + Cp[4]*gy_rel + Cp[5]*gz_rel) * h;
                vz_aff += (Cp[6]*gx_rel + Cp[7]*gy_rel + Cp[8]*gz_rel) * h;
            }
            if (!deterministic) {
                p_grid[3*gi+0] += w * mp * vx_aff;
                p_grid[3*gi+1] += w * mp * vy_aff;
                p_grid[3*gi+2] += w * mp * vz_aff;
            } else {
                const float mm = w * mp;
                contribs.push_back(Contrib{gi, (int)i, mm, mm * vx_aff, mm * vy_aff, mm * vz_aff});
            }
        }
    }
    if (deterministic) {
        std::stable_sort(contribs.begin(), contribs.end(), [](const Contrib& a, const Contrib& b){
            if (a.gi != b.gi) return a.gi < b.gi;
            return a.pid < b.pid;
        });
        int cur = -1;
        float msum = 0.0f, mxsum = 0.0f, mysum = 0.0f, mzsum = 0.0f;
        auto flush = [&](){ if (cur >= 0) { m_grid[cur] += msum; p_grid[3*cur+0] += mxsum; p_grid[3*cur+1] += mysum; p_grid[3*cur+2] += mzsum; } };
        for (const auto& c : contribs) {
            if (c.gi != cur) { flush(); cur = c.gi; msum = mxsum = mysum = mzsum = 0.0f; }
            msum += c.m; mxsum += c.mx; mysum += c.my; mzsum += c.mz;
        }
        flush();
    }
}

void zx_apic_g2p_ref(
    size_t N,
    const float* pos,
    float*       out_vel,
    float*       out_C,
    const float origin[3],
    float h,
    int nx, int ny, int nz,
    const float* ZX_RESTRICT m_grid,
    const float* ZX_RESTRICT v_grid)
{
    const float invh = 1.0f / h;
    for (size_t i = 0; i < N; ++i) {
        const float px = pos[3*i+0];
        const float py = pos[3*i+1];
        const float pz = pos[3*i+2];

        const float gx = (px - origin[0]) * invh;
        const float gy = (py - origin[1]) * invh;
        const float gz = (pz - origin[2]) * invh;

        const int base_x = (int)std::floor(gx - 0.5f);
        const int base_y = (int)std::floor(gy - 0.5f);
        const int base_z = (int)std::floor(gz - 0.5f);

        float wx[3], wy[3], wz[3];
        zx_bspline_w3(gx - (base_x + 1), wx);
        zx_bspline_w3(gy - (base_y + 1), wy);
        zx_bspline_w3(gz - (base_z + 1), wz);

        float vx=0, vy=0, vz=0;
        float Cx[3]={0,0,0}, Cy[3]={0,0,0}, Cz[3]={0,0,0};

        for (int dz = 0; dz < 3; ++dz) for (int dy = 0; dy < 3; ++dy) for (int dx = 0; dx < 3; ++dx) {
            const int ix = base_x + dx;
            const int iy = base_y + dy;
            const int iz = base_z + dz;
            if ((unsigned)ix >= (unsigned)nx || (unsigned)iy >= (unsigned)ny || (unsigned)iz >= (unsigned)nz) continue;
            const int gi = idx3(ix,iy,iz,nx,ny);
            const float m = m_grid[gi];
            if (m <= 0.0f) continue;
            const float w = wx[dx] * wy[dy] * wz[dz];
            const float vx_i = v_grid[3*gi+0];
            const float vy_i = v_grid[3*gi+1];
            const float vz_i = v_grid[3*gi+2];
            vx += w * vx_i; vy += w * vy_i; vz += w * vz_i;

            if (out_C) {
                const float gx_rel = (ix + 0.0f - gx);
                const float gy_rel = (iy + 0.0f - gy);
                const float gz_rel = (iz + 0.0f - gz);
                Cx[0] += w * vx_i * gx_rel; Cx[1] += w * vx_i * gy_rel; Cx[2] += w * vx_i * gz_rel;
                Cy[0] += w * vy_i * gx_rel; Cy[1] += w * vy_i * gy_rel; Cy[2] += w * vy_i * gz_rel;
                Cz[0] += w * vz_i * gx_rel; Cz[1] += w * vz_i * gy_rel; Cz[2] += w * vz_i * gz_rel;
            }
        }
        out_vel[3*i+0] = vx;
        out_vel[3*i+1] = vy;
        out_vel[3*i+2] = vz;
        if (out_C) {
            const float s = 4.0f * invh * invh; // common APIC scaling
            out_C[9*i+0] = s * Cx[0]; out_C[9*i+1] = s * Cx[1]; out_C[9*i+2] = s * Cx[2];
            out_C[9*i+3] = s * Cy[0]; out_C[9*i+4] = s * Cy[1]; out_C[9*i+5] = s * Cy[2];
            out_C[9*i+6] = s * Cz[0]; out_C[9*i+7] = s * Cz[1]; out_C[9*i+8] = s * Cz[2];
        }
    }
}


