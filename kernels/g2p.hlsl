//!
// \file g2p.hlsl
// \brief APIC/MLS G2P gather: grid velocities to particles, update C.
// \author Colin Macritchie (Ripple Group, LLC)
// \credits See docs: ALGORITHMS.md (P2G/G2P Transfers)

#include "apic_basis.hlsl"

struct ParticlesRW
{
  RWStructuredBuffer<float> vel_x;
  RWStructuredBuffer<float> vel_y;
  RWStructuredBuffer<float> vel_z;
  RWStructuredBuffer<float> C;  // 9*N
};

StructuredBuffer<float> g_mass;  // nx*ny*nz
StructuredBuffer<float> g_vel;   // 3*nx*ny*nz
StructuredBuffer<float> pos_x, pos_y, pos_z;

cbuffer Params : register(b0)
{
  uint num_particles;
  uint nx;
  uint ny;
  uint nz;
  float3 origin;
  float h;
}

[numthreads(256, 1, 1)] void cs_g2p(uint3 tid : SV_DispatchThreadID)
{
  uint i = tid.x;
  if (i >= num_particles)
    return;
  float3 xp = float3(pos_x[i], pos_y[i], pos_z[i]);
  float3 gp = (xp - origin) / h;
  int3 base = (int3) floor(gp - 0.5);
  float3 wx = wx3(gp.x - (base.x + 1));
  float3 wy = wx3(gp.y - (base.y + 1));
  float3 wz = wx3(gp.z - (base.z + 1));

  float3 v      = 0.0;
  float3x3 Cacc = 0.0;
  [unroll] for (int dz = 0; dz < 3; ++dz)[unroll] for (int dy = 0; dy < 3;
                                                       ++dy)[unroll] for (int dx = 0; dx < 3; ++dx)
  {
    int3 gi = base + int3(dx, dy, dz);
    if (gi.x < 0 || gi.y < 0 || gi.z < 0 || gi.x >= (int) nx || gi.y >= (int) ny ||
        gi.z >= (int) nz)
      continue;
    uint idx = (gi.z * ny + gi.y) * nx + gi.x;
    float m  = g_mass[idx];
    if (m <= 0.0)
      continue;
    float w   = wx[dx] * wy[dy] * wz[dz];
    float3 vi = float3(g_vel[3 * idx + 0], g_vel[3 * idx + 1], g_vel[3 * idx + 2]);
    v += w * vi;
    float3 xrel = (float3(gi) - gp);
    Cacc[0][0] += w * vi.x * xrel.x;
    Cacc[0][1] += w * vi.x * xrel.y;
    Cacc[0][2] += w * vi.x * xrel.z;
    Cacc[1][0] += w * vi.y * xrel.x;
    Cacc[1][1] += w * vi.y * xrel.y;
    Cacc[1][2] += w * vi.y * xrel.z;
    Cacc[2][0] += w * vi.z * xrel.x;
    Cacc[2][1] += w * vi.z * xrel.y;
    Cacc[2][2] += w * vi.z * xrel.z;
  }
  float s               = 4.0 / (h * h);
  ParticlesRW.vel_x[i]  = v.x;
  ParticlesRW.vel_y[i]  = v.y;
  ParticlesRW.vel_z[i]  = v.z;
  uint ci               = i * 9;
  ParticlesRW.C[ci + 0] = s * Cacc[0][0];
  ParticlesRW.C[ci + 1] = s * Cacc[0][1];
  ParticlesRW.C[ci + 2] = s * Cacc[0][2];
  ParticlesRW.C[ci + 3] = s * Cacc[1][0];
  ParticlesRW.C[ci + 4] = s * Cacc[1][1];
  ParticlesRW.C[ci + 5] = s * Cacc[1][2];
  ParticlesRW.C[ci + 6] = s * Cacc[2][0];
  ParticlesRW.C[ci + 7] = s * Cacc[2][1];
  ParticlesRW.C[ci + 8] = s * Cacc[2][2];
}
