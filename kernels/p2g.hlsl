//!
// \file p2g.hlsl
// \brief APIC/MLS P2G scatter: particle mass/momentum to grid.
// \author Colin Macritchie (Ripple Group, LLC)
// \credits See docs: ALGORITHMS.md (P2G/G2P Transfers, APIC/MLS) and ARCHITECTURE.md (Execution
// Graph)

#include "apic_basis.hlsl"

struct Particles
{
  StructuredBuffer<float> pos_x;
  StructuredBuffer<float> pos_y;
  StructuredBuffer<float> pos_z;
  StructuredBuffer<float> vel_x;
  StructuredBuffer<float> vel_y;
  StructuredBuffer<float> vel_z;
  StructuredBuffer<float> mass;
  StructuredBuffer<float> C;  // 9*N, row-major
};

RWStructuredBuffer<float> g_mass;  // nx*ny*nz
RWStructuredBuffer<float> g_mom;   // 3*nx*ny*nz

cbuffer Params : register(b0)
{
  uint num_particles;
  uint nx;
  uint ny;
  uint nz;
  float3 origin;
  float h;
}

[numthreads(256, 1, 1)] void cs_p2g(uint3 tid : SV_DispatchThreadID)
{
  uint i = tid.x;
  if (i >= num_particles)
  {
    return;
  }

  float3 xp = float3(Particles.pos_x[i], Particles.pos_y[i], Particles.pos_z[i]);
  float3 vp = float3(Particles.vel_x[i], Particles.vel_y[i], Particles.vel_z[i]);
  float mp  = Particles.mass[i];

  float3 gp = (xp - origin) / h;
  int3 base = (int3) floor(gp - 0.5);
  float3 wx = wx3(gp.x - (base.x + 1));
  float3 wy = wx3(gp.y - (base.y + 1));
  float3 wz = wx3(gp.z - (base.z + 1));

  float3x3 Cp = 0;
  if (Particles.C.GetDimensions() > 0)
  {
    uint ci  = i * 9;
    Cp[0][0] = Particles.C[ci + 0];
    Cp[0][1] = Particles.C[ci + 1];
    Cp[0][2] = Particles.C[ci + 2];
    Cp[1][0] = Particles.C[ci + 3];
    Cp[1][1] = Particles.C[ci + 4];
    Cp[1][2] = Particles.C[ci + 5];
    Cp[2][0] = Particles.C[ci + 6];
    Cp[2][1] = Particles.C[ci + 7];
    Cp[2][2] = Particles.C[ci + 8];
  }

  [unroll] for (int dz = 0; dz < 3; ++dz)
  {
    [unroll] for (int dy = 0; dy < 3; ++dy)
    {
      [unroll] for (int dx = 0; dx < 3; ++dx)
      {
        int3 gi = base + int3(dx, dy, dz);
        if (gi.x < 0 || gi.y < 0 || gi.z < 0 || gi.x >= (int) nx || gi.y >= (int) ny ||
            gi.z >= (int) nz)
        {
          continue;
        }
        uint idx     = (gi.z * ny + gi.y) * nx + gi.x;
        float w      = wx[dx] * wy[dy] * wz[dz];
        float3 xrel  = (float3(gi) - gp) * h;
        float3 v_aff = vp + mul(Cp, xrel);
        InterlockedAdd(g_mass[idx], w * mp);
        InterlockedAdd(g_mom[3 * idx + 0], w * mp * v_aff.x);
        InterlockedAdd(g_mom[3 * idx + 1], w * mp * v_aff.y);
        InterlockedAdd(g_mom[3 * idx + 2], w * mp * v_aff.z);
      }
    }
  }
}
