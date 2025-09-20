//!
// \file contact_project.hlsl
// \brief Grid-node contact projection: non-penetration + Coulomb friction with compliance.
// \credits See ALGORITHMS.md Contact section for method.

StructuredBuffer<float> g_mass;   // nx*ny*nz
RWStructuredBuffer<float> g_vel;  // 3*nx*ny*nz
StructuredBuffer<float> g_norm;   // 3*nx*ny*nz (surface normal per node)
StructuredBuffer<float> g_phi;    // nx*ny*nz (signed distance)

cbuffer Params : register(b0)
{
  uint nx;
  uint ny;
  uint nz;
  float mu;
  float kappa_n;  // friction and compliance
}

[numthreads(256, 1, 1)] void cs_contact(uint3 tid : SV_DispatchThreadID)
{
  uint idx   = tid.x;
  uint total = nx * ny * nz;
  if (idx >= total)
  {
    return;
  }
  float m = g_mass[idx];
  if (m <= 0.0)
  {
    return;
  }
  float3 v  = float3(g_vel[3 * idx + 0], g_vel[3 * idx + 1], g_vel[3 * idx + 2]);
  float3 n  = float3(g_norm[3 * idx + 0], g_norm[3 * idx + 1], g_norm[3 * idx + 2]);
  float phi = g_phi[idx];

  float vn  = dot(v, n);
  float3 vt = v - vn * n;
  if (phi < 0.0)
  {
    float allowance = -phi * kappa_n;
    if (vn < -allowance)
      vn = -allowance;
  }
  float vt_len = length(vt);
  if (vt_len > 1e-12)
  {
    float max_t = mu * max(0.0, vn);
    if (vt_len > max_t)
      vt *= (max_t / vt_len);
  }
  v                  = vn * n + vt;
  g_vel[3 * idx + 0] = v.x;
  g_vel[3 * idx + 1] = v.y;
  g_vel[3 * idx + 2] = v.z;
}
