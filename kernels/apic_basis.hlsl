//!
// \file apic_basis.hlsl
// \brief Quadratic B-spline weights and gradients for APIC/MLS transfers.
// \author Colin Macritchie (Ripple Group, LLC)
//

static float3 wx3(float x) {
    float x0 = 0.5 - x;
    float x1 = 0.5 + x;
    float w0 = 0.5 * x0 * x0;
    float w1 = 0.75 - x * x;
    float w2 = 0.5 * x1 * x1;
    return float3(w0, w1, w2);
}

static float3 dwx3(float x) {
    float x0 = 0.5 - x;
    float x1 = 0.5 + x;
    float g0 = -x0;      // d(0.5*x0^2)/dx = -x0
    float g1 = -2.0*x;   // d(0.75 - x^2)/dx = -2x
    float g2 =  x1;      // d(0.5*x1^2)/dx = x1
    return float3(g0, g1, g2);
}


