/*!
 * \file apic_basis_test.cpp
 * \brief CPU-side verification of quadratic B-spline weights: partition of unity.
 */

#include <cassert>
#include <cmath>

static void cpu_wx3(float x, float w[3]) {
    float x0 = 0.5f - x;
    float x1 = 0.5f + x;
    w[0] = 0.5f * x0 * x0;
    w[1] = 0.75f - x * x;
    w[2] = 0.5f * x1 * x1;
}

int main() {
    for (int i = -50; i <= 50; ++i) {
        float x = i / 50.0f * 0.5f; // sample in [-0.5, 0.5]
        float w[3]; cpu_wx3(x, w);
        float s = w[0] + w[1] + w[2];
        assert(std::fabs(s - 1.0f) < 1e-5f);
    }
    return 0;
}


