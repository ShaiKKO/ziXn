/*!
 * \file zx_cli.cpp
 * \brief Minimal CLI harness to query validation proxies (inclined plane, column collapse).
 */

#include "zx/zx_validation.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

static void print_usage() {
    std::printf("ziXn CLI\n");
    std::printf("Usage:\n");
    std::printf("  zx_cli --mode inclined --phi <deg>\n");
    std::printf("  zx_cli --mode collapse --phi <deg> --ar <H_over_L0>\n");
}

int main(int argc, char** argv) {
    const char* mode = nullptr;
    float phi_deg = 30.0f;
    float ar = 1.0f;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--mode") == 0 && i+1 < argc) { mode = argv[++i]; }
        else if (std::strcmp(argv[i], "--phi") == 0 && i+1 < argc) { phi_deg = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--ar") == 0 && i+1 < argc) { ar = std::atof(argv[++i]); }
        else if (std::strcmp(argv[i], "--help") == 0) { print_usage(); return 0; }
    }
    if (!mode) { print_usage(); return 1; }

    if (std::strcmp(mode, "inclined") == 0) {
        zx_mc_params mc{ phi_deg, 0.0f };
        float theta_c = zx_validation_inclined_plane_theta_c(&mc);
        std::printf("inclined_plane: phi=%.3f deg -> theta_c=%.3f deg\n", phi_deg, theta_c*180.0/3.14159265358979323846);
        return 0;
    } else if (std::strcmp(mode, "collapse") == 0) {
        float ratio = zx_validation_column_collapse_runout_ratio(phi_deg, ar);
        std::printf("column_collapse: phi=%.3f deg, ar=%.3f -> L/H=%.3f\n", phi_deg, ar, ratio);
        return 0;
    } else {
        print_usage();
        return 1;
    }
}


