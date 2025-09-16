/*!
 * \file presets_test.cpp
 * \brief Verify preset list and parameter retrieval.
 */

#include "zx/zx_presets.h"
#include <cassert>

int main(){
    int n = zx_preset_count();
    assert(n > 0);
    const char* name0 = zx_preset_name(0);
    assert(name0 && name0[0] != '\0');
    zx_elastic_params ep{}; zx_mc_params mc{}; zx_norsand_params ns{}; zx_norsand_state st{};
    int ok = zx_preset_get(name0, &ep, &mc, &ns, &st);
    assert(ok == 1);
    assert(ep.young_E > 0 && mc.friction_deg > 0);
    return 0;
}


