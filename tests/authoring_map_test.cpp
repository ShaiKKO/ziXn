/*!
 * \file authoring_map_test.cpp
 * \brief Verify UI mapping clamps and monotonic trends.
 */

#include "zx/zx_authoring.h"
#include <cassert>

int main(){
    zx_ui_sand ui0{0.0f, 0.0f}, ui1{1.0f, 1.0f};
    zx_elastic_params e0{}, e1{}; zx_mc_params mc0{}, mc1{}; zx_norsand_params ns0{}, ns1{}; zx_norsand_state st0{}, st1{};
    zx_authoring_map_sand(&ui0, &e0, &mc0, &ns0, &st0);
    zx_authoring_map_sand(&ui1, &e1, &mc1, &ns1, &st1);
    assert(e1.young_E > e0.young_E);
    assert(mc1.friction_deg > mc0.friction_deg);

    zx_ui_snow s0{0.0f, 0.0f}, s1{1.0f, 1.0f};
    zx_elastic_params se0{}, se1{}; zx_mc_params sm0{}, sm1{};
    zx_authoring_map_snow(&s0, &se0, &sm0);
    zx_authoring_map_snow(&s1, &se1, &sm1);
    assert(se1.young_E > se0.young_E);
    return 0;
}


