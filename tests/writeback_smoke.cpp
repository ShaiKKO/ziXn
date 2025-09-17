/*!
 * \file writeback_smoke.cpp
 * \brief Validate displacement writeback clipping.
 */

#include "zx/zx_writeback.h"
#include <vector>
#include <cassert>

int main(){
    const uint32_t W=8,H=8; std::vector<float> disp(W*H, 0.0f);
    zx_write_target t{ disp.data(), W, H, (uint32_t)(W*sizeof(float)), nullptr, nullptr, nullptr };
    std::vector<float> src(4*4, 1.0f);
    zx_rect rect{ 6,6,4,4 };
    zx_writeback_copy_displacement(src.data(), 4,4, &t, rect);
    // Expect corners within [6,7] to be 1.0
    assert(disp[6 + 6*W] == 1.0f);
    assert(disp[7 + 7*W] == 1.0f);
    return 0;
}


