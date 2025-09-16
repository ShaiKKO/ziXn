/*!
 * \file abi_smoke.cpp
 * \brief Smoke test for zxGetProcTable and basic size/version gates.
 * \author Colin Macritchie (Ripple Group, LLC)
 */

#include "zx/zx_abi.h"
#include <cassert>
#include <cstdio>

int main() {
    zx_procs P{};
    // Expect unsupported version to fail
    zx_status st = zxGetProcTable(0xDEADBEEFu, &P);
    assert(st == ZX_E_UNSUPPORTED);

    // Now query with expected version
    st = zxGetProcTable(ZX_ABI_VERSION, &P);
    assert(st == ZX_OK);
    assert(P.size == sizeof(zx_procs));
    assert(P.create_context && P.error_string);

    zx_context C{};
    zx_context_desc cd{}; cd.size = sizeof(cd);
    st = P.create_context(&cd, &C);
    assert(st == ZX_OK && C != 0);

    const char* msg = P.error_string(ZX_OK);
    assert(msg && msg[0] != '\0');

    std::puts("abi_smoke: OK");
    return 0;
}


