/*!
 * \file abi_conformance_test.cpp
 * \brief ABI conformance tests: size/version gates, nulls, zero-on-failure, idempotent destroys.
 * \author Colin Macritchie (Ripple Group, LLC)
 * \license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
 */

#include "zx/zx_abi.h"
#include <cassert>
#include <cstring>
#include <thread>
#include <vector>

static void assert_zeroed(const void* p, size_t n){
    const unsigned char* b = static_cast<const unsigned char*>(p);
    for (size_t i=0;i<n;++i) assert(b[i] == 0u);
}

int main(){
    // Case 1: invalid version must return ZX_E_UNSUPPORTED and zero out_procs if non-null
    zx_procs P1; std::memset(&P1, 0xAB, sizeof(P1));
    zx_status st = zxGetProcTable(0xDEADBEEFu, &P1);
    assert(st == ZX_E_UNSUPPORTED);
    assert_zeroed(&P1, sizeof(P1));

    // Case 2: null out_procs must return ZX_E_INVALID
    st = zxGetProcTable(ZX_ABI_VERSION, nullptr);
    assert(st == ZX_E_INVALID);

    // Case 3: valid path populates proc table and non-null function pointers
    zx_procs P{};
    st = zxGetProcTable(ZX_ABI_VERSION, &P);
    assert(st == ZX_OK);
    assert(P.size == sizeof(zx_procs));
    assert(P.create_context && P.destroy_context);
    assert(P.error_string);

    // Case 4: double-destroy safety for representative handles
    zx_context ctx{};
    zx_context_desc cd{}; cd.size = sizeof(cd);
    st = P.create_context(&cd, &ctx);
    assert(st == ZX_OK && ctx != 0);
    P.destroy_context(ctx);
    P.destroy_context(ctx); // idempotent no-op

    // For device unbind: unbind returns zx_status; call twice on a dummy handle
    zx_device_desc dd{}; dd.size = sizeof(dd);
    zx_device dev{};
    st = P.bind_device(ctx, &dd, &dev);
    assert(st == ZX_OK && dev != 0);
    st = P.unbind_device(dev);
    assert(st == ZX_OK);
    st = P.unbind_device(dev);
    assert(st == ZX_OK);

    // Case 5: concurrent zxGetProcTable calls return stable identical pointers
    zx_procs Pref{}; assert(zxGetProcTable(ZX_ABI_VERSION, &Pref) == ZX_OK);
    auto same = true;
    auto worker = [&](){
        zx_procs Px{}; zx_status s = zxGetProcTable(ZX_ABI_VERSION, &Px);
        if (s != ZX_OK) { same = false; return; }
        if (Px.create_context != Pref.create_context || Px.error_string != Pref.error_string) same = false;
    };
    std::vector<std::thread> threads;
    for (int i=0;i<8;++i) threads.emplace_back(worker);
    for (auto& t: threads) t.join();
    assert(same);

    return 0;
}


