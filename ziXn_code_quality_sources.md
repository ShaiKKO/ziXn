# ziXn Code Quality & Architecture — Curated Source Links (C/C++)

> A focused, file‑level reading list of exemplary code and docs you can **lift patterns from directly** for ziXn’s GPU/CPU simulation core, memory management, tasking, and engine integration.

---

## 1) GPU Compute & Real‑Time Simulation Patterns

- **DirectX 12 compute: N‑Body Gravity sample (dispatch, barriers, root signature)**  
  `D3D12nBodyGravity.cpp` — Microsoft DirectX‑Graphics‑Samples  
  https://github.com/microsoft/DirectX-Graphics-Samples/blob/master/Samples/D3D12nBodyGravity/src/D3D12nBodyGravity.cpp  
  `NBodyGravityCS.hlsl` (compute kernel) — Microsoft DirectX‑Graphics‑Samples  
  https://github.com/microsoft/DirectX-Graphics-Samples/blob/master/Samples/D3D12nBodyGravity/src/NBodyGravityCS.hlsl

- **Vulkan compute: particle system with descriptor & pipeline setup**  
  `computeparticles.cpp` — Sascha Willems Vulkan examples  
  https://github.com/SaschaWillems/Vulkan/blob/master/examples/computeparticles/computeparticles.cpp

- **Bindless‑style resource indexing (Vulkan): production allocator header**  
  `vk_mem_alloc.h` — AMD **Vulkan Memory Allocator (VMA)**  
  https://skia.googlesource.com/external/github.com/GPUOpen-LibrariesAndSDKs/VulkanMemoryAllocator/+/HEAD/include/vk_mem_alloc.h

- **D3D12 resource suballocation (buddy/ring allocators, residency): production header**  
  `src/D3D12MemAlloc.h` — AMD **D3D12 Memory Allocator (D3D12MA)**  
  https://github.com/GPUOpen-LibrariesAndSDKs/D3D12MemoryAllocator/blob/master/src/D3D12MemAlloc.h

---

## 2) Task Systems, CPU Parallelism & Lock‑Free Queues

- **Work‑stealing task scheduler (lightweight, production‑grade)**  
  `src/TaskScheduler.h` — **enkiTS**  
  https://github.com/dougbinks/enkiTS/blob/master/src/TaskScheduler.h

- **SPSC/MPMC queue patterns (SoA‑friendly producer/consumer hot path)**  
  `concurrentqueue.h` — **moodycamel** lock‑free MPMC queue (single header)  
  https://github.com/cameron314/concurrentqueue/blob/master/concurrentqueue.h

- **Cache‑efficient ECS registry (entity iteration, sparse sets, data‑oriented design)**  
  `src/entt/entity/registry.hpp` — **EnTT**  
  https://github.com/skypjack/entt/blob/master/src/entt/entity/registry.hpp

---

## 3) High‑Performance Containers & Utilities

- **Flat Hash Map (SwissTable) — predictable probes & cache locality**  
  `absl/container/flat_hash_map.h` — **abseil‑cpp**  
  https://github.com/abseil/abseil-cpp/blob/master/absl/container/flat_hash_map.h

- **Robin Hood hashing (excellent for tile/particle spatial maps)**  
  `src/include/robin_hood.h` — **martinus/robin-hood-hashing**  
  https://github.com/martinus/robin-hood-hashing/blob/master/src/include/robin_hood.h

- **Fast logging (header & sinks) for frame‑time instrumentation**  
  `include/spdlog/spdlog.h` — **spdlog**  
  https://github.com/gabime/spdlog/blob/v1.x/include/spdlog/spdlog.h

- **Frame profiler you can embed (zones, GPU/CPU timings, locks)**  
  `public/tracy/Tracy.hpp` — **Tracy**  
  https://github.com/wolfpld/tracy/blob/master/public/tracy/Tracy.hpp

---

## 4) Engine Integration (Unity / Unreal) — Plugin Patterns

- **Unity native plugin: graphics interop, fences, render event hook**  
  `PluginSource/source/RenderingPlugin.cpp` — Unity **NativeRenderingPlugin**  
  https://github.com/Unity-Technologies/NativeRenderingPlugin/blob/master/PluginSource/source/RenderingPlugin.cpp

- **Unreal coding standard & engine‑side patterns (module style, headers, naming)**  
  Epic’s official **C++ coding standard** (useful when writing a UE plugin around ziXn)  
  https://dev.epicgames.com/documentation/en-us/unreal-engine/epic-cplusplus-coding-standard-for-unreal-engine

---

## 5) Core Guidelines & Code Health (to bake into CI)

- **C++ Core Guidelines (Stroustrup/Sutter)** — Adopt via clang‑tidy `cppcoreguidelines-*` checks  
  Main index: https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines  
  Repo: https://github.com/isocpp/CppCoreGuidelines

- **Google‑style flat hash map, abseil style/utility headers (as used by many engines)**  
  See `absl/container/flat_hash_map.h` (above) and companion docs/readme in abseil‑cpp.

---

## 6) Nice‑to‑Have Testing & Benchmarking

- **Google Benchmark (microbench harness for kernel hot‑loops)**  
  `include/benchmark/benchmark.h` — **google/benchmark**  
  https://github.com/google/benchmark/blob/main/include/benchmark/benchmark.h

- **Catch2 single‑header test framework**  
  `single_include/catch2/catch.hpp` — **catchorg/Catch2**  
  https://github.com/catchorg/Catch2/blob/devel/single_include/catch2/catch.hpp

---

### Notes on Fit for ziXn

- **GPU compute samples (DX12/Vulkan)** above mirror ziXn’s dispatch model (bindless descriptors, UAV barriers, timestamp queries).  
- **Allocators (VMA/D3D12MA)** match ziXn’s tile‑pool + transient scratch needs.  
- **Tasking (enkiTS) + lock‑free queues** map to CPU fallback & background staging.  
- **ECS/SoA (EnTT)** helps with particle/tile ownership, iteration order, and cache‑friendly updates.  
- **Containers (abseil/robin_hood)** are ideal for the active‑tile hash and per‑tile registries.  
- **Unity/Unreal links** show the exact places to hang GPU interop & fences for writeback.  
- **Guidelines/tests/benchmarks** help keep the codebase maintainable while optimizing.

