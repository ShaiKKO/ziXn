# ziXn — Architecture & ABI
**Granular, Slurry, and Snow Terrain Solver for Real‑Time Engines**  
_Target platforms: PC/Console (DX12/Vulkan/Metal) and scalable CPU fallback._

> This document defines the C‑ABI, object model, memory layout, execution graph, backend abstractions, and engine interop surfaces for **ziXn**. The goal is **deterministic‑enough**, **SIMT‑friendly**, and **engine‑agnostic** simulation of granular and slurry terrain at 60–120 Hz.

---

## Table of Contents
1. [Design Tenets](#design-tenets)
2. [External Requirements & Budgets](#external-requirements--budgets)
3. [C‑ABI Overview](#c-abi-overview)
   - [Versioning & Proc Table](#versioning--proc-table)
   - [Opaque Handles & Object Lifetime](#opaque-handles--object-lifetime)
   - [Error Model](#error-model)
   - [Threading & Reentrancy](#threading--reentrancy)
4. [Object Model](#object-model)
   - [Context, Device, Scene, Material, Rig, TileSet](#context-device-scene-material-rig-tileset)
   - [Resource Lifetimes & Pools](#resource-lifetimes--pools)
   - [Scene Graph & Instances](#scene-graph--instances)
5. [Data Layout & Memory](#data-layout--memory)
   - [Particle SoA, Grid Tiles, and Sparse Hash](#particle-soa-grid-tiles-and-sparse-hash)
   - [Resident Sets, Heaps, and Staging](#resident-sets-heaps-and-staging)
   - [Alignment & Packing](#alignment--packing)
   - [Quantization & Mixed Precision](#quantization--mixed-precision)
6. [Execution Graph](#execution-graph)
   - [Frame Phases](#frame-phases)
   - [Task Graph Nodes](#task-graph-nodes)
   - [Synchronization Strategy](#synchronization-strategy)
7. [GPU Backends](#gpu-backends)
   - [Direct3D 12](#direct3d-12)
   - [Vulkan](#vulkan)
   - [Metal](#metal)
   - [Shader Source Strategy](#shader-source-strategy)
8. [CPU Backend](#cpu-backend)
9. [Interop Surfaces](#interop-surfaces)
   - [Unreal (RDG/RHI) Surface](#unreal-rdgrhi-surface)
   - [Unity (Native Plugin) Surface](#unity-native-plugin-surface)
   - [Shared Terrain Writeback](#shared-terrain-writeback)
10. [Determinism Strategy](#determinism-strategy)
11. [Instrumentation & Telemetry](#instrumentation--telemetry)
12. [Build, Packaging, & Platform Notes](#build-packaging--platform-notes)
13. [Security & Robustness](#security--robustness)
14. [Do’s & Don’ts](#dos--donts)
15. [Appendix: ABI Types](#appendix-abi-types)
16. [Appendix: Execution Timings](#appendix-execution-timings)

---

## Design Tenets
- **C‑ABI first**: single dynamic library per platform; zero STL in API; stable POD structs.
- **Opaque handles**: outside world never sees internal pointers or layouts.
- **Tile‑sparse**: compute only where matter exists; keep active set bounded.
- **SIMT‑friendly**: minimize divergence; prefer structure‑of‑arrays; use shared memory reductions.
- **Deterministic‑enough**: stable reductions, quantized I/O, reproducible seeds.
- **Engine‑agnostic**: interop is narrow—buffers in/out, command callback hooks, no hard engine deps.
- **Fail soft**: graceful degradation (presentation‑only writeback) when budgets are exceeded.
- **Observability**: built‑in counters, GPU/CPU timestamps, checkers for invariants.

## External Requirements & Budgets
- **Targets**: 60–120 Hz; **simulation budget** 1.0–3.5 ms on high‑end, 3–6 ms on mid‑tier, 6–8 ms on consoles’ heavy scenes.
- **Memory**: ≤ 512 MB GPU budget for ziXn on high‑end scenes; ≤ 256 MB on mid; ≤ 192 MB on last‑gen.
- **Scalability knobs**: particle density per tile, grid h, solver precision, model LOD (dry granular ⇄ slurry ⇄ two‑phase).

---

## C‑ABI Overview

### Versioning & Proc Table
- Single export: `zxGetProcTable(uint32_t abi_version, zx_procs* out) -> zx_status`.
- **Semantic versioning** in the table; strict **size checks** on all input structs (`sizeof` gates).
- All functions are **C** (`extern "C"`) and **stdcall/cdecl** per platform (resolved at build).

```c
typedef struct zx_procs {
    uint32_t size;      /* must be set by caller to sizeof(zx_procs) */
    uint32_t version;   /* semantic ABI version, e.g., 0x00010002 */
    /* creation */
    zx_status (*create_context)(const zx_context_desc*, zx_context*);
    void      (*destroy_context)(zx_context);
    /* device/backends */
    zx_status (*bind_device)(zx_context, const zx_device_desc*, zx_device*);
    zx_status (*unbind_device)(zx_device);
    /* scene/material/rig lifecycle */
    zx_status (*create_scene)(zx_context, const zx_scene_desc*, zx_scene*);
    void      (*destroy_scene)(zx_scene);
    zx_status (*create_material)(zx_context, const zx_material_desc*, zx_material*);
    void      (*destroy_material)(zx_material);
    zx_status (*create_rig)(zx_context, const zx_rig_desc*, zx_rig*);
    void      (*destroy_rig)(zx_rig);
    /* frame */
    zx_status (*begin_frame)(zx_scene, const zx_frame_begin*);
    zx_status (*simulate)(zx_scene, const zx_sim_params*);
    zx_status (*writeback)(zx_scene, const zx_writeback_desc*);
    zx_status (*end_frame)(zx_scene, const zx_frame_end*);
    /* telemetry */
    zx_status (*get_counters)(zx_context, zx_counters*, uint32_t* count);
    /* errors/logging */
    const char* (*error_string)(zx_status);
} zx_procs;
```

### Opaque Handles & Object Lifetime
- Handles are **64‑bit non‑zero tokens**; validated by epoch and pool index.
- Creation requires **descriptor** structs; destruction is idempotent (safe to call twice).
- **Thread‑safe** for distinct handles; per‑scene functions require external synchronization unless marked `*_ts`.

### Error Model
- `zx_status` codes are negative for errors (`ZX_E_*`), zero for OK, positive for soft warnings.
- On error, functions **never partially write** output objects.

### Threading & Reentrancy
- ABI calls are **non‑blocking** unless documented. GPU work is enqueued; completion is observed via fences/timestamps returned by `end_frame` or engine fences passed in `zx_device_desc`.

---

## Object Model

### Context, Device, Scene, Material, Rig, TileSet
- **Context**: host allocators, logging sinks, feature bits, and global pools.
- **Device**: binding to one GPU backend or CPU fallback; owns heaps, queues, pipelines.
- **Scene**: simulation world (tile sparse grid, particle sets, constraints, emitters).
- **Material**: immutable physical parameters + authoring name and UI ranges.
- **Rig**: kinematic actors (wheels, feet, tracks) and colliders; can be instanced.
- **TileSet**: internal object (not exposed); manages active grid tiles and page tables.

### Resource Lifetimes & Pools
- All resources are **pooled** with **generation counters** to avoid ABA issues.
- GPU memory uses **sub‑allocation** from large heaps (buddy allocator + ring for transient).

### Scene Graph & Instances
- Scenes can hold many **instances** of rigs and material regions. Instances reference shared immutable data (materials, geometry) and per‑instance state (poses, masks, densities).

---

## Data Layout & Memory

### Particle SoA, Grid Tiles, and Sparse Hash
- **Particles (SoA)**: positions (float3), velocities (float3), mass (f32), volume (f32), deformation gradient **F** (9 floats), affine velocity **C** (9 floats, APIC), plastic state (varies), material id (u16), flags (u16).
- **Grid**: uniform spacing `h` per **tile**. Tiles are `BxBxB` nodes (e.g., 8³ or 16³) with **ghost** margins. Active tiles tracked in a **hash map** keyed by tile coords; each tile has compact arrays for node mass, momentum, and aux fields.
- **Neighbor lists**: per‑tile neighbor indices for stencil transfers (P2G/G2P).

### Resident Sets, Heaps, and Staging
- **GPU**: one **persistent heap** for long‑lived buffers (particles, tiles, materials), one or more **transient heaps** (P2G staging, sort keys, scans). Bindless descriptors for SRV/UAV tables.
- **CPU**: NUMA‑aware pools; large pages when available; **pinned** staging for fast upload/download.
- **Streaming**: tiles can be promoted/demoted each frame; hysteresis avoids thrash.

### Alignment & Packing
- All GPU buffers aligned to **256 bytes** minimum; per‑structure alignment is documented in [Appendix](#appendix-abi-types). Arrays are padded to warp/wave multiples to avoid divergence at tails.

### Quantization & Mixed Precision
- Support **F16** for many per‑particle scalars (mass, volume, plastic metrics) with **F32** accumulators on grid. Optionally quantize deformations (F) to 20‑bit pack where LOD permits. Determinism path uses full **F32** everywhere.

---

## Execution Graph

### Frame Phases
1. **Begin**: read engine inputs (rig transforms, wheel velocities, emitters), update material/scene params, allocate transient workspaces.
2. **Activation**: stamp particles and colliders into tile hash; (de)activate tiles based on support.
3. **P2G**: scatter particle mass/momentum/affine to grid nodes (tile‑local reduction, then global merge).
4. **Grid Solve**: (a) external & internal forces (constitutive update), (b) contact & constraints, (c) optional pressure solve for slurries/two‑phase.
5. **G2P**: gather nodal velocities back to particles (APIC/MLS), update positions and history.
6. **Writeback**: export height/displacement/masks to engine terrain targets; optional decals/FX maps.
7. **End**: publish counters, timestamps, and fences; compact inactive tiles if under pressure.

### Task Graph Nodes
- Nodes are coarse **compute passes** with explicit resources: e.g., `p2g_scatter`, `tile_merge`, `constitutive`, `contact_project`, `pressure_pcg`, `g2p_update`, `writeback_encode`.
- Dependencies are tracked as **barriers** (UAV/aliasing) and **timeline semaphores** across queues.

### Synchronization Strategy
- Prefer **asynchronous compute** for most passes; graphics only for visualization/debug.
- **Tile‑local shared memory** reductions eliminate global atomics during P2G/G2P where possible; fall back to atomics on tile merge.
- Cross‑queue sync is minimized: `pressure_pcg` may run concurrently with **presentation LOD generation**, merging before `g2p_update`.

---

## GPU Backends

### Direct3D 12
- **Descriptor heaps**: global SRV/UAV/CBV heap; per‑frame dynamic table for transient resources.
- **Root signatures**: minimal—one SRV table, one UAV table, constants via root CBV for tiny kernels.
- **Heaps**: large **DEFAULT** heaps, **UPLOAD** rings for staging.
- **Timestamp queries** for pass‑level timings; PIX markers for captures.

### Vulkan
- **Descriptor sets**: one persistent set for bindless SRV/UAV and one transient set per frame.
- **Memory**: VMA or custom allocator partitions; prefer large dedicated allocations for tile pools.
- **Queues**: compute + transfer; use **timeline semaphores** for coarse sync.
- **Spir‑V** generated via DXC; specialization constants choose B (tile size), kernel variants.

### Metal
- **Argument buffers** emulate bindless; use **MTLHeaps** for residency control.
- **Command queues**: separate compute for simulation; **blit** for compaction and copies.
- Shader translation via **MSL** backend; care for threadgroup memory limits on older iOS‑class HW (if used).

### Shader Source Strategy
- Single source of truth in **HLSL** templates; compiled to DXIL/SPIR‑V/MSL.  
- Kernel **variants** generated via macros: dimensionality (2D/3D), model (dry/slurry/two‑phase), precision (F32/F16), and debug (checksums).

---

## CPU Backend
- Parallel for and tasks via **OpenMP tasks** (or TBB fallback); NUMA pinning and memory interleave.
- Vectorization: **AVX2/AVX‑512/SVE** kernels for P2G/G2P and constitutive loops; gather/scatter optimized with SOA‑friendly layouts and precomputed neighbor indices.
- CPU backend is **feature‑equivalent** but designed for low particle densities or headless servers.

---

## Interop Surfaces

### Unreal (RDG/RHI) Surface
- The core ABI does **not** depend on UE. A thin UE module:
  - Injects an **RDG pass** that calls `simulate` with RHI buffer aliases.
  - Exposes **Niagara data interfaces** to sample displacement/masks (separate doc covers details).
  - Provides **RVT/VHM writeback** helpers (format conversion, tiling, clipping).

### Unity (Native Plugin) Surface
- Unity native plugin:
  - Manages `GraphicsBuffer/ComputeBuffer` objects bound to ziXn buffers.
  - C# interop wraps handles (safe `IntPtr`) and passes fences via `GraphicsFence` when available.
  - TerrainData writeback helpers (heightmap, holes, detail masks) with region tiling.

### Shared Terrain Writeback
- A minimal cross‑engine spec for writeback targets:
  - **Displacement** (R32F), **Compression/Moisture** (R16F), **Contact Mask** (R8U), **Debug** (R8U flags).
  - Tiled, clip‑safe, and **mip‑consistent** so engines can stream LODs without cracks.

---

## Determinism Strategy
- **Math**: clamp to F32; disable fused‑multiply‑add variance; optional **fixed‑point** path for reductions.
- **Reductions**: perform **tree‑based** tile reductions with deterministic order.
- **Randomness**: seeded per‑scene/particle via **counter‑based RNG** (Philox‑like) used only for emission/jitter.
- **I/O Quantization**: exporter quantizes to grid; consistent on all backends.
- **Cross‑platform drift tests**: checksum particle/world state post‑frame.

---

## Instrumentation & Telemetry
- **Counters**: particle count, active tiles, atomics, compactions, PCG iters, dt clamps, substeps.
- **Timings**: start/end GPU timestamps per pass; CPU wall times.
- **Invariants**: mass/momentum conservation, contact gap, bounds checks; failure escalates with context tag and frame index.
- **Capture hooks**: emit PIX/RenderDoc/Nsight markers; optional auto‑capture on regression.

---

## Build, Packaging, & Platform Notes
- Toolchain: MSVC/Clang; **DXC** for HLSL→DXIL/SPIR‑V; **SPIRV‑Cross** or custom MSL backend.
- Continuous integration builds **fat binaries** per platform with only the narrow C‑ABI exposed.
- **Static analyzers**: AddressSanitizer/UBSan on CPU path; GPU validation toggles for debug.
- Packaging: shared library + headers + JSON capabilities manifest (features, limits, shader variants).

---

## Security & Robustness
- Validate all counts/strides; refuse oversized inputs.
- Guard against **NaNs/Inf** on inputs; sanitize before compute.
- Avoid unbounded growth of tiles—enforce **hard caps** with graceful degrade.
- Harden destructors and out‑of‑order destruction via reference counts & generation checks.

---

## Do’s & Don’ts
**Do**
- Keep per‑frame allocations off the CPU heap; use rings/arenas.
- Batch small scenes to keep GPU busy; let scheduler over‑subscribe waves.
- Use presentation LOD writeback when solver substeps exceed the frame budget.

**Don’t**
- Expose raw pointers in the ABI.
- Depend on engine internals in the core library.
- Rely on unordered global atomics for correctness.

---

## Appendix: ABI Types
> Sizes and alignments are fixed across platforms for the ABI boundary.

```c
typedef uint64_t zx_handle;
typedef zx_handle zx_context;
typedef zx_handle zx_device;
typedef zx_handle zx_scene;
typedef zx_handle zx_material;
typedef zx_handle zx_rig;

typedef enum zx_status {
    ZX_OK = 0,
    ZX_E_INVALID = -1,
    ZX_E_UNSUPPORTED = -2,
    ZX_E_OOM = -3,
    ZX_E_DEVICE = -4,
    ZX_E_TIMEOUT = -5,
    ZX_W_SOFT = 1
} zx_status;

typedef struct zx_vec3 { float x,y,z; } zx_vec3;
typedef struct zx_mat3 { float m[9]; } zx_mat3;

typedef struct zx_context_desc {
    uint32_t size;
    void*    (*host_alloc)(size_t, void*);
    void     (*host_free)(void*, void*);
    void*    user;
    uint32_t feature_bits; /* e.g., ZX_FEATURE_F16, ZX_FEATURE_DET */
} zx_context_desc;

/* ... (trimmed for brevity within doc) ... */
```

### Binary Compatibility Notes
- All public structs start with `uint32_t size` for forward/back compat.
- New fields are appended; old clients still pass smaller `size` and get defaults.

---

## Appendix: Execution Timings
A representative 120 Hz profile (PC high‑end, ~2.0 M particles, 16³ tiles):
- Activation: 0.10 ms  
- P2G scatter: 0.55 ms  
- Tile merge: 0.25 ms  
- Constitutive: 0.60 ms  
- Contact project: 0.35 ms  
- Pressure PCG (slurry): 0.80 ms (8–12 iters)  
- G2P update: 0.55 ms  
- Writeback encode: 0.25 ms  
- **Total**: ~3.45 ms  
Use these numbers as **targets**, not guarantees; scale via particle density and tile LOD.
