# ziXn — Forge the Ground
_A high‑performance C‑ABI terrain continuum library for **granular**, **snow**, and **slurry** materials, built for real‑time engines (Unreal, Unity, custom)._

---

## 1) What ziXn is
**ziXn** is a portable simulation SDK that makes **deformable terrain** a first‑class runtime system: footprints in powder, ruts in wet sand, wheel trenching and berms, bogging in mud, slush and landslide‑like flows, and layered beds (e.g., sand over clay). It provides:
- A **stable C ABI** with opaque handles; deterministic stepping; zero‑copy interop.
- **GPU‑first** execution (D3D12/Vulkan/Metal) with a robust **CPU fallback**.
- A **tile‑sparse** world, **APIC/MLS‑MPM** solvers, frictional/viscoplastic/two‑phase models, convex contact, and game‑ready writeback to engine terrain systems. fileciteturn0file1 fileciteturn0file8

**Outcome:** believable, controllable, _deterministic‑enough_ ground deformation at 60–120 Hz, packaged as a small native library you can ship across platforms. fileciteturn0file1

---

## 2) Why it exists (the problem)
Games typically fake deformable ground with decals, heightmap nudges, or VFX particles; these don’t affect locomotion or networked gameplay consistently. Studios need a **portable** and **deterministic** way to have terrain push back: sinkage, slip, ruts that persist, puddles that creep, and materials that transition with wetness—without being locked to one engine or GPU vendor. ziXn fills this production gap with a compact, cross‑engine SDK. fileciteturn0file1

---

## 3) What ziXn does (capabilities)
- **Materials:** dry/wet granular (sand, gravel), **snow** (powder→packed→slush), **mud/slurries**, layered beds, optional soil–water two‑phase coupling. fileciteturn0file4
- **Gameplay physics:** wheel/track/foot interaction, sinkage & traction metrics, berm/heave formation, bogging thresholds, avalanching cues. fileciteturn0file8
- **Runtime surfaces:** height/displacement & masks for wetness/debris to drive engine materials & effects. fileciteturn0file1
- **Performance:** 1–4 ms budgets on modern GPUs via tile sparsity, shared‑mem reductions, coarse pressure solves, and persistent scheduling. CPU fallback for servers/tools. fileciteturn0file1
- **Determinism:** fixed‑step, stable reductions, quantized I/O; strict and practical modes for netcode. fileciteturn0file1 fileciteturn0file6

---

## 4) How ziXn works (system overview)
**C‑ABI & object model.** A minimal C89 ABI exposes context/device/scene/material/rig/solver; handles are opaque; all functions return explicit status codes; no hidden threads. fileciteturn0file1

**Execution graph.** Each frame compiles a DAG of passes: **Activate tiles → P2G → Constitutive update → Contact projection → (Optional) pore‑pressure/viscous solve → G2P → Writeback → Diagnostics**. GPU implementations fuse tile‑local work and overlap independent passes; CPU mirrors this with OpenMP tasks. fileciteturn0file1

**Backends.** D3D12 (descriptor heaps, timeline fences), Vulkan (descriptor indexing, timeline semaphores), Metal (argument buffers, heaps); CPU fallback with SIMD loops. Shader source (HLSL) cross‑compiled for each backend. fileciteturn0file1

**Engine interop.**
- **Unreal:** ziXn emits RDG passes that composite **Runtime Virtual Texture (RVT)** height/masks; **Virtual Heightfield Mesh (VHM)** renders displacement; Niagara DI exposes fields for VFX. fileciteturn0file3
- **Unity:** native plugin + C# bindings; GPU→GPU height updates via **TerrainData.CopyActiveRenderTextureToHeightmap**; masks via SRP; `GraphicsBuffer/ComputeBuffer` interop for zero‑copy. fileciteturn0file2

---

## 5) Core algorithms (the “physics”)
- **APIC/MLS‑MPM transfers** to reduce dissipation and preserve angular momentum; terrain‑aligned slabs with thin vertical extent. fileciteturn0file8 fileciteturn0file0
- **Frictional plasticity:** Mohr–Coulomb and Drucker–Prager(+cap); **NorSand** for density‑dependent dilatancy. **Snow** uses an elasto‑plastic backbone with densification. Return‑mapping is robust, with clamps and line search. fileciteturn0file8 fileciteturn0file0
- **Mud/slurry:** **Herschel–Bulkley** with **Papanastasiou** regularization for yield‑stress flow; semi‑implicit viscous update on a coarse grid. fileciteturn0file0
- **Two‑phase (soil–water):** mixture drag `β≈μ/k`, effective stress `σ′=σ−αp_f I`, **pore pressure** diffusion solved (PCG/AMG) on coarsened tiles; saturation/infiltration optional. fileciteturn0file0
- **Contact:** grid‑space non‑penetration + Coulomb friction via projected/convex solvers; anisotropic ellipses for tracks/boots; compliance layer for stability. fileciteturn0file0

**Stability & LOD.** CFL and viscous/drag bounds determine substeps; model LOD swaps (granular ↔ HBP ↔ two‑phase) by moisture/gameplay need; tile/presentation LOD manage far‑field cost. fileciteturn0file8

---

## 6) Materials & calibration (authoring)
ziXn includes **exhaustive presets** with pragmatic ranges: sands (dry/damp), gravel, silt/loam/peat, clays, mud (HBP), and snow (powder→packed→slush). Each defines ρ, φ, c, ψ_d, E/ν or cap, HBP (τ_y,K,n,μ,m), and two‑phase (k,n,α, S). Calibrate with **column collapse** (fit φ), **inclined plane** (φ,c), **wheel trenching** (ψ_d, τ_y, k), **footprints** (cap, cohesion), **dam‑break over porous bed** (HBP+k). Metrics are concrete (runout error, berm volume, PCG residuals, sinkage RMS). fileciteturn0file4

---

## 7) Performance, determinism, and testing
- **Budgets & scaling.** Tile activation, shared‑mem reductions, coarse pressure solves, and async compute hit 1–4 ms step budgets on modern GPUs; CPU fallback mirrors the graph with fewer particles. fileciteturn0file1
- **Determinism.** Fixed dt, ordered reductions, quantized exports; strict mode disables fast‑math and uses fixed‑point paths where needed; practical mode allows tiny drift with correction. fileciteturn0file1
- **Validation & CI.** Invariants (mass/momentum bounds, non‑negative plastic work), golden scenarios (column collapse, inclined plane, footprints, wheel trenching, dam‑break), automated PIX/RenderDoc/Nsight captures, perf gates, nightly soaks. fileciteturn0file6

---

## 8) Roadmap (where this is going)
- **Phase 0:** APIC/MLS‑MPM slice + UE/Unity writeback + telemetry.  
- **Phase 1:** Plasticity suite, contact options, authoring presets.  
- **Phase 2:** HBP slurry + infiltration + coarse two‑phase.  
- **Phase 3:** Determinism & server build, tile streaming.  
- **Phase 4:** Scale & consoles, partner integrations. fileciteturn0file5

The **Agentic Dev Plan** breaks these into testable tasks with acceptance criteria so teams (or agents) can run workstreams in parallel. fileciteturn0file7

---

## 9) Getting started (orientation)
1. **Link the C ABI**: load the proc table (`zxGetProcTable`), create a context/device/scene, and choose a backend (D3D12/Vulkan/Metal or CPU). fileciteturn0file1  
2. **Create materials**: select a preset (e.g., _Dense Sand_, _Thick Mud_, _Wind‑Packed Snow_) and adjust φ/c/ψ_d or HBP/`k`. fileciteturn0file4  
3. **Bind engine targets**: Unreal → RVT/VHM; Unity → TerrainData + SRP; pass native resource handles for zero‑copy. fileciteturn0file3 fileciteturn0file2  
4. **Sim loop**: for each frame → activate tiles, `simulate(dt)`, `writeback()`, then draw. Inspect counters and timings; tune LOD if budgets exceed targets. fileciteturn0file1

**Minimal call flow (illustrative):**
```c
zx_procs P; zxGetProcTable(ZX_ABI_VERSION, &P);
zx_context C; P.create_context(&ctx_desc, &C);
zx_device  D; P.bind_device(C, &dev_desc, &D);
zx_scene   S; P.create_scene(C, &scene_desc, &S);
zx_material M_sand; P.create_material(C, &sand_desc, &M_sand);
for (;;) {
  P.begin_frame(S, &begin);
  P.simulate(S, &sim_params);   // fixed dt, optional substeps
  P.writeback(S, &wb);          // height/masks → engine
  P.end_frame(S, &end);
}
```

---

## 10) Repository layout (suggested)
```
/core        # C ABI, object model, allocators, logging
/kernels     # HLSL kernels (APIC/MLS‑MPM, constitutive, contact, pressure)
/backends    # d3d12 / vk / metal / cpu
/integrations/ue   # RDG/RVT/VHM + Niagara DI
/integrations/unity# Native plugin + C# + SRP
/tools       # material fitter, profiler, capture/replay
/docs        # architecture, algorithms, materials, testing, roadmap
/samples     # footprints, wheel trenching, bogging, dam‑break
```

---

## 11) License, support, and contributions
- **License**: choose MIT/BSD‑2/Commercial hybrid as suits your business model.
- **Support**: issue tracker + Discord/Slack; paid support available for source licensees.
- **Contributions**: CLA + style guide; CI must pass invariants/perf gates.

---

## 12) Sources (internal design docs)
- Architecture overview. fileciteturn0file1
- Algorithms detail and derivations. fileciteturn0file8 fileciteturn0file0
- Materials & calibration library. fileciteturn0file4
- Unreal integration (RDG/RVT/VHM/NI). fileciteturn0file3
- Unity integration (Native plugin/TerrainData/SRP). fileciteturn0file2
- Tooling & testing (CI/perf/invariants). fileciteturn0file6
- Roadmap & milestones. fileciteturn0file5
- Agentic development plan. fileciteturn0file7
