# ziXn

Note: Invalid Unicode symbols replaced with ASCII citation/marker text.

A high-performance C-ABI terrain continuum library for **granular**, **snow**, and **slurry** materials, built for real-time engines (Unreal, Unity, custom).

---

## 1) What ziXn is

**ziXn** is a portable simulation SDK that makes **deformable terrain** a first-class runtime system: footprints in powder, ruts in wet sand, wheel trenching and berms, bogging in mud, slush and landslide-like flows, and layered beds (e.g., sand over clay). It provides:

- A **stable C ABI** with opaque handles; deterministic stepping; zero-copy interop.
- **GPU-first** execution (D3D12/Vulkan/Metal) with a robust **CPU fallback**.
- A **tile-sparse** world, **APIC/MLS-MPM** solvers, frictional/viscoplastic/two-phase models, convex contact, and game-ready writeback to engine terrain systems.

**Outcome:** believable, controllable, *deterministic-enough* ground deformation at 60–120 Hz, packaged as a small native library you can ship across platforms.

---

## 2) Why it exists (the problem)

Games typically fake deformable ground with decals, heightmap nudges, or VFX particles; these don't affect locomotion or networked gameplay consistently. Studios need a **portable** and **deterministic** way to have terrain push back: sinkage, slip, ruts that persist, puddles that creep, and materials that transition with wetness—without being locked to one engine or GPU vendor. ziXn fills this production gap with a compact, cross-engine SDK.

---

## 3) What ziXn does (capabilities)

- **Materials:** dry/wet granular (sand, gravel), **snow** (powder→packed→slush), **mud/slurries**, layered beds, optional soil–water two-phase coupling.
- **Gameplay physics:** wheel/track/foot interaction, sinkage & traction metrics, berm/heave formation, bogging thresholds, avalanching cues.
- **Runtime surfaces:** height/displacement & masks for wetness/debris to drive engine materials & effects.
- **Performance:** 1–4 ms budgets on modern GPUs via tile sparsity, shared-mem reductions, coarse pressure solves, and persistent scheduling. CPU fallback for servers/tools.
- **Determinism:** fixed-step, stable reductions, quantized I/O; strict and practical modes for netcode.

---

## 4) How ziXn works (system overview)

**C-ABI & object model.** A minimal C89 ABI exposes context/device/scene/material/rig/solver; handles are opaque; all functions return explicit status codes; no hidden threads.

**Execution graph.** Each frame compiles a DAG of passes: **Activate tiles → P2G → Constitutive update → Contact projection → (Optional) pore-pressure/viscous solve → G2P → Writeback → Diagnostics**. GPU implementations fuse tile-local work and overlap independent passes; CPU mirrors this with OpenMP tasks.

**Backends.** D3D12 (descriptor heaps, timeline fences), Vulkan (descriptor indexing, timeline semaphores), Metal (argument buffers, heaps); CPU fallback with SIMD loops. Shader source (HLSL) cross-compiled for each backend.

**Engine interop.**

- **Unreal:** ziXn emits RDG passes that composite **Runtime Virtual Texture (RVT)** height/masks; **Virtual Heightfield Mesh (VHM)** renders displacement; Niagara DI exposes fields for VFX.
- **Unity:** native plugin + C# bindings; GPU→GPU height updates via **TerrainData.CopyActiveRenderTextureToHeightmap**; masks via SRP; `GraphicsBuffer/ComputeBuffer` interop for zero-copy.

---

## 5) Core algorithms (the "physics")

- **APIC/MLS-MPM transfers** to reduce dissipation and preserve angular momentum; terrain-aligned slabs with thin vertical extent.
- **Frictional plasticity:** Mohr–Coulomb and Drucker–Prager(+cap); **NorSand** for density-dependent dilatancy. **Snow** uses an elasto-plastic backbone with densification. Return-mapping is robust, with clamps and line search.
- **Mud/slurry:** **Herschel–Bulkley** with **Papanastasiou** regularization for yield-stress flow; semi-implicit viscous update on a coarse grid.
