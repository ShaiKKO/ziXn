
# ziXn Roadmap — From Prototype to Full Shipping Library
_Granular / slurry / snow terrain for real‑time engines (C‑ABI; GPU‑first; deterministic‑enough)_

**Scope of this roadmap.** Synthesizes all decisions and targets from our conversation into a single execution plan that carries **ziXn** from early prototype through a **1.0 production release**, with parallel **BeamNG test harness** work and subsequent **Unreal/Unity adapters**. Each phase lists **goals, deliverables, acceptance/exit criteria, resourcing, and risks**.

---

## Phase 1 — Core Spine (Foundations) — _Weeks 0–6_
**Goal.** Stand up the minimal vertical slice: C‑ABI, opaque handles, fixed‑step simulate, tile‑sparse world, APIC/MLS‑MPM transfers, dry granular MC/DP, basic frictional grid contact, CPU fallback, GPU backend skeleton, and counters/telemetry.

**Deliverables**
- `libzixn` C‑ABI: context/device/scene/material/rig; error model; logging/telemetry callbacks.
- Data model: particle **SoA**, tile‑sparse grid (16×16×Nz), Morton sort, per‑tile neighbor lists.
- Kernels: **P2G**, elastic predictor, MC/DP return mapping, grid update, **contact projection (PGS)**, **G2P**; determinism scaffolding (ordered reductions).
- CPU backend (OpenMP + AVX2/SVE2) parity; D3D12 backend skeleton (descriptor heap, root sig, timestamp queries).
- Tooling: minimal capture (particle header dump), perf counters, assertions.

**Exit criteria**
- 1M particles @ h≈8–12 cm runs **≤3.5 ms** on mid‑range GPU; CPU fallback runs ≤30 ms.
- Mass conservation and non‑negative plastic work invariants hold over 60 s soak.
- Determinism smoke test: same checksum across two runs, same machine, strict mode.

**Risks & notes**
- Kernel divergence; mitigate with SoA and shared‑mem tile reductions.
- Early over‑engineering; keep models minimal (MC or DP without cap).

---

## Phase 2 — Materials Suite (Granular & Snow) — _Weeks 6–12_
**Goal.** Expand constitutive coverage and stability for gameplay fidelity; add authorable presets + calibration flows.

**Deliverables**
- Constitutive: **DP+Cap**, **NorSand** (critical‑state), **Snow** (elasto‑plastic with densification proxy).
- Return‑mapping with line search and clamping; cap evolution.
- Authoring presets for: loose/dense sand, gravel, silty/loamy soil, soft clay, powder/wind‑packed/wet snow.
- Calibration scenes: column collapse, inclined plane, wheel trenching, footprints.
- Presentation LOD: displacement + wetness/debris mask export (engine‑agnostic targets).

**Exit criteria**
- Column collapse runout error <10% vs baseline; slope angle within ±2°.
- Wheel trenching test shows berm/heave qualitatively correct; reproducible across strict mode.
- 60–120 Hz stable at target particle densities.

**Risks**
- NorSand stability at large Δt → bound plastic increment; add substep triggers.

---

## Phase 3 — Slurry & Two‑Phase Tier — _Weeks 12–20_
**Goal.** Introduce mud/slurry with **Herschel–Bulkley + Papanastasiou (HBP)** and a **coarsened pore‑pressure** two‑phase option.

**Deliverables**
- HBP viscous term (semi‑implicit, coarse grid PCG) with μ_eff clamps.
- Two‑phase mixture: drag `β≈μ/k`, effective stress `σ′=σ−αp_f I`, coarse PCG for `p_f` with AMG/Jacobi preconditioner.
- Model LOD switching: dry ↔ HBP ↔ two‑phase by saturation/pressure heuristics.
- Presets: thin mud, thick mud, saturated sands; permeability and porosity tables.
- Counters: PCG iterations, residuals, drag work, pressure cache hit rate.

**Exit criteria**
- Coarse pressure solve ≤0.9 ms on target GPU for typical active footprint; PCG residual <1e‑3.
- Bogging test: correct stuck/free threshold tunable by τ_y, k; repeatable sinkage curves.

**Risks**
- Over‑tight coupling costs; overlap pressure with presentation LOD; allow async degradation.

---

## Phase 4 — Contact, Rigs, and Anisotropy — _Weeks 20–26_
**Goal.** Robust terrain–rig coupling and traction realism.

**Deliverables**
- Grid‑space **convex** contact option (variational) for harsh stacks; default PGS remains.
- **Anisotropic friction ellipses** and **compliance layers** for wheels/tracks/boots.
- Rig API: wheels/feet/track bands authoring; per‑rig material overrides (e.g., studs/chains).
- Determinism mode hardening: fixed‑point reductions option; stable partitioning.

**Exit criteria**
- Vehicle traction tests reproduce expected slip/μ curves across materials.
- Convex mode survives aggressive Δt and stack scenarios without chatter.

**Risks**
- Convex solver cost; keep iteration caps and fall back to PGS on budget pressure.

---

## Phase 5 — Tooling, Determinism & CI — _Weeks 26–32_
**Goal.** Make correctness and performance observable and enforceable.

**Deliverables**
- Invariants: mass/momentum checks; penetration caps; NaN/Inf guards.
- GPU/CPU telemetry: per‑pass timestamps; stall breakdown; atomics pressure.
- Golden scenes & CI: column collapse, inclined plane, wheel trenching, footprints, dam‑break; perf gates per backend.
- Snapshot/Replay: delta‑compressed tile+particle snapshots; deterministic replays across strict mode.

**Exit criteria**
- Nightly CI runs all golden scenes under time + invariant budgets.
- Cross‑device reproducibility (strict determinism path) within checksum tolerance.

---

## Phase 6 — BeamNG Test Harness (Adapter) — _Weeks 16–34 (overlaps P3–P5)_
**Goal.** Validate ziXn in a high‑fidelity driving sandbox without touching core.

**Deliverables**
- **Out‑of‑process native adapter** (`zixn-beamng-adapter.exe`) + Lua bridges (`zixn_bridge` GELUA, `zixn_probe` VLUA).
- Telemetry path (wheels/contact → contact tiles) and fixed‑rate ziXn stepping.
- Physics feedback (sinkage/drag) via VLUA `applyForce`; visual ruts via decals; optional bounded terrain edits.
- Debug UI overlays (contact heatmap, rut depth, perf HUD); recording/replay across adapter boundary.

**Exit criteria**
- 2 vehicles @ 60 Hz with adapter: ≤1 ms GELUA and ≤0.1 ms VLUA/veh overhead; end‑to‑end feedback latency ≤2 frames.
- Clear traction differences across presets; visible ruts where enabled; stable long‑run soak (30 min).

**Risks**
- Terrain edit costs → keep visual edits bounded; prefer decals; inject physics via forces primarily.

---

## Phase 7 — Engine Adapters (Unreal & Unity) — _Weeks 30–44_
**Goal.** Ship reference adapters while keeping core engine‑agnostic.

**Deliverables (Unreal)**
- RDG pass integration; **RVT** writeback (height/masks); **VHM** consumption; Niagara Data Interface for flow/wetness.
- Zero‑copy buffers via RHI handles; async compute scheduling; sample map and materials.

**Deliverables (Unity)**
- Native plugin (IUnityGraphics); `GraphicsBuffer/ComputeBuffer` interop; `TerrainData` writeback helpers.
- SRP example (HDRP/URP) with shader graph for masks/displacement.

**Exit criteria**
- Unreal & Unity samples run @ 60 Hz with visible deformation + physics feedback hooks (force/friction scalars).

---

## Phase 8 — Scale, Consoles & Optimization — _Weeks 44–56_
**Goal.** Production hardening on target hardware and content scales.

**Deliverables**
- Tile streaming & residency policies; LOD tuning; async overlap (pressure ↔ presentation).
- GPU optimizations: subgroup ops, LDS tiling, atomics minimization; pipeline cache prewarm.
- Console bring‑up checklists; memory budgets & allocator telemetry.

**Exit criteria**
- Large scene (≥2 km² with active 256–512 m² patches) holds 60 Hz on target SKUs within memory budgets.
- Heatmaps show no pathological stalls; PCG solves within iteration caps.

---

## Phase 9 — 1.0 Release & SDK — _Weeks 56–64_
**Goal.** Ship a supported SDK with docs, samples, and versioned ABI.

**Deliverables**
- Stable **C‑ABI** versioned; header docs; samples for Unreal/Unity; BeamNG adapter optional showcase.
- Material library & fitter; profiling guide; determinism/netcode guide.
- Licensing & packaging; crash‑safe release builds; symbol server.

**Exit criteria**
- External team (pilot studio) integrates ziXn using only public docs + samples.
- Support processes in place (issue triage, crash intake, nightly artifacts).

---

## Cross‑cutting Concerns (all phases)
- **Determinism modes:** strict vs practical; document flags and costs.
- **Error handling:** explicit status codes; never partial writes on failure.
- **Security/robustness:** bounds checks; watchdog clamp on Δt; graceful degrade.
- **Docs:** living specs (`ARCHITECTURE.md`, `ALGORITHMS.md`, `MATERIALS_AND_CALIBRATION.md`, `INTEGRATION_*`).
- **Agentic work breakdown:** each phase decomposed to ≤3‑day tasks with artifacts & tests.

---

## Staffing & Roles (suggested)
- **Core physics (2):** kernels, constitutive models, stability.
- **GPU engineer (1):** backends, async scheduling, profiling.
- **Systems/interop (1):** C‑ABI, allocators, determinism, snapshots.
- **Adapter engineer (1):** BeamNG, Unreal/Unity bridges, tooling.
- **QA/TD (0.5):** validation scenes, perf gates, doc polish.

---

## KPIs & Telemetry
- Frame time per pass (GPU/CPU), active tiles, particles/cell, PCG iters, contact saturation %.
- Validation metrics: runout %, berm volume %, sinkage RMS, determinism checksum drift.
- Adapter metrics: end‑to‑end latency, VLUA/GELUA budget, packet loss.

---

## Risk Register (top items)
- Stability under high cohesion/viscosity → clamp μ_eff, substep triggers, line search.
- Terrain writeback perf (BeamNG) → prefer decals/forces; bound edits.
- Cross‑platform determinism → fixed partitioning and ordered reductions; disable fast‑math in strict mode.

---

## Versioning & Release Cadence
- **0.x** feature releases bi‑weekly; ABI may change.  
- **1.0‑rc** freeze ABI; fix only defects and perf regressions.  
- **1.y** minor: additive API; **1.y.z** patch: bug‑fix only.

---

## Artifacts per Phase
- Binary builds (core + backends), kernel cache, symbol files.
- Samples (footprints, wheel trenching, dam‑break), adapter demos.
- Docs & change log; perf dashboards; replay bundles for golden scenes.
