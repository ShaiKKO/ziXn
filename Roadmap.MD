/*!
\file Roadmap.MD
\brief Program roadmap and delivery plan for ziXn.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
*/

<!--
Proprietary Notice:
This document contains confidential and proprietary information of Colin Macritchie / Ripple Group, LLC. 
No part of this work may be disclosed, reproduced, or distributed without prior written consent.
-->

## ziXn — Program Roadmap (2025–2026)

### Executive Summary
ziXn is a portable, engine‑agnostic terrain continuum SDK delivering deformable granular, snow, and slurry materials in real time. This roadmap aligns algorithmic capabilities (APIC/MLS‑MPM; Mohr–Coulomb/Drucker–Prager+cap; NorSand; HBP slurries; two‑phase soil–water), the C‑ABI and backend architecture (D3D12/Vulkan/Metal; CPU fallback), and authoring/calibration workflows into phased milestones with measurable KPIs, rigorous acceptance criteria, and production hardening. Targets are 60–120 Hz with 1–4 ms simulation budgets on modern GPUs, deterministic‑enough stepping, and robust engine writeback for Unreal and Unity.

### Program Themes
- **Physics depth with scalability**: particle–grid hybrids, frictional/viscoplastic/two‑phase models with LOD across model, spatial, and presentation layers.
- **Stable C‑ABI**: opaque handles, strict size/versioning, deterministic reductions, engine‑agnostic interop.
- **GPU‑first with CPU parity**: SIMT‑friendly SoA, shared‑memory reductions, coarse pressure solves; OpenMP/SIMD for CPU.
- **Authoring & calibration**: comprehensive presets and repeatable fits (column collapse, inclined plane, wheel trenching, footprints, dam‑break).
- **Observability, security, and CI**: counters, traces, invariants, captures; input validation; hard caps; binary/package hygiene.

### Success Criteria (Top‑Level KPIs)
- **Performance**: ≤3.5 ms/frame @ 2.0 M particles (PC high‑end) with pressure solve; ≤6–8 ms on consoles’ heavy scenes.
- **Determinism**: identical determinism checksum across supported backends within tolerance; strict mode disables fast‑math as needed.
- **Stability**: mass drift <0.2% over 60 s; PCG residual <1e‑3 within 10–16 iterations (coarse grid).
- **Interop**: Unreal RDG/RVT/VHM and Unity TerrainData/SRP writeback with zero‑copy handles and sample scenes.
- **Calibration**: presets validated to within target error bands (runout <10%, berm volume <15%).

---

## Phases, Milestones, and Deliverables

### Phase 0 — Foundations (M0–M1)
- **Scope**: C‑ABI skeleton, context/device/scene lifecycles; tile‑sparse grid; APIC/MLS transfers; basic contact; engine writeback; counters & timings.
- **Deliverables**:
  - `zx_procs` table with versioning/size checks; opaque handles and idempotent destroy.
  - Tile manager (hash, activation, neighbor lists), particle/grid SoA, APIC P2G/G2P.
  - Deterministic tile‑local reductions; initial contact projection; boundary clamps.
  - Unreal and Unity integration stubs with displacement/mask writeback; sample scenes (footprints, basic ruts).
  - Telemetry: counters and GPU/CPU timestamps; invariants framework (mass, energy, contact gap).
- **KPIs**: 120 Hz dry granular slice ≤2.0 ms @ 1.2 M particles; determinism checksum stable across CPU/GPU on small scenes.
- **Acceptance**: API conformance tests; binary compatibility checks; capture scripts; docs baseline.

### Phase 1 — Plasticity Suite & Contact (M2–M3)
- **Scope**: robust return mapping for MC/DP(+cap), NorSand state, anisotropic/convex contact, presentation LOD; authoring presets v1.
- **Deliverables**:
  - Elastic predictor with polar decomposition; MC, DP(+cap) return mapping; NorSand state parameter flows.
  - Contact projection with anisotropic tangential cone; rolling resistance proxy; compliance layer.
  - Authoring presets for sands, gravel, soils, clays, snow; presentation LOD filters.
- **KPIs**: heave fidelity (dilatancy) matched in wheel/foot tests; footprint crispness via cohesion/cap; ≤3.0 ms at parity particle count.
- **Acceptance**: column collapse and inclined plane validation scenes pass targets; contact gap metrics within bounds.

### Phase 2 — Slurry & Two‑Phase (M4–M5)
- **Scope**: HBP with Papanastasiou regularization, semi‑implicit viscous update on coarse grid, two‑phase mixture pressure solve.
- **Deliverables**:
  - Effective viscosity clamp; Jacobi/PCG with AMG/multigrid smoothing on coarsened tiles.
  - Mixture coupling (drag β≈μ/k), effective stress, pore‑pressure diffusion; saturation hooks.
  - Integration tests: dam‑break over porous bed, bog memory tuning.
- **KPIs**: pressure residual <1e‑3 in ≤16 iters; stable at dt targets via substeps; ≤3.5 ms total at 2.0 M particles.
- **Acceptance**: slurry/mud presets calibrated; validation scenes pass; telemetry publishes PCG traces.

### Phase 3 — Determinism, Streaming, and Server (M6–M7)
- **Scope**: strict determinism path; tile streaming and residency; CPU parity for headless/server.
- **Deliverables**:
  - Fixed‑point or F32‑only strict reductions; FMAD variance controls; cross‑platform checksum harness.
  - Resident set manager (heaps/MTLHeaps/VMA), streaming hysteresis; tile compaction.
  - CPU backend with SIMD kernels (AVX2/AVX‑512/SVE) and NUMA pinning.
- **KPIs**: checksum stable across D3D12/Vulkan/Metal/CPU; server tick within budget for low particle densities.
- **Acceptance**: network determinism soak tests; capture/replay parity; streaming thrash below threshold.

### Phase 4 — Consoles & Partner Integrations (M8–M9)
- **Scope**: platform hardening, memory pressure tuning, partner engine integrations.
- **Deliverables**:
  - Console backends with tuned threadgroup limits; memory budgets met.
  - Partnered UE/Unity samples and docs polish; authoring UX refinement.
- **KPIs**: consoles heavy scene ≤8 ms worst‑case; no validation regressions.
- **Acceptance**: partner sign‑offs; CI green with perf gates on all platforms.

### Phase 5 — Tooling, Calibration, and Commercialization (M10–M12)
- **Scope**: material fitter tool, advanced authoring, hardening, and packaging.
- **Deliverables**:
  - Material fitter minimizing composite loss (runout, rut depth, slip curve RMS).
  - Docs: algorithms, architecture, materials, testing, integration; samples suite complete.
  - Packaging: shared library + headers + JSON capabilities manifest; support/process docs.
- **KPIs**: fitter converges on target scenes; onboarding ≤1 day for new teams.
- **Acceptance**: end‑to‑end validation battery passes; readiness review complete.

---

## Cross‑Cutting Engineering Tracks
- **Data layout & precision**: F32 determinism path; optional F16 for selected fields with F32 accumulators; 256‑byte alignment.
- **Execution graph**: asynchronous compute, explicit barriers, shared‑mem tile reductions; co‑execution of pressure and presentation LOD where possible.
- **Telemetry**: counters (particles, tiles, atomics, PCG iters, dt clamps), timings, invariant failures with frame context; auto‑capture hooks.
- **Security & robustness**: strict input validation, hard caps on tiles, NaN/Inf sanitization, ref‑counted lifetimes, ABA‑safe pools.
- **CI/CD**: fat binaries per platform; GPU validation toggles; CPU sanitizers; perf gates; nightly soak tests and drift audits.

---

## Risks and Mitigations
- **Solver instability at extreme stickiness** → widen μ_eff clamp, adaptive substeps, contact compliance tuning.
- **Tile thrash under streaming** → activation hysteresis, region pinning, prefetch hints, present‑LOD fallback.
- **Cross‑backend numeric drift** → deterministic reductions, fixed‑point paths where needed, checksum gates in CI.
- **Memory pressure on consoles** → resident set compaction, quantization toggles, model/presentation LOD trade‑offs.
- **Interop fragility** → narrow surfaces (buffers + fences), exhaustive conformance tests, sample parity.

---

## Versioning, Release Cadence, and Governance
- **Semantic ABI versioning** in `zx_procs` table; structs begin with `size` for forward/back compatibility.
- **Minor releases** after each phase with patch releases for perf/robustness; long‑term support branch after Phase 4.
- **Change control** via RFCs; design reviews for ABI changes; performance regression triage every sprint.

---

## Acceptance Test Matrix (excerpt)
- Dry granular column collapse: runout error <10%, mass drift <0.2%/60 s.
- Inclined plane: fitted φ (±2°), c within targeted range.
- Wheel trenching: berm volume error <15%; heave profile within tolerance; slip curve RMS minimized.
- Dam‑break (slurry/two‑phase): PCG residual and iters within bounds; visual creep consistent.
- Cross‑backend checksum: identical within tolerance for reference seeds.

---

## Staffing & Roles (indicative)
- **Core/ABI/Backends**: 1–2 senior GPU/low‑level engineers.
- **Simulation/Models**: 1 computational physicist/graphics engineer.
- **Interop/Integrations**: 1 engine integration engineer.
- **Tooling/CI**: 0.5–1 tools engineer.
- **QA/Validation**: 0.5 SDET with physics background.

---

## Artifacts per Phase
- Updated headers and shared library; validated kernels; backend binaries; integration plugins; samples; captures; performance reports; validation dossiers; updated docs.


