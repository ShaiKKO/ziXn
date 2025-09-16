/*!
\file Todo.md
\brief Execution backlog for ziXn with prioritized tasks and acceptance criteria.
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC. All rights reserved.
*/

<!-- Proprietary and Confidential. Distribution prohibited without written consent. -->

## Conventions
- **Priority**: P0 (must), P1 (should), P2 (nice).
- **Type**: FEAT, TECH, DOC, TEST, PERF, SEC.
- **AC**: Acceptance Criteria (measurable, testable).
- **Traceability**: Each item maps to Roadmap phase and document anchors (`ALGORITHMS.md`, `ARCHITECTURE.md`, `MATERIALS_AND_CALIBRATION.md`).

---

## Phase 0 — Foundations

1. P0 • FEAT • C‑ABI skeleton and proc table
   - Desc: Implement `zxGetProcTable`, opaque handles, version/size checks, status codes, idempotent destroy.
   - AC: ABI conformance tests pass; size‑gate rejections covered; error strings mapped; thread‑safety doc published.

2. P0 • FEAT • Tile‑sparse grid and particle SoA
   - Desc: Tile hash, activation/deactivation with hysteresis; neighbor lists; particle/grid SoA per docs.
   - AC: Activation ≤0.15 ms @ 2.0 M particles; active tile count matches reference in sample scenes.

3. P0 • FEAT • APIC/MLS P2G/G2P transfers
   - Desc: Quadratic B‑splines; affine velocity; internal force scatter; deterministic tile‑local reduction.
   - AC: Energy budget sane; determinism checksum stable across 3 consecutive runs; unit tests for transfer identities.

4. P0 • FEAT • Basic contact and boundary
   - Desc: Grid‑node contact projection (normal clamp + Coulomb tangent); boundary clamps.
   - AC: Max penetration below threshold; friction cone saturation metrics emitted; footprint sample works.

5. P0 • FEAT • Engine writeback (UE/Unity stubs)
   - Desc: Displacement + masks writeback surfaces; zero‑copy handles; sample scenes boot.
   - AC: UE RVT/VHM renders displacement; Unity TerrainData path updates height; frame budgets respected.

6. P0 • TEST • Telemetry & invariants
   - Desc: Counters, GPU/CPU timings; mass/momentum checks; capture hooks (PIX/Nsight/RenderDoc).
   - AC: Counters surfaced; invariant failures escalate with context; auto‑capture on failure.

---

## Phase 1 — Plasticity Suite & Contact

7. P0 • FEAT • Elastic predictor + MC/DP(+cap) return mapping
   - AC: Column collapse/inclined plane tests within targets; plastic work non‑negative; line search stabilizes tough cases.

8. P0 • FEAT • NorSand state and dilatancy
   - AC: Wheel/foot heave profiles match references; state parameter ψ tracked and persisted.

9. P1 • FEAT • Anisotropic/convex contact
   - AC: Track/boot anisotropy ellipse validated; rolling resistance proxy responsive to slip.

10. P1 • DOC • Authoring presets v1
    - AC: Presets for sands/gravel/soils/clays/snow documented with safe ranges and UI clamps.

11. P1 • PERF • Presentation LOD filters
    - AC: Far‑field cost reduced with no visible popping; height/mask mip consistency verified.

---

## Phase 2 — Slurry & Two‑Phase

12. P0 • FEAT • HBP viscous update with Papanastasiou
    - AC: Effective viscosity clamped; dam‑break scene stable; ≤0.9 ms on coarse grid at targets.

13. P0 • FEAT • Two‑phase mixture pressure solve
    - AC: PCG residual <1e‑3 in ≤16 iters; drag β≈μ/k tunable; bog memory knob documented.

14. P1 • TEST • Integration scenes and metrics
    - AC: Dam‑break over porous bed, wheel bogging, puddle creep scenes pass metrics; telemetry logs exported.

---

## Phase 3 — Determinism, Streaming, Server

15. P0 • TECH • Strict determinism mode
    - AC: Cross‑backend checksum parity; FMAD controls toggled; fixed‑point or F32‑only reductions documented.

16. P0 • FEAT • Tile streaming and residency
    - AC: Thrash below threshold; page‑in/out traces stable; memory caps enforced with graceful degrade.

17. P1 • FEAT • CPU backend parity
    - AC: OpenMP/SIMD kernels for P2G/G2P/constitutive pass unit tests and timing thresholds.

---

## Phase 4 — Consoles & Partner Integrations

18. P0 • PERF • Console backends hardening
    - AC: Threadgroup memory tuned; worst‑case heavy scene ≤8 ms; validation parity maintained.

19. P1 • FEAT • Partner UE/Unity samples & UX
    - AC: Partner sign‑offs; onboarding guide polished; sample parity across engines.

---

## Phase 5 — Tooling, Calibration, Commercialization

20. P0 • FEAT • Material fitter tool
    - AC: Minimizes composite loss (runout, rut depth, slip RMS); converges on preset scenes.

21. P1 • DOC • Documentation completion
    - AC: Algorithms/Architecture/Materials/Testing/Integration docs finalized; API reference generated.

22. P1 • TECH • Packaging & manifest
    - AC: Headers + shared lib + JSON capabilities manifest; CI publishes per‑platform artifacts.

---

## Continuous Tracks

23. P0 • SEC • Input validation & hard caps
    - AC: Fuzz tests; oversize inputs rejected; no partial writes on failure.

24. P0 • TEST • CI gates and perf budgets
    - AC: Invariants, perf gates, nightly soaks, drift audits; auto‑capture on regression.

25. P1 • OBS • Telemetry depth
    - AC: Counters expanded (atomics, compactions, dt clamps, substeps); error triage with frame context.

---

## Detailed Backlog (Actionable Checklists)

### 1. C‑ABI skeleton and proc table (P0 • FEAT)
- Owner: Core/ABI
- Estimate: 2–3 weeks
- Dependencies: None
- Artifacts: `zx_procs.h`, shared library exports, ABI conformance tests
- Checklist:
  - [ ] Define `zx_status`, opaque `zx_handle` types, and base structs with `size` prefix
  - [ ] Implement `zxGetProcTable(abi_version, zx_procs*)` with semantic versioning
  - [ ] Populate proc table: context/device/scene/material/rig/frame/telemetry
  - [ ] Enforce size/version gates and zero‑on‑failure out‑params
  - [ ] Error string table and logging sink integration
  - [ ] Threading and reentrancy policy documentation
  - [ ] ABI conformance test suite (size mismatch, nulls, double destroy)
  - [ ] CI job to run ABI tests on Windows/Linux/macOS

### 2. Tile‑sparse grid and particle SoA (P0 • FEAT)
- Owner: Core/Backend GPU
- Estimate: 2 weeks
- Dependencies: 1
- Artifacts: `tiles.hlsl`, `tile_hash.cpp`, SoA buffers
- Checklist:
  - [ ] Design tile coordinates, `B` size, and ghost margins
  - [ ] Implement GPU/CPU tile hash map with generation counters
  - [ ] Activation/deactivation with hysteresis and residency caps
  - [ ] Neighbor list precompute per tile
  - [ ] Particle SoA layout (positions, velocities, mass, volume, F, C, plastic state)
  - [ ] Bounds checks and NaN/Inf sanitization paths
  - [ ] Activation timing benchmark ≤0.15 ms @ spec scene

### 3. APIC/MLS P2G/G2P transfers (P0 • FEAT)
- Owner: Simulation/Backend GPU
- Estimate: 2 weeks
- Dependencies: 2
- Artifacts: `p2g.hlsl`, `g2p.hlsl`, unit tests
- Checklist:
  - [ ] Quadratic B‑spline weights and gradients implementation
  - [ ] APIC affine velocity scatter/gather with clamping
  - [ ] Internal force scatter using MLS‑consistent stress divergence
  - [ ] Tile‑local shared memory reduction; deterministic merge
  - [ ] Unit tests for partition of unity and momentum consistency
  - [ ] Determinism checksum harness for transfers

### 4. Basic contact and boundary (P0 • FEAT)
- Owner: Simulation
- Estimate: 1–2 weeks
- Dependencies: 3
- Artifacts: `contact.hlsl`, SDF interface, tests
- Checklist:
  - [ ] Grid‑node penetration detection via SDF sampling
  - [ ] Normal projection with compliance `κ_n`
  - [ ] Tangential Coulomb cone projection; saturation metrics
  - [ ] Boundary clamps and sticky flags where needed
  - [ ] Footprint sample validation and metrics logging

### 5. Engine writeback (UE/Unity stubs) (P0 • FEAT)
- Owner: Interop
- Estimate: 2 weeks
- Dependencies: 3,4
- Artifacts: UE module, Unity plugin, sample scenes
- Checklist:
  - [ ] Unreal RDG pass that binds ziXn buffers and triggers `simulate`
  - [ ] RVT/VHM displacement/mask writeback with tiling and clipping
  - [ ] Unity native plugin with `GraphicsBuffer/ComputeBuffer` interop
  - [ ] TerrainData height/mask update path and SRP hookup
  - [ ] Sample scenes boot and render deformations at target frame times

### 6. Telemetry & invariants (P0 • TEST)
- Owner: Tools/CI + Core
- Estimate: 1 week
- Dependencies: 3–5
- Artifacts: counters API, timings, invariant checks, capture hooks
- Checklist:
  - [ ] Counters: particles, active tiles, atomics, compactions, substeps
  - [ ] GPU/CPU timestamps per pass; frame budget summary
  - [ ] Invariants: mass/momentum conservation, contact gap bounds
  - [ ] Auto‑capture hooks for PIX/Nsight/RenderDoc on failure
  - [ ] CI export of telemetry logs and threshold gates

### 7. Elastic predictor + MC/DP(+cap) (P0 • FEAT)
- Owner: Simulation
- Estimate: 3 weeks
- Dependencies: 3,6
- Artifacts: `constitutive.hlsl`, unit tests, validation scenes
- Checklist:
  - [ ] Polar decomposition and elastic trial stress (Neo‑Hookean/corotated)
  - [ ] MC/DP yield surfaces; parameter mapping MC→DP
  - [ ] Cap surface with plastic volumetric strain evolution
  - [ ] Robust return mapping with line search and clamps
  - [ ] Column collapse and inclined plane validations within targets

### 8. NorSand state and dilatancy (P0 • FEAT)
- Owner: Simulation
- Estimate: 2 weeks
- Dependencies: 7
- Artifacts: NorSand state tracking, tests
- Checklist:
  - [ ] State parameter `ψ` from density/void ratio; critical state curve
  - [ ] Evolving `M(ψ)` if enabled; non‑associated flow tuning
  - [ ] Wheel/foot heave profile validations

### 9. Anisotropic/convex contact (P1 • FEAT)
- Owner: Simulation
- Estimate: 2 weeks
- Dependencies: 4,7
- Artifacts: Anisotropic tangent cone, tests
- Checklist:
  - [ ] Elliptical tangential limits per node
  - [ ] Rolling resistance proxy linked to slip ratio
  - [ ] Track/boot validation scenes and metrics

### 10. Authoring presets v1 (P1 • DOC)
- Owner: Simulation + Docs
- Estimate: 1 week
- Dependencies: 7,8
- Artifacts: Preset JSON/TOML, docs, UI ranges
- Checklist:
  - [ ] Presets for sands, gravel, soils, clays, snow with safe ranges
  - [ ] Authoring UI mappings and clamps
  - [ ] Examples and validation numbers

### 11. Presentation LOD filters (P1 • PERF)
- Owner: Backend GPU + Interop
- Estimate: 1 week
- Dependencies: 5
- Artifacts: LOD kernels, mip rules
- Checklist:
  - [ ] Far‑field coarsening and displacement filtering
  - [ ] Mip consistency and crack avoidance
  - [ ] Visual tests without popping artifacts

### 12. HBP viscous update (P0 • FEAT)
- Owner: Simulation
- Estimate: 2 weeks
- Dependencies: 7
- Artifacts: `viscous.hlsl`, coarse grid solver
- Checklist:
  - [ ] Herschel–Bulkley with Papanastasiou regularization
  - [ ] Semi‑implicit update on coarsened grid (Jacobi/PCG)
  - [ ] μ_eff clamp policy; stability at target dt
  - [ ] Dam‑break scene stability and timing

### 13. Two‑phase mixture pressure solve (P0 • FEAT)
- Owner: Simulation
- Estimate: 3 weeks
- Dependencies: 12
- Artifacts: `pressure_pcg.hlsl`, mixture fields
- Checklist:
  - [ ] Effective stress, drag coupling β≈μ/k, porosity/saturation fields
  - [ ] Diffusion/Poisson assembly on coarse tiles
  - [ ] PCG with AMG/multigrid smoothing
  - [ ] Bog memory knob and tests

### 14. Integration scenes and metrics (P1 • TEST)
- Owner: QA/SDET + Simulation
- Estimate: 1 week
- Dependencies: 12,13
- Artifacts: Automated scenes, metric extractors
- Checklist:
  - [ ] Dam‑break over porous bed with target residuals
  - [ ] Wheel bogging and puddle creep metrics
  - [ ] Telemetry export and dashboards

### 15. Strict determinism mode (P0 • TECH)
- Owner: Core/Backend + Tools/CI
- Estimate: 2 weeks
- Dependencies: 3,6
- Artifacts: Build flags, checksum harness, CI jobs
- Checklist:
  - [ ] Disable FMAD variance; fixed reduction order paths
  - [ ] Optional fixed‑point reductions where required
  - [ ] Cross‑backend checksum parity tests

### 16. Tile streaming and residency (P0 • FEAT)
- Owner: Core/Backend GPU
- Estimate: 2 weeks
- Dependencies: 2,5
- Artifacts: Residency manager, compaction tools
- Checklist:
  - [ ] Resident set manager for heaps/MTLHeaps/VMA
  - [ ] Activation hysteresis, region pinning, prefetch hints
  - [ ] Present‑LOD fallback under pressure

### 17. CPU backend parity (P1 • FEAT)
- Owner: Core/CPU
- Estimate: 3 weeks
- Dependencies: 3,7
- Artifacts: OpenMP/SIMD kernels, parity tests
- Checklist:
  - [ ] AVX2/AVX‑512/SVE kernels for P2G/G2P/constitutive loops
  - [ ] NUMA pinning and memory interleave
  - [ ] Parity tests with GPU small scenes

### 18. Console backends hardening (P0 • PERF)
- Owner: Backend GPU
- Estimate: 3–4 weeks
- Dependencies: 11,16
- Artifacts: Console builds, perf profiles
- Checklist:
  - [ ] Threadgroup memory limits tuned per platform
  - [ ] Heavy scene timing ≤8 ms; capture reviews
  - [ ] Validation parity against PC baseline

### 19. Partner UE/Unity samples & UX (P1 • FEAT)
- Owner: Interop
- Estimate: 2 weeks
- Dependencies: 5,11
- Artifacts: Sample projects, onboarding guide
- Checklist:
  - [ ] UE and Unity sample parity scenarios
  - [ ] Onboarding doc with step‑by‑step setup
  - [ ] Partner feedback incorporation and sign‑off

### 20. Material fitter tool (P0 • FEAT)
- Owner: Tools/CI + Simulation
- Estimate: 3 weeks
- Dependencies: 7–14
- Artifacts: CLI/GUI fitter, scene runners
- Checklist:
  - [ ] Loss function across runout, rut depth, slip RMS
  - [ ] Batch runner and convergence criteria
  - [ ] Export/import of presets

### 21. Documentation completion (P1 • DOC)
- Owner: Docs + All
- Estimate: 1–2 weeks
- Dependencies: 1–20
- Artifacts: Finalized docs, API reference
- Checklist:
  - [ ] Algorithms/Architecture/Materials/Testing/Integration finalized
  - [ ] C‑ABI reference generated and linted
  - [ ] Release notes and migration guide

### 22. Packaging & manifest (P1 • TECH)
- Owner: Core/Build + Tools/CI
- Estimate: 1 week
- Dependencies: 1–21
- Artifacts: Headers, shared libs, JSON capabilities manifest
- Checklist:
  - [ ] Multi‑platform builds with symbol/version checks
  - [ ] Manifest enumerating features, limits, shader variants
  - [ ] Publish pipeline to artifacts storage

### 23. Input validation & hard caps (P0 • SEC)
- Owner: Core
- Estimate: 1 week (ongoing)
- Dependencies: 1–22
- Artifacts: Validation layer, fuzz corpus
- Checklist:
  - [ ] Validate counts/strides; reject oversize inputs
  - [ ] No partial writes on failure; transactional patterns
  - [ ] Fuzzers for ABI calls and buffer parameters

### 24. CI gates and perf budgets (P0 • TEST)
- Owner: Tools/CI
- Estimate: 1 week (ongoing)
- Dependencies: 6
- Artifacts: CI pipelines, perf dashboards
- Checklist:
  - [ ] Invariants and perf gates with thresholds
  - [ ] Nightly soaks and drift audits
  - [ ] Auto‑capture on regression with attachments

### 25. Telemetry depth (P1 • OBS)
- Owner: Tools/CI + Core
- Estimate: 1 week
- Dependencies: 6
- Artifacts: Extended counters, error triage
- Checklist:
  - [ ] Counters for atomics, compactions, dt clamps, substeps
-  [ ] Error triage with frame context and recent telemetry window
  - [ ] Export to JSON/CSV and visual dashboards

