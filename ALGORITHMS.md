# ziXn — Algorithms
**Numerical models for granular, snow, and slurry terrain suitable for real‑time execution**

> This document specifies the numerical schemes implemented in ziXn: particle–grid hybrids (APIC/MLS‑MPM), constitutive updates (Mohr–Coulomb, Drucker–Prager + cap, NorSand), snow densification, viscoplastic slurries (Herschel–Bulkley with Papanastasiou regularization), two‑phase mixture coupling (soil–water), contact projection, time integration, and stability/LOD heuristics.

---

## Table of Contents
1. [Discretization](#discretization)
2. [State Carried by Particles and Grid](#state-carried-by-particles-and-grid)
3. [P2G/G2P Transfers](#p2gg2p-transfers)
4. [Constitutive Models](#constitutive-models)
   - [Elastic Predictor](#elastic-predictor)
   - [Mohr–Coulomb (MC)](#mohrcoulomb-mc)
   - [Drucker–Prager (DP) and DP+Cap](#druckerprager-dp-and-dpcap)
   - [NorSand (Critical‑State Sand)](#norsand-criticalstate-sand)
   - [Snow (Elasto‑Plastic with Densification)](#snow-elastoplastic-with-densification)
5. [Viscoplastic Slurries: HBP](#viscoplastic-slurries-hbp)
6. [Two‑Phase Mixture (Soil–Water)](#twophase-mixture-soilwater)
7. [Contact: Frictional, Anisotropic, and Convex Projection](#contact-frictional-anisotropic-and-convex-projection)
8. [Time Integration & Substepping](#time-integration--substepping)
9. [Stability, CFL, and Regularization](#stability-cfl-and-regularization)
10. [LOD Strategy](#lod-strategy)
11. [Reference Pipeline (Per Frame)](#reference-pipeline-per-frame)
12. [Numerical Diagnostics](#numerical-diagnostics)
13. [Parameter Sensitivities](#parameter-sensitivities)

---

## Discretization
We use **particle–grid hybrids** (MPM family) to handle large deformation and history‑dependent materials. Particles carry mass and history; a background grid supports local momentum updates and contact. We prefer **APIC** transfers with an **MLS‑consistent** stress divergence for reduced dissipation and better angular momentum conservation.

- **Grid**: regular lattice with spacing `h` (per tile). Nodes store mass `m_i`, momentum `p_i = m_i v_i`, and auxiliaries (e.g., flags, contact normals, pore pressure terms).
- **Shape functions**: quadratic or cubic B‑splines. We default to quadratic for performance; cubic is available where aliasing demands smoother transfers.

## State Carried by Particles and Grid

**Particles (per material point p):**  
- Position `x_p`, velocity `v_p`, mass `m_p`, volume `V_p^0` (rest) and current `V_p`.  
- Deformation gradient `F_p` (3×3), affine velocity `C_p` (3×3) for APIC.  
- Plastic state (accumulated plastic strain `ε_p^pl`, density or packing `ρ_p`, state parameter `ψ` for NorSand, densification `η` for snow).  
- Material id `mat_p`.

**Grid nodes (per i):**  
- Mass `m_i`, velocity `v_i`, external force `f_i^ext`, internal force `f_i^int`, contact data (normal `n_i`, tangential basis), optional pore pressure `p_i` and drag data for two‑phase.

---

## P2G/G2P Transfers

### P2G (Particle → Grid)
For each particle `p`, for neighbors `i` with weight `w_{ip}` and gradient `∇w_{ip}`:
- **Mass**: `m_i += w_{ip} m_p`
- **Momentum**: `p_i += w_{ip} m_p (v_p + C_p (x_i - x_p))` (APIC affine term)
- **Internal force** (stress divergence): `f_i^int -= V_p σ_p ∇w_{ip}` where `σ_p` is Cauchy stress from constitutive update.
- Accumulate auxiliary quantities (for contact / pressure) tile‑locally; merge with deterministic order.

### Grid Update
- Update velocities with forces: `v_i ← (p_i + Δt (f_i^ext + f_i^int)) / m_i` for active nodes where `m_i>ε`.
- Apply **constraints**: contact projection (see below), boundary clamps, and optional pressure solve.

### G2P (Grid → Particle)
- Interpolate velocity and its gradient: `v_p ← Σ_i w_{ip} v_i`; `C_p ← 4 Σ_i w_{ip} v_i (x_i - x_p)^T / h^2` (one APIC form).
- Update position: `x_p ← x_p + Δt v_p`.
- Update deformation: `F_p ← (I + Δt ∇v_p) F_p` (with `∇v_p` from grid interpolation).
- Apply **plasticity** post‑update via return mapping (see models).

---

## Constitutive Models

### Elastic Predictor
- We use **Neo‑Hookean** or **corotated linear** elasticity for the **trial** stress.  
- Compute the elastic deformation `F_e` via polar decomposition: `F = R S`; take `F_e = R S` or a clamped variant to keep singular values within `[s_min,s_max]`.  
- Trial Kirchhoff stress `τ_trial = μ (F_e F_e^T - I) + λ ln(J_e) I` with `J_e = det(F_e)`; convert to Cauchy `σ_trial = τ_trial / J_e`.

### Mohr–Coulomb (MC)
- Yield function `f = τ_max + α p - c ≤ 0` with friction angle `φ` (through `α = sinφ / (1 - sinφ)`), cohesion `c`, and pressure `p = -tr(σ)/3`.
- **Return mapping**: project `σ_trial` onto MC cone; compute plastic multiplier `γ` using radial return in principal stress space; adjust deviatoric part `s` and pressure as needed; enforce **non‑negative plastic work**.
- **Dilatancy**: optional non‑associated flow with dilatancy angle `ψ_d` (`ψ_d ≤ φ`), controlling volumetric plastic strain.

### Drucker–Prager (DP) and DP+Cap
- Smooth yield `f = α I_1 + √J_2 - k ≤ 0` with `I_1 = tr(σ)` and `J_2 = 1/2 s:s`. Map MC parameters to DP (`α, k`) for stability.  
- **Cap** adds an elliptical surface in `p–ε_v^pl` to capture compaction; state evolves with plastic volumetric strain; suitable for snowpack compaction and wet sands.

### NorSand (Critical‑State Sand)
- Captures **dilatancy vs contraction** based on **state parameter** `ψ = e - e_cs(p)` (void ratio vs critical).  
- Yield and plastic potential with parameters `M` (critical stress ratio), `λ, κ` (compression lines), and reference void ratio curve `e_cs(p)`.  
- Implementation: track `e` via `ρ` (bulk density) and solids density; update `ψ` each step; perform return mapping with evolving `M(ψ)` if desired.

### Snow (Elasto‑Plastic with Densification)
- Elasto‑plastic with cohesion and evolving stiffness modulated by **density/temperature surrogate** `η ∈ [0,1]` (authoring proxy).  
- Use DP+cap with small cohesion and rate‑dependent flow to achieve powder → packed transitions; add **fracture‑like softening** via cap evolution when shear surpasses threshold.

---

## Viscoplastic Slurries: HBP
For mud‑like flows we use **Herschel–Bulkley** with **Papanastasiou regularization** to avoid `1/γ̇` singularities:
- Effective viscosity:
\[
\mu_\mathrm{eff} = \mu_0 + \frac{\tau_y}{\dot\gamma} \left(1 - e^{-m \dot\gamma}\right) + K \dot\gamma^{n-1},
\]
where `μ0` is base viscosity of the carrier fluid (e.g., water), `τ_y` yield stress, `K` consistency, `n` shear exponent, and `m` is a large regularization factor (e.g., 100–500 s).
- Implement on grid as a **viscous force**: `f_i^visc = ∇·(2 μ_eff D)` where `D = sym(∇v)`; treat semi‑implicitly using Jacobi/PCG on a **coarsened grid** to keep frame time bounded.
- Clamp `μ_eff` to `[μ_min, μ_max]` to avoid timestep collapse.

---

## Two‑Phase Mixture (Soil–Water)
To capture **pore pressure** and **relative motion** between grains and water (bogging, slush, quick‑sand behavior), we solve a reduced **mixture** model.

### Variables
- Solid phase (s) and fluid phase (f): velocities `v_s, v_f`, volume fractions `φ_s, φ_f = 1-φ_s`, densities `ρ_s, ρ_f`, pore pressure `p_f`.
- **Momentum balance (mixture form)**:
\[
ρ_s φ_s \dot v_s = ∇·σ' + ρ_s φ_s g - β (v_s - v_f),
\]
\[
ρ_f φ_f \dot v_f = -∇ p_f + ρ_f φ_f g + β (v_s - v_f) - viscous_f,
\]
where `σ' = σ - α p_f I` is **effective stress** (α≈φ_s), and `β ≈ μ / k` couples phases using **permeability** `k` and fluid viscosity `μ`.

### Pore Pressure
- **Darcy + mass balance** leads to a diffusion/Poisson‑like equation in `p_f`:
\[
∂(φ_f ρ_f)/∂t + ∇·( - \frac{k}{μ} ∇ p_f + φ_f v_s ) = sources.
\]
- Discretization: coarsened grid (e.g., 2× tile spacing), semi‑implicit solve via **PCG** with **AMG** or **multigrid smoothing**; boundary conditions from open boundaries or impermeable interfaces.

### Implementation Choice
- Use **double particle sets** (s,f) sharing the grid, or a **single particle set** with additional two‑phase fields; ziXn uses **single particle set** with fluid fields stored per tile to reduce bandwidth. Relative slip `v_rel = v_f - v_s` is reconstructed on the grid for drag.

---

## Contact: Frictional, Anisotropic, and Convex Projection
- Grid‑node **contact detection** uses signed distance fields (SDF) for rigs/terrain features. Nodes with penetration perform a **velocity projection**:
  - Normal component: clamp to non‑penetration with compliance `κ_n` for stability.
  - Tangential: Coulomb cone (coefficient `μ_t`), optionally **anisotropic** (ellipse) to represent tread/track patterns; projection solved via small per‑node QP or projected Gauss–Seidel with fixed iteration count.
- **Rolling resistance** proxy: scale tangential limit based on local slip ratio to emulate track grousers without meshing.

---

## Time Integration & Substepping
- Base integrator is **explicit/semi‑implicit** Euler on the grid with plasticity applied at particles.  
- Use **substeps** when any of the stability limits are exceeded (contact, viscous, pressure). Substeps scale `Δt_sub = Δt / N_sub` and repeat P2G–G2P.  
- Two‑phase pressure solve runs **once per frame** on a coarse grid; velocities are then corrected via gradient of `p_f`.

---

## Stability, CFL, and Regularization
Let `c_s` be elastic wave speed (from `λ, μ`), `c_d` characteristic drift speed, and `ν_eff` viscosity proxy.
- **Elastic CFL**: `Δt ≤ C1 h / (c_s + |v|)` with `C1≈0.4` for quadratic B‑splines.
- **Plasticity limiter**: cap plastic increment per step `||Δε_pl|| ≤ ε_max`.
- **Viscous limit** (slurry/HBP): `Δt ≤ C2 h^2 / ν_eff` with `C2≈0.25`.
- **Drag limit** (two‑phase): `Δt ≤ C3 / (β / ρ_mix)`.
- **Smoothing**: strain‑rate smoothing or **APIC** affine clamping reduces checkerboarding at low particle count.

---

## LOD Strategy
- **Model LOD**: dry granular ↔ HBP slurry ↔ two‑phase mixture; auto‑select by local moisture and gameplay need.  
- **Spatial LOD**: coarser `h` and fewer particles outside interaction frustum; keep activation hysteresis.  
- **Presentation LOD**: height/displacement filtering and decoupled visual foam/slush layers.

---

## Reference Pipeline (Per Frame)
1. Emit/retire particles; update materials and moisture sources/sinks.
2. Activate tiles; rebuild neighbor lists (lazy when topology stable).
3. **P2G**: scatter mass/momentum/affine and internal forces.
4. **Grid update**: external forces → constitutive update → contact projection.
5. If slurry/two‑phase: assemble pressure/viscous system and run **PCG** (coarse grid).
6. **G2P**: update velocities/positions and deformation; do plastic return mapping.
7. Fix‑up: boundary clamping, invariant checks, LOD transitions.
8. **Writeback**: terrain displacement/masks; optional particles to FX systems.
9. Publish counters and timings.

---

## Numerical Diagnostics
- **Mass drift** (per material and total).  
- **Energy budget**: ensure plastic work non‑negative; viscous dissipation matches expectations.  
- **Contact gap**: max penetration and number of nodes saturating friction cone.  
- **Pressure residual**: PCG iterations and final residual norm.  
- **Determinism checksum**: hash of particle headers after `G2P`.

---

## Parameter Sensitivities
- **Friction angle φ**: drives runout length; ±2° noticeably changes dune/rut stability.  
- **Cohesion c**: affects footprint crispness; too high produces unrealistic arches.  
- **Dilatancy ψ_d**: controls heave around tracks/tires; tune carefully with wheel tests.  
- **HBP τ_y, K, n**: yield stress sets stuck/free threshold; `n<1` (shear‑thinning) improves flow under stirring.  
- **Permeability k**: smaller k → longer pore pressure memory (more “bogging”).  
- **Cap parameters**: too aggressive → volume locking; monitor volumetric strain.  
