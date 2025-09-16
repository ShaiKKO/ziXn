# ziXn — Materials & Calibration
**Exhaustive practical parameter ranges and fitting workflows for granular, snow, and slurry terrain.**

> This document translates physics into **authorable presets** and **calibration workflows**. Ranges reflect typical game‑facing behavior rather than lab exactness. Use them as **starting points**; finalize with validation scenes (column collapse, slope runout, wheel trenching, footprints).

---

## Table of Contents
1. [Authoring Model](#authoring-model)
2. [Common Symbols & Units](#common-symbols--units)
3. [Global Controls](#global-controls)
4. [Material Families & Presets](#material-families--presets)
   - [Dry Sands](#dry-sands)
   - [Damp/Wet Sands](#dampwet-sands)
   - [Gravel & Crushed Rock](#gravel--crushed-rock)
   - [Soils: Silt, Loam, Peat](#soils-silt-loam-peat)
   - [Clays](#clays)
   - [Mud (Slurry, HBP)](#mud-slurry-hbp)
   - [Snow (Powder → Packed → Slush)](#snow-powder--packed--slush)
   - [Mixed/Layered Beds](#mixedlayered-beds)
5. [Two‑Phase Parameters](#twophase-parameters)
6. [Contact & Interaction](#contact--interaction)
7. [Calibration Workflows](#calibration-workflows)
8. [Validation Scenes & Metrics](#validation-scenes--metrics)
9. [Authoring UI & LOD Mappings](#authoring-ui--lod-mappings)
10. [Troubleshooting & Smell Tests](#troubleshooting--smell-tests)

---

## Authoring Model
Each **material preset** ultimately maps to solver parameters:
- **Frictional plasticity** (MC/DP): friction angle `φ`, cohesion `c`, dilatancy `ψ_d`, elastic `E, ν` (or `λ, μ`).
- **Snow cap**: cap size/evolution controlling compaction.  
- **Viscoplastic slurry (HBP)**: `τ_y, K, n, μ0, m`.
- **Two‑phase**: permeability `k`, porosity `n`, saturation `S`, Biot `α`.
- **Density**: bulk or solids density as appropriate.  
- **Presentation coupling**: displacement weight, compression/mask ranges.

We provide **author‑facing sliders** that are linearized and clamped to safe simulation bounds.

---

## Common Symbols & Units
- Density `ρ` [kg·m⁻³] (bulk unless specified), Solids density `ρ_s`.
- Friction angle `φ` [deg], Cohesion `c` [kPa].
- Elastic: Young’s `E` [MPa], Poisson `ν` [–].
- HBP: yield stress `τ_y` [Pa], consistency `K` [Pa·sⁿ], exponent `n` [–], base viscosity `μ0` [Pa·s], regularization `m` [s].
- Two‑phase: permeability `k` [m²] or hydraulic conductivity `K_h` [m·s⁻¹], porosity `n` [–], saturation `S` [–].
- Moisture content `w` [%], Temperature proxy `η` [0–1] for snow strength.

**Conversions**: if author prefers **hydraulic conductivity** `K_h`, approximate **permeability** by `k ≈ K_h μ/ (ρ g)` with `μ≈1.0e‑3 Pa·s` (water at ~20 °C).

---

## Global Controls
- **Grid spacing h**: 5–12 cm typical for gameplay ruts/footprints; smaller h improves detail but increases cost ∝ h⁻³.
- **Particles per cell**: 4–16; fewer for dry sand (APIC preserves detail), more for snow/slurry.
- **Cap on plastic increment** and **μ_eff clamp** (HBP) to ensure stability at target `Δt`.

---

## Material Families & Presets

### Dry Sands
**Use cases**: dunes, beaches (dry), desert roads, sandbox piles.  
**Model**: MC or DP (no cohesion), optional small dilatancy.

| Parameter | Typical Range | Notes |
|---|---:|---|
| Bulk density `ρ` | 1400–1700 | Depends on packing; solids ~2650 |
| Friction angle `φ` | 30–38° | Dune sands near 32–34° |
| Cohesion `c` | 0–2 kPa | ≈0 for very dry |
| Dilatancy `ψ_d` | 2–10° | Heave around tires/feet |
| Young’s `E` | 2–20 MPa | For elastic trial; tune feel |
| Poisson `ν` | 0.2–0.35 |  |
| Permeability `k` | 1e‑11–1e‑9 m² | Very permeable (high) |
| Porosity `n` | 0.35–0.45 |  |

**Authoring presets**: _Fine Sand_, _Coarse Sand_, _Dune Crest_.

### Damp/Wet Sands
**Use cases**: beach near waterline, sand castles, wet footprints.  
**Model**: MC/DP with small **cohesion**; if soupy, switch to **HBP**.

| Parameter | Typical Range | Notes |
|---|---:|---|
| `ρ` | 1600–1900 | Heavier due to water |
| `φ` | 28–35° | Friction reduced with water |
| `c` | 0.5–5 kPa | Capillary cohesion |
| `ψ_d` | 0–6° | Lower than dry |
| HBP `τ_y` | 10–80 Pa | For soupy patches |
| HBP `K, n` | 0.5–5 Pa·sⁿ, 0.5–1.0 |  |

### Gravel & Crushed Rock
**Use cases**: rally stages, construction yards.  
**Model**: MC with large φ, no cohesion; stiffer elastic trial.

| Parameter | Typical Range | Notes |
|---|---:|---|
| `ρ` | 1700–2200 |  |
| `φ` | 38–48° | Angular particles |
| `c` | 0–2 kPa |  |
| `ψ_d` | 6–14° | Strong heave |
| `E` | 10–80 MPa | Stiffer |
| `k` | 1e‑12–1e‑10 m² |  |

### Soils: Silt, Loam, Peat
**Use cases**: soft forest beds, agricultural fields.  
**Model**: DP with small cohesion; consider two‑phase for bogging (peat).

| Parameter | Typical Range | Notes |
|---|---:|---|
| `ρ` | 1200–1700 | Lower for peat |
| `φ` | 20–32° |  |
| `c` | 1–15 kPa | Organic content ↑ cohesion |
| `ψ_d` | 0–4° |  |
| `k` | 1e‑14–1e‑11 m² | Lower permeability |
| `n` | 0.45–0.85 | High porosity for peat |

### Clays
**Use cases**: sticky mud ruts, slide‑prone cuts.  
**Model**: **Undrained** behavior (two‑phase) at short times; **drained** DP+cap otherwise.

| Parameter | Typical Range | Notes |
|---|---:|---|
| `ρ` | 1500–1900 |  |
| `φ` | 10–25° | Low friction |
| `c` | 15–100+ kPa | Wide range |
| `E` | 5–50 MPa |  |
| `k` | 1e‑17–1e‑14 m² | Very low |
| HBP `τ_y` | 50–300 Pa | For muddy behavior |
| HBP `n` | 0.3–0.8 | Shear‑thinning |

### Mud (Slurry, HBP)
**Use cases**: bogs, puddles, wheel churning.  
**Model**: **HBP**; optionally coupled with two‑phase for pore memory.

| Parameter | Typical Range | Notes |
|---|---:|---|
| Base viscosity `μ0` | 1e‑3–1e‑2 Pa·s | Water → dirty water |
| Yield `τ_y` | 20–300 Pa | Stickiness |
| Consistency `K` | 0.2–10 Pa·sⁿ | |
| Exponent `n` | 0.3–1.0 | <1 shear‑thinning |
| Regularization `m` | 50–500 s |  |
| Density `ρ` | 1100–1600 |  |

### Snow (Powder → Packed → Slush)
**Use cases**: ski slopes, tundra, roadsides.  
**Model**: DP+cap with evolving strength `η` (temperature/packing proxy). Switch to HBP‑like slush near melt.

| Parameter | Powder | Packed | Wet/Slush |
|---|---:|---:|---:|
| `ρ` | 80–200 | 200–400 | 400–700 |
| `φ` | 10–25° | 20–30° | 15–25° |
| `c` | 0.5–3 kPa | 2–6 kPa | 1–5 kPa |
| `E` | 0.5–5 MPa | 2–15 MPa | 0.5–4 MPa |
| HBP `τ_y` | – | – | 10–80 Pa |
| Notes | Frangible | Crisp ruts | Slushy smear |

### Mixed/Layered Beds
Stack layers (e.g., sand over clay). Provide per‑layer parameters + thickness; ziXn blends via **tile fractions** and enforces contact based on the **top layer**; infiltration changes local saturation and shifts model LOD.

---

## Two‑Phase Parameters
- **Porosity `n`**: sands 0.35–0.45; clays 0.3–0.6; peat up to 0.85.  
- **Permeability `k`**: sands `1e‑11–1e‑9 m²`, silts `1e‑14–1e‑12`, clays `1e‑17–1e‑14`.  
- **Biot α**: 0.6–1.0 (use 0.9 default).  
- **Van Genuchten**: `α_vg = 0.5–10 m⁻¹`, `n_vg = 1.2–2.0` for unsaturated curves; Mualem `l ≈ 0.5`.

**Gameplay knobs**: “Bog memory” scales `β ∝ μ/k` and pressure relaxation; high memory → vehicles remain stuck longer after load removal.

---

## Contact & Interaction
- **Friction coefficients** for tire/boot vs material (tangential cone): sand 0.5–0.8, gravel 0.7–1.0, mud 0.3–0.6, snow 0.2–0.5.  
- **Anisotropy** for tracks: ellipse ratio 1.0–3.0 (longitudinal vs lateral).  
- **Compliance**: small normal compliance `κ_n` prevents chatter; 1e‑4–1e‑3 s²·m⁻¹ recommended.

---

## Calibration Workflows

### 1) Column Collapse (Dry Granular)
- Build column aspect ratio `H/L`; remove gate; measure **runout** and final slope angle.
- Fit **φ** to match runout; adjust **ψ_d** for heave around base.

### 2) Inclined Plane (Friction & Cohesion)
- Increase slope until steady sliding; fit **φ** (dry) or **φ, c** (damp).

### 3) Wheel Trenching & Slip Curves
- Drive at constant torque; record rut depth, berm volume, and slip ratio.
- Tune **ψ_d** (heave), **c** (rut crispness), **τ_y** (mud stickiness), then **k** (pore memory).

### 4) Footprints Under Cyclic Loads
- Apply periodic load; evaluate permanent set and rebound. Calibrate **cap** parameters and **E**.

### 5) Dam‑Break Over Porous Bed (Slurry/Infiltration)
- Fit **HBP** and **k** together; ensure pressure residuals are controlled (PCG iters).

**Optimization Loop**: Use ziXn’s material fitter to minimize multi‑metric loss (runout error + rut depth + slip curve RMS).

---

## Validation Scenes & Metrics
- **Mass conservation** (<0.2% per 60 s).  
- **Runout error** (<10%).  
- **Ridge/berm volume** error (<15%).  
- **PCG residual** (<1e‑3 in 10–16 iters on coarse grid).  
- **Determinism checksum** identical across target backends within tolerance.

---

## Authoring UI & LOD Mappings
- Sliders: _Dryness_ ↔ adjusts `c` and `k`; _Grain size_ ↔ `φ` and dilatancy; _Stickiness_ ↔ `τ_y`.
- **LOD rules**: as `k` drops and `S` rises, switch to two‑phase (coarse pressure). In far‑field, reduce particles per cell and clamp `μ_eff` to widen stability.

---

## Troubleshooting & Smell Tests
- **Sand flows too far** → raise `φ`, reduce `ψ_d`.  
- **Mud looks watery** → increase `τ_y` or `K`, raise `μ0`, check `m` regularization.  
- **Vehicle never gets unstuck** → reduce `β` (raise `k`) or shorten pressure memory.  
- **Snow lacks crisp edges** → increase cohesion `c` and cap stiffness; ensure APIC not over‑clamped.  
- **Solver unstable at strong stickiness** → increase `μ_eff` clamp, reduce Δt or enable substeps; verify contact compliance.
