
**Key design choices**
- **Out-of-process native C++** gives full control of threading, allocators, and portability.  
- **Lua bridges** are minimal, versionable, and disposable.  
- **VLUA → GELUA**: VLUA collects **high-rate** wheel+contact signals; GELUA **coalesces and ships** them at a lower rate to the adapter; adapter steps ziXn at a **fixed rate (e.g., 60–120 Hz)**.  
- **Backpressure**: GELUA throttles upload if adapter cannot keep up; VLUA down‑samples at source.

---

## 3. Update order & rates

- **VLUA (per vehicle)**: `onPhysicsStep(dt_phys)` @ up to ~2000 Hz → sample wheels (normal force, slip ratio/angle, ground normal, contact pos), vehicle pose/vel, and (optional) per‑wheel node id. Ship to mailbox (`be:sendToMailbox`) or queue to GELUA (`obj:queueGameEngineLua(...)`).  
- **GELUA (global)**: `onUpdate(dt_gfx)` @ framerate → drains per‑vehicle mailboxes, merges sensor data into a **compact per‑tile contact field**, sends a **single message** to the adapter containing: vehicle summaries, contact patches, terrain tile IDs, and frame timestamps.  
- **Adapter (native)**: receives payload, updates ziXn scene/rigs, calls `simulate(dt_fixed)`, returns **(a)** height/mask deltas (visual), **(b)** per‑wheel feedback (additional sinkage/drag forces or friction scalars), **(c)** material state (wetness/porosity).  
- **GELUA writeback**: applies **visual** updates (terrain patch edits or decal/mesh overlays) and **physics** feedback routing commands to each VLUA. VLUA applies forces (`obj:applyForce`) or sets wheel friction scalars via available controllers.

---

## 4. BeamNG integration points (concrete)

### 4.1 Virtual Machines & hooks
- **Files:**  
  - `lua/ge/extensions/zixn_bridge.lua` (loaded via `extensions.load('zixn_bridge')`)  
  - `lua/vehicle/extensions/zixn_probe.lua` (added to vehicles through BeamNGpy or by `extensions.load` per-vehicle)
- **Hooks:**  
  - `M.onExtensionLoaded()` / `M.onUpdate(dt)` in GELUA.  
  - `onPhysicsStep(dt)` / `onGraphicsStep(dt)` in VLUA.  
- **Inter-VM calls** (official examples):  
  - Enqueue into VLUA: `be:getPlayerVehicle(0):queueLuaCommand('myLua()')` (from GELUA).  
  - Enqueue into GELUA: `obj:queueGameEngineLua('myLua()')` (from VLUA).  
  - Mailboxes: `be:sendToMailbox(addr, data)` (GELUA) and `obj:getLastMailbox(addr)` (VLUA).

### 4.2 Remote control / IPC (choose one; both supported)
- **TCP (preferred)** — implement the **BeamNGpy Control API** protocol subset:  
  - Connect to **port 25252** (default since v1.31; legacy **64256** supported).  
  - Use `control.queue_lua_command('...')` semantics to issue Lua; register a **zixn channel** in Lua that streams/receives JSON or binary blobs (see §4.3).  
- **UDP (custom)** — implement a **custom Protocol**:  
  - Place `USER_FOLDER/mods/unpacked/ziXnProto/lua/vehicle/protocols/ziXnProto.lua`.  
  - Use a unique header; encode per‑wheel telemetry and receive feedback packets.

### 4.3 GELUA bridge responsibilities
- Lifecycle: `onExtensionLoaded()` starts TCP/UDP client, loads `zixn_probe` for existing vehicles, subscribes to spawn/despawn.  
- Per‑frame: read mailboxes; build **contact tiles** `{tileId, avg_n, avg_p, slip, load, materialTag}`; send to adapter.  
- Visual writeback:  
  - **Option A (terrain edits):** call terrain editor Lua to modify **heightmap** in bounded patches (brush operations that set/raise/smooth), update material masks (e.g., **mud track** layer).  
  - **Option B (runtime mesh/decal):** render a dynamic displacement mesh or decals under the vehicle; cheaper and reversible.  
- Diagnostics: draw overlays (contact normals, rut depth, tile occupancy), and per‑frame timings.

### 4.4 VLUA probe responsibilities
- Sample: wheel contact point, contact normal, normal force/load, slip ratio/angle, wheel angular speed, vehicle linear speed; optionally **node ids** under wheels.  
- Ship: pack into a compact array per tick; on interval (e.g., every N physics ticks or on change), push to GELUA via **mailbox**, or use `obj:queueGameEngineLua()` to append to a GELUA aggregator.  
- Feedback: accept per‑wheel commands from GELUA to apply `obj:applyForce(Fx,Fy,Fz,nodeId)` or adjust available friction scalars (where supported by controllers).

> **Notes:** exact wheel data fields vary per vehicle; leverage existing **Vehicle Controller** and **wheels** data structures. Start with load, slip, angular/linear speeds, and contact position/normal; extend as needed.

---

## 5. Terrain writeback strategies

### 5.1 Visual fidelity targets
- **Near-field**: ruts, berms/heave, splash/sling trails, footprint detail.  
- **Mid-field**: smoothed displacement and wetness/dirt masks.  
- **Far-field**: static decals only (no edits).

### 5.2 Terrain physics coupling (test harness)
- **Primary (robust):** leave BeamNG terrain collision untouched; **inject forces** to vehicles to simulate sinkage/rolling resistance and traction loss. Guarantees stability across maps/mods.  
- **Secondary (optional):** for small patches, **edit heightmap** to form shallow ruts (visual+collision). Gate behind safety budget (tile size, brush radius, write frequency).

### 5.3 Programmatic terrain edits (GELUA)
- Use the **World Editor’s terrain tools** programmatically from GELUA to modify heightmap tiles (set/raise/smooth) and **paint material layers** (groundmodels map physics). This mirrors what users do interactively, but driven by our bridge in bounded areas around contact tiles.  
- **Budgeting:** patch size ≤ 32×32 cells; amortize to ≤ 1–2 updates/frame; throttle under low FPS; undo/restore buffers per tile.

> **Caveat:** Run-time terrain edits are expensive; prefer visual decals/meshes for heavy interactions and restrict height edits to slow, dramatic effects (e.g., mud bogging berms).

---

## 6. Data contracts

### 6.1 From BeamNG → adapter
```text
FrameHeader { sim_time, frame_idx, vehicles[n] }
VehicleSummary {
  vid, pose (x,y,z,q), v_lin (x,y,z), v_ang (x,y,z),
  wheels[m] { id, pos, normal, load_N, slip_ratio, slip_angle, omega, radius }
}
ContactTile {
  tile_id, center, half_extent, material_tag, avg_normal, avg_load_N,
  avg_slip_ratio, avg_slip_angle, moisture_hint, temperature_hint
}
6.2 From adapter → BeamNG
text
Always show details

Copy code
Writeback {
  tiles[k] { tile_id, height_delta_grid[r×r], material_mask_delta[r×r] },
  per_wheel_feedback[m] { vid, wheel_id, sinkage_m, Fx,Fy,Fz, friction_scale }
}
Binary encoding: flatbuffers/Cap’n Proto or length‑prefixed JSON; gzip for large tiles. All timestamps monotonic to debug latency.

7. Scheduling, determinism & performance
Simulation rate: fixed 60–120 Hz in adapter; substep internally as ziXn requires.

GELUA rate: at graphics rate; coalesce VLUA samples to last Δt window.

VLUA rate: sample at 200–500 Hz (down from 2000 Hz) to bound CPU.

Determinism: queue ordering and reductions serialized in GELUA; quantize tile outputs; replay files stored in adapter.

Performance gates: adapter budgets (CPU/GPU) + GELUA/VLUA budgets (≤ 1 ms GELUA, ≤ 0.1 ms VLUA/veh).

8. Security & failure handling
Sandboxing: ziXn runs out‑of‑process; a crash cannot take down BeamNG.

Graceful degradation: if adapter disconnects, GELUA stops edits and disables feedback; a hotkey toggles overlay and connection.

Debugging: start BeamNG with -luadebug and attach VS Code Lua Debugger; adapter logs human‑readable telemetry; capture JSON packets for offline replay.

9. Build & deployment
Adapter: single C++17 exe; Windows x64; statically link runtime; config via zixn.toml.

Lua bridges: ship as a mod ZIP (automatic install to USER_FOLDER/mods/).

Startup: via BeamNGpy‑style launcher or command line (-lua to autoload our bridges).

10. Milestones & exit criteria
M0 — Bootstrap (1–2 wks)

GELUA bridge loads/unloads; VLUA probe attaches to vehicles; TCP loopback; UI toggles.

Exit: round‑trip ping/pong at ≥ 200 Hz; no hitches.

M1 — Telemetry (1–2 wks)

Wheels/contact stream; per‑tile coalescing; perf counters.

Exit: correct contact heatmaps vs. debug rays; ≤ 1 ms GELUA cost for 2 vehicles.

M2 — ziXn step & feedback (2–3 wks)

Adapter drives ziXn; returns per‑wheel sinkage/drag; VLUA applies forces.

Exit: measurable slowdown in sand/mud; controlled slip; stable frames.

M3 — Visual writeback (2–3 wks)

Option A decals/mesh; Option B bounded terrain edits; groundmodel mask painting.

Exit: visible ruts/berms with bounded perf cost; undo works.

M4 — Scale & soak (2–4 wks)

4–8 vehicles; long‑run stability; replay capture+replay.

Exit: 30‑minute soak, no leaks; consistent timings.

11. Risks & mitigations
In‑process C++ plugins are not a public mod surface → use external C++ + Lua bridges.

Runtime terrain edits are heavy → minimize area/frequency; prefer decals/meshes; inject physics via forces.

API evolution (BeamNG updates) → keep bridges minimal, behind feature flags; CI against latest stable.

Vehicle diversity → probe only stable fields (load, slip, contact pos/normal); add per‑vehicle fallbacks.

12. Documentation & API pointers (curated)
Virtual Machines (VLUA/GELUA/UI), rates, and inter‑VM APIs — queues/mailboxes and example calls (obj:queueGameEngineLua, be:queueAllObjectLua, be:sendToMailbox, obj:getLastMailbox).

Virtual Machines (update rates, 2000 Hz physics, comms examples)

Extensions (how to author M.onUpdate, load/unload, dependency graph)

Remote control

BeamNGpy Control API — queue_lua_command, port default 25252 (v1.31+), backwards 64256; per‑vehicle extensions loading parameter.

Arguments & Settings — -lua, -luastdin to inject Lua; -luadebug port 21110.

BeamNG.tech support — core C/C++ with Lua exposure (confirms extension surface).

Protocols (UDP) — OutGauge, MotionSim, and custom protocols via lua/vehicle/protocols/*.lua (since v0.32), shipped under permissive bCDDL; packaging & hot‑reload.

Terrain editing & groundmodels — Terrain Editor tools (raise/smooth/set height), material‑to‑groundmodel mapping and file locations.

For deep vehicle data structures and function names in VLUA/GE Lua, see the community‑maintained List of BeamNG Functions and Fields (tree of modules/fields) as a practical complement to the official docs.

13. Minimal “shape” of the Lua bridges (non‑code)
GELUA zixn_bridge

Exposes: onExtensionLoaded, onUpdate, onVehicleSpawned, onVehicleDestroyed

Starts TCP/UDP; keeps vehState[vid] = { lastMailbox, lastPose, … }

Aggregates contact tiles; visual writeback via terrain functions or decals

VLUA zixn_probe

Hooks onPhysicsStep; reads wheel load/slip/pos/normal; obj:getVelocity(); obj:applyForce() for feedback

Sends snapshots to GELUA mailbox at interval

14. Acceptance tests (adapter ready)
Telemetry sanity: contact heatmap aligns with tire tracks in dust; overlay dot cloud follows wheels.

Feedback effect: same vehicle, same throttle on sand vs asphalt shows higher slip and lower peak speed when adapter feedback is enabled.

Visual tracks: shallow ruts appear under slow bogging; undo on reset.

Performance: ≤ 1 ms GELUA, ≤ 0.1 ms VLUA/veh; end‑to‑end latency ≤ 2 frames.

Resilience: disconnect/reconnect; reload VLUA/GELUA (Ctrl‑R/Ctrl‑L) without crash.

Appendix A — Tick diagram
ruby
Always show details

Copy code
VLUA:onPhysicsStep ──► mailbox += contact sample (N Hz)
          ▲                         │
          │                         ▼
    per-vehicle cmds       GELUA:onUpdate (fps)
          │                   drain mailboxes → compact contact tiles
          │                   TCP/UDP send → adapter
          │                         │
          │                         ▼
          ◄──────────── adapter: ziXn.simulate(dt_fixed) ── writeback (tiles, feedback)
Appendix B — BeamNG terminology crosswalk
VLUA/GELUA/UI — Vehicle, GameEngine, and UI virtual machines (Lua sandboxes).

Extensions — Lua modules returning a table (M) with special hooks (onUpdate, etc.).

Mailboxes/Queues — official inter‑VM data and code pathways.

Protocols — Lua modules emitting/consuming UDP telemetry.

BeamNGpy protocol — TCP control channel used by the official Python client (queue_lua_command, etc.).

References
Virtual Machines (update rates, comms, examples) — https://documentation.beamng.com/modding/programming/virtualmachines/

Extensions (hooks, load/unload) — https://documentation.beamng.com/modding/programming/extensions/

BeamNGpy Control API & queue_lua_command (TCP) — https://documentation.beamng.com/api/beamngpy/ (see ControlApi, vehicle extensions parameter), plus changelog noting default port 25252 since v1.31

BeamNG.tech Support — core C/C++ with Lua exposure — https://beamng.tech/support/

Protocols (OutGauge, MotionSim, custom UDP via lua/vehicle/protocols) — https://documentation.beamng.com/modding/protocols/

Arguments & Settings (-lua, -luastdin, -luadebug) — https://documentation.beamng.com/beamng_tech/arguments_and_settings/

Terrain Editor & Groundmodels mapping — https://documentation.beamng.com/world_editor/tools/terrain_editor/

BeamNGpy docs/changelog (port change) — https://beamngpy.readthedocs.io / https://documentation.beamng.com/api/beamngpy/*/changelog.html

Community function tree (practical discovery) — https://angelo234.github.io/List-of-BeamNG-Functions-And-Fields/