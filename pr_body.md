Implements present-LOD fallback (hysteresis + blend) with telemetry and CLI wiring.

Highlights
- C-ABI: policy/state in core/include/zx/zx_lod.h; impl in core/src/zx_lod.cpp
- CLI: --fallback on|off|auto, --fallback-thresholds; deterministic path suppresses timers
- Telemetry: fallback_activations, fallback_active_frames, fallback_blend_frames
- CI: Windows auto-format, bash format check, ctest logs, BugBot context

Risks
- Noisy timing signals -> EWMA/median in auto; determinism unaffected when fallback off.
