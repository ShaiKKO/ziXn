# Backends

This folder will contain platform backends:

- `d3d12/`
- `vk/`
- `metal/`
- `cpu/`

Each backend adheres to the narrow C-ABI surfaces defined in `core/` and must not depend on engine internals.
