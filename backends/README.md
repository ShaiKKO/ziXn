/*!
\file README.md
\brief Backend directory overview (D3D12/Vulkan/Metal/CPU).
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.
*/

This folder will contain platform backends:
- `d3d12/`
- `vk/`
- `metal/`
- `cpu/`

Each backend adheres to the narrow C-ABI surfaces defined in `core/` and must not depend on engine internals.

