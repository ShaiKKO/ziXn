/*!
\file zx_tiles.h
\brief Tile-sparse grid and particle SoA declarations (scaffolding).
\author Colin Macritchie (Ripple Group, LLC)
\version 1.0.0
\date 2025-09-16
\license Proprietary — Copyright (c) 2025 Colin Macritchie / Ripple Group, LLC.

Units:
- Tile dimension B (nodes per side) is given by ZX_TILE_B.
- Node momentum mom_* uses unit density proxy; velocity = mom/mass.
*/

#ifndef ZX_TILES_H
#define ZX_TILES_H

#include <stdint.h>

enum
{
  ZX_TILE_B = 16
};
/* nodes per dimension per tile */  // NOLINT(cppcoreguidelines-use-enum-class,performance-enum-size)

typedef struct zx_particle_soa
{
  /* Pointers to arrays owned by core; lengths tracked separately */
  float* pos_x;
  float* pos_y;
  float* pos_z;
  float* vel_x;
  float* vel_y;
  float* vel_z;
  float* mass;
  float* volume;
  float* F; /* 9*N */
  float* C; /* 9*N */
  uint16_t* mat_id;
  uint16_t* flags;
} zx_particle_soa;

typedef struct zx_tile_node
{
  float mass;
  float mom_x, mom_y, mom_z;
} zx_tile_node;

typedef struct zx_tile
{
  int32_t coord_x, coord_y, coord_z;
  uint32_t generation;
  zx_tile_node
      nodes[ZX_TILE_B * ZX_TILE_B * ZX_TILE_B];  // NOLINT(cppcoreguidelines-avoid-c-arrays,
                                                 // modernize-avoid-c-arrays)
} zx_tile;

#endif /* ZX_TILES_H */
