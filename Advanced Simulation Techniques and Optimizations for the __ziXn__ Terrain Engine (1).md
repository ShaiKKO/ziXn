# Advanced Simulation Techniques and Optimizations for the **ziXn** Terrain Engine

This document outlines the **Phase 2 development tasks** for enhancing **ziXn** – focusing on multi-tile simulation, GPU/CPU performance, and engine integration – and provides in-depth research context for each area. The goal is to achieve real-time, high-fidelity terrain deformation with network-deterministic behavior and cross-platform support. We detail algorithmic approaches for each task and cite relevant techniques from literature to ensure the implementation is both robust and state-of-the-art.

## Phase 2 Development & Optimization Tasks for **ziXn**

### Multi-Tile Simulation Loop & Data Management

1. **Tile Activation & Deactivation Mechanism:** Implement logic to dynamically spawn or remove simulation tiles based on particle presence. A **TileSet** manager will track active tile coordinates in a spatial hash map[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=assemble%20the%20correct%20vertex%20data,be%20analyzed%20for%20both%20approaches). As particles move or new ones emit into empty regions, allocate new tiles from the pool (using `zxCreateTilePool` to get a free tile index[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=2%20PARALLELIZATION%20STRATEGIES%202,illustrates%20the%20communication%20steps%20for)) and mark them active. Conversely, if a tile’s support (mass or particle count) drops to zero, schedule it for deactivation – using hysteresis (e.g. a frame delay or threshold) to avoid thrashing tiles on/off each frame[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y). This mechanism touches the **Scene** object’s update loop and a tile management module (likely in `zx_scene.cpp`).
    
     
    
    **Algorithm Steps:**
    
    - For each active particle, compute its spatial tile index (e.g. floor of position / tile size). Mark that tile as **required** for this frame.
        
    - If a required tile is not currently active, allocate a tile from the pool and initialize it (grid nodes reset, neighbor links set)[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=2%20PARALLELIZATION%20STRATEGIES%202,illustrates%20the%20communication%20steps%20for). Add it to the active TileSet.
        
    - If an active tile has no particles (and negligible mass) for several frames, mark it as **inactive** and eligible for removal. Only actually free the tile after a cooldown period to prevent rapid toggle. (In adaptive schemes, “ghost” tiles that lose all neighbors are also removed to shrink the domain[wang-mengdi.github.io](https://wang-mengdi.github.io/proj/25-cirrus/cirrus-preprint.pdf#:~:text=%5BPDF%5D%20Cirrus%3A%20Adaptive%20Hybrid%20Particle,bottom%20line%29%20only).)
        
    - Update the scene’s spatial hash map of active tiles each frame, so lookups (tile by coordinate) are O(1). This map is used for quick existence checks when particles move.
        
2. **Neighbor List Construction & Ghost Region Updates:** Maintain connectivity info for each active tile[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=assemble%20the%20correct%20vertex%20data,be%20analyzed%20for%20both%20approaches). Each tile needs to know its 26 neighbors in 3D (including diagonal adjacencies) to facilitate cross-tile interactions. Implement a function to rebuild neighbor lists whenever tiles are activated or deactivated[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=assemble%20the%20correct%20vertex%20data,be%20analyzed%20for%20both%20approaches). Along with this, handle **ghost regions** at tile boundaries: either extend each tile’s grid to include an overlap layer, or explicitly exchange boundary data with neighbors each step. After **Particle-to-Grid (P2G)** transfers on a tile, share the edge node masses/momenta with adjacent tiles (e.g. copy values or perform summed reductions on the overlapping one-voxel border)[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=cell%20approach%2C%20where%20the%20summation,be%20analyzed%20for%20both%20approaches). This ensures continuity across tile boundaries in the multi-tile simulation loop. The approach is analogous to domain decomposition in parallel MPM: one method is to mirror “ghost particles” onto neighbor domains to assemble boundary data[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y), and another (chosen here) is to exchange summed node data at boundaries[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=assemble%20the%20correct%20vertex%20data,be%20analyzed%20for%20both%20approaches). The code for this will reside in the core simulation step, likely a tile boundary merge function called after each tile’s local grid update.
    
     
    
    **Algorithm Steps:**
    
    - When a tile is activated, also link or create any missing neighbor tiles if they become required (for example, if a particle is near a tile’s boundary, ensure the neighbor tile exists or is created so ghost nodes can be updated).
        
    - Represent neighbor connectivity in a structure (e.g. bitmask or array of 26 entries). Recompute this whenever a new tile is added or removed.
        
    - After each tile computes its grid quantities (mass, momentum, etc.), perform a **ghost exchange**: for each active tile, identify all neighbor tiles and copy the border node data to the neighbor’s ghost nodes. Overlapping node contributions are summed deterministically – e.g. if two tiles share a face, each computes that face’s nodes and we sum them once[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=cell%20approach%2C%20where%20the%20summation,be%20analyzed%20for%20both%20approaches).
        
    - Use synchronized or atomic operations for these ghost updates to avoid races (since two tiles will attempt to write the same ghost node). Alternatively, designate an order (like always add neighbor tile’s data into the tile with smaller index) to maintain determinism.
        
    - If a neighbor is inactive (missing), treat it as a solid boundary or simply ignore out-of-bounds ghost data.
        
3. **Particle Migration Between Tiles:** Ensure particles move seamlessly across tile boundaries. When updating a particle’s position (during **Grid-to-Particle, G2P**), detect if it has left its current tile’s bounds. If so, reassign it to the appropriate neighboring tile. This requires updating the particle’s tile index (or spatial hash) and possibly activating a new tile if it stepped into an inactive region[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y). Data structures like the particle SoA (`zx_particle_soa`) should support efficient removal and insertion of particles. Use spatial hashing or a grid-to-tile lookup to find the new tile index based on particle coordinates (since the world is sparsely tiled). The migration should be O(1) per particle on average.
    
     
    
    **Algorithm Steps:**
    
    - Calculate particle’s new tile coordinate = `(floor(new_position.x/B), floor(new_position.y/B), floor(new_position.z/B))` where B is tile length.
        
    - If this coordinate differs from the particle’s current tile, mark the particle for transfer. Remove it from the old tile’s list and add it to the new tile’s list.
        
    - If the new tile is not active, invoke the tile activation mechanism (Task 1) to allocate it[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=2%20PARALLELIZATION%20STRATEGIES%202,illustrates%20the%20communication%20steps%20for). Conversely, if the old tile becomes empty after removing this particle (and others), it may trigger deactivation (Task 1).
        
    - Handle edge cases: a particle might cross multiple tiles in one large time-step (though unlikely if time-step is small relative to tile size). If it does, ensure the algorithm accounts for that by possibly iterating migration until particle’s tile is correct.
        
    - This approach avoids needing “ghost particles”. (In distributed MPM, ghost particles are used to keep copies in adjacent domains[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y), but here we simply move the particle to the single authoritative tile that covers its new position.)
        
4. **Deterministic Tile Merge Reduction:** After each tile computes local grid contributions (mass, momentum, forces) in P2G, implement a **tile-merge** step to combine data on overlapping tile borders in a deterministic order[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=assemble%20the%20correct%20vertex%20data,be%20analyzed%20for%20both%20approaches). The goal is to sum contributions of shared nodes exactly once regardless of thread scheduling. One strategy is to assign a canonical ordering to tiles (e.g. by tile ID or coordinates) and always have the tile with the lowest ID perform the sum for a shared face or edge. Alternatively, gather all overlapping node values into a temporary array and sum them in a fixed sequence (ensuring floating-point sums occur in the same order every frame). This avoids race conditions and floating-point nondeterminism across tiles. On GPU, this can be done with atomic adds to a unified grid buffer, but to keep determinism, we can partition the work such that each overlap region is summed by a single thread or warp in a controlled manner. For example, Gao _et al._ (2018) avoided nondeterministic atomics by using warp-synchronous programming: threads in a warp first collectively accumulate values for the same cell using shuffle operations, then write once[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=the%20same%20location%20are%20conventionally,cell%20idpr%20ev%20then)[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=iteration%2C%20and%20the%20total%20number,adds%20to%20the). A similar technique can be applied at tile boundaries.
    
     
    
    **Algorithm Steps:**
    
    - Identify all duplicate grid nodes that belong to multiple tiles (e.g. a node on a shared boundary face belongs to 2 tiles; on an edge to 4 tiles; on a corner to 8 tiles). Create a list of these “overlap nodes” via neighbor links (Task 2).
        
    - For each overlap node, choose one tile (e.g. the one with lowest ID or a specific positional criterion) as the **owner** that will compute the merged value. The other tiles will not write that node to the final buffer.
        
    - During the P2G accumulation on GPU, use shared memory or warp-local reduction to sum contributions within each tile without atomics[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=the%20same%20location%20are%20conventionally,cell%20idpr%20ev%20then). Then, for overlap nodes, use an atomic add only when the owner thread writes to global memory, which happens exactly once per overlap[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=iteration%2C%20and%20the%20total%20number,adds%20to%20the). Since only the owner performs the global write, the result is order-independent and deterministic.
        
    - On CPU, simply perform the merges in a single-threaded loop after multithreaded P2G: iterate over each overlap node and sum values from each contributing tile in a fixed order (e.g. sorted list of contributing tiles). This ensures reproducibility run-to-run.
        
    - Verify determinism by running two simulations with different thread counts or GPU execution orders and computing checksums of grid state. They should match if this merge is correctly implemented[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=is%20widely%20believed%20that%20highly,from%20scattering%20without%20atomic%20op%02erations) (Gao’s approach shows that eliminating most atomic conflicts yields identical results to sequential summation).
        

### GPU Pipeline Integration

5. **P2G Compute Shader Dispatch:** Implement the **Particle-to-Grid (P2G)** transfer on the GPU. Each frame (or substep), dispatch a compute shader (HLSL) that processes particles for each active tile. A straightforward mapping is one thread per particle; however, to minimize divergence and memory contention, group threads by tile or by grid cell. For example, launch one thread group per tile, and within it use shared memory to accumulate particle data into that tile’s node array. Each thread can load one particle’s data (position, velocity, volume, etc.) from the SoA buffers (bound as SRV) and scatter its mass and momentum to the 3×3×3 neighboring nodes of that particle’s cell[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=1,the%20same%20used%20in%20FEM). Use **shared memory** to perform these accumulations locally: threads writing to the same node first add their contributions in shared memory, and only after that does one thread write the final sum to the global grid buffer (avoiding many global atomics). This technique is inspired by Gao _et al._ (2018), who sorted particles by cell to achieve coalesced access and used warp-level intrinsics to handle write conflicts without atomic operations[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=is%20widely%20believed%20that%20highly,from%20scattering%20without%20atomic%20op%02erations)[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=operations.%20We%20exploit%20the%20warp,5%3A%20boundary%20%E2%86%90%20t%20rue). Before dispatching P2G, ensure any newly activated tiles have their grid arrays zeroed for the new frame.
    
     
    
    **Algorithm Steps:**
    
    - **Buffer setup:** Bind the particle attribute buffers as read-only (SRV) and the grid node buffer (or per-tile buffers) as unordered-write (UAV). Ensure grid buffers are in UAV state. Zero out the grid buffers for new frame (either by compute shader or during tile activation).
        
    - **Dispatch geometry:** Launch `num_groups = active_tiles` (one group per tile) or another partition strategy. Each group knows which tile (by index) it handles, either via a constant or by reading a list of active tiles.
        
    - **Within shader (per group/tile):** Use group-shared memory for an array of B×B×B nodes (with padding for ghost). Each thread loads one particle from that tile (using tile’s particle range)[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=1,the%20same%20used%20in%20FEM), computes its base cell index (i,j,k within tile), and its weight contributions to the 8 (or 27 for quadratic) surrounding nodes. Accumulate its mass and momentum to the shared memory node array using atomic adds or, better, thread synchronization with local reduction if multiple threads target the same node.
        
    - After all particles in the tile group are processed, have one thread (or a few) write the shared memory node data to the global grid buffer (which might be a large array indexed by tile and node index). Thanks to grouping, global memory writes are minimized and mostly coalesced.
        
    - Insert a GPU memory barrier (UAV barrier) after P2G dispatch to ensure all writes to grid buffers are complete before the next phase[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=CUDA%20Driver%20checks%2C%20if%20the,page%20copy). This prepares for the grid update shader.
        
6. **Constitutive & Internal Forces Kernel:** Implement the GPU pass that updates each grid node’s state based on material constitutive laws (this corresponds to the **grid physics** phase: computing internal forces, applying elasticity/plasticity, and external forces like gravity). This kernel will read the intermediate grid (post-P2G mass and momentum) and compute forces and new velocities for each grid node[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=2,extrapolated%20back%20to%20material%20points). Key tasks within this kernel:
    
    - Compute each node’s velocity $\mathbf{v} = \mathbf{p}/m$ (momentum over mass) and acceleration $\mathbf{a} = \mathbf{F}_{\text{net}}/m$[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=2,extrapolated%20back%20to%20material%20points).
        
    - Compute stress divergence for solid materials: Each tile stores particle deformation or volume data; interpolate those to nodes or compute node-wise strain rates by comparing neighboring velocities[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=extrapolated%20from%20the%20surrounding%20nodes,depending%20on%20%2067%20integration). Then apply the constitutive model (e.g. linear elasticity or hyperelasticity with a yield criterion).
        
    - For plasticity: perform a **yield check** at each node. For example, use a Drucker–Prager yield surface for soils and snow, which relates shear stress and pressure. Compute trial stresses and if the yield criterion is violated (node’s stress state outside allowed domain), project the stress back onto the yield surfaceweb.cs.ucla.edu by reducing the deviatoric (shear) part and possibly volume change according to the flow rule. This is essentially a per-node return-mapping algorithm that ensures no node’s stress exceeds the material’s yield limit.
        
    - Apply the internal force (which is divergence of stress tensor) as well as gravity (a constant force) to update node momenta.
        
    
    The framework might separate these into multiple shaders for clarity: e.g. one to compute strain and trial stress, one to do yield correction. However, it can be done in one compute pass for efficiency. Use UAVs for writing updated node velocities. Each node’s update reads its current momentum, mass, and the neighboring nodes’ momentum (for spatial gradients).
    
     
    
    **Algorithm Steps:**
    
    - Launch a compute shader with one thread per grid node (e.g. grouped by tiles). Each thread loads the mass and momentum of one node (and perhaps an average deformation gradient or volume from nearby particles, if needed for stress).
        
    - Compute velocity $\mathbf{v}_i = \mathbf{p}_i / m_i$ for the node. If $m_i=0`, skip the node (no material).
        
    - Compute strain rate or deformation gradient update at the node by looking at neighboring node velocities[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=extrapolated%20from%20the%20surrounding%20nodes,depending%20on%20%2067%20integration). (If using APIC, an affine velocity field per particle might be used instead to compute stress directly from particle data; but on grid we might reconstruct an approximate strain.)
        
    - Calculate the trial Cauchy stress $\boldsymbol{\sigma}_\text{trial}$ using the material’s constitutive model (e.g. linear elastic: $\boldsymbol{\sigma} = \mathbf{C} : \boldsymbol{\epsilon}$). For fluids, the “stress” is pressure $p = -K \Delta V$ etc.
        
    - Check yield criteria: e.g. **Mohr–Coulomb** or **Drucker–Prager** for soil. Mohr–Coulomb defines yield when $\tau = c + \sigma_n \tan\phi$ (shear vs normal stress), whereas Drucker–Prager gives an approximation in invariant formweb.cs.ucla.edu. We use Drucker–Prager because it’s easier to implement in 3D (smooth cone) and it approximates Mohr–Coulomb’s frictional behaviorweb.cs.ucla.edu. Compute the yield function $F(\sigma)$ (for D-P, something like $F = \alpha, \text{tr}(\sigma) + \sqrt{J_2} - k$) and if $F>0$, scale back the deviatoric stress to satisfy $F=0$ (this is the plastic correction).
        
    - Apply the internal force: $\mathbf{f}_{\text{int}} = \nabla \cdot \boldsymbol{\sigma}$ for that node (discretized by differences with neighbor node stresses). Subtract it from node momentum (because it’s internal resistance). Add external forces like gravity $\mathbf{f}_\text{ext} = m \mathbf{g}$ (straight to momentum).
        
    - Write out the updated momentum (or velocity) for each node to the grid buffer. At this point, the grid nodes have forces applied and are ready for contact and integration.
        
    
    _Implementation note:_ The plastic yield and flow should preserve volume for incompressible plastic flows (non-associative flow rule can be usedweb.cs.ucla.edu). The algorithm essentially does an elastic predictor – plastic corrector (return mapping) per node, which can be done analytically for simple yield surfaces or via a few iterations for complex ones[sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0266352X2400836X#:~:text=The%20goal%20of%20the%20return,Let%20us). The computed plastic strain increment could be stored back to particles later (during G2P) for accumulated deformation.
    
7. **Contact Projection & Constraints (GPU):** Integrate a GPU kernel to enforce **contact constraints** on grid nodes after internal forces are applied. In ziXn, contact refers to interactions between the deformable ground and rigid bodies (wheels, tracks, feet, etc.) including friction limits at the interface[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=develop%20the%20MPM%2C%20with%20funding,known%20as%20CRAMP). We implement a projection of node velocities to satisfy two conditions at contact points: (a) no penetration (normal velocity ≥ 0 relative to the collider surface), and (b) Coulomb friction limit (tangential shear force ≤ μ * normal force). This is similar to applying a frictional boundary condition as in Daviet and Bertails-Descoubes (2016)[yzhu.io](https://yzhu.io/publication/mpmmls2018siggraph/paper.pdf#:~:text=,to%20the%20adoption%20of) or Bardenhagen’s multi-material contact in MPM[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=develop%20the%20MPM%2C%20with%20funding,known%20as%20CRAMP). The algorithm operates on grid nodes that lie inside or near a collider (e.g. within a wheel’s contact patch or below ground level):
    
    - Identify contact-active nodes: those which are in or very near a rigid colliding geometry. For each such node, compute the contact normal (usually provided by collider shape).
        
    - Project the node’s post-force velocity $\mathbf{v}$: if $\mathbf{v}$ has a negative component along the contact normal (meaning penetration velocity), set that normal component to zero (or to a small positive value if using bouncy contact)[cof.orst.edu](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=Given%20normal%20vector%2C%20early%20MPM,condition%20and%20the%20definition%20of)[cof.orst.edu](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=compression,1A%20shows%20a%201D%20material). This ensures no further interpenetration.
        
    - Compute the maximum allowed tangential velocity magnitude by $v_t^{\max} = \mu, v_n$ (where $v_n$ is the remaining normal velocity, typically 0 for sticking contact). If the node’s tangential velocity $|\mathbf{v}_t|$ exceeds this, scale $\mathbf{v}_t$ down to $v_t^{\max}$ in the direction of motion (enforcing the friction cone).
        
    - If the tangential velocity is below the threshold, we leave it (sticking contact).
        
    - Update the node’s momentum accordingly (since $p = m v$).
        
    
    This can be done in one shader pass over all grid nodes, with each thread checking if its node is in contact with any rig. The **zx_rig** data (which contains positions, orientations, wheel radii, etc.) should be accessible (perhaps via a constant buffer or structured buffer). For example, for ground contact, we can mark any node with world height < ground collider height as contact node. For wheel contact, any node inside a wheel’s volume might be a contact node with normal pointing upward.
    
     
    
    **Algorithm Steps:**
    
    - Determine contact geometry each frame (from the physics engine or game): e.g. positions of wheels, terrain level for ground plane, etc. Prepare this as input (list of colliders).
        
    - Dispatch a kernel with one thread per grid node. Each thread finds if its node is within a collider’s influence (this could be a simple check like if node is below ground plane or inside a wheel radius in 2D projection).
        
    - If node i is a contact candidate, compute penetration depth or relative velocity. For a ground plane, if node_y < ground_y, it’s penetrating. For a wheel, compute distance from wheel center in horizontal plane.
        
    - Compute contact normal (for ground it’s (0,1,0) upward; for wheel, normal is radial from wheel center at node position).
        
    - Project velocity: let $v_n = \mathbf{v}\cdot \mathbf{n}$. If $v_n < 0$ (moving into surface), set $v_n=0$. Recompose $\mathbf{v}$.
        
    - Compute $v_t$ (tangential component of old velocity) = $\mathbf{v} - v_n \mathbf{n}$. If $|\mathbf{v}_t| > \mu, v_n$ (with $v_n$ now nonnegative), clamp it: $\mathbf{v}_t := \mathbf{v}_t \min(1, \frac{\mu v_n}{|\mathbf{v}_t|})$[cof.orst.edu](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=,however%2C%20the%20Nc%20definition%20still).
        
    - Reassign $\mathbf{v} = v_n \mathbf{n} + \mathbf{v}_t$. Update node momentum = $m \mathbf{v}$ in the grid buffer.
        
    - (Optionally, count how many nodes were clamped for diagnostics – this can increment `contact_saturated` in `zx_counters`.)
        
    
    After this pass, grid node velocities obey non-penetration and friction limits. This approach is essentially a one-iteration solve of the contact complementarity problem, which is typically sufficient for small time-steps. It mirrors the approach of enforcing Coulomb friction constraints on MPM nodes as done by Bardenhagen et al. (2001) and later work: multiple velocity fields or velocity splits can also be used[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=develop%20the%20MPM%2C%20with%20funding,known%20as%20CRAMP), but the projection method is simpler for our real-time use.
    
8. **Optional Slurry Pressure Solve (PCG):** For two-phase materials like mud or wet soil, incorporate an **iterative pressure solve** to enforce incompressibility of the fluid phase. This is analogous to a fluid simulation’s pressure projection step – often done with a Preconditioned Conjugate Gradient (PCG) solver on a Poisson equation. In ziXn, this would apply if a material is flagged as nearly incompressible (liquid or fully saturated soil). The algorithm:
    
    - Identify nodes or cells that contain fluid (or pore fluid in soil). Compute the divergence of the velocity field at those nodes (after the contact step, before G2P).
        
    - Form a linear system $A p = b$, where $p$ is the pressure at each fluid node and $b$ is the negative divergence (mass change)[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/MPM-particle-laden-flow.pdf#:~:text=Flows%20pages,Tampubolon%20et%20al). The matrix $A$ comes from the discrete Laplacian on the fluid region (taking into account open boundaries or boundary conditions at solid walls).
        
    - Use PCG to solve for pressure $p$ that makes the divergence nearly zero. Typically, this involves a few iterations (since real-time constraints may limit full convergence). A Jacobi or multigrid preconditioner can accelerate convergence.
        
    - Adjust node velocities: subtract the gradient of the solved pressure from each node’s velocity (this is the pressure impulse that removes divergence). Specifically, for each neighbor pair of fluid nodes, apply $\Delta \mathbf{v} \propto -(p_{\text{neighbor}} - p_{\text{node}})$ along the connecting axis.
        
    - This will ensure that fluid nodes satisfy continuity (no net volume change) at the end of the step.
        
    
    To implement efficiently, note that fluid nodes might span many tiles – it’s effectively a global solve. We can still perform it on GPU: assemble the RHS $b$ in a compute shader, then run a loop of PCG iterations. Each PCG iteration can be split into kernels for sparse matrix-vector multiply (Laplacian * pressure), vector dot product, and update steps. Modern GPUs can handle moderate system sizes (tens of thousands of unknowns) in sub-millisecond times with Jacobi or diagonal preconditioning.
    
     
    
    **Algorithm Steps:**
    
    - Mark fluid-active grid nodes (e.g. those belonging to slurry particles). Also determine connectivity – each node connects to its 6 neighbors (in a regular grid) that are also fluid.
        
    - Compute divergence at each fluid node: $b_i = -\sum_{\text{neighbors } j} (\mathbf{v}_j - \mathbf{v}_i)\cdot \mathbf{n}_{ij}$, which for a regular grid simplifies to finite differences of velocity. This is done in a shader after contact (using the velocities computed in Task 6-7).
        
    - Initialize pressure array $p_i = 0$. Set tolerance or iteration count limit.
        
    - Iterate PCG (on GPU or CPU):
        
        1. **Matrix-vector multiply:** For each fluid node, compute $(A p)_i = \sum_{j \in \text{nbrs}(i)} (p_i - p_j)$ (Laplacian action). This can be one compute kernel.
            
        2. **Compute residual:** $r_i = b_i - (A p)_i$. Compute global residual norm (sum of $r_i^2$) using a reduction (another kernel). If below tolerance, break.
            
        3. **Precondition (optional):** e.g. Jacobi preconditioner would divide $r$ by diagonal of $A` (which is number of neighbors, but we can just incorporate that into relaxation parameter).
            
        4. **Update search direction:** standard PCG uses $d^{(k+1)} = r^{(k+1)} + \beta d^{(k)}$ with $\beta = (r^{(k+1)}\cdot z^{(k+1)})/(r^{(k)}\cdot z^{(k)})$. Implement dot products with parallel reductions.
            
        5. **Update pressure:** $p := p + \alpha d$, with $\alpha = \frac{r\cdot z}{d\cdot A d}$ computed via reductions and a SpMV from step (i).
            
    - After iterations, use the final pressure field to adjust velocities: for each fluid node, for each neighbor that is fluid, apply $v_{i,\text{new}} = v_i - \frac{\Delta t}{m_i} (p_j - p_i) A_{ij}$, which effectively is $-\nabla p$. In practice, since we use equal volumes, this simplifies to distributing pressure difference equally into velocity.
        
    - Insert a barrier to ensure all pressure updates complete.
        
    
    This advanced feature would run after contact (Task 7) and before writeback. For performance, we might only do a few iterations each frame – enough to remove the bulk of compressibility (making the material behave like mud rather than water). Research has shown multi-phase MPM can use a split solver like this to simulate saturated sand[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/MPM-particle-laden-flow.pdf#:~:text=Flows%20pages,Tampubolon%20et%20al). If computation is too heavy, we can run this on a subset of tiles (e.g. only where slurry exists) or at a coarser grid resolution (half resolution pressure solve) to save time. Overlap it with other GPU work if possible (Task 12 suggests async compute streams).
    
9. **G2P Kernel & Particle Update:** Once grid node velocities are finalized (after internal forces, contact, and pressure solve if any), perform the **Grid-to-Particle (G2P)** transfer to update particle velocities and positions[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=extrapolated%20from%20the%20surrounding%20nodes,depending%20on%20%2067%20integration). This is implemented as a GPU compute shader (`kernels/g2p.hlsl`). Each particle will gather velocity from its surrounding grid nodes (the same ones used in P2G) using the interpolation weights. Update the particle’s velocity and position: $\mathbf{v}_p^{\text{new}} = \sum_{i \in \text{nodes}} N_i(\mathbf{x}_p), \mathbf{v}_i$ and $\mathbf{x}_p^{\text{new}} = \mathbf{x}_p + \mathbf{v}_p^{\text{new}} \Delta t$. Also update the particle’s affine momentum matrix **C** (for APIC) or deformation gradient for FLIP/MLSPIC[arxiv.org](https://arxiv.org/html/2308.01060v2#:~:text=Jiang%20et%20al%20,numerical%20dissipation%20and%20improving%20preservation). In APIC, we maintain an affine velocity matrix per particle that enhances rotational motion. After G2P, update this via $ \mathbf{C}_p^{\text{new}} = \mathbf{B}_p = \sum_{i} N_i(\mathbf{x}_p) (\mathbf{v}_i - \mathbf{v}_p)\otimes (\mathbf{x}_i - \mathbf{x}_p)$[arxiv.org](https://arxiv.org/html/2308.01060v2#:~:text=Jiang%20et%20al%20,numerical%20dissipation%20and%20improving%20preservation), which effectively stores the local velocity gradient[arxiv.org](https://arxiv.org/html/2308.01060v2#:~:text=reducing%20the%20numerical%20dissipation%20of,numerical%20dissipation%20and%20improving%20preservation). This matrix $\mathbf{C}_p$ will be used in the next P2G to carry forward sub-grid velocity information (reducing dissipation).
    
     
    
    Additionally, apply any particle-level plasticity updates. For example, for solid particles (snow or soil), we may have stored an elastic deformation gradient $F_p$. After updating particle velocity (which gives an elastic predictor), we can do a **particle “return mapping”**: ensure $F_p$ is clamped such that the stress lies on the yield surface[sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0266352X2400836X#:~:text=The%20goal%20of%20the%20return,Let%20us). In practice, this might mean decomposing $F_p$ into elastic and plastic parts and capping the amount of elastic strain (for snow, one might reduce $F_p$ volume if it exceeded compaction limits, etc. as per ALGORITHMS doc). The specifics depend on material model – for now, we ensure consistency by possibly adjusting particle volume or deformation based on the node pressures computed.
    
     
    
    **Algorithm Steps:**
    
    - Dispatch one thread per particle. Each thread reads the particle’s current position $\mathbf{x}_p$ and old velocity $\mathbf{v}_p$ from the SoA.
        
    - Determine the grid cell index and the shape function weights $N_i$ for neighboring nodes (same as P2G, but now used for interpolation). Fetch each neighbor node’s velocity $\mathbf{v}_i$ from the grid buffer (SRV).
        
    - Compute new particle velocity: $\mathbf{v}_p^{\text{new}} = \sum_i N_i, \mathbf{v}_i$. Also compute the APIC affine matrix: initialize it to zero, then for each neighbor node do $\mathbf{C}_p += N_i, (\mathbf{v}_i - \mathbf{v}_p^{\text{new}})\otimes (\mathbf{x}_i - \mathbf{x}_p)$[arxiv.org](https://arxiv.org/html/2308.01060v2#:~:text=Jiang%20et%20al%20,numerical%20dissipation%20and%20improving%20preservation). Store $\mathbf{C}_p$ in the particle buffer.
        
    - Update position: $\mathbf{x}_p^{\text{new}} = \mathbf{x}_p + \mathbf{v}_p^{\text{new}} \Delta t$ (with $\Delta t$ from sim params). If using substeps, this might be a fractional $\Delta t$ per step.
        
    - If the material is elastic/plastic, update the particle’s deformation: e.g. compute new $F_p = (I + \nabla \mathbf{v}_p \Delta t), F_p$ if we track deformation gradient. Then enforce any yielding: if $\det(F_p)$ exceeds a compaction limit (snow densification) or shear part exceeds yield, scale it back. This is analogous to a return mapping in continuum mechanics[sciencedirect.com](https://www.sciencedirect.com/science/article/pii/S0266352X2400836X#:~:text=The%20goal%20of%20the%20return,Let%20us) but done on the particle. For example, for snow, one might reduce the trace of $F_p$ if pressure > yield, mimicking the Drucker–Prager cap model[math.ucdavis.edu](https://math.ucdavis.edu/~jteran/papers/SSCTS13.pdf#:~:text=Mathematics%20math,Drucker%20and).
        
    - Write back updated particle velocity, position (and possibly deformation $F_p$, volume, etc.) to the particle buffers.
        
    - If a particle’s new position is outside the bounds of its current tile, mark it (e.g. set a flag in an array) for tile migration. (The actual migration will be handled by the CPU or a subsequent GPU stream, as implemented in Task 3.)
        
    
    After G2P, all particle states are updated to the end of the time step. The APIC method ensures that angular momentum is better preserved than FLIP/PIC by virtue of the affine velocity storage[arxiv.org](https://arxiv.org/html/2308.01060v2#:~:text=reducing%20the%20numerical%20dissipation%20of,numerical%20dissipation%20and%20improving%20preservation). This step, along with P2G, should be verified to conserve momentum and energy (except where plasticity intentionally dissipates energy). We can compare the total kinetic energy before and after as a check (small loss can occur due to interpolation, but it should be consistent each frame for determinism).
    
10. **Terrain Writeback Encoding:** Develop the final GPU pass to translate simulation results into engine-readable output layers (heightmap and material masks). **Writeback** involves sampling the simulation state (which lives on particles or grid) and writing it into textures that the game engine (Unreal/Unity) uses for rendering and gameplay. ZiXn defines several output layers: e.g. **heightfield displacement**, **compression/wetness mask**, **contact friction mask**, etc., as described in the architecture[dev.epicgames.com](https://dev.epicgames.com/documentation/en-us/unreal-engine/runtime-virtual-texturing-in-unreal-engine#:~:text=A%20Runtime%20Virtual%20Texture%20,to%20conform%20to%20the%20terrain). Implementation details:
    
    - **Heightfield Displacement:** We maintain a base terrain and output a texture of height displacements (deltas). We can generate this by taking the vertical position of the simulation’s “ground surface” and subtracting the base terrain height. One way is to rasterize particles onto a height grid: for each tile’s particles, take their height and write to a global heightmap at the corresponding x-y location (possibly using maximum or average if multiple particles overlap a cell). Alternatively, use the grid: after G2P, the grid node positions or an averaged particle height per cell can give terrain height.
        
    - **Compression/Wetness Mask:** If simulating snow or mud, output a mask indicating compression (how dense/compressed the material is) or wetness (saturated vs dry). This could be based on particle volume change ($J = \det(F)$) or a material state variable. We can accumulate this per cell and output as a normalized value 0-1.
        
    - **Contact Mask:** Areas where a collider touched the ground can be marked (for footprints, track marks). We can output 1 where a contact was active this frame (e.g. any node had contact force > 0) and 0 elsewhere, possibly blurred or accumulated over time for persistent effects.
        
    
    These outputs are typically low-frequency (heightfield might be 256×256 or 512×512 for the terrain region). We can dedicate a compute shader to generate each layer. For height, each thread can correspond to one texel, sample nearby particles or grid points, and write the height value (or height delta) to an R32F UAV texture. Masks can be R8 or R16 textures.
    
     
    
    For integration with Unreal’s **Runtime Virtual Texture (RVT)** system, we might directly write into a buffer that Unreal’s rendering thread will blit into the RVT. Unreal expects height in world units (cm by default) encoded in the virtual texture[forums.unrealengine.com](https://forums.unrealengine.com/t/rvt-virtual-heightfield-mesh-multiple-questions-issues/548820#:~:text=RVT%20%2B%20Virtual%20Heightfield%20mesh,according%20to%20youtube%20tutorials%2C). We ensure that our output height map aligns and scales properly with the engine’s coordinate system (e.g. 1 unit = 1 cm). Unity’s terrain uses a heightmap (often normalized 0–1); we’ll need to scale our output to that range and update via their API.
    
     
    
    **Algorithm Steps:**
    
    - Allocate output textures (or use those provided by engine via interop) for displacement and masks. These could be created as UAV resources at the needed resolution.
        
    - **Height pass:** For each output texel (x,y), determine the corresponding simulation world position. Sample the simulation data to get terrain height. If using particles, we could take the highest particle below that (x,y) location (for ground surface) or interpolate from nearest particles. If using the grid, we might take the highest solid node or do a tri-linear interp of node heights. Write the height minus base_terrain_height into the texture.[dev.epicgames.com](https://dev.epicgames.com/documentation/en-us/unreal-engine/runtime-virtual-texturing-in-unreal-engine#:~:text=A%20Runtime%20Virtual%20Texture%20,to%20conform%20to%20the%20terrain)
        
    - **Compression mask pass:** For each texel, find nearby particles and average their compression (e.g. $J_p$ or plastic strain). Write out a value (0 = no compression, 1 = fully compacted or wet) to a channel (e.g. R16).
        
    - **Contact mask pass:** For each texel, determine if a contact happened near that location (we could mark texels touched by a wheel or foot). This can be done by rasterizing contact points: during the contact solve (Task 7), we can record positions of nodes where projection occurred. Blur these slightly and output an R8 mask.
        
    - Synchronize and transition these textures for use on the CPU (or rendering). For Unreal, if we have a pointer to a GPU texture that backs a Runtime Virtual Texture, we might not need to copy – just signal the fence (next task) for the engine to know it’s updated. In Unity, we may copy this to a Texture2D via `GraphicsBuffer.GetNativeBufferPtr` and a plugin callback.
        
    
    The result of writeback is that the game engine’s terrain is updated **the same frame** as the simulation. For example, in Unreal Engine, we render the Virtual Heightfield Mesh using the updated RVT which now contains the displacement map we wrote[forums.unrealengine.com](https://forums.unrealengine.com/t/rvt-virtual-heightfield-mesh-multiple-questions-issues/548820#:~:text=RVT%20%2B%20Virtual%20Heightfield%20mesh,according%20to%20youtube%20tutorials%2C). These output layers can also drive particle effects or audio (e.g. a wetness mask making puddle particles, or compression affecting footstep sound). We ensure the outputs are properly double-buffered if needed (so we’re not writing to a texture that’s being read for rendering simultaneously).
    
11. **GPU Synchronization & Fencing:** Implement proper synchronization around GPU tasks to maintain execution order and allow the CPU to know when the simulation work is complete. In Direct3D12, use **UAV barriers** between unordered-access phases that need ordering (for example, after P2G and before the constitutive kernel, insert `UAVBarrier()` to ensure all thread groups’ writes to the grid are visible)[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=CUDA%20Driver%20checks%2C%20if%20the,page%20copy). Also manage resource state transitions: e.g. transition the grid buffer from UAV (writeable) to SRV (readable) when moving from P2G to the grid update shader, and back to UAV if writing again. We will accumulate a command list that encodes the entire simulation pipeline each frame:
    
    1. Record transitions for used resources (particle buffers SRV, grid buffer UAV, output textures UAV).
        
    2. Dispatch P2G shader. Insert `UAVBarrier()`.
        
    3. Dispatch grid update shader (internal forces). `UAVBarrier()`.
        
    4. Dispatch contact shader. `UAVBarrier()`.
        
    5. If pressure solve, dispatch multiple compute shaders for PCG iterations (with UAV barriers as needed, or use separate command list/queue as in Task 12).
        
    6. Dispatch G2P shader (reading grid, writing particle data). `UAVBarrier()`.
        
    7. Dispatch writeback shaders (writing textures).
        
    8. Transition output textures to copy or present state if CPU needs them.
        
    
    After enqueueing all these, signal a GPU **Fence** (ID3D12Fence) to mark completion. The `zx_frame_end.fence` value (as defined in `zx_abi.h`) should carry the fence value the GPU will signal when the command list finishes[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31). The engine (Unreal/Unity) will wait for this fence before using the output (for example, Unreal might wait before reading the updated height texture). We create this fence at device initialization and use `CommandQueue->Signal(fence, value)` at the end of each frame’s list.
    
     
    
    Also handle multiple frames in flight: if the engine allows two frames of simulation to overlap, we might use two fence values and double-buffered resources. For simplicity, we can start with a single frame in flight (i.e. simulate only after the previous is done). The timeline can be improved later.
    
     
    
    In Vulkan, use a similar approach with a **VkFence** or timeline semaphore. Since Vulkan operations are recorded and submitted, we’d use `vkQueueSubmit` with a fence that signals when GPU work is done. Alternatively, use a timeline semaphore and signal a value corresponding to frame count[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31).
    
     
    
    **Algorithm Steps:**
    
    - During `zxBindDevice`, create a fence object (D3D12) or timeline semaphore (Vulkan). Initialize fence value to 0.
        
    - Each `zxBeginFrame`, record commands as above. At `zxEndFrame`, increment the fence value and signal it on the command queue.
        
    - Store the fence value in `zx_frame_end.fence` so the engine can wait on it[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31). In D3D12, we also provide an event or use `WaitForFence` if the engine requests a CPU wait (for dedicated server mode, perhaps).
        
    - Insert necessary barriers: D3D12 UAV barriers after each unordered phase, and resource transitions when switching a buffer from read to write. For example, after writing the height texture in compute, transition it to PIXEL_SHADER_RESOURCE for the rendering.
        
    - Test by running a frame and then immediately reading back a small portion of the output (or using PIX GPU captures) to ensure all work completed by the fence.
        
12. **Asynchronous Compute Overlap:** Optimize GPU utilization by overlapping independent parts of the simulation using multiple command queues. Modern GPUs (and especially consoles) can run compute tasks concurrently if they use different units or have free resources. According to the design, certain passes can run in parallel[gpuopen.com](https://gpuopen.com/download/GDC2017-Asynchronous-Compute-Deep-Dive.pdf#:~:text=Async%20Compute%20%E2%86%92%20More%20Performance,CPU). For example, while the heavy grid update (internal forces and pressure solve) is running on one compute queue, the writeback of the previous frame’s data (which might be a texture copy or compression) could run on another queue. Another scenario: if the pressure solve (Task 8) is on a large 3D grid, it might benefit from running on a compute queue while, in parallel, the graphics queue handles rendering of the previous frame’s terrain.
    
     
    
    Specifically, we can designate:
    
    - **Compute Queue 0:** Main simulation tasks (P2G, Grid, Contact, G2P).
        
    - **Compute/Copy/Graphics Queue 1:** Writeback tasks (which might just be copying data to textures or updating rendering structures).
        
    
    Use timeline semaphores (or fences) to synchronize between queues. For example, the G2P must finish (on Queue 0) before we start using those particle results for rendering, but the writeback of frame N and simulation of frame N+1 could potentially overlap. In D3D12, one can use two command lists recorded in parallel threads and then submit to two different queues. Use `ID3D12CommandQueue::Wait` to make one queue wait for an event from the other when needed.
    
     
    
    As a concrete plan: issue the long-running pressure solve on Compute Queue 0. Immediately after kicking it, we can start the Writeback passes on Compute Queue 1 (which only depend on final grid velocities). However, the pressure solve modifies grid velocities (for fluid nodes), so actually we need to wait until it’s done to do final G2P for those particles. If we assume the pressure solve is a small part or not always used, another overlap opportunity is running **graphics rendering** (like the terrain drawing) in parallel with the next frame’s simulation. In Unreal Engine, for instance, we could overlap ziXn simulation (compute) with the rendering of the frame if the GPU has spare compute units.
    
     
    
    We will implement basic support: create two D3D12 command queues (one for Direct compute, one for copy or compute), and route certain workloads to the second queue. Insert fence waits accordingly. For example:
    
    - Submit P2G, grid, contact, G2P on Compute0. After G2P, signal a fence F0.
        
    - Meanwhile on Compute1, once F0 is signaled, do Writeback. (If Writeback doesn’t depend on fluid solve, it could run slightly earlier, but generally it depends on final particle positions, so after G2P anyway.)
        
    - Submit rendering (on graphics queue) that uses the updated textures with a wait on the same fence F0 to ensure physics is done.
        
    
    On Vulkan, we can create two `VkQueue`s (from same family if needed) and use `vkQueueSubmit` with semaphores to synchronize tasks. For instance, a semaphore signals after physics dispatches, and rendering queue waits on it before sampling the updated height texture.
    
     
    
    **Algorithm Steps:**
    
    - Expose multiple command queues in `zx_context` (e.g. `context->gpu_compute_queue` and maybe reuse the graphics queue from the engine for certain tasks).
        
    - Determine which parts can overlap. A safe starting point is to overlap _terrain rendering_ with _particle simulation for next frame_. This requires splitting `zxEndFrame` such that it signals the fence as soon as writeback is done, and the engine can start rendering, while we begin next frame’s simulation concurrently.
        
    - Use fine-grained timing: measure each stage’s GPU time (Task 21) to identify bottlenecks. If, say, pressure solve takes 0.2 ms and writeback 0.3 ms, overlapping them could save ~0.2 ms.
        
    - Implement queue synchronization: in D3D12, use `ComputeQueue1->Wait(ComputeQueue0, fenceValue)` to ensure ordering when needed. In Vulkan, use `vkCmdPipelineBarrier` with `VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT` waits or submit-level semaphores.
        
    - Test with scenarios where GPU is heavily loaded to see if overlapping improves frame rate (in some cases, it might not if one queue saturates the GPU). Use Nvidia Nsight or similar to verify concurrency.
        
    
    The goal is to reduce idle time on the GPU. By overlapping independent tasks (or tasks from consecutive frames), we strive to approach full GPU utilization[developer.nvidia.com](https://developer.nvidia.com/blog/advanced-api-performance-async-compute-and-overlap/#:~:text=This%20post%20covers%20best%20practices,applications%2C%20see%20all%20Advanced). This can be especially beneficial on AMD GPUs which have separate async compute engines – e.g. physics can run on ACE while graphics runs on GFX, improving throughput.
    
13. **Direct3D12 Backend Finalization:** Complete the implementation of the D3D12 backend (in `backends/d3d12/`). This includes creating and managing:
    
    - A **descriptor heap** for all SRV/UAV/CBV resources. We will use a shader-visible descriptor heap with enough capacity for all tiles’ buffers, particle buffers, and output textures. For simplicity, we can allocate descriptors statically (e.g. a heap with 1000 entries) and store indices in `zx_gpu_buffer` objects.
        
    - A root signature that likely has a descriptor table for SRVs (particle data), a descriptor table for UAVs (grid, outputs), and root constants (like `zx_frame_params` with dt, etc.)[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=,locked%20through%20cudaMallocHost%20or%20cudaHostAlloc).
        
    - Pipeline State Objects (PSOs) for each compute shader (P2G, grid update, contact, G2P, writeback). Compile HLSL to DXIL and create PSOs at init or on-demand.
        
    - Memory allocators: use `ID3D12Heap` to allocate large chunks for tile grids and particle arrays. We might allocate one large default heap for all tile grids and suballocate (as per Task 19, using a buddy allocator). Same for upload heaps for data transfer.
        
    - **Timestamp queries**: Create a query heap for timestamps (D3D12_QUERY_HEAP_TYPE_TIMESTAMP) and a readback buffer. Before and after each major dispatch, insert `EndQuery`. After execution, resolve the query heap to the readback buffer. Map it on CPU and compute timing differences, store in `zx_counters` (which has fields for P2G, simulation, etc.)[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=CUDA%20Driver%20checks%2C%20if%20the,page%20copy)[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31).
        
    - Error handling: wrap D3D calls with checks. If a device is lost or out-of-memory (D3D debug layer will output messages), translate to `ZX_E_DEVICE` or `ZX_E_OOM`. Also ensure to release all D3D resources in `zxDestroyContext`.
        
    
    We also need to implement `zxBindDevice()` for D3D12, which creates the D3D12Device, command queues, command allocators, and lists, as well as descriptor heaps and fence objects. **Bindless** resource management (as mentioned in architecture) can be achieved by the single large descriptor heap where each GPU buffer knows its descriptor index (like a pointer). We might not implement fully variable indexing in shaders (which requires DX12 Tier 2), but rather bind the heaps and use indices as needed.
    
     
    
    **Algorithm Steps:**
    
    - **Initialization:** Create D3D12 device via D3D12CreateDevice (with adapter selection if needed). Create a command queue (and a second for async if using Task 12). Create command allocators (one per frame in flight) and a command list.
        
    - Create descriptor heap: e.g. D3D12_DESCRIPTOR_HEAP_DESC with Type = CBV_SRV_UAV, NumDescriptors = N, Flags = SHADER_VISIBLE. Keep CPU handles for writing descriptors.
        
    - Create root signature: define ranges for SRV (particle buffers, read-only), UAV (grid buffers, writeback textures), and perhaps CBV for constants. Use D3D12_ROOT_SIGNATURE_FLAG_NONE (allowing compute). Serialize and create it.
        
    - Create compute shaders: compile HLSL files (or load precompiled DXIL) for each kernel. For each, create D3D12_COMPUTE_PIPELINE_STATE_DESC linking the root signature and the CS bytecode.
        
    - Setup memory allocator for tile buffers: use `ID3D12Device->CreateHeap` for a large chunk of DEFAULT heap memory (GPU local). Use alignment = 256 bytes for UAV. Implement suballoc (buddy system as in Task 19) – e.g. maintain a free list of blocks. When a tile is created, allocate B^3 * node_struct_size bytes from the heap (or use an existing buffer from a pool) and create a ID3D12Resource buffer placed in that heap at the offset. This yields efficient memory reuse and minimal fragmentation[bruce-lee-ly.medium.com](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=)[bruce-lee-ly.medium.com](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=The%20advantage%20of%20the%20buddy,be%20a%20waste%20of%20memory).
        
    - Implement `zxMapBuffer`/`zxUnmapBuffer` for CPU readback if needed (for debugging or sync).
        
    - During simulation, record commands in the command list as per tasks above. Close and execute the command list on the queue, then wait on the fence (or let engine wait).
        
    - Use PIX or Nsight to capture a frame and ensure all resources are bound correctly (e.g. check that descriptor indices in shaders match the layout we set in the heap).
        
    - Fill out `zx_counters`: after each frame, when the query results are available (which might be next frame), calculate time differences in nanoseconds for each region and provide in the struct for the engine to read[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=CUDA%20Driver%20checks%2C%20if%20the,page%20copy)[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31).
        
    
    After finalizing, test on an NVIDIA GPU first (since these are quite forgiving with descriptor heap sizes and minor errors). Then test on an AMD GPU, which will ensure we handle things like alignment and barrier rules strictly. We should observe the simulation running and profile the timing to ensure we meet the ~1–4 ms frame budget on a high-end GPU for moderate tile counts.
    
14. **Vulkan Backend & Cross-Vendor Support:** Begin implementing the Vulkan backend (`backends/vk/`). Many structures mirror D3D12, but the API details differ:
    
    - Use **Vulkan-HPP** or C API to create a `VkInstance`, select a `VkPhysicalDevice` (discrete GPU preferred), create `VkDevice` with a compute queue. Since we also want to support rendering interop, we may create a device with graphics queue as well (or accept one from engine – but better to create our own for now).
        
    - Set up descriptor pools and descriptor set layouts analogous to the D3D root signature. In Vulkan, we can use a single large descriptor set with bindless indexing by enabling the extension `VK_EXT_descriptor_indexing` (and the feature for runtime descriptor array). This allows us to create a descriptor set with, say, 1000 storage buffers and sample them with an index in shader[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=,locked%20through%20cudaMallocHost%20or%20cudaHostAlloc). Alternatively, since our number of resource types is limited, we could also use push descriptors or a combination of sets (but bindless is cleaner).
        
    - Memory allocation: use Vulkan Memory Allocator (VMA) or roll our own similar to D3D12. Vulkan requires explicit memory type selection. We will allocate large `VkDeviceMemory` blocks for buffers and textures. Using VMA is convenient as it will handle a buddy system internally; we just provide size and usage.
        
    - Pipeline creation: similar to D3D12 PSO, create `VkShaderModule` for each SPIR-V shader and a `VkPipeline` with `vkCreateComputePipelines`. The descriptor set layout must match shader bindings.
        
    - Command buffer recording: for each frame, record commands for each pipeline dispatch, with barriers (Vulkan memory barriers for UAV accesses). Use `vkCmdPipelineBarrier` with `VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT` and `VK_ACCESS_SHADER_WRITE_BIT` etc., to enforce ordering between P2G and grid update[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=CUDA%20Driver%20checks%2C%20if%20the,page%20copy).
        
    - Fences and synchronization: use `vkQueueSubmit` with a `VkFence` for simulation completion, analogous to D3D12 fence. Use timeline semaphores if multiple queues (VK_KHR_timeline_semaphore).
        
    - Debugging: enable validation layers in debug mode to catch issues (they will flag any misuse of descriptor indexing or barriers).
        
    
    Cross-vendor considerations:
    
    - Ensure we test on AMD hardware (which in Vulkan may use subgroup size 64 vs NVIDIA’s 32). Avoid assumptions about warp size in shaders; if we use subgroup operations (like ballot), query `SubgroupSize` or write shaders that adapt (HLSL’s WaveGetLaneCount or GLSL’s `gl_SubgroupSize`).
        
    - Memory alignment: Vulkan’s `minUniformBufferOffsetAlignment` and others need to be respected. For storage buffers, typically 16-byte alignment is fine, but we align to 256 as in D3D12.
        
    - Different vendors’ compilers might be sensitive: e.g. avoid using `groupshared` memory beyond 32KB as older AMD GCN have 32KB LDS limit per workgroup. Limit workgroup size and shared memory per the lowest common denominator (we target 256 threads per group and ~16KB shared memory, which is safe).
        
    - Test on an Intel Arc as well if available, since Arc supports Vulkan and has different architecture (to ensure no vendor-specific issue).
        
    
    **Algorithm Steps:**
    
    - Abstract common parts of GPU backend in a way that most simulation code (command recording logic) can be shared, with API-specific calls hidden behind function pointers or macros. This avoids duplicating all 30k lines of code for each API.
        
    - Implement `zxCreateContextVk` analogous to D3D12: pick physical device, create logical device (with `VkQueue`), create descriptor set layout and allocate a big descriptor pool.
        
    - Translate shaders to SPIR-V: either use DXC to compile HLSL to SPIR-V offline, or use GLSL/GLSLang. We might use `#ifdef VK` in shaders to handle slight differences (like descriptor binding syntax). The architecture mentioned cross-compiling HLSL to SPIR-V[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=,locked%20through%20cudaMallocHost%20or%20cudaHostAlloc)[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=developer,driver%20can%20detect%20pinned%20memory), so we’ll likely use DXC with `-spirv`.
        
    - Allocate buffers and descriptor sets for particle data, grid, etc. Update descriptors with `vkUpdateDescriptorSets`.
        
    - Record command buffers for one frame’s simulation. Submit and wait on fence.
        
    - Compare results between Vulkan and D3D12 on the same GPU. They should produce identical outcomes (this is a good test of determinism as well – floating-point differences between APIs should be none if we use same math and no undefined behavior).
        
    - Optimize if needed: e.g. Vulkan might incur more CPU overhead per submission; we can mitigate by reusing command buffers and just updating descriptors via `vkCmdBindDescriptorSets` with dynamic offsets for different tiles if needed.
        
    
    By having both backends, we ensure ziXn runs on a wide range of hardware: D3D12 covers Windows (NVIDIA/AMD), Vulkan covers both Windows (especially AMD, Intel) and potential future Linux/Android platforms. Down the line, we could also consider using DirectX12 on Xbox and Vulkan or a Metal backend for other platforms. The design remains API-agnostic, and our shader code is portable through common HLSL. We’ll keep an eye on **wave32 vs wave64** issues (for example, if we use Wave Intrinsics in HLSL, on AMD’s wave64 hardware the behavior might differ – thus for determinism, it might be better not to use them or to test thoroughly).
    

### CPU Fallback & Multithreading

15. **CPU Simulation Loop Implementation:** Develop the CPU backend so that `zxSimulate` can run without a GPU (for dedicated servers or when GPU is not available)[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=,explicit%20formulation)[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=2,extrapolated%20back%20to%20material%20points). This involves executing each phase of the MPM simulation on the CPU using the same data structures. We should aim for the CPU results to closely match the GPU results (within floating-point rounding differences). The tasks:
    
    - **P2G (Particle to Grid):** Iterate over all particles, for each particle compute its weighting to 8 (or 27) surrounding grid nodes[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=1,the%20same%20used%20in%20FEM). Accumulate particle mass and momentum to those nodes. We can use the same code as in `zx_apic_ref.cpp` (which presumably already has a reference implementation[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=1,the%20same%20used%20in%20FEM)). The difference is that now we support multi-tile: so for each tile’s particle list, we offset the node indices by the tile’s origin in the global grid.
        
    - **Grid update:** For each active tile’s grid nodes, compute forces and update velocities. On CPU, we can mirror what the GPU does: loop over nodes, compute stress from particle data or from stored elastic state, apply gravity, etc. We might reuse some of the `zx_apic_ref` functions or the ALGORITHMS formulations (ensuring to include yield criteria). Because the CPU might be used in a server authoritative simulation, determinism is crucial here as well – we will avoid any parallel reduction that isn’t deterministic.
        
    - **Contact:** Check each contact node (e.g. ground nodes at y=0, or near rig contacts) – adjust velocities as per Task 7. CPU can handle this with straightforward conditionals and math.
        
    - **Pressure solve (if enabled):** Solve the linear system on CPU. We can reuse a CG implementation from a library or write a simple one. For moderate grid sizes, a CPU can do a few iterations quickly. If the server runs at a fixed lower tick (say 30 Hz), it might manage a larger iteration count.
        
    - **G2P (Grid to Particle):** Loop over all particles, interpolate velocity from grid[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=extrapolated%20from%20the%20surrounding%20nodes,depending%20on%20%2067%20integration), update positions, and update the APIC matrix. This again can follow the reference APIC code[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=extrapolated%20from%20the%20surrounding%20nodes,depending%20on%20%2067%20integration).
        
    
    We will integrate this into `zxSimulate` such that if no GPU context is bound, it automatically uses the CPU path. We might introduce a flag or simply detect `ctx->gpu_device == NULL`.
    
     
    
    **Algorithm Steps:**
    
    - Ensure data structures are accessible on CPU: the tile grid data might currently be GPU-only (allocated in GPU memory). We need a CPU mirror or to allocate CPU memory for grid nodes. One approach: always allocate core simulation data in host memory (possibly aligned for SIMD) and have the GPU backend upload/download as needed. Since we want minimal overhead, we might maintain two copies when GPU is used, but for CPU it would just use the main copy.
        
    - P2G: Initialize a global (or per-tile) grid array to zero. For each particle: get tile index, then the 8 neighbor node indices (this requires knowing the particle’s position within the tile and the tile’s origin in the global grid). Compute weights (N_x, N_y, N_z) in each dimension based on particle’s fractional position[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=1,the%20same%20used%20in%20FEM). Distribute mass and momentum to the 8 nodes. This involves adding to arrays, which is straightforward single-threaded; for multi-thread, see Task 16.
        
    - Grid update: For each tile’s nodes:
        
        - Compute velocity = momentum/mass for each node (guard against division by zero).
            
        - Compute internal forces: e.g. using finite differences, for each node get its neighbors’ velocities to estimate strain rate, then stress. If material info is particle-based, we might loop particles to update stress (this is what APIC does: update particle’s C or F, but some models accumulate to grid).
            
        - Update momentum: $\mathbf{p}_i += (\mathbf{f}_{\text{int}} + \mathbf{f}_{\text{ext}}) \Delta t$.
            
        - Apply yield: limit each node’s stress or velocity change if above yield (the CPU can directly do the yield check similar to GPU).
            
    - Contact: For each node that is below ground (y < 0) or inside an object, project velocity: if v_y < 0, set v_y = 0; clamp horizontal v as $\le \mu v_y$[cof.orst.edu](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=%29%20%C2%B7%20n%CB%86%20,In%20other%20words%2C%20this%20condition). Update momentum.
        
    - Pressure solve (if applicable): Build a sparse matrix of size = (# fluid nodes). Use an iterative solver (Jacobi or PCG). This could be time-consuming on CPU for large grids, so possibly skip or do few iterations.
        
    - G2P: For each particle, get its 8 neighbor nodes, interpolate new velocity[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=extrapolated%20from%20the%20surrounding%20nodes,depending%20on%20%2067%20integration). Update pos: `x += v * dt`. Update APIC C: similar to P2G but using node velocities minus particle velocity.
        
    - Particle state update: If using deformation gradient, update that here for each particle.
        
    - Tile changes: Check each particle’s new tile; if moved, mark for migration (can actually move it in the particle array immediately).
        
    - End of frame: free any tiles marked for deactivation, etc.
        
    
    Once this is implemented, we test with a simple scenario (e.g. one tile with a few particles) to ensure CPU and GPU produce similar results for one step. This single-thread implementation will likely be slow for large particle counts, but it sets the baseline.
    
16. **OpenMP Parallelization of CPU Path:** Optimize the CPU backend using multithreading, leveraging multi-core processors[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y). We can use OpenMP (since it’s cross-platform and easy to integrate) or C++17 `<thread>` with thread pools. The objective is to divide the work of the simulation loop across threads to achieve near-linear scaling with core count. Key parallel regions:
    
    - **Particle loops (P2G and G2P):** Particles can be processed in parallel because each particle updates different grid nodes (except at boundaries which we handle with atomics or local sums). We can do `#pragma omp parallel for` over particles. However, care must be taken for particles that lie on tile boundaries to avoid race conditions on grid nodes. A solution: partition by tiles first (since tiles are mostly independent). For each tile, handle its particles in one thread – this avoids writing to the same node from two threads in most cases (except at tile boundaries, where one tile’s thread might write a ghost node that another thread also writes). To handle that, we use the deterministic merge (Task 4) after the parallel loop.
        
    - **Grid node loops:** Updating nodes can also be parallel per tile or even across tiles. Since each tile’s interior nodes are independent, we can parallelize by tile easily. For boundary nodes that might interact, that’s handled either in the merge step or by using atomic operations for writing forces.
        
    - **Contact and Pressure:** These can also be parallelized by splitting the domain (e.g. each thread handles a subset of nodes for contact checks, since each node’s projection is independent). Pressure solve (if using PCG) can parallelize matrix-vector multiply and vector updates by splitting the index range.
        
    
    OpenMP tasks or sections can also mirror the GPU execution graph[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y). For instance, one could launch an OpenMP task for P2G, another for contact, etc., and use `#pragma omp taskwait` for synchronization points analogous to GPU barriers. However, a simpler approach is to parallelize each major loop and run them sequentially on CPU.
    
     
    
    **Algorithm Steps:**
    
    - Compile with OpenMP support (ensure the build uses `-fopenmp` or MSVC equivalent). Determine an optimal scheduling (static scheduling should be fine if work per particle is uniform).
        
    - **Parallel P2G:** Use `#pragma omp parallel for` over all particles (or over tiles, then inner loop over particles in tile). Within the loop, accumulate to grid using thread-local buffers for each tile to avoid atomic clashes. E.g., allocate an array of grid node accumulators per thread (or use reduction). However, per-thread grid arrays would be large. Instead, use fine-grained locking or atomic adds for the rare case of simultaneous write. We could also do: each thread accumulates into a private array of size (B+1)^3 for the tile it’s processing (including ghost area), then after the loop, merge these arrays into the global grid.
        
    - **Parallel grid update:** For each tile’s nodes, run in parallel. Each thread processes one tile’s node range, computing forces. Since different tiles might update shared boundary nodes (e.g. computing force needs neighbor info), ensure ghost exchange was done or simply let each tile compute internal forces only for its interior, and maybe have a second pass for boundaries.
        
    - Alternatively, parallelize by node index globally (if indexing grid globally). This might cause false sharing if adjacent nodes are updated by different threads. But if we use a coloring approach (like checkerboard update), it could be done.
        
    - **Parallel G2P:** Use `omp parallel for` over particles to interpolate velocities and update positions. This is straightforward since each particle reads multiple nodes – no writes to shared data except particle arrays (each thread writes its own particle, which is fine).
        
    - Use atomic operations or critical sections for any accumulation races. For example, if two threads might add to the same grid node, use `#pragma omp atomic` for the add. Modern CPUs handle millions of atomics per second reasonably well, but we prefer to minimize them (hence tile grouping).
        
    - Balance load: If particle counts vary by tile, static scheduling by tile might leave some threads idle while others work. We could use dynamic scheduling or break large tiles’ particles into chunks.
        
    - Test scaling: run a fixed scenario on 1, 2, 4,... threads (set `omp_set_num_threads`) and measure `zx_counters.frame_time`. We expect roughly inverse scaling until memory bandwidth becomes the bottleneck.
        
    
    Careful attention to memory layout (Task 18) will enhance scaling by reducing cache misses. Also, be mindful of **NUMA** on multi-socket servers – pin threads to sockets and allocate data per socket if needed (though likely overkill unless running >64 threads).
    
17. **SIMD Vectorization on CPU:** Utilize SIMD instructions (AVX2 or AVX-512) to speed up math-heavy loops. The inner computations of MPM – weight calculations, vector operations, etc. – are good candidates for vectorization. Since we use structure-of-arrays (SoA) for particles, we can load multiple particles’ data into SIMD registers easily. For example, process 4 particles at once:
    
    - Load 4 x positions (x4 floats) into an AVX register, do operations, store 4 results.
        
    - Similarly for velocities and other attributes.  
        Many compilers can auto-vectorize simple loops, but for complicated indexing (like scatter/gather to grid) we may need to use intrinsics or reorganize computations.
        
    
    One challenging part is the scatter of particle contributions to grid – AVX2 has gather instructions but scatter only in AVX-512. Even gather is somewhat slow when memory is non-contiguous, often hitting memory bandwidth limits[johnnysswlab.com](https://johnnysswlab.com/when-vectorization-hits-the-memory-wall-investigating-the-avx2-memory-gather-instruction/#:~:text=Vectorization%20puts%20a%20large%20additional,vectorized%20way). To mitigate this:
    
    - We can vectorize the computation of weights and node indices, then perform 4 separate scatter operations (since atomic scatter via SIMD is not straightforward on AVX2). Alternatively, use **structure of arrays** for grid too (separate arrays for mass, momentum x, y, z) and use `_mm256_i32gather_ps` to fetch values, then `_mm256_i32scatter_ps` on AVX-512 or do scalar loops for the scatter part.
        
    - Another approach: sort particles by cell (like on GPU) so that memory accesses to grid are localized. Then process one cell’s particles sequentially (which may not vectorize well, but improves cache).
        
    - Compute multiple particles’ weight contributions concurrently: since each particle affects 8 nodes, we could compute the weight for each of the 8 corners in SIMD across 4 particles at once (this leads to 4x8 operations stored). But scattering 4x8 values might negate the gain.
        
    
    However, vectorization can still be applied to other steps:
    
    - Grid node updates: operations on each node (like computing stress from velocity gradient) can be vectorized across multiple nodes. E.g. load 8 node masses in an AVX register, do a fused multiply-add on them, etc.
        
    - Particle advection: updating positions `x += v*dt` for many particles is easily vectorized by the compiler if memory is aligned.
        
    
    We should use compiler hints like `#pragma omp simd` or eigen library for small matrix ops. For instance, computing the deformation gradient update could benefit from SIMD if we unroll matrix operations.
    
     
    
    **Algorithm Steps:**
    
    - Align data: ensure particle arrays (positions, velocities) are aligned to 32 or 64-byte boundaries (see Task 18) so that `_mm256_load_ps` can be used.
        
    - Enable AVX2: use compiler flags (`/arch:AVX2` on MSVC, `-mavx2` on GCC) and ensure the CPU supports it.
        
    - Identify hotspots by profiling the CPU version: likely P2G and G2P loops are hotspots (lots of arithmetic).
        
    - Attempt auto-vectorization: write loops in a way that is friendly to the compiler (avoid data-dependent conditions inside loops, use arrays instead of pointer arithmetic, etc.). Check compiler output or reports to see if vectorized.
        
    - For manual vectorization: use intrinsics. For example, in P2G weight computation:
        
        `__m256 px = _mm256_load_ps(&particle_x[i]); // load 8 floats // Compute fx = (px - floor(px)) for 8 particles at once: __m256 cellX = _mm256_floor_ps(px); __m256 fx = _mm256_sub_ps(px, cellX); __m256 w0 = _mm256_sub_ps(_mm256_set1_ps(1.0f), fx); __m256 w1 = fx;`
        
        This yields two weight vectors (for simplicity 1D example). For 3D, we’d do similarly for y and z to get 8 weight combinations.  
        We then multiply appropriate combinations for each of 8 nodes in 8-wide vectors. We might then have to scatter to 8 node locations – if those are contiguous in memory (which they might not be), gather/scatter is needed. If the grid indexing for those 8 nodes can be computed (we can compute their indices in an `__m256i` vector), AVX2 gather can load their current mass values. Adding and scattering back might require AVX512. Without AVX512, we might fall back to storing the results in an array and then doing scalar adds in a loop.
        
    - Another micro-optimization: use **SIMD horizontal sum** for reductions. For example, summing partial results in stable summation can use `_mm256_hadd_ps`.
        
    - Also, ensure math operations (like sqrt for stress invariants) use SIMD versions (`_mm256_sqrt_ps`).
        
    - Keep an eye on the memory wall: Vectorizing can quadruple memory throughput demands. For scatter-heavy patterns, it might not help if memory becomes the bottleneck[johnnysswlab.com](https://johnnysswlab.com/when-vectorization-hits-the-memory-wall-investigating-the-avx2-memory-gather-instruction/#:~:text=Vectorization%20puts%20a%20large%20additional,vectorized%20way). Profile with CPU counters (cache misses, etc.). Johnny’s Software Lab analysis indicates that non-contiguous gathers can saturate memory bandwidth, yielding diminishing returns beyond a certain vector width[johnnysswlab.com](https://johnnysswlab.com/when-vectorization-hits-the-memory-wall-investigating-the-avx2-memory-gather-instruction/#:~:text=Vectorization%20puts%20a%20large%20additional,vectorized%20way). If we see that, we might restrict vectorization to the contiguous parts of the algorithm.
        
    
    Ultimately, the CPU path may not reach the performance of the GPU for large problems, but these optimizations ensure that for moderate sizes (or on powerful servers with many cores and AVX-512), the simulation runs fast enough (potentially 30Hz for gameplay). The deterministic nature of a single thread must be preserved in parallel runs – we will have to confirm bitwise identical results with fixed thread counts (which may require turning off some reordering optimizations or summing in double precision to reduce order dependency).
    
18. **Memory Alignment & Cache Optimization:** Audit and adjust memory layouts to optimize cache usage and avoid false sharing. Key points:
    
    - Align frequently accessed structures to cache line boundaries (64 bytes on most CPUs). For example, the `zx_tile` struct could be padded so that each begins at a 64-byte boundary. This ensures that two different tiles’ data don’t share a cache line, which avoids false sharing when threads update different tiles in parallel[stackoverflow.com](https://stackoverflow.com/questions/22766191/what-is-false-sharing-how-to-reproduce-avoid-it#:~:text=Emma%20Jane%20Bonestell%20Over%20a,year%20ago). We can use `alignas(64)` in C++ on such structs.
        
    - Particle SoA arrays should be aligned to at least 32 bytes for AVX. We can allocate with `_aligned_malloc` or `std::aligned_alloc` and also pad the array length to a multiple of the SIMD width (so that vector loops need not handle remainders).
        
    - Avoid placing unrelated data next to each other if they are updated by different threads. For example, if we maintain an array of `tile_active_flags` (bool per tile), each might be just 1 byte – 64 of them could fit in one cache line, meaning 64 threads toggling different flags could thrash. Instead, pad that array so each flag is on a separate cache line (this is excessive for memory, but ensures no interference). This can be done by making it an array of struct with `char flag; char pad[63];`.
        
    - Revisit the spatial and temporal locality: the pattern of particle access to grid is memory-intensive. We could arrange memory so that each tile’s grid nodes are contiguous (already likely true by design). If a tile is BxBxB, then nodes are stored in an array of that length. We should allocate that array contiguously and aligned. When iterating over nodes, ensure the loop order corresponds to memory order (usually z inside y inside x).
        
    - Use **cache blocking** if necessary: for example, if updating grid and then particles, data might ping-pong in and out of cache. We can consider processing smaller chunks that fit in L2 cache. However, given the moderate size of a tile (maybe 16^3 nodes, which is 4096 nodes; times data size ~ maybe 64 bytes each = 256KB), a tile’s grid might fit in L2 of modern CPUs. If tile is bigger, consider splitting computation in slices.
        
    
    Also take note of hardware prefetchers: On Intel, adjacent cache line prefetcher might fetch the next line whenever you access one line. If a structure spans 2 lines (which it will if >64 bytes), sequential access is fine. But if we have an array of structs and we often use only one part of each struct, it might be better to split the struct (SoA layout already does this for particles vs grid).  
    The stackoverflow comment by Emma (2022) points out that even aligned data can suffer false sharing because hardware prefetch might pull in the next line as well[stackoverflow.com](https://stackoverflow.com/questions/22766191/what-is-false-sharing-how-to-reproduce-avoid-it#:~:text=Emma%20Jane%20Bonestell%20Over%20a,year%20ago). Therefore, to truly avoid contention, sometimes aligning on 128 bytes (two cache lines) is used. We might consider aligning critical arrays on 128-byte boundaries to thwart spatial prefetch interference[stackoverflow.com](https://stackoverflow.com/questions/22766191/what-is-false-sharing-how-to-reproduce-avoid-it#:~:text=Emma%20Jane%20Bonestell%20Over%20a,year%20ago).
    
     
    
    **Implementation Steps:**
    
    - Define `ZX_CACHELINE_SIZE` (64) and use it to align global or static structures. For dynamic arrays (like particle arrays), allocate memory with alignment. If using C++17, `std::aligned_alloc(64, size)` can be used.
        
    - In OpenMP loops, use the `schedule(static, chunk_size)` to ensure each thread works on a contiguous chunk of data, reducing jumping around memory.
        
    - Place padding in structs that are shared among threads. E.g. if we have a global stats struct updated by threads, align it or give each thread a separate instance to aggregate into.
        
    - Test false sharing scenarios: a quick test is to run with threads and see if throughput scales linearly. If not, and if CPU is not fully utilized, possibly cache line contention is happening. We can also intentionally add padding and see if it improves (thereby confirming false sharing was an issue).
        
    - Use tools like Intel VTune to check cache miss rates. Ensure that memory bandwidth usage is within reason. If we saturate memory (which can be ~50 GB/s on a CPU), that becomes our limit. E.g., if each particle updates 8 nodes with 32 bytes each (mass+momentum) = 256 bytes, 1e6 particles -> 256 MB of data, which at 60 Hz is 15 GB/s, within limits but significant. Optimizing memory access (by tiling, sorting, or SIMD gather) is critical to keep under bandwidth.
        
    
    These optimizations are low-level but crucial for making the CPU solver viable for real-time or large-scale uses. We prioritize avoiding _false sharing_ – for instance, by aligning tile data so that two threads writing different tiles never contend on the same cache line[stackoverflow.com](https://stackoverflow.com/questions/22766191/what-is-false-sharing-how-to-reproduce-avoid-it#:~:text=Emma%20Jane%20Bonestell%20Over%20a,year%20ago) – and improving _spatial locality_, so that when one particle’s data is loaded, adjacent memory likely needed next is already in cache.
    
19. **Memory Allocator for GPU Buffers:** Implement a suballocation strategy for GPU memory to efficiently manage buffers for particles and tiles[bruce-lee-ly.medium.com](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=)[bruce-lee-ly.medium.com](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=The%20advantage%20of%20the%20buddy,be%20a%20waste%20of%20memory). Rather than creating a separate `ID3D12Resource` or `VkDeviceMemory` for each tile or each buffer, we allocate large memory blocks and carve them up. This avoids overhead and fragmentation over time:
    
    - Use a **Buddy Allocator** for tile grids: for example, allocate a 256 MB heap upfront. When a tile of size X is needed, round up to the next power of two (buddy system) and allocate that chunk. Keep track of free blocks. When a tile is destroyed, free its block and coalesce buddies if both free[bruce-lee-ly.medium.com](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=continue%20checking%20upward%20until%20it,).
        
    - Use linear or **ring allocators** for transient buffers that are reused every frame, such as those for copy or staging. For instance, if we need a staging buffer for readback of a small amount of data each frame, we can reuse a single buffer by rotating the region we write to (so GPU writes to region for frame N while CPU reads back region from frame N-1).
        
    - Also manage descriptor allocations similarly: instead of creating descriptor heaps frequently, create one big heap and allocate descriptor handles as needed for each resource. Our descriptor indices can be offset into this heap.
        
    
    For Vulkan, using **VMA (VulkanMemoryAllocator)** can simplify this: we can instruct it to allocate from large blocks and it will do fragmentation management. For D3D12, we implement our own (since D3D12 is lower-level in resource allocation).
    
    - We might maintain separate pools for buffers of different types (GPU local vs upload vs readback).
        
    - For instance, particle buffers are GPU-read/CPU-write (upload) or GPU-read/CPU-read (if we read back states). Grid buffers are GPU-only (DEFAULT heap).
        
    - The buddy allocator works well for variable sizes; our tile grids might all be the same size if each tile is same dimension, so a simpler fixed-size pool could work too (like keep an array of free indices for tile buffers, since each tile’s buffer is B^3 nodes * nodeSize bytes – if B is fixed, tile buffer size is fixed, making it trivial).
        
    
    **Algorithm Steps:**
    
    - On context creation, allocate a big chunk of GPU memory (e.g. `ID3D12Heap` of size say 256 MB for grid data, another for particle data if needed). Or allocate a big `ID3D12Resource` buffer that is large enough to hold max particles and tile nodes, and we will suballocate within it (this can be simpler: a large buffer resource that we treat like memory).
        
    - Implement functions: `Allocate(size, alignment)` and `Free(offset, size)`. For buddy system, maintain an array of free lists by size power-of-two. When allocating, find the smallest free block ≥ size, split until block size ~ size, return offset. Mark blocks as used. For free, mark and try to merge with buddy.
        
    - If memory runs out, we can allocate a second heap (buddy can also be extended by treating it as separate pool).
        
    - Use this for tile grid buffers: each tile activation calls Allocate(tileBufferSize). Deactivation calls Free.
        
    - For particle buffers, since max particle count is likely fixed or slowly changing, we can allocate one big buffer for all particles (structured buffer of capacity N_max). We then don’t need to allocate per particle; we just manage an integer `particle_count`. If particles are added/removed, we recycle indices (like free list for dead particles). This is more efficient than frequent alloc/free. So for particles, design the data structure to handle holes or out-of-order data (maybe swap removed particle with last particle, etc. to keep array dense).
        
    - Test memory usage patterns: allocate and free many tiles in random order to see if buddy system fragments (it shouldn’t badly, by design). Also, ensure alignment requirements (D3D12 often requires 256-byte alignment for buffers to be bound as UAV). Our buddy allocator should always return aligned addresses (buddy system inherently works with power of 2 sizes, thus naturally aligned to that size, satisfying alignment).
        
    - On Vulkan, if not using VMA, implement similarly with `VkDeviceMemory`. Or just use VMA’s pooling flag for simplicity, which will allocate big blocks and suballocate.
        
    
    This custom allocator prevents bogging down the driver with hundreds of small allocations each time a tile is created. It also means less OS-level memory management overhead. By reducing fragmentation, we ensure that even after many activate/deactivate cycles, we don’t run out of memory due to fragmentation. Essentially, we trade a bit of memory slack (due to buddy rounding up) for speed and predictability.
    
20. **Pinned Memory & Async Transfers:** Optimize CPU-GPU data transfers by using **pinned (page-locked) memory** for any staging buffers. Pinned memory (also called page-locked host memory) allows direct DMA transfers by the GPU without an extra copy by the OS[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31). In CUDA terms, `cudaHostAlloc` yields pinned memory which drastically increases transfer speed by avoiding paging[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31). In our context:
    
    - Use an UPLOAD heap in D3D12 for buffers that CPU writes into for GPU to read (e.g. initial particle data upload, or constant buffers each frame). The UPLOAD heap is backed by pinned system memory. This way, when the GPU reads from it, it can do so via PCIe directly.
        
    - Use a READBACK heap for buffers that GPU writes and CPU reads (like our timestamp readback or any debug data). These are also page-locked.
        
    - Avoid using default pageable memory for intermediate buffers. For example, if we need to update 100 KB of constants per frame, allocate a constant buffer in an upload heap at init and reuse it, rather than allocating new memory each time.
        
    - When reading back terrain data to CPU (if ever needed), do so via a buffer in a readback heap to skip the extra copy. In D3D12, copying from default heap to a readback heap buffer then mapping it gives pinned memory access on CPU.
        
    
    Another benefit: using pinned memory enables truly asynchronous transfers overlapping with computation. For instance, we can initiate a `ID3D12CommandQueue->CopyResource` from a default heap to an upload (or readback) heap resource. The CPU can map the upload resource even while GPU writes to it (with proper fencing). The transfer will use DMA and not stall the GPU if bandwidth is available.
    
     
    
    We should note that excessive pinned memory can reduce overall system performance (it reduces available pageable RAM for other apps). So we allocate only what we need. But e.g. a 256 MB upload heap is fine on a system with many GB of RAM.
    
     
    
    **Implementation Steps:**
    
    - Allocate upload buffers with `D3D12_HEAP_TYPE_UPLOAD` (in D3D12, resources in this heap are automatically considered pinned). In Vulkan, use memory with `HOST_VISIBLE | HOST_COHERENT` flags (which typically ends up pinned).
        
    - Examples: the large particle buffer could be in an UPLOAD heap if we plan to update it from CPU often. However, usually we update only small constants, not all particles. So we might keep particle buffer in GPU memory and only occasionally map it (inefficient to map GPU memory because it requires the driver to do staging anyway if not pinned).
        
    - Instead, for dynamic data we plan to stream, use upload heaps. For static data, use default heaps (GPU-only).
        
    - Ensure to align constant buffers to 256 bytes as D3D12 requires. Possibly allocate a single big upload buffer and suballocate for various constant data each frame (similar to how we did for GPU memory, but CPU accessible).
        
    - Asynchronous copies: if we have large data to transfer (e.g. initial terrain height map into particles or so), initiate it as a non-blocking copy and use a fence to sync. For instance, to upload initial particles, one could map an upload buffer, fill it, then issue a GPU copy from upload to default buffer. The CPU can continue to do other work (like initializing other structures) while the DMA happens.
        
    - Document the performance difference: pinned host memory enables the GPU to use DMA directly[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31), thus a copy that would normally require CPU involvement (pageable) now is handled by the GPU’s bus master. This can roughly double throughput for large transfers and avoids stalling the CPU.
        
    
    By leveraging pinned memory for the few places we do host-device transfers (such as texture output readback or initial data loads), we keep the pipeline smooth. In a networking scenario, one could even map an upload buffer and directly write network-received data into it for the GPU to consume, avoiding an extra copy. Overall, pinned memory usage aligns with best practices for high-performance GPU apps: use it for frequent transfers to save the extra copy through a temporary, since non-pinned (pageable) memory would force the driver to copy into a page-locked buffer internally[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31).
    
21. **Profiling & Instrumentation:** Integrate fine-grained profiling to guide optimization. We want to measure how long each phase of the simulation takes on both GPU and CPU:
    
    - On GPU, we utilize timestamp queries around each major dispatch. In D3D12, create a timestamp query heap (with, say, 16 slots). Place `EndQuery` before and after P2G, grid, contact, G2P, etc.[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=CUDA%20Driver%20checks%2C%20if%20the,page%20copy). Resolve them to a readback buffer. The CPU can then retrieve the timestamp values (in ticks) next frame and convert to milliseconds using the GPU timestamp frequency (queried via `ID3D12CommandQueue::GetTimestampFrequency`). Fill these durations into `zx_counters.P2G_time`, `grid_time`, etc.[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31).
        
    - On CPU, we can use high-resolution timers (std::chrono or platform-specific like `QueryPerformanceCounter`). Surround sections of code and accumulate times. Since the simulation loop is short, ensure to use a high-res clock. Store these in `zx_counters` as well (there might be fields for `cpu_p2g_time`, etc., or we repurpose the same).
        
    - Also track number of particles, active tiles each frame in counters (these are cheap to compute). This helps us see scaling behavior in logs.
        
    - Use conditional compilation or a runtime flag for profiling so we can disable it in shipping (though the overhead of a few queries is small).
        
    - Instrumentation: use markers for external tools. In D3D12, we can use `PIXBeginEvent`/`PIXEndEvent` around groups of commands[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=,locked%20through%20cudaMallocHost%20or%20cudaHostAlloc). In Vulkan, use `vkCmdBeginDebugUtilsLabelEXT`. This allows GPU debugging tools (PIX, RenderDoc) to show our passes with friendly names (“P2G”, “G2P”, etc.).
        
    - Logging invariants: at end of each frame on CPU, compute things like total mass or energy if possible and log if it drifts significantly (mass should stay roughly constant minus outflow). We can compute total mass by summing particle masses (should remain constant) and store it in `zx_counters.mass` for debugging[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=,explicit%20formulation). Similarly, track max penetration depth of any particle below ground, or number of contact nodes saturated (as we counted in Task 7).
        
    
    This instrumentation will be invaluable during development:
    
    - We can confirm that GPU timings meet targets (e.g. P2G 0.5 ms, Grid 0.4 ms, etc., totaling ~2 ms)[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=CHENFANFU%20JIANG%2C%20University%20of%20Pennsylvania,importance%20of%20regularity%2C%20MPM%20is). If one phase is slower, we know where to optimize.
        
    - On CPU, profiling might show which loop is bottlenecking (maybe memory copy in P2G). Then we can focus optimization efforts there.
        
    - The counters can be exposed to the game (e.g. Unreal plugin might draw them on screen or log them) to monitor performance in real-time.
        
    
    **Implementation Steps:**
    
    - In `zx_context`, have a `struct zx_counters` that contains timing fields (float milliseconds or long long ticks)[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31). Already defined fields like `p2g_time`, etc., can be filled.
        
    - After `zxEndFrame`, when GPU has finished, map the query readback buffer and calculate each duration: e.g. `p2g_time = (query[1] - query[0]) / freq`.
        
    - For CPU, use `auto t0 = Clock::now()` and `auto t1 = Clock::now()` around sections, accumulate or assign to counters.
        
    - Provide a function `zxGetCounters()` that returns these values for external use, or include them in the `zx_frame_end` struct.
        
    - Include PIX markers: around each dispatch, call `PIXBeginEvent(cmdList, 0, L"P2G")` (the first param is a color, say 0 for auto) and `PIXEndEvent` after. This way, if a developer runs PIX on the game, they see labeled GPU events[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=,locked%20through%20cudaMallocHost%20or%20cudaHostAlloc).
        
    - Similarly, for CPU debugging, one could integrate with Chrome Trace or similar by writing out events to a JSON file (optional).
        
    
    With these profiling tools, we will validate that our Phase 2 improvements achieve the expected performance: e.g., see that overlapping (Task 12) actually reduced frame time, or that SIMD (Task 17) lowered the CPU loop time by X%. We’ll also catch any unexpected stalls (if, say, a GPU stage is waiting on something, we’d see a gap in timeline). The invariant checks ensure the simulation stays healthy: if mass or energy blows up, we can catch it in logs rather than in uncontrolled behavior[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=develop%20the%20MPM%2C%20with%20funding,known%20as%20CRAMP).
    
22. **Cross-GPU and Vendor Optimization:** Ensure the system runs efficiently on different GPU brands (NVIDIA vs AMD vs Intel). Different vendors have varying strengths:
    
    - **NVIDIA (CUDA architecture):** Warps of 32, very fast context switching, good at handling many small thread groups. Shared memory is 96KB per block on latest (Ampere) for 32 threads per warp. Atomics are generally fast on Nvidia, but still can serialize.
        
    - **AMD (GCN/RDNA architecture):** Wavefronts of 64 threads, potentially higher ALU throughput but can suffer if control divergence is high (as it affects 64 threads). Shared memory (LDS) typically 64KB. AMD cards often rely on occupancy to hide latency; large thread groups (e.g. 256-512) are common. Also, memory access patterns that are coalesced in 128-byte segments are ideal.
        
    - **Intel GPUs:** (Xe-HPG) have 8 or 16-wide subgroups and a focus on simd8 and simd16 compilation. They benefit from coherent memory access and have good atomics too.
        
    
    To optimize across:
    
    - Use **subgroup operations** in a cross-vendor way: HLSL’s Wave ops will map to NVIDIA warps (size 32) and AMD wavefronts (size 64). We must ensure any algorithm using these (like warp-sum in P2G) still works if the wave size differs. For instance, Gao’s method that assumed warp=32 might need slight modification for wave64. One can query `WaveGetLaneCount()` in HLSL to adapt. Alternatively, limit ourselves to algorithms that work regardless of exact wave size (like prefix sums which work in multiples).
        
    - Avoid using `GroupShared` memory size beyond AMD’s limit. We have B^3 node values in shared mem for P2G; if B=16, that’s 4096 nodes * maybe 16 bytes each = 64KB, which is exactly AMD’s typical max. We should be careful if we consider larger tiles – maybe stick to B=16 or less for GPU.
        
    - Test on AMD: They often have lower clock but more cores, so if our code is heavy on single-thread work (like one thread per tile), it might underutilize an AMD GPU. Instead, launching more threads (even if some do redundant work) could help. We may experiment with splitting work differently for AMD if needed. For example, if tile count is low but particle count high, launch by particle or by smaller tile subdivisions to get more work groups.
        
    - Shader compiler differences: run with Vulkan on AMD and D3D12 on NVIDIA and compare performance counters. We might find, for instance, that AMD benefits from using matrix operations differently or that memory bound sections need reordering.
        
    - Pay attention to atomic throughput: a research noted that on heavy loads, even though AMD can match NVIDIA at light loads, at heavier contention NVIDIA performed better[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=Comparison%20with%20GPU%20implementation%20with,that%20we%20did%20not%20find). Our algorithm tries to minimize atomics, which should benefit both. If any remain (like a few atomic adds to shared memory), they should be okay.
        
    - Use vendor-specific profilers: Nsight for NVIDIA, Radeon GPU Profiler for AMD. Check for issues like divergent branches (should be minimal in our compute kernels as they’re mostly straight-line per node/particle).
        
    - Check half-precision or tensor core usage: Possibly not relevant to our workloads (we need precision for physics, so float32 is minimum). But if, say, we detect that using FP16 for some computations (like weight factors) is acceptable, we could consider it. However, determinism and precision might suffer, so likely skip.
        
    - Keep an eye on occupancy: for wave64, our thread group size of 256 means 4 wavefronts per group, which is okay. On wave32 (NVIDIA), 256 means 8 warps. These are reasonable. We might even increase group size to 512 on AMD if it helps occupancy, but it could hurt NVIDIA (which likes <=1024 threads per SM but splitting them in warps of 32).
        
    - If implementing multi-queue, note that NVIDIA and AMD handle async differently (NVIDIA can overlap but gains are sometimes limited, AMD async often shows more benefit when graphics pipeline is busy[forums.anandtech.com](https://forums.anandtech.com/threads/amd-vs-nvidia-asynchronous-compute-performance.2504980/#:~:text=AnandTech%20forums,at%20least%20on)).
        
    
    In summary, maintain a balance and test. We aim that our kernels are not tuned with assumptions that break on other vendors (like not assuming warp = 32 explicitly). Use HLSL intrinsics that are well-supported (avoid NV-specific ones like `WaveMatch` that may not exist on AMD’s older DX12). If necessary, have fallback paths (e.g. an alternate method for reduction if `WaveActiveSum` isn’t available or behaves differently).
    
     
    
    **Implementation Steps:**
    
    - Incorporate a GPU info query in `zxBindDevice` to detect vendor/architecture (if available via DX12 ID3D12Device->GetAdapterLuid and querying the adapter, or in Vulkan via `vkGetPhysicalDeviceProperties`). This can be used to toggle certain shader defines (like `#define WAVE_SIZE 64` for AMD, 32 for NV).
        
    - Test on an AMD GPU: run sample scenarios, profile, ensure no validation layer warnings. Pay attention to VRAM usage as well – AMD’s drivers might be more sensitive to how memory is allocated and aliased. Use our counters to see if any phase is unusually slow on AMD compared to NVIDIA for similar work – that could indicate a need to adjust (e.g. maybe AMD benefits from combining certain passes to reduce memory traffic).
        
    - Do the same on Intel if possible (Intel’s GPU architecture is quite different; if something runs slow there, it might hint at memory access issues).
        
    - Document any workarounds: e.g., “On AMD GCN cards, using wave32 instead of wave64 by forcing `WaveOpsIncludeHalfLanes` improved performance by X.” or “On NVIDIA, the contact kernel was memory bound, so we did Y.”
        
    
    By proactively addressing cross-vendor performance, we ensure ziXn is not locked to a single GPU type and behaves consistently. This also reinforces determinism: if we get identical results on different GPUs, it’s a good sign we have no hidden race conditions (because GPU float differences are usually minor if the algorithm is order-consistent).
    
23. **Scalability & Frame-Budget Tuning:** Provide mechanisms to scale the simulation to meet performance targets on various hardware. Not all users will have a high-end GPU or CPU, so we need quality levels:
    
    - **Particle density LOD:** We can emit fewer particles in distant or less important areas. For example, near the camera or player path, use full resolution material (e.g. 4 particles per cell), but far away, reduce to 1 per cell or even use a simplified heightmap deformation only. We can implement this by varying the emission in `zx_emit_particles` based on distance, or culling particles beyond a range (or compressing them into fewer “super-particles”).
        
    - **Tile activation limit:** If too many tiles become active (exceeding budget), decide on a strategy: we could stop activating new tiles (cap the simulation area) and treat new areas as rigid/un-deformable until others deactivate[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=develop%20the%20MPM%2C%20with%20funding,known%20as%20CRAMP). Or we could remove the furthest tile (least important) to free budget for a new one near the player (this is a form of cache).
        
    - **Simplified physics mode:** Provide a toggle to use a cheaper solver for less critical cases. For example, if the ground is very hard, we could skip the pressure solve or limit substeps. Or use an alternative solver like Projective Dynamics with fewer iterations for large time steps (though not in current scope).
        
    - **Adaptive substepping:** If frame time is too high, dynamically reduce the number of substeps. ZiXn might normally use e.g. 4 substeps for stability; on low-end, we might drop to 2 and accept some loss in accuracy. Monitor the actual simulation time via counters; if it consistently exceeds a threshold (say 10 ms), lower the quality.
        
    - **Multi-threading scaling:** On CPU, if the engine indicates only 4 cores available (maybe console or reserved cores), we might limit our thread usage. Or allow the engine to pass a `feature_bits` flag to disable CPU physics entirely and rely on GPU.
        
    - **Memory vs performance trade:** On consoles or lower GPUs, memory might be constrained. We could allow a smaller max particle count or tile count on those platforms (freeing memory and also speeding up simulation inherently). Expose such constants via config (ini or engine integration).
        
    
    Exposing these controls:
    
    - In Unreal, perhaps tie into the scalability settings (like “Effects Quality” could adjust ziXn particle density).
        
    - In Unity, provide a scriptable setting for quality level.
        
    
    Also implement failsafes: if at any frame we detect we exceeded budget (e.g. GPU time > target), print a warning and automatically reduce something. For instance, if active_tiles > 1000 (way too many), start deactivating far tiles or merge them:
    
    - Merging tiles: if adjacent tiles are both sparsely filled, we could combine them into one larger tile to reduce overhead (though our architecture isn’t built for variable tile sizes, so this is complex).
        
    - Hard limits: define MAX_PARTICLES (maybe 1e6). If an emitter tries to add beyond that, either reject new particles or remove oldest. This prevents runaway memory usage.
        
    - Provide feedback to game: maybe a callback when quality is reduced (so design can notify player or just log).
        
    
    **Algorithm Steps:**
    
    - Add fields in `zx_simulation_params` or context for max tiles, max particles, quality level, etc.
        
    - During tile activation (Task 1), if we’re at the limit (e.g. active_tiles == max_tiles), decide which tile to drop: possibly find the tile farthest from any player or camera that is currently not heavily deformed, write its deformation into the heightmap (so it ‘freezes’ in place), then remove it. This prevents new tiles from spawning uncontrolled[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=develop%20the%20MPM%2C%20with%20funding,known%20as%20CRAMP).
        
    - During particle emission or motion, if particle count > max_particles, either stop emitting or cull some particles (maybe merge small clusters into one with bigger mass).
        
    - Allow the engine to request a specific quality: e.g. a function `zx_set_quality(ctx, level)` that sets internal factors like substeps (lower quality -> 1 substep), pressure iterations (maybe skip on low), particle multiplier (maybe 50% of normal on medium).
        
    - Use profiling data (Task 21) at runtime to adjust: e.g., if last 10 frames average > target ms, reduce quality by one notch until within target (dynamic scalability).
        
    - Ensure determinism when scaling changes: ideally, don’t change quality during a critical network session unless all clients do the same. So dynamic adjustments might be single-player only or require host authority to broadcast.
        
    
    This ensures ziXn can run on a spectrum of hardware: high-end PC can use ultra settings (dense particles, multiple substeps for high precision, full fluid solve), whereas a low-end GPU or a switch to CPU fallback can use coarse settings (sparser particles, maybe no fluid solve, etc.) while still providing deformation albeit less detailed. By designing these now, we avoid scrambling later if performance is an issue – we have knobs to turn.
    

### Engine Integration & Application

24. **Unreal Engine 5 Plugin Integration:** Develop an Unreal Engine 5 plugin module that integrates ziXn’s C API into Unreal’s ecosystem[dev.epicgames.com](https://dev.epicgames.com/documentation/en-us/unreal-engine/runtime-virtual-texturing-in-unreal-engine#:~:text=A%20Runtime%20Virtual%20Texture%20,to%20conform%20to%20the%20terrain). The tasks include:
    
    - Creating a UZiXnSystem UObject or UActorComponent that will manage a ziXn context.
        
    - On game start or level load, initialize ziXn (e.g. `zxCreateContext` and `zxBindDevice` using D3D12 with the `ID3D12Device` from Unreal’s RHI). Unreal’s rendering RHI (DirectX12) can be accessed via `GRHIDevice` or similar, but integration may allow us to supply our own D3D12 device. However, better is to use Unreal’s device to avoid duplication. There is an Unreal extension to get the native ID3D12Device pointer.
        
    - Create a ziXn scene (`zxCreateScene`) for the landscape. Set the scene’s material properties corresponding to the Unreal material (friction, etc.).
        
    - Each tick (or rather, each physics substep), gather inputs: e.g. vehicle wheel transforms, foot locations. Unreal’s physics or animation can provide these. Feed them to ziXn as rigids: e.g. for each wheel, call `zxRigSetTransform(rig, position, orientation)` and perhaps `zxRigApplyForce` if needed.
        
    - Call `zxBeginFrame` with delta time, then `zxSimulate`, then `zxEndFrame`. This will queue GPU work.
        
    - After simulation, retrieve the output terrain data. For Unreal’s **Virtual Heightfield Mesh (VHFM)** system, we likely have a runtime virtual texture that stores height and possibly a normal map or mask[forums.unrealengine.com](https://forums.unrealengine.com/t/rvt-virtual-heightfield-mesh-multiple-questions-issues/548820#:~:text=RVT%20%2B%20Virtual%20Heightfield%20mesh,according%20to%20youtube%20tutorials%2C). We want ziXn to update this texture. Ideally, we do this on the GPU: if we can get the pointer or UAV of the RVT texture, we can direct ziXn’s writeback to it. In Unreal, one can create a UAV for an RVT by using the `FRHITexture` pointer.
        
    - If direct writing is complicated, an alternative is to copy ziXn’s output to a system memory buffer and then use `VirtualHeightfieldMesh->EditHeight()` function to update it. But that may be slow, better to stay GPU.
        
    - Unreal’s RVT is updated by marking pages dirty. We might call an Unreal API to notify that certain regions changed, so it refreshes those tiles in the virtual texture.
        
    - Additionally, integrate with **Niagara** (Unreal’s particle system) if needed: e.g. expose the ground velocity field to VFX, so dust can be emitted where ground is disturbed. This can be done by writing the velocity field to a texture or buffer that Niagara reads.
        
    - Networking: if Unreal is running in network mode, ensure to replicate necessary data (foot impacts, etc.) to remote machines because ziXn has to run deterministically. Possibly run ziXn only on server and broadcast resulting deformation to clients (cheaper, since clients would just get the heightmap).
        
    - Provide editor integration: perhaps a custom detail panel for the component to set material properties (snow vs mud), and tools to paint initial material (like where is mud vs rock).
        
    - Performance: Use Unreal’s RHI thread or GPU compute queues properly. We might register a callback with the renderer to execute our compute after the main frame’s rendering of the previous frame (to hide latency). Unreal’s `ENQUEUE_RENDER_COMMAND` can be used to push a lambda that calls ziXn simulation on the render thread, ensuring it overlaps with game thread.
        
    
    **Algorithm Steps (per frame in-game):**
    
    - On game tick (e.g. UZiXnSystem::TickComponent):
        
        1. Update `zx_rig` inputs from Unreal actors (wheels, etc.).
            
        2. Call `zxBeginFrame(scene, simParams)`.
            
        3. Dispatch ziXn simulation – possibly on Unreal’s RHI thread via `ENQUEUE_RENDER_COMMAND` so it runs in parallel with game thread. This will call `zxSimulate` which kicks GPU work.
            
        4. After GPU fence signals (we wait or use Unreal’s `FRHICommandListImmediate::WriteGPUFence`), copy or use output textures. If using RVT, ensure the GPU work is done before rendering that mesh.
            
        5. Mark the VirtualHeightfieldMesh to update. Possibly call `UVirtualHeightfieldMeshComponent::ForceRefresh()` or update its virtual texture.
            
    - The deformed terrain will then be rendered by the VHFM. It uses a material that samples the RVT for height. We should see tracks and craters appear.
        
    - Also update any gameplay things: e.g. if we want physical foot sinking, we might adjust character position based on the current heightfield.
        
    
    Test in Unreal Editor: create a landscape or VHFM, attach our ZiXn component, set a material (snow, etc.), add a vehicle or character with a rig component. As the vehicle moves, we expect the ground to deform. Use Unreal’s profiling (stat GPU) to ensure our compute doesn’t break budgets. We also ensure that in PIE (Play In Editor) with multiple clients, either only server runs physics (and replicates heightmap) or if each client runs, results stay in sync (this is tricky – likely authoritative on server).
    
25. **Unity Engine Plugin & C# API:** Create a Unity plugin (probably as a native DLL) to allow using ziXn in Unity projects[reddit.com](https://www.reddit.com/r/unrealengine/comments/qewicz/runtime_virtual_heightfield_mesh_help_desperately/#:~:text=Runtime%20virtual%20heightfield%20mesh%20help,so%20much%20in%20practical%20use). Steps:
    
    - Unity uses C# scripts, so we expose C API via PInvoke. Write a C# class `ZiXnPlugin` with `[DllImport]` for functions like `zxCreateContext`, `zxCreateScene`, etc.
        
    - Unity’s rendering is often DirectX11/12 on Windows. We can attempt to get the D3D11 device pointer via `UnityEngine.SystemInfo` – but Unity’s low-level Native Plugin interface might be needed to tie into GPU command execution. Alternatively, we run ziXn in its own D3D12 device and then share textures.
        
    - The simplest path: do CPU simulation only in Unity (not ideal performance) by using the CPU fallback. But that loses the visual fidelity.
        
    - For GPU, Unity offers a `GraphicsBuffer` or `ComputeBuffer` that can be shared with native code. We could create a ComputeBuffer in Unity for the heightmap (as RWStructuredBuffer) and pass its pointer to the plugin (Unity has functions to get native buffer pointer in DX11).
        
    - Another approach: use Unity’s command buffer API to inject our own compute shaders. However, converting all ziXn HLSL to Unity compute shaders would be major work and lose determinism.
        
    - Instead, maybe run ziXn entirely on CPU in Unity for now (depending on target usage).
        
    - Focus on output: Unity’s Terrain system expects a heightmap array. We can get the terrain’s `TerrainData.heightmapResolution` and directly modify the height array via `TerrainData.SetHeights` (this is heavy if done every frame for large terrains). A better approach is Unity’s new Mesh-based terrain or VFX with compute. Possibly Unity’s Visual Effect Graph could ingest a height texture.
        
    - Another possibility: use a RenderTexture for height, update it via native plugin. Unity allows native code to register a callback at the end of frame to write into a RenderTexture using low-level graphics API (e.g. DirectX12, one could open a handle to the D3D12 resource).
        
    - For simplicity, a route is: have ziXn write the heightfield to a CPU array each frame (since Unity Terrain’s API is CPU-based). If resolution is not too high (say 256x256), this might be acceptable. For better performance, one might use `ComputeBuffer.SetData` to push a float array to a compute buffer and a Unity shader to apply it to terrain mesh.
        
    
    **Algorithm Steps:**
    
    - Write a C++ DLL that wraps ziXn and provides functions Unity can call: e.g. `InitZiXn(width, height)` to create context and scene, `SimulateZiXn(dt)` to step, `GetHeightData(Float[] heights)` to copy out the height field.
        
    - In Unity C#, create a MonoBehaviour script that holds a reference to a Terrain or mesh. On Start, call plugin init (with appropriate terrain size).
        
    - Each frame, gather rig inputs: Unity’s wheel colliders or character positions. Pass them to plugin (e.g. call `ApplyRigForce(int id, float x, y, z)` for each).
        
    - Call `SimulateZiXn(Time.deltaTime)`. This could run GPU internally; if so, we might need to wait for the GPU fence internally and copy the height data to CPU (the plugin can do that).
        
    - After simulate, call `GetHeightData(array)`. Apply this to Unity’s terrain: either `terrainData.SetHeights(0,0,array)` or update a mesh’s vertices.
        
    - Alternatively, instead of copying heights CPU, we could create a Unity `Texture2D` asset for height and use `Texture2D.LoadRawTextureData` with the data, then `Apply()` to update it.
        
    
    Test in a Unity scene: a simple car driving on a terrain. We expect the terrain to deform where wheels go. Unity’s Terrain system might not update collider immediately when heights change each frame, so we might need to force collider recompute (could be expensive). Possibly use an alternative approach for collision (like have wheels always at terrain height via raycast into our height array).
    
     
    
    Performance might be a challenge with naive approach (since Unity Terrain wasn’t designed for per-frame edits across large areas). But as a prototype, it would demonstrate integration. For a more optimized Unity integration, one might bypass Terrain and use a custom mesh or shader-based displacement using the texture.
    
26. **Consistent Terrain Writeback & Scaling:** As part of engine integration, ensure that ziXn’s coordinate system and output format align with the engine’s terrain coordinate system[dev.epicgames.com](https://dev.epicgames.com/documentation/en-us/unreal-engine/runtime-virtual-texturing-in-unreal-engine#:~:text=A%20Runtime%20Virtual%20Texture%20,to%20conform%20to%20the%20terrain). Steps:
    
    - If Unreal’s landscape is scaled or offset, we must apply the same transforms to ziXn’s output. For example, if the landscape component has a scale of (100,100,1) (common if each vertex = 100 cm spacing), then our height displacements need to be scaled accordingly. We can store a conversion factor (Unreal units to ziXn units).
        
    - Set ziXn scene origin to match engine world origin (or landscape actor origin). ZiXn outputs world positions for particles; to get terrain height in engine coordinates, we ensure we started with same reference.
        
    - Output format: Unreal expects height in a texture, likely as a normalized value or a real float (for RVT, we can use R32F which matches our float). Unity expects either a float array [0..1] or actual meters. We must adjust accordingly. E.g., if ziXn says ground moved -0.2 m, in Unity TerrainData you might set a new height value = baseHeight + (-0.2) * (1/terrainHeight) in [0..1] normalized.
        
    - If multiple tiles in ziXn correspond to one continuous terrain, ensure tile seams do not produce cracks in the engine. The ghost region blending should guarantee continuity in height. We should output the final unified heightfield (if using a single texture). If engine’s terrain LOD or streaming could cause misalignment, we might need to update multiple mips or lock the highest detail.
        
    - Mipmaps: For Unreal RVT, we might generate coarse mips automatically by writing to the texture. The VHFM might sample lower mips at distance. Check that the deformation appears at distance without gaps – if needed, we might have to recalc those mips or rely on RVT’s built-in.
        
    - Tiling: If the terrain is larger than simulation domain, consider how to handle edges. Possibly beyond active tiles, we leave height unchanged (so heightmap is static outside).
        
    - Ensure that height updates don’t accumulate error: e.g., if we keep applying small changes, floating-point error could accumulate. It might be better to track an absolute heightfield state rather than incremental delta. ZiXn provides displacement per frame, but we can maintain a persistent heightmap that we update additively. Or if ziXn provides absolute height (if we initialised it with the base terrain shape), then we directly set that each frame which is more stable.
        
    
    For scaling:
    
    - If the engine’s frame rate is dynamic (e.g. Unity might drop to 30 FPS on slow device), ziXn’s `dt` should correspond, which it does by design.
        
    - If running physics at a fixed rate (like Unreal often runs physics at 60 Hz independent of frame), we might call ziXn multiple times per frame if needed (substepping).
        
    - Provide a user adjustment for strength of deformation. E.g. if one wants exaggerated tracks, multiply the output displacement by a factor before applying to engine terrain. This is not physically accurate but an artistic control.
        
    
    After integration, validate with known shapes: if we drop a known object that creates a 10 cm deep imprint in ziXn, measure in engine that it is indeed 10 cm. If not, apply a scale factor until it matches (and figure out why: likely unit conversion issues).
    
27. **Robust API Error Handling:** Improve robustness by validating inputs and states in the ziXn API[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=,locked%20through%20cudaMallocHost%20or%20cudaHostAlloc). This involves:
    
    - Checking all pointers passed to API calls: if a `zxScene` or `zxContext` handle is null or invalid (we can maintain an internal registry of created scenes to verify), return `ZX_ERR_INVALID`.
        
    - Verify parameters: e.g. if user calls `zxCreateMaterial` with nonsensical parameters (negative density, etc.), print a warning and clamp or return error. We can define sane ranges for material constants (like Young’s modulus > 0, friction 0-∞).
        
    - During simulation, detect NaNs or infinities: after grid update, we can scan a few random nodes for NaN velocity (or check if any force became NaN by adding an assertion). If detected, maybe call `zxResetScene` or at least stop simulation to avoid wild behavior.
        
    - Use `ZX_ASSERT` internally for conditions that should never happen (like memory allocation failing to find a free block when we think it should).
        
    - Provide error messages: The `zx_status` type could be accompanied by a function `zxGetLastError(context)` that returns a string or code for the last error. Or simply log to stderr (in dev mode).
        
    - Ensure that destroying a scene or context multiple times doesn’t crash. The API should gracefully ignore double-destroy or return an error.
        
    - Resource release: after `zxDestroyContext`, check that all internal allocations (tiles, buffers) have been freed (maybe track a counter, and if not zero, log a leak warning).
        
    - Out-of-memory handling: If our memory allocator (Task 19) cannot fulfill a request (like creating too many tiles), have it return a null or error, and propagate that to the user (maybe through `zxSimulate` returning ZX_ERR_OOM or similar). This way, the game can know to reduce load or at least not crash. Similarly, GPU out-of-memory from D3D/Vulkan should be caught (D3D12 will fail resource creation calls; we catch and translate).
        
    - Thread-safety: if the API is called from multiple threads (not typical for our design, but just in case), document that it’s not thread-safe except certain functions like reading outputs which could be.
        
    - Testing: Write a suite of tests (maybe in `tests/` directory) that deliberately misuse the API (null inputs, large values) to ensure the error handling works and doesn’t crash.
        
    
    This robustification will make the engine integration more stable. For example, if the game accidentally sets dt=0 or extremely large (which could cause numerical blow-up), we detect and clamp dt to a reasonable max or return an error[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=developer,driver%20can%20detect%20pinned%20memory). We should also guard against extremely small tile sizes or extremely large particle counts that could overflow internal arrays (if a user set tile size such that B^3 > int max, obviously not feasible; clamp B or error out).
    
     
    
    Provide clear documentation on these checks: e.g., in the README or header, say “zxSimulate returns ZX_ERR_INVALID if dt<=0 or scene not ready” so users know.
    
28. **Invariants Checking & Determinism Validation:** Add end-of-frame checks to verify the simulation’s physical invariants and reproducibility. Some invariants:
    
    - **Mass conservation:** Sum of particle masses should remain constant (assuming no particles added/removed). After P2G, we can sum grid mass and compare to initial total mass (they should match, aside from floating error). After G2P, sum particle mass again. If discrepancy > epsilon, log it[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=,explicit%20formulation).
        
    - **Energy (optional):** Compute total energy (kinetic + potential + maybe strain energy) at start and end of frame. Energy might not be strictly conserved due to plasticity (dissipation) or friction (energy lost) – so this is tricky. But we expect non-increasing total energy in a system with dissipation. If we ever see a big increase in energy (e.g. simulation exploded numerically), that’s a red flag.
        
    - **No NaNs:** After each major step, scan a few arrays to ensure no `NaN` or `Inf`. If found, abort simulation or correct it. A common cause could be 0 mass leading to division by zero. We protect against that by small epsilon in divisions, but if a NaN slips through, catch it.
        
    - **Determinism:** To validate, run the simulation twice with the same inputs and ensure outputs are identical. We can implement a debug mode where we store a checksum of state each frame and compare to a reference. For example, compute a 64-bit CRC of all particle positions and velocities at end of frame[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y). If the CRC deviates between runs or between host and client, determinism was broken. We could expose a function `zxCalculateChecksum(scene)` that returns a hash of state (maybe just summing up some values to reduce overhead)[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y). This is useful in multiplayer: each client could compute this and compare to ensure sync (though typically, one wouldn’t trust the clients to do it).
        
    - Use high precision for checksum to avoid floating rounding differences affecting it when they are trivial. One approach: quantize particle positions to a fixed grid (e.g. multiply by 1000 and cast to int) then sum or hash, so that tiny FP differences don’t flip the checksum if within tolerance. Or use an tolerance-based comparison rather than bitwise exact.
        
    - Provide a debug mode where every arithmetic step uses double precision or a fixed seed for any random choices, to test determinism across hardware.
        
    
    We’ll include these checks under a compile flag (so release builds don’t pay the cost). In debug or a special build, one can enable “paranoid checks” which will degrade performance but catch issues quickly.
    
     
    
    For network determinism specifically, we might allow the engine to get an array of important state to send over network occasionally to resync. Or at least log if a divergence is detected.
    
     
    
    Example: at frame end, compute `checksum = f(sum_i (p_i.x + 3*p_i.y + 7*p_i.z))` (some linear combination to mix values)[congress.cimne.com](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y). Print it. Run two identical runs – if any difference in checksums appears, investigate.
    
     
    
    These measures ensure that as we optimize, we don’t break the core stability and determinism. If a certain optimization (like using FMA or changing thread ordering) introduces nondeterminism, we’ll catch it by these tests. Similarly, if the physics goes off (mass loss, etc.), we’ll know immediately rather than noticing only in results.
    
29. **Resource Limits & Graceful Degradation:** Implement hard limits on simulation size and a strategy to degrade gracefully when limits are hit. This complements Task 23, but here we enforce absolute boundaries to avoid crashing or huge slowdowns:
    
    - **Max tiles:** define a constant or configurable `ZX_MAX_TILES`. If the simulation tries to activate beyond this, we simply refuse to activate more. Possibly print a warning “Max tiles reached, further deformation outside area will be ignored.” This ensures memory usage stays bounded (each tile uses memory).
        
    - **Max particles:** similarly `ZX_MAX_PARTICLES`. If an emitter tries to create more, either stop emitting new particles or implement particle pooling (e.g. recycle oldest inactive particles).
        
    - **Max substeps:** if user sets an absurd number of substeps (like 100), maybe clamp it to a safe max (maybe 10) to avoid extremely long frames.
        
    - If memory allocator fails (as per Task 19, which can happen if we genuinely run out of GPU mem or exceed our reserved heap), catch that and deactivate some content. For instance, if we cannot allocate a new tile, we could trigger the removal of the farthest tile (like a cache eviction) to free space. This is a more advanced policy akin to virtual memory: rarely used tiles get swapped out (we could even save their deformation to a heightmap and restore it later if needed).
        
    - Provide debug stats so a user can see if they are hitting limits (e.g. a counter of “tiles skipped due to cap”).
        
    
    The idea is that instead of a sudden crash or stall, the system just stops adding new detail when overloaded. The already active simulation continues, so the player might notice edges of deformable area but the game continues to run.
    
     
    
    Implementation:
    
    - In `zxActivateTile` (or wherever new tile is spawned), check current `scene->active_tile_count`. If >= max, decide a victim tile to remove or simply abort activation. Possibly remove the tile farthest from camera (if we have that info via engine) or oldest used. One can maintain a simple FIFO or LRU for tiles.
        
    - When removing a tile for this reason, we might want to preserve its deformation: e.g. bake its deformed state into a coarse heightfield that remains as static terrain. That way, if player goes back, it’s not pristine ground. This could be complex, but maybe treat it as out-of-scope for now, instead we warn the content designer that the area is too large.
        
    - On particle emission, if > max, either don’t spawn new ones or cull some. If it’s an effect like debris, culling is okay. If it’s soil particles critical for simulation, better to not spawn than overload. Document the cap so designers can adjust emission rates.
        
    
    Testing: deliberately set very low max values and ensure engine doesn’t crash: e.g. max_tiles=10, drive out of area, check that after 10 tiles, the 11th doesn’t allocate (we might see a slight missing deformation beyond there but no crash).
    
     
    
    These limits should rarely be hit if content is designed within bounds, but they are important safety nets, especially in open-world scenarios. By logging when they happen, we can inform level designers to perhaps confine deformable regions or up the budget if hardware allows.
    
30. **Forward-Looking Extensions:** Note down tasks for future phases and ensure current design can accommodate them:
    
    - **Multi-GPU support:** Possibly running different regions of the terrain on different GPUs or splitting particle sets[yzhu.io](https://yzhu.io/publication/mpmgpu2020siggraph/paper.pdf#:~:text=Method%20yzhu,reorders%20pipeline%20by%20moving). Our tile-based approach could allow domain decomposition: assign half the tiles to GPU0, half to GPU1. Communication would happen at tile borders (similar to CPU domain decomposition in HPC). We should keep the tile data structures flexible enough to support multiple devices or contexts. This might mean abstracting device pointers or having an API to migrate tile data between devices for overlaps. For now, just bear it in mind (no hardcoding single-GPU assumptions too deep).
        
    - **Mobile/Console support:** E.g. a Vulkan backend paves way for Linux/Android. For consoles (PS5, Xbox), their APIs are similar (PS5 uses GNMX or a thin wrapper, Xbox uses DX12). If our D3D12 backend is robust, porting to Xbox is trivial. For PlayStation, we might eventually need a separate backend.
        
    - **Regional time stepping (adaptive substeps):** Research by Fu et al. (2018) implemented regionally varying time steps[cg.cs.tsinghua.edu.cn](https://cg.cs.tsinghua.edu.cn/papers/CGF-2018-mpm.pdf#:~:text=,hybrid%20particle%2Fgrid%20nature%20of%20MPM). This is complex but could allow, for example, high-rate simulation near a fast-moving object and low-rate far away to save CPU/GPU. Our engine would need the ability to have asynchronous tile updates (some tiles update more frequently). The current design is synchronous for simplicity. But we can mention this in docs as a potential area (maybe if doing a large-scale sim, one could manually subdivide scenes).
        
    - **New material models:** e.g. coupling to vegetation (grass that bends), or adding water flow (could integrate SPH or FLIP fluid on top of MPM soil). The architecture can likely extend by adding new particle types or additional field variables. We should ensure adding an extra field (like pore pressure or temperature) is not too difficult. Possibly leave room in `zxParticle` struct for extra data or design functions to handle per-particle custom data. Forward-looking: maybe snow melting (phase change) or drying of mud (progressive hardening) – these would require time-dependent material properties.
        
    - **Editor tools:** In Unreal, one may want to paint initial deformation or material distribution. Possibly integrate with Unreal’s editor to allow stamping shapes or defining which areas are deformable.
        
    - **Machine learning integration:** An idea (for way future) could be using an ML model to emulate the physics in far distance to reduce cost, or to up-res the deformation detail beyond what physics sim did. This is speculative, but ensuring data can be extracted (for training ML) or injected (from ML output) without rewriting the core would be nice.
        
    - Keep track of latest research: e.g. 2025 papers on MPM might introduce better solvers or better GPU techniques (like the Chen 2025 method we cited[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=is%20widely%20believed%20that%20highly,from%20scattering%20without%20atomic%20op%02erations)). Maintaining a modular structure will help replace components (for instance, if a new P2G method is developed that is even faster, we can swap our P2G implementation).
        
    
    Summarize these ideas in the documentation so stakeholders know what could be added and that the current design isn’t a dead end for them. For example, multi-GPU: “The tiling approach lends itself to parallelization across GPUs[dl.acm.org](https://dl.acm.org/doi/10.1145/3386569.3392442#:~:text=A%20massively%20parallel%20and%20scalable,structure%20that%20promotes%20coalesced%20memory); a future update could assign tiles to different GPU nodes, communicating at tile borders with minimal overhead.”
    

In conclusion, Phase 2 sets up ziXn with a strong foundation of optimized, multi-platform code ready for real-world use. The tasks covered above ensure that by the end of Phase 2, ziXn can simulate rich deformable terrain in real-time within game engines, and we have a clear roadmap for further enhancements such as larger scale distribution and feature additions.

 

**Sources:** The above tasks and designs were informed by the ziXn architecture documentation[dev.epicgames.com](https://dev.epicgames.com/documentation/en-us/unreal-engine/runtime-virtual-texturing-in-unreal-engine#:~:text=A%20Runtime%20Virtual%20Texture%20,to%20conform%20to%20the%20terrain)[en.wikipedia.org](https://en.wikipedia.org/wiki/Material_point_method#:~:text=develop%20the%20MPM%2C%20with%20funding,known%20as%20CRAMP) and by state-of-the-art techniques in physics simulation and GPU computing. We referenced research on anisotropic soil modelingweb.cs.ucla.edu, frictional contact algorithms[cof.orst.edu](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=,however%2C%20the%20Nc%20definition%20still), GPU MPM optimizations[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=the%20same%20location%20are%20conventionally,cell%20idpr%20ev%20then)[pages.cs.wisc.edu](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=is%20widely%20believed%20that%20highly,from%20scattering%20without%20atomic%20op%02erations), and best practices for parallel computation (avoiding false sharing[stackoverflow.com](https://stackoverflow.com/questions/22766191/what-is-false-sharing-how-to-reproduce-avoid-it#:~:text=Emma%20Jane%20Bonestell%20Over%20a,year%20ago), using pinned memory for DMA[stackoverflow.com](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31), etc.). This ensures our implementation is grounded in proven methods and is ready for the demands of modern games and simulations.

Citations

[

https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf

](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=assemble%20the%20correct%20vertex%20data,be%20analyzed%20for%20both%20approaches)[

https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf

](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=2%20PARALLELIZATION%20STRATEGIES%202,illustrates%20the%20communication%20steps%20for)[

https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf

](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=The%20ghost%20particle%20approach%20was,broken%20into%20separate%20X%2C%20Y)[

[PDF] Cirrus: Adaptive Hybrid Particle-Grid Flow Maps on GPU

https://wang-mengdi.github.io/proj/25-cirrus/cirrus-preprint.pdf

](https://wang-mengdi.github.io/proj/25-cirrus/cirrus-preprint.pdf#:~:text=%5BPDF%5D%20Cirrus%3A%20Adaptive%20Hybrid%20Particle,bottom%20line%29%20only)[

https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf

](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=cell%20approach%2C%20where%20the%20summation,be%20analyzed%20for%20both%20approaches)[

GPU Optimization of Material Point Methods

https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf

](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=the%20same%20location%20are%20conventionally,cell%20idpr%20ev%20then)[

GPU Optimization of Material Point Methods

https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf

](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=iteration%2C%20and%20the%20total%20number,adds%20to%20the)[

GPU Optimization of Material Point Methods

https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf

](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=is%20widely%20believed%20that%20highly,from%20scattering%20without%20atomic%20op%02erations)[

![](https://www.google.com/s2/favicons?domain=https://en.wikipedia.org&sz=32)

Material point method - Wikipedia

https://en.wikipedia.org/wiki/Material_point_method

](https://en.wikipedia.org/wiki/Material_point_method#:~:text=1,the%20same%20used%20in%20FEM)[

GPU Optimization of Material Point Methods

https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf

](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=operations.%20We%20exploit%20the%20warp,5%3A%20boundary%20%E2%86%90%20t%20rue)[

![](https://www.google.com/s2/favicons?domain=https://stackoverflow.com&sz=32)

c++ - Why is CUDA pinned memory so fast? - Stack Overflow

https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast

](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=CUDA%20Driver%20checks%2C%20if%20the,page%20copy)[

![](https://www.google.com/s2/favicons?domain=https://en.wikipedia.org&sz=32)

Material point method - Wikipedia

https://en.wikipedia.org/wiki/Material_point_method

](https://en.wikipedia.org/wiki/Material_point_method#:~:text=2,extrapolated%20back%20to%20material%20points)[

![](https://www.google.com/s2/favicons?domain=https://en.wikipedia.org&sz=32)

Material point method - Wikipedia

https://en.wikipedia.org/wiki/Material_point_method

](https://en.wikipedia.org/wiki/Material_point_method#:~:text=extrapolated%20from%20the%20surrounding%20nodes,depending%20on%20%2067%20integration)[

G. Klar – Dissertation

http://web.cs.ucla.edu/~dt/theses/klar-thesis.pdf

](http://web.cs.ucla.edu/~dt/theses/klar-thesis.pdf#:~:text=We%20model%20the%20dynamics%20with,Our%20approach%20is%20able%20to)[

G. Klar – Dissertation

http://web.cs.ucla.edu/~dt/theses/klar-thesis.pdf

](http://web.cs.ucla.edu/~dt/theses/klar-thesis.pdf#:~:text=of%20material%3A%20the%20von%20Mises,Coulomb%20yield)[

![](https://www.google.com/s2/favicons?domain=https://www.sciencedirect.com&sz=32)

Nonlinear Mohr–Coulomb yield criterion: Integration in ADMM ...

https://www.sciencedirect.com/science/article/pii/S0266352X2400836X

](https://www.sciencedirect.com/science/article/pii/S0266352X2400836X#:~:text=The%20goal%20of%20the%20return,Let%20us)[

![](https://www.google.com/s2/favicons?domain=https://en.wikipedia.org&sz=32)

Material point method - Wikipedia

https://en.wikipedia.org/wiki/Material_point_method

](https://en.wikipedia.org/wiki/Material_point_method#:~:text=develop%20the%20MPM%2C%20with%20funding,known%20as%20CRAMP)[

![](https://www.google.com/s2/favicons?domain=https://yzhu.io&sz=32)

[PDF] A Moving Least Squares Material Point Method with Displacement ...

https://yzhu.io/publication/mpmmls2018siggraph/paper.pdf

](https://yzhu.io/publication/mpmmls2018siggraph/paper.pdf#:~:text=,to%20the%20adoption%20of)[

![](https://www.google.com/s2/favicons?domain=https://www.cof.orst.edu&sz=32)

https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf

](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=Given%20normal%20vector%2C%20early%20MPM,condition%20and%20the%20definition%20of)[

![](https://www.google.com/s2/favicons?domain=https://www.cof.orst.edu&sz=32)

https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf

](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=compression,1A%20shows%20a%201D%20material)[

![](https://www.google.com/s2/favicons?domain=https://www.cof.orst.edu&sz=32)

https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf

](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=,however%2C%20the%20Nc%20definition%20still)[

[PDF] Animating Fluid Sediment Mixture in Particle-Laden Flows

https://pages.cs.wisc.edu/~sifakis/papers/MPM-particle-laden-flow.pdf

](https://pages.cs.wisc.edu/~sifakis/papers/MPM-particle-laden-flow.pdf#:~:text=Flows%20pages,Tampubolon%20et%20al)[

![](https://www.google.com/s2/favicons?domain=https://arxiv.org&sz=32)

Using The Polynomial Particle-in-Cell Method for Liquid-Fabric Interaction

https://arxiv.org/html/2308.01060v2

](https://arxiv.org/html/2308.01060v2#:~:text=Jiang%20et%20al%20,numerical%20dissipation%20and%20improving%20preservation)[

![](https://www.google.com/s2/favicons?domain=https://arxiv.org&sz=32)

Using The Polynomial Particle-in-Cell Method for Liquid-Fabric Interaction

https://arxiv.org/html/2308.01060v2

](https://arxiv.org/html/2308.01060v2#:~:text=reducing%20the%20numerical%20dissipation%20of,numerical%20dissipation%20and%20improving%20preservation)[

[PDF] A material point method for snow simulation - UC Davis Mathematics

https://math.ucdavis.edu/~jteran/papers/SSCTS13.pdf

](https://math.ucdavis.edu/~jteran/papers/SSCTS13.pdf#:~:text=Mathematics%20math,Drucker%20and)[

![](https://www.google.com/s2/favicons?domain=https://dev.epicgames.com&sz=32)

Runtime Virtual Texturing in Unreal Engine | Unreal Engine 5.6 Documentation | Epic Developer Community

https://dev.epicgames.com/documentation/en-us/unreal-engine/runtime-virtual-texturing-in-unreal-engine

](https://dev.epicgames.com/documentation/en-us/unreal-engine/runtime-virtual-texturing-in-unreal-engine#:~:text=A%20Runtime%20Virtual%20Texture%20,to%20conform%20to%20the%20terrain)[

![](https://www.google.com/s2/favicons?domain=https://forums.unrealengine.com&sz=32)

RVT + Virtual Heightfield mesh - multiple questions/issues

https://forums.unrealengine.com/t/rvt-virtual-heightfield-mesh-multiple-questions-issues/548820

](https://forums.unrealengine.com/t/rvt-virtual-heightfield-mesh-multiple-questions-issues/548820#:~:text=RVT%20%2B%20Virtual%20Heightfield%20mesh,according%20to%20youtube%20tutorials%2C)[

![](https://www.google.com/s2/favicons?domain=https://stackoverflow.com&sz=32)

c++ - Why is CUDA pinned memory so fast? - Stack Overflow

https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast

](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=31)[

![](https://www.google.com/s2/favicons?domain=https://gpuopen.com&sz=32)

[PDF] Deep Dive: Asynchronous Compute - AMD GPUOpen

https://gpuopen.com/download/GDC2017-Asynchronous-Compute-Deep-Dive.pdf

](https://gpuopen.com/download/GDC2017-Asynchronous-Compute-Deep-Dive.pdf#:~:text=Async%20Compute%20%E2%86%92%20More%20Performance,CPU)[

![](https://www.google.com/s2/favicons?domain=https://developer.nvidia.com&sz=32)

Advanced API Performance: Async Compute and Overlap

https://developer.nvidia.com/blog/advanced-api-performance-async-compute-and-overlap/

](https://developer.nvidia.com/blog/advanced-api-performance-async-compute-and-overlap/#:~:text=This%20post%20covers%20best%20practices,applications%2C%20see%20all%20Advanced)[

![](https://www.google.com/s2/favicons?domain=https://stackoverflow.com&sz=32)

c++ - Why is CUDA pinned memory so fast? - Stack Overflow

https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast

](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=,locked%20through%20cudaMallocHost%20or%20cudaHostAlloc)[

![](https://www.google.com/s2/favicons?domain=https://bruce-lee-ly.medium.com&sz=32)

Nvidia GPU Memory Pool-BFC. How to pool gpu memory? | by Bruce-Lee-LY | Medium

https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82

](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=)[

![](https://www.google.com/s2/favicons?domain=https://bruce-lee-ly.medium.com&sz=32)

Nvidia GPU Memory Pool-BFC. How to pool gpu memory? | by Bruce-Lee-LY | Medium

https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82

](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=The%20advantage%20of%20the%20buddy,be%20a%20waste%20of%20memory)[

![](https://www.google.com/s2/favicons?domain=https://stackoverflow.com&sz=32)

c++ - Why is CUDA pinned memory so fast? - Stack Overflow

https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast

](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=developer,driver%20can%20detect%20pinned%20memory)[

![](https://www.google.com/s2/favicons?domain=https://en.wikipedia.org&sz=32)

Material point method - Wikipedia

https://en.wikipedia.org/wiki/Material_point_method

](https://en.wikipedia.org/wiki/Material_point_method#:~:text=,explicit%20formulation)[

![](https://www.google.com/s2/favicons?domain=https://www.cof.orst.edu&sz=32)

https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf

](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=%29%20%C2%B7%20n%CB%86%20,In%20other%20words%2C%20this%20condition)[

![](https://www.google.com/s2/favicons?domain=https://johnnysswlab.com&sz=32)

When vectorization hits the memory wall: investigating the AVX2 memory gather instruction - Johnny's Software Lab

https://johnnysswlab.com/when-vectorization-hits-the-memory-wall-investigating-the-avx2-memory-gather-instruction/

](https://johnnysswlab.com/when-vectorization-hits-the-memory-wall-investigating-the-avx2-memory-gather-instruction/#:~:text=Vectorization%20puts%20a%20large%20additional,vectorized%20way)[

![](https://www.google.com/s2/favicons?domain=https://stackoverflow.com&sz=32)

caching - What is "false sharing"? How to reproduce / avoid it? - Stack Overflow

https://stackoverflow.com/questions/22766191/what-is-false-sharing-how-to-reproduce-avoid-it

](https://stackoverflow.com/questions/22766191/what-is-false-sharing-how-to-reproduce-avoid-it#:~:text=Emma%20Jane%20Bonestell%20Over%20a,year%20ago)[

![](https://www.google.com/s2/favicons?domain=https://bruce-lee-ly.medium.com&sz=32)

Nvidia GPU Memory Pool-BFC. How to pool gpu memory? | by Bruce-Lee-LY | Medium

https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82

](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=continue%20checking%20upward%20until%20it,)[

GPU Optimization of Material Point Methods

https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf

](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=CHENFANFU%20JIANG%2C%20University%20of%20Pennsylvania,importance%20of%20regularity%2C%20MPM%20is)[

GPU Optimization of Material Point Methods

https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf

](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=Comparison%20with%20GPU%20implementation%20with,that%20we%20did%20not%20find)[

![](https://www.google.com/s2/favicons?domain=https://forums.anandtech.com&sz=32)

AMD vs NVidia asynchronous compute performance - AnandTech

https://forums.anandtech.com/threads/amd-vs-nvidia-asynchronous-compute-performance.2504980/

](https://forums.anandtech.com/threads/amd-vs-nvidia-asynchronous-compute-performance.2504980/#:~:text=AnandTech%20forums,at%20least%20on)[

![](https://www.google.com/s2/favicons?domain=https://www.reddit.com&sz=32)

Runtime virtual heightfield mesh help desperately required - Reddit

https://www.reddit.com/r/unrealengine/comments/qewicz/runtime_virtual_heightfield_mesh_help_desperately/

](https://www.reddit.com/r/unrealengine/comments/qewicz/runtime_virtual_heightfield_mesh_help_desperately/#:~:text=Runtime%20virtual%20heightfield%20mesh%20help,so%20much%20in%20practical%20use)[

![](https://www.google.com/s2/favicons?domain=https://stackoverflow.com&sz=32)

c++ - Why is CUDA pinned memory so fast? - Stack Overflow

https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast

](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=developer,driver%20can%20detect%20pinned%20memory)[

![](https://www.google.com/s2/favicons?domain=https://yzhu.io&sz=32)

[PDF] A Massively Parallel and Scalable Multi-GPU Material Point Method

https://yzhu.io/publication/mpmgpu2020siggraph/paper.pdf

](https://yzhu.io/publication/mpmgpu2020siggraph/paper.pdf#:~:text=Method%20yzhu,reorders%20pipeline%20by%20moving)[

[PDF] A Temporally Adaptive Material Point Method with Regional Time ...

https://cg.cs.tsinghua.edu.cn/papers/CGF-2018-mpm.pdf

](https://cg.cs.tsinghua.edu.cn/papers/CGF-2018-mpm.pdf#:~:text=,hybrid%20particle%2Fgrid%20nature%20of%20MPM)[

![](https://www.google.com/s2/favicons?domain=https://dl.acm.org&sz=32)

A massively parallel and scalable multi-GPU material point method

https://dl.acm.org/doi/10.1145/3386569.3392442

](https://dl.acm.org/doi/10.1145/3386569.3392442#:~:text=A%20massively%20parallel%20and%20scalable,structure%20that%20promotes%20coalesced%20memory)[

G. Klar – Dissertation

http://web.cs.ucla.edu/~dt/theses/klar-thesis.pdf

](http://web.cs.ucla.edu/~dt/theses/klar-thesis.pdf#:~:text=of%20the%20granular%20media%20is,Our%20approach%20is%20able%20to)

All Sources

[

congress.cimne

](https://congress.cimne.com/iacm-eccomas2014/admin/files/filePaper/p701.pdf#:~:text=assemble%20the%20correct%20vertex%20data,be%20analyzed%20for%20both%20approaches)[

wang-mengdi.github

](https://wang-mengdi.github.io/proj/25-cirrus/cirrus-preprint.pdf#:~:text=%5BPDF%5D%20Cirrus%3A%20Adaptive%20Hybrid%20Particle,bottom%20line%29%20only)[

pages.cs.wisc

](https://pages.cs.wisc.edu/~sifakis/papers/GPU_MPM.pdf#:~:text=the%20same%20location%20are%20conventionally,cell%20idpr%20ev%20then)[

![](https://www.google.com/s2/favicons?domain=https://en.wikipedia.org&sz=32)

en.wikipedia

](https://en.wikipedia.org/wiki/Material_point_method#:~:text=1,the%20same%20used%20in%20FEM)[

![](https://www.google.com/s2/favicons?domain=https://stackoverflow.com&sz=32)

stackoverflow

](https://stackoverflow.com/questions/5736968/why-is-cuda-pinned-memory-so-fast#:~:text=CUDA%20Driver%20checks%2C%20if%20the,page%20copy)[

web.cs.ucla

](http://web.cs.ucla.edu/~dt/theses/klar-thesis.pdf#:~:text=We%20model%20the%20dynamics%20with,Our%20approach%20is%20able%20to)[

![](https://www.google.com/s2/favicons?domain=https://www.sciencedirect.com&sz=32)

sciencedirect

](https://www.sciencedirect.com/science/article/pii/S0266352X2400836X#:~:text=The%20goal%20of%20the%20return,Let%20us)[

![](https://www.google.com/s2/favicons?domain=https://yzhu.io&sz=32)

yzhu

](https://yzhu.io/publication/mpmmls2018siggraph/paper.pdf#:~:text=,to%20the%20adoption%20of)[

![](https://www.google.com/s2/favicons?domain=https://www.cof.orst.edu&sz=32)

cof.orst

](https://www.cof.orst.edu/cof/wse/faculty/Nairn/papers/MPMContactRevisited.pdf#:~:text=Given%20normal%20vector%2C%20early%20MPM,condition%20and%20the%20definition%20of)[

![](https://www.google.com/s2/favicons?domain=https://arxiv.org&sz=32)

arxiv

](https://arxiv.org/html/2308.01060v2#:~:text=Jiang%20et%20al%20,numerical%20dissipation%20and%20improving%20preservation)[

math.ucdavis

](https://math.ucdavis.edu/~jteran/papers/SSCTS13.pdf#:~:text=Mathematics%20math,Drucker%20and)[

![](https://www.google.com/s2/favicons?domain=https://dev.epicgames.com&sz=32)

dev.epicgames

](https://dev.epicgames.com/documentation/en-us/unreal-engine/runtime-virtual-texturing-in-unreal-engine#:~:text=A%20Runtime%20Virtual%20Texture%20,to%20conform%20to%20the%20terrain)[

![](https://www.google.com/s2/favicons?domain=https://forums.unrealengine.com&sz=32)

forums.unrealengine

](https://forums.unrealengine.com/t/rvt-virtual-heightfield-mesh-multiple-questions-issues/548820#:~:text=RVT%20%2B%20Virtual%20Heightfield%20mesh,according%20to%20youtube%20tutorials%2C)[

![](https://www.google.com/s2/favicons?domain=https://gpuopen.com&sz=32)

gpuopen

](https://gpuopen.com/download/GDC2017-Asynchronous-Compute-Deep-Dive.pdf#:~:text=Async%20Compute%20%E2%86%92%20More%20Performance,CPU)[

![](https://www.google.com/s2/favicons?domain=https://developer.nvidia.com&sz=32)

developer.nvidia

](https://developer.nvidia.com/blog/advanced-api-performance-async-compute-and-overlap/#:~:text=This%20post%20covers%20best%20practices,applications%2C%20see%20all%20Advanced)[

![](https://www.google.com/s2/favicons?domain=https://bruce-lee-ly.medium.com&sz=32)

bruce-lee-ly.medium

](https://bruce-lee-ly.medium.com/nvidia-gpu-memory-pool-bfc-d3502b355a82#:~:text=)[

![](https://www.google.com/s2/favicons?domain=https://johnnysswlab.com&sz=32)

johnnysswlab

](https://johnnysswlab.com/when-vectorization-hits-the-memory-wall-investigating-the-avx2-memory-gather-instruction/#:~:text=Vectorization%20puts%20a%20large%20additional,vectorized%20way)[

![](https://www.google.com/s2/favicons?domain=https://forums.anandtech.com&sz=32)

forums.anandtech

](https://forums.anandtech.com/threads/amd-vs-nvidia-asynchronous-compute-performance.2504980/#:~:text=AnandTech%20forums,at%20least%20on)[

![](https://www.google.com/s2/favicons?domain=https://www.reddit.com&sz=32)

reddit

](https://www.reddit.com/r/unrealengine/comments/qewicz/runtime_virtual_heightfield_mesh_help_desperately/#:~:text=Runtime%20virtual%20heightfield%20mesh%20help,so%20much%20in%20practical%20use)[

cg.cs.tsinghua.edu

](https://cg.cs.tsinghua.edu.cn/papers/CGF-2018-mpm.pdf#:~:text=,hybrid%20particle%2Fgrid%20nature%20of%20MPM)[

![](https://www.google.com/s2/favicons?domain=https://dl.acm.org&sz=32)

dl.acm

](https://dl.acm.org/doi/10.1145/3386569.3392442#:~:text=A%20massively%20parallel%20and%20scalable,structure%20that%20promotes%20coalesced%20memory)