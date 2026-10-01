# Migration gaps

Work list for the next phase: every legacy mathematical function (commit `84f3ee89`) whose row in `docs/MIGRATION_LEDGER.md` is still `pending` or `partial`, grouped by domain, largest gap first. Each line says what the legacy implementation did (read from its source) and, for `partial` rows, what the new home lacks. Sections owned by the other migration agent (cryptography, error correction, finite fields, graphs, topology, fractals, computer graphics, groups) are listed from the ledger as they stood when this file was written.

Legacy dropped-by-design areas (Expr/DAG/handle/FFI/plugin/JIT/config plumbing) are not listed.

| domain | pending | partial |
|---|---|---|
| Physics formula libraries (symbolic) | 83 | 9 |
| Graph theory (other agent) | 52 | 0 |
| Complex analysis and multi-valued functions | 30 | 16 |
| Coding theory and finite fields (other agent) | 46 | 0 |
| Differential geometry, tensors and coordinates | 33 | 9 |
| Rewriting, units, radicals, optimisation, functional analysis, geometric algebra | 16 | 24 |
| Partial differential equations (symbolic) | 37 | 0 |
| Integral transforms and convolution | 26 | 2 |
| Topology (other agent) | 24 | 0 |
| Group theory and Lie algebras | 24 | 0 |
| Computer graphics (other agent) | 22 | 0 |
| Special functions | 14 | 4 |
| Calculus, integration and series | 7 | 10 |
| Polynomials, factorisation and algebraic geometry | 7 | 10 |
| Cryptography (other agent) | 15 | 0 |
| Integral equations and calculus of variations | 14 | 0 |
| Output formats and plotting | 13 | 0 |
| Fractals and chaos (other agent) | 12 | 0 |
| Indefinite sums and products | 11 | 0 |
| Symbolic ODEs | 6 | 2 |
| Simulation scenarios and numeric solvers (test gaps) | 0 | 7 |
| Verification (proof) | 7 | 0 |
| Linear algebra, statistics and combinatorics | 0 | 5 |

## Physics formula libraries (symbolic) (92)

### `src/symbolic/classical_mechanics.rs`

- `new`: Creates a new `Kinematics` state from a given position expression. **pending**
- `newtons_second_law`: Calculates the force `F` using Newton's second law, `F = m * a`. **pending**
- `momentum`: Calculates the momentum `p` of an object, `p = m * v`. **partial**: `sim::classical::momentum` (f64 only; tests/sim/classical.rs test_particle3d_momentum)
- `kinetic_energy`: Calculates the kinetic energy `T` of an object, `T = 1/2 * m * v^2`. **partial**: `sim::classical::kinetic_energy` (f64 only; test_particle3d_kinetic_energy)
- `potential_energy_gravity_uniform`: Calculates the gravitational potential energy near Earth's surface, `V = m * g * h`. **pending**
- `potential_energy_gravity_universal`: Calculates the universal gravitational potential energy, `V = -G * m1 * m2 / r`. **pending**
- `potential_energy_spring`: Calculates the potential energy of a spring, `V = 1/2 * k * x^2`. **pending**
- `work_constant_force`: Calculates the mechanical work done by a constant force, `W = F · d`. **pending**
- `work_line_integral`: Calculates the work done by a variable force field along a path. **pending**
- `power`: Calculates the power delivered by a force, `P = F · v`. **pending**
- `torque`: Calculates the torque, `τ = r × F`. **pending**
- `angular_momentum`: Calculates the angular momentum, `L = r × p`. **pending**
- `centripetal_acceleration`: Calculates the centripetal acceleration, `a_c = v^2 / r`. **pending**
- `moment_of_inertia_point_mass`: Calculates the moment of inertia for a point mass, `I = m * r^2`. **pending**
- `rotational_kinetic_energy`: Calculates the rotational kinetic energy, `T_rot = 1/2 * I * ω^2`. **pending**
- `lagrangian`: Calculates the Lagrangian `L = T - V`. **pending**
- `hamiltonian`: Calculates the Hamiltonian `H = T + V`. **pending**
- `euler_lagrange_equation`: Computes the left-hand side of the Euler-Lagrange equation. **pending**
- `poisson_bracket`: Calculates the Poisson bracket `{f, g} = (∂f/∂q)(∂g/∂p) - (∂f/∂p)(∂g/∂q)`. **pending**

### `src/symbolic/electromagnetism.rs`

- `new`: Constructs Maxwell's equations from the given fields and sources. **pending**
- `lorentz_force`: Calculates the Lorentz force on a particle with charge $q$. **partial**: `sim::classical::lorentz_force` (f64 only; test_lorentz_force)
- `electric_field_from_potentials`: Calculates the electric field $\mathbf{E}$ from scalar potential $V$ and vector potential $\mathbf{A}$. **pending**
- `electric_field_from_potential`: Calculates the electric field $\mathbf{E}$ from the scalar electric potential $V$. **pending**
- `magnetic_field_from_vector_potential`: Calculates the magnetic field $\mathbf{B}$ from the vector potential $\mathbf{A}$. **pending**
- `poynting_vector`: Calculates the Poynting vector, representing the directional energy flux density. **pending**
- `energy_density`: Calculates the electromagnetic energy density $u$. **pending**
- `coulombs_law`: Represents Coulomb's Law for the electric field produced by a point charge. **partial**: `sim::classical::electric_field_point_charge` / `coulomb_force` (f64 only)

### `src/symbolic/quantum_field_theory.rs`

- `dirac_adjoint`: Computes the Dirac adjoint of a fermion field: **pending**
- `feynman_slash`: Computes the Feynman slash notation: **pending**
- `scalar_field_lagrangian`: Lagrangian density for a free real scalar field (Klein-Gordon): **pending**
- `qed_lagrangian`: Lagrangian density for Quantum Electrodynamics (QED): **pending**
- `qcd_lagrangian`: Lagrangian density for Quantum Chromodynamics (QCD): **pending**
- `propagator`: Represents a propagator for a particle in QFT. **pending**
- `scattering_cross_section`: Scattering cross-section: **pending**
- `feynman_propagator_position_space`: Feynman propagator in position space (symbolic integral representation). **pending**

### `src/symbolic/quantum_mechanics.rs`

- `bra_ket`: Computes the inner product of a Bra and a Ket, `<Bra|Ket>`. **pending**
- `bra_ket_internal`: internal solver of `bra_ket` **pending**
- `new`: Creates a new operator from an expression. **pending**
- `apply`: Applies an operator to a Ket, `O|Ket>`. **pending**
- `commutator`: Computes the commutator of two operators: **pending**
- `commutator_internal`: internal solver of `commutator` **pending**
- `expectation_value`: Computes the expectation value of an operator: **pending**
- `expectation_value_internal`: internal solver of `expectation_value` **pending**
- `uncertainty`: Computes the uncertainty (standard deviation) of an operator: **pending**
- `probability_density`: Computes the probability density at a point: **pending**
- `hamiltonian_free_particle`: Hamiltonian for a free particle: **pending**
- `hamiltonian_harmonic_oscillator`: Hamiltonian for a harmonic oscillator: **pending**
- `angular_momentum_z`: Angular momentum operator `L_z`: **pending**
- `pauli_matrices`: Returns the Pauli matrices: **pending**
- `spin_operator`: Spin operator: **pending**
- `solve_time_independent_schrodinger`: Solves the time-independent Schrödinger equation: **pending**
- `time_dependent_schrodinger_equation`: Time-dependent Schrödinger equation: **pending**
- `dirac_equation`: Dirac equation for a free particle: **pending**
- `klein_gordon_equation`: Klein-Gordon equation: **pending**
- `first_order_energy_correction`: Computes the first-order energy correction in perturbation theory: **pending**
- `scattering_amplitude`: Scattering amplitude in quantum mechanics: **pending**

### `src/symbolic/relativity.rs`

- `lorentz_factor`: Calculates the Lorentz factor: **partial**: `sim::classical::lorentz_factor` (f64 only; test_lorentz_factor_high_speed)
- `lorentz_transformation_x`: Performs a Lorentz transformation in the x-direction. **pending**
- `velocity_addition`: Calculates relativistic velocity addition. **partial**: `sim::classical::relativistic_velocity_addition` (f64 only)
- `mass_energy_equivalence`: Calculates mass-energy equivalence: **partial**: `sim::classical::mass_energy` (f64 only)
- `relativistic_momentum`: Calculates relativistic momentum: **partial**: `sim::classical::relativistic_momentum` (f64 only)
- `doppler_effect`: Calculates the Relativistic Doppler Effect for source and observer moving apart. **pending**
- `schwarzschild_radius`: Calculates the Schwarzschild Radius: **pending**
- `gravitational_time_dilation`: Calculates gravitational time dilation in the Schwarzschild metric. **pending**
- `einstein_tensor`: Represents the simplified Einstein Field Equations (LHS). **pending**
- `geodesic_acceleration`: Represents the Geodesic Equation term: **pending**

### `src/symbolic/solid_state_physics.rs`

- `new`: Creates a new crystal lattice. **pending**
- `volume`: Computes the volume of the unit cell: **pending**
- `reciprocal_lattice_vectors`: Computes the reciprocal lattice vectors: **pending**
- `bloch_theorem`: Represents Bloch's Theorem: **pending**
- `energy_band`: Represents a simple energy band model using the parabolic band approximation. **pending**
- `density_of_states_3d`: Computes the Density of States (DOS) for a 3D electron gas. **pending**
- `fermi_energy_3d`: Fermi Energy for a 3D electron gas: **pending**
- `drude_conductivity`: Drude model electrical conductivity: **pending**
- `hall_coefficient`: Hall Coefficient: **pending**
- `debye_frequency`: Debye Frequency: **pending**
- `einstein_heat_capacity`: Einstein Heat Capacity: **pending**
- `plasma_frequency`: Plasma Frequency: **pending**
- `london_penetration_depth`: London penetration depth: **pending**

### `src/symbolic/thermodynamics.rs`

- `first_law_thermodynamics`: Represents the First Law of Thermodynamics: **pending**
- `ideal_gas_law`: Represents the Ideal Gas Law: **partial**: `sim::classical::{ideal_gas_pressure,ideal_gas_volume,ideal_gas_temperature}` (f64 only)
- `enthalpy`: Calculates Enthalpy: **pending**
- `helmholtz_free_energy`: Calculates Helmholtz Free Energy: **pending**
- `gibbs_free_energy`: Calculates Gibbs Free Energy: **pending**
- `boltzmann_entropy`: Calculates Entropy via Boltzmann's formula: **pending**
- `carnot_efficiency`: Calculates the efficiency of a Carnot engine: **pending**
- `boltzmann_distribution`: Represents the Boltzmann Distribution: **pending**
- `partition_function`: Calculates the Partition Function: **pending**
- `fermi_dirac_distribution`: Represents the Fermi-Dirac Distribution for fermions. **pending**
- `bose_einstein_distribution`: Represents the Bose-Einstein Distribution for bosons. **pending**
- `work_isothermal_expansion`: Calculates the work done during an isothermal expansion: **pending**
- `verify_maxwell_relation_helmholtz`: Maxwell Relation examples: **pending**

## Graph theory (other agent) (52)

### `src/symbolic/graph.rs`

- `new`: Creates a new graph. **pending**
- `nodes`: Returns a reference to the nodes in the graph. **pending**
- `node_count`: Returns the number of nodes in the graph. **pending**
- `is_directed`: Returns true if the graph is directed. **pending**
- `add_node`: Adds a node with a given label to the graph. **pending**
- `add_edge`: Adds an edge between two nodes. **pending**
- `get_node_id`: Gets the internal ID of a node given its label. **pending**
- `neighbors`: Gets the neighbors of a node. **pending**
- `out_degree`: Gets the out-degree of a node. **pending**
- `in_degree`: Gets the in-degree of a node. **pending**
- `get_edges`: Returns a list of all edges in the graph. **pending**
- `add_hyperedge`: Adds a hyperedge that connects a set of vertices. **pending**
- `to_adjacency_matrix`: Returns the adjacency matrix of the graph. **pending**
- `to_incidence_matrix`: Returns the incidence matrix of the graph. **pending**
- `to_laplacian_matrix`: Returns the Laplacian matrix of the graph (`L = D - A`). **pending**

### `src/symbolic/graph_algorithms.rs`

- `dfs`: Performs a Depth-First Search (DFS) traversal on a graph. **pending**
- `bfs`: Performs a Breadth-First Search (BFS) traversal on a graph. **pending**
- `connected_components`: Finds all connected components of an undirected graph. **pending**
- `is_connected`: Checks if the graph is connected. **pending**
- `strongly_connected_components`: Finds all strongly connected components (SCCs) of a directed graph using Tarjan's algorithm. **pending**
- `has_cycle`: Detects if a cycle exists in the graph. **pending**
- `find_bridges_and_articulation_points`: Finds all bridges and articulation points (cut vertices) in a graph using Tarjan's algorithm. **pending**
- `kruskal_mst`: Finds the Minimum Spanning Tree (MST) of a graph using Kruskal's algorithm. **pending**
- `edmonds_karp_max_flow`: Finds the maximum flow from a source `s` to a sink `t` in a flow network using the Edmonds-Karp algorithm. **pending**
- `dinic_max_flow`: Finds the maximum flow from a source `s` to a sink `t` in a flow network using Dinic's algorithm. **pending**
- `bellman_ford`: Finds the shortest paths from a single source in a graph with possible negative edge weights. **pending**
- `min_cost_max_flow`: Solves the Minimum-Cost Maximum-Flow problem using the successive shortest path algorithm with Bellman-Ford. **pending**
- `is_bipartite`: Checks if a graph is bipartite using BFS-based 2-coloring. **pending**
- `bipartite_maximum_matching`: Finds the maximum cardinality matching in a bipartite graph by reducing it to a max-flow problem. **pending**
- `prim_mst`: Finds the Minimum Spanning Tree (MST) of a graph using Prim's algorithm. **pending**
- `topological_sort_kahn`: Performs a topological sort on a directed acyclic graph (DAG) using Kahn's algorithm (BFS-based). **pending**
- `topological_sort_dfs`: Performs a topological sort on a directed acyclic graph (DAG) using a DFS-based algorithm. **pending**
- `topological_sort`: Performs a topological sort on a directed acyclic graph (DAG). **pending**
- `bipartite_minimum_vertex_cover`: Finds the minimum vertex cover of a bipartite graph using Kőnig's theorem. **pending**
- `hopcroft_karp_bipartite_matching`: Finds the maximum cardinality matching in a bipartite graph using the Hopcroft-Karp algorithm. **pending**
- `blossom_algorithm`: Finds the maximum cardinality matching in a general graph using Edmonds's Blossom Algorithm. **pending**
- `shortest_path_unweighted`: Finds the shortest path in an unweighted graph from a source node using BFS. **pending**
- `dijkstra`: Finds the shortest paths from a single source using Dijkstra's algorithm. **pending**
- `floyd_warshall`: Finds all-pairs shortest paths using the Floyd-Warshall algorithm. **pending**
- `spectral_analysis`: Performs spectral analysis on a graph matrix (e.g., Adjacency or Laplacian). **pending**
- `algebraic_connectivity`: Computes the algebraic connectivity of a graph. **pending**

### `src/symbolic/graph_isomorphism_and_coloring.rs`

- `are_isomorphic_heuristic`: Checks if two graphs are potentially isomorphic using the Weisfeiler-Lehman test (Color Refinement). **pending**
- `greedy_coloring`: Finds a valid vertex coloring using a greedy heuristic (Welsh-Powell algorithm). **pending**
- `chromatic_number_exact`: Finds the chromatic number of a graph using exhaustive backtracking. **pending**

### `src/symbolic/graph_operations.rs`

- `induced_subgraph`: Creates an induced subgraph from a given set of node labels. **pending**
- `union`: Computes the union of two graphs. **pending**
- `intersection`: Computes the intersection of two graphs. **pending**
- `cartesian_product`: Computes the Cartesian product of two graphs. **pending**
- `tensor_product`: Computes the Tensor product of two graphs. **pending**
- `complement`: Computes the complement of a graph. **pending**
- `disjoint_union`: Computes the disjoint union of two graphs. **pending**
- `join`: Computes the join of two graphs (Zykov sum). **pending**

## Complex analysis and multi-valued functions (46)

There is no complex evaluation of terms; `kernels::complex` offers numeric closure-based contour integrals, residues, zero/pole counting and Mobius maps.

### `src/numerical/complex_analysis.rs`

- `contour_integral_expr`: Performs numerical contour integration of a symbolic expression. **partial**: `kernels::complex::{contour_integral,residue}` take closures; no term-level operator (needs complex evaluation of terms)
- `residue_expr`: Calculates the residue of a symbolic expression at a point. **partial**: `kernels::complex::{contour_integral,residue}` take closures; no term-level operator (needs complex evaluation of terms)
- `eval_complex_expr`: Evaluates a symbolic expression to a numerical `Complex<f64>` value. **pending**

### `src/numerical/multi_valued.rs`

- `newton_method_complex`: Finds a root of a complex function `f(z) = 0` using Newton's method. **pending**
- `complex_log_k`: Computes the k-th branch of the complex logarithm. **pending**
- `complex_sqrt_k`: Computes the k-th branch of the complex square root. **pending**
- `complex_pow_k`: Computes the k-th branch of the complex power z^w. **pending**
- `complex_nth_root_k`: Computes the k-th branch of the complex n-th root. **pending**
- `complex_arcsin_k`: Computes the k-th branch of the complex arcsine. **pending**
- `complex_arccos_k`: Computes the k-th branch of the complex arccosine. **pending**
- `complex_arctan_k`: Computes the k-th branch of the complex arctangent. **pending**

### `src/symbolic/complex_analysis.rs`

- `new`: Creates a new analytic continuation starting with a Taylor series for `func` centered at `start_point`. **partial**: `kernels::complex::Mobius::new` (numeric coefficients; the analytic-continuation `new` has no home)
- `continue_along_path`: Continues the function along a given path. **pending**
- `get_final_expression`: Returns the final expression (Taylor series) after continuation to the last point. **pending**
- `estimate_radius_of_convergence`: Estimates the radius of convergence for a Taylor series using the ratio test. **pending**
- `complex_distance`: Calculates the Euclidean distance between two complex points. **pending**
- `classify_singularity`: Classifies the type of singularity at a point. **pending**
- `laurent_series`: Computes the Laurent series expansion around a point. **partial**: `rules::calculus` `laurent` (real expansion variable; test taylor_and_laurent)
- `calculate_residue`: Calculates the residue of a function at a singularity. **partial**: `kernels::complex::residue` numeric only
- `calculate_residue_internal`: internal solver of `calculate_residue` **partial**: `kernels::complex::residue` numeric only
- `contour_integral_residue_theorem`: Evaluates a contour integral using the residue theorem. **partial**: `kernels::complex::contour_integral` numeric (no residue sum)
- `contour_integral_residue_theorem_internal`: internal solver of `contour_integral_residue_theorem` **partial**: `kernels::complex::contour_integral` numeric (no residue sum)
- `identity`: Creates the identity transformation. **partial**: `kernels::complex::Mobius::identity` (numeric)
- `apply`: Applies the transformation to a point z. **partial**: `kernels::complex::Mobius::apply` (numeric)
- `compose`: Composes two Möbius transformations. **partial**: `kernels::complex::Mobius::compose` (numeric)
- `inverse`: Computes the inverse transformation. **partial**: `kernels::complex::Mobius::inverse` (numeric)
- `cauchy_integral_formula`: Evaluates f(z0) using Cauchy's integral formula. **partial**: `kernels::complex::circle_integral` numeric
- `cauchy_integral_formula_internal`: internal solver of `cauchy_integral_formula` **partial**: `kernels::complex::circle_integral` numeric
- `cauchy_derivative_formula`: Computes the n-th derivative using Cauchy's formula for derivatives. **partial**: `kernels::complex::cauchy_derivative` numeric (unit test derivatives)
- `cauchy_derivative_formula_internal`: internal solver of `cauchy_derivative_formula` **partial**: `kernels::complex::cauchy_derivative` numeric (unit test derivatives)
- `complex_exp`: Computes the complex exponential e^z = e^(x+iy) = e^x(cos(y) + i*sin(y)) **pending**
- `complex_log`: Computes the principal branch of complex logarithm. **pending**
- `complex_arg`: Computes the argument (angle) of a complex number. **pending**
- `complex_modulus`: Computes the modulus (absolute value) of a complex number. **pending**

### `src/symbolic/multi_valued.rs`

- `arg`: Returns the principal argument of a complex expression `z`. **pending**
- `abs`: Returns the absolute value (magnitude) of a complex expression `z`. **pending**
- `general_log`: Computes the general multi-valued logarithm of a complex expression `z`. **pending**
- `general_sqrt`: Computes the general multi-valued square root of a complex expression `z`. **pending**
- `general_power`: Computes the general multi-valued power `z^w`. **pending**
- `general_nth_root`: Computes the general multi-valued n-th root of a complex expression `z`. **pending**
- `general_arcsin`: Computes the general multi-valued arcsin of a complex expression `z`. **pending**
- `general_arccos`: Computes the general multi-valued arccos of a complex expression `z`. **pending**
- `general_arctan`: Computes the general multi-valued arctan of a complex expression `z`. **pending**
- `general_arcsinh`: Computes the general multi-valued inverse hyperbolic sine (arcsinh). **pending**
- `general_arccosh`: Computes the general multi-valued inverse hyperbolic cosine (arccosh). **pending**
- `general_arctanh`: Computes the general multi-valued inverse hyperbolic tangent (arctanh). **pending**

## Coding theory and finite fields (other agent) (46)

### `src/symbolic/error_correction.rs`

- `hamming_distance`: Computes the Hamming distance between two byte slices. **pending**
- `hamming_weight`: Computes the Hamming weight (number of 1s) of a byte slice. **pending**
- `hamming_encode`: Encodes a 4-bit data block into a 7-bit Hamming(7,4) codeword. **pending**
- `hamming_check`: Checks if a Hamming(7,4) codeword is valid without correcting. **pending**
- `hamming_decode`: Decodes a 7-bit Hamming(7,4) codeword, correcting a single-bit error if found. **pending**
- `rs_encode`: Encodes a data message using a Reed-Solomon code, adding `n_sym` error correction symbols. **pending**
- `rs_check`: Checks if a Reed-Solomon codeword is valid without attempting correction. **pending**
- `rs_error_count`: Estimates the number of errors in a Reed-Solomon codeword. **pending**
- `rs_decode`: Decodes a Reed-Solomon codeword, correcting errors if found. **pending**
- `crc32_compute`: Computes CRC-32 checksum of data. **pending**
- `crc32_verify`: Verifies CRC-32 checksum of data. **pending**
- `crc32_update`: Updates an existing CRC-32 with additional data. **pending**
- `crc32_finalize`: Finalizes a CRC-32 computation started with `crc32_update`. **pending**

### `src/symbolic/error_correction_helper.rs`

- `new`: Creates a new finite field `GF(modulus)`. **pending**
- `from_bigint`: Creates a finite field from a `BigInt` modulus. **pending**
- `is_zero`: Returns true if this element is zero. **pending**
- `is_one`: Returns true if this element is one. **pending**
- `inverse`: Computes the multiplicative inverse of the element in the finite field. **pending**
- `pow`: Computes a^exp mod p using binary exponentiation. **pending**
- `gf256_exp`: Computes the exponentiation (anti-logarithm) in GF(2^8). **pending**
- `gf256_log`: Computes the discrete logarithm in GF(2^8). **pending**
- `gf256_add`: Performs addition in the finite field GF(2^8). **pending**
- `gf256_mul`: Performs multiplication in GF(2^8) using precomputed lookup tables. **pending**
- `gf256_inv`: Computes the multiplicative inverse of an element in GF(2^8). **pending**
- `gf256_div`: Performs division in GF(2^8). **pending**
- `gf256_pow`: Computes a^exp in GF(2^8). **pending**
- `poly_eval_gf256`: Evaluates a polynomial over GF(2^8) at a given point `x`. **pending**
- `poly_add_gf256`: Adds two polynomials over GF(2^8). **pending**
- `poly_mul_gf256`: Multiplies two polynomials over GF(2^8). **pending**
- `poly_scale_gf256`: Scales a polynomial by a constant in GF(2^8). **pending**
- `poly_derivative_gf256`: Computes the formal derivative of a polynomial in GF(2^8). **pending**
- `poly_gcd_gf256`: Computes the GCD of two polynomials over GF(2^8) using Euclidean algorithm. **pending**
- `poly_div_gf256`: Divides two polynomials over GF(2^8). **pending**
- `poly_add_gf`: Adds two polynomials whose coefficients are `FieldElement`s from a given finite field. **pending**
- `poly_mul_gf`: Multiplies two polynomials whose coefficients are `FieldElement`s from a given finite field. **pending**
- `poly_div_gf`: Divides two polynomials whose coefficients are `FieldElement`s from a given finite field. **pending**

### `src/symbolic/finite_field.rs`

- `serialize`:  **pending**
- `deserialize`:  **pending**
- `new`: Creates a new prime field `GF(p)` with the given modulus. **pending**
- `inverse`: Computes the multiplicative inverse of the element in the prime field. **pending**
- `degree`: Returns the degree of the polynomial. **pending**
- `long_division`: Performs polynomial long division over the prime field. **pending**
- `add`: Adds two extension field elements. **pending**
- `sub`: Subtracts one extension field element from another. **pending**
- `mul`: Multiplies two extension field elements. **pending**
- `div`: Divides one extension field element by another. **pending**

## Differential geometry, tensors and coordinates (42)

### `src/numerical/coordinates.rs`

- `transform_point`: Transforms a numerical point from one coordinate system to another. **pending**
- `numerical_jacobian`: Computes the numerical Jacobian matrix of a coordinate transformation at a specific point. **pending**
- `transform_point_pure`: Transforms a numerical point using direct `f64` calculations for high performance. **pending**

### `src/numerical/differential_geometry.rs`

- `metric_tensor_at_point`: Evaluates the metric tensor at a given point for a coordinate system. **pending**
- `christoffel_symbols`: Computes the Christoffel symbols of the second kind at a given point. **pending**
- `riemann_tensor`: Computes the Riemann curvature tensor at a given point. **pending**
- `ricci_tensor`: Computes the Ricci tensor at a given point. **pending**
- `ricci_scalar`: Computes the Ricci scalar at a given point. **pending**

### `src/symbolic/coordinates.rs`

- `transform_point`: Transforms a point from one coordinate system to another. **pending**
- `transform_expression`: Transforms a symbolic expression from one coordinate system to another. **pending**
- `get_transform_rules`: Helper function to get the variables and transformation rules between two coordinate systems. **pending**
- `get_to_cartesian_rules`: Provides the transformation rules from a given coordinate system to Cartesian coordinates. **pending**
- `transform_contravariant_vector`: Transforms a contravariant vector field (e.g., velocity) from one coordinate system to another. **pending**
- `transform_covariant_vector`: Transforms a covariant vector field (e.g., gradient) from one coordinate system to another. **pending**
- `transform_tensor2`: Transforms a rank-2 tensor field from one coordinate system to another. **pending**
- `symbolic_mat_mat_mul`: Performs symbolic matrix-matrix multiplication. **pending**
- `get_metric_tensor`: Computes and returns the metric tensor for a given orthogonal coordinate system. **pending**
- `transform_divergence`: Computes the divergence of a contravariant vector field in any orthogonal coordinate system. **pending**
- `transform_curl`: Computes the curl of a covariant vector field in any orthogonal coordinate system. **pending**
- `transform_gradient`: Transforms the gradient of a scalar field from one coordinate system to another. **pending**

### `src/symbolic/differential_geometry.rs`

- `exterior_derivative`: Computes the exterior derivative of a k-form, resulting in a (k+1)-form. **pending**
- `wedge_product`: Computes the wedge product (exterior product) of two differential forms. **pending**
- `boundary`: Represents the boundary of a domain, denoted as `∂M` for a manifold `M`. **pending**
- `generalized_stokes_theorem`: Represents the generalized Stokes' Theorem. **pending**
- `gauss_theorem`: Represents Gauss's Theorem (Divergence Theorem) as a special case of Stokes' Theorem. **pending**
- `stokes_theorem`: Represents the classical Stokes' Theorem as a special case of the generalized theorem. **pending**
- `greens_theorem`: Represents Green's Theorem as a 2D special case of Stokes' Theorem. **pending**

### `src/symbolic/tensor.rs`

- `new`: Creates a new `Tensor` with the given components and shape. **partial**: `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components)
- `rank`: Returns the rank (order) of the tensor. **partial**: `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components)
- `get`: Returns an immutable reference to the component at the specified indices. **partial**: `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components)
- `add`: Performs tensor addition with another tensor. **partial**: `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components)
- `sub`: Performs tensor subtraction with another tensor. **partial**: `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components)
- `scalar_mul`: Multiplies the tensor by a scalar expression. **partial**: `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components)
- `outer_product`: Computes the outer product of this tensor with another tensor. **partial**: `kernels::tensor::outer_product` (f64 ndarray only)
- `contract`: Contracts two specified axes of the tensor. **partial**: `kernels::tensor::contract` (f64 ndarray only)
- `to_matrix_expr`: Converts a rank-2 tensor into an `Expr::Matrix`. **partial**: `kernels::tensor::TensorData` / `ndarray` arrays (f64 only; no symbolic components)
- `raise_index`: Raises an index of a covector (rank-1 tensor with lower index) to a vector (upper index). **pending**
- `lower_index`: Lowers an index of a vector (rank-1 tensor with upper index) to a covector (lower index). **pending**
- `christoffel_symbols_first_kind`: Computes the Christoffel symbols of the first kind `Γ_{ijk}`. **pending**
- `christoffel_symbols_second_kind`: Computes the Christoffel symbols of the second kind `Γ^i_{jk}`. **pending**
- `riemann_curvature_tensor`: Computes the Riemann curvature tensor `R^i_{jkl}`. **pending**
- `covariant_derivative_vector`: Computes the covariant derivative of a vector field `V^i` with respect to a coordinate `x^k`. **pending**

## Rewriting, units, radicals, optimisation, functional analysis, geometric algebra (40)

### `src/numerical/elementary.rs`

- `floor`: Floor rounding. **pending**
- `ceil`: Ceil rounding. **pending**
- `round`: Round to nearest integer. **pending**

### `src/numerical/optimize.rs`

- `auto_solve_conjugate_gradient`: Automatically configures and solves a problem using the Conjugate Gradient method. **partial**: `kernels::optimize::EquationOptimizer::auto_solve_conjugate_gradient` (no test: needs a custom argmin `Operator` impl, more than a few lines)
- `auto_solve`: Automatically select solver and solve Returns an error if the optimization process fails. **partial**: `kernels::optimize::EquationOptimizer::auto_solve` (no test: `problem: P` is the `Array1<f64>` alias, which is never a `CostFunction`, so it cannot be called)

### `src/symbolic/elementary.rs`

- `expand`: Expands a symbolic expression by applying distributive, power, and trigonometric identities. **partial**: `rules::poly` `expand` (polynomial expansion); sum-angle expansion only in the Explore tier of `rules::elementary`
- `expand_internal`: internal solver of `expand` **partial**: `rules::poly` `expand` (polynomial expansion); sum-angle expansion only in the Explore tier of `rules::elementary`

### `src/symbolic/functional_analysis.rs`

- `new`: Creates a new L^2 space on the interval `[a, b]`. **pending**
- `apply`: Applies the operator to a given expression (function). **pending**
- `inner_product`: Computes the inner product of two functions, `f` and `g`, in a given Hilbert space. **partial**: `kernels::functional_analysis::inner_product` (discrete samples only, no symbolic integrals)
- `inner_product_internal`: internal solver of `inner_product` **partial**: `kernels::functional_analysis::inner_product` (discrete samples only, no symbolic integrals)
- `norm`: Computes the norm of a function `f` in a given Hilbert space. **partial**: `kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}` (discrete samples; no general L^p)
- `norm_internal`: internal solver of `norm` **partial**: `kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}` (discrete samples; no general L^p)
- `banach_norm`: Computes the L^p norm of a function `f` in a given Banach space. **partial**: `kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}` (discrete samples; no general L^p)
- `banach_norm_internal`: internal solver of `banach_norm` **partial**: `kernels::functional_analysis::{l1_norm,l2_norm,infinity_norm}` (discrete samples; no general L^p)
- `are_orthogonal`: Checks if two functions are orthogonal in a given Hilbert space. **pending**
- `project`: Computes the projection of function `f` onto function `g` in a given Hilbert space. **partial**: `kernels::functional_analysis::project` (discrete samples)
- `project_internal`: internal solver of `project` **partial**: `kernels::functional_analysis::project` (discrete samples)
- `gram_schmidt`: Performs the Gram-Schmidt process to orthogonalize a set of functions. **partial**: `kernels::functional_analysis::gram_schmidt` (discrete samples)
- `gram_schmidt_orthonormal`: Performs the Gram-Schmidt process to orthonormalize a set of functions. **partial**: `kernels::functional_analysis::gram_schmidt_orthonormal` (discrete samples)

### `src/symbolic/geometric_algebra.rs`

- `new`: Creates a new, empty multivector for a given algebra signature. **partial**: `Multivector3D` fields (numeric G3 only)
- `scalar`: Creates a new multivector representing a scalar value. **partial**: `Multivector3D` fields (numeric G3 only)
- `vector`: Creates a new multivector representing a vector (grade-1 element). **partial**: `Multivector3D` fields (numeric G3 only)
- `geometric_product`: Computes the geometric product of this multivector with another. **partial**: `kernels::geometric_algebra::Multivector3D` `Mul` (numeric G3 only)
- `grade_projection`: Extracts all terms of a specific grade from the multivector. **pending**
- `outer_product`: Computes the outer (or wedge) product of this multivector with another. **partial**: `Multivector3D::wedge` (numeric G3 only)
- `inner_product`: Computes the inner (or left contraction) product of this multivector with another. **partial**: `Multivector3D::dot` (numeric G3 only)
- `reverse`: Computes the reverse of the multivector. **partial**: `Multivector3D::reverse` (numeric G3 only)
- `magnitude`: Computes the magnitude (norm) of the multivector. **partial**: `Multivector3D::norm` (numeric G3 only)
- `dual`: Computes the dual of the multivector with respect to the pseudoscalar. **pending**
- `normalize`: Normalizes the multivector to unit magnitude. **pending**

### `src/symbolic/numeric.rs`

- `evaluate_complex`: Evaluates a symbolic expression to a complex numerical `Complex64` value, if possible. **pending**

### `src/symbolic/optimize.rs`

- `find_extrema`: Finds and classifies the critical points of a multivariate function. **pending**
- `find_constrained_extrema`: Finds the extrema of a function subject to equality constraints using the method of Lagrange Multipliers. **pending**

### `src/symbolic/radicals.rs`

- `simplify_radicals`: Recursively simplifies radical expressions in the given expression tree. **pending**
- `denest_sqrt`: Attempts to denest a nested square root of the form `sqrt(A ± B*sqrt(C))`. **pending**

### `src/symbolic/real_roots.rs`

- `count_real_roots_in_interval`: Counts the number of distinct real roots of a polynomial in an interval `(a, b]`. **partial**: `kernels::real_roots::isolate_real_roots` returns intervals (count = length); no interval-count function

### `src/symbolic/rewriting.rs`

- `apply_rules_to_normal_form`: Applies a set of rewrite rules repeatedly to an expression until a fixed point is reached. **partial**: `Session::with_rules` + `Engine` saturation (rules fixed at session build, not per call)
- `knuth_bendix`: Attempts to produce a complete term-rewriting system from a set of equations using the Knuth-Bendix completion algorithm. **pending**

### `src/symbolic/unit_unification.rs`

- `unify_expression`: Unifies an expression containing quantities with units. **pending**

## Partial differential equations (symbolic) (37)

Legacy `pde.rs` was a dispatcher (`solve_pde`) over closed-form templates: it matched the equation structure and returned the textbook general solution (D'Alembert `F(x-ct)+G(x+ct)`, separation of variables as a Fourier series with unknown coefficients, Green-function integrals, plane-wave decompositions). Only numeric finite-difference solvers exist now (`kernels::pde`, `sim::*`).

### `src/symbolic/pde.rs`

- `solve_pde`: dispatcher trying separation of variables, characteristics, Green, second order, Fourier **pending**
- `solve_pde_internal`: internal solver of `solve_pde` **pending**
- `solve_pde_by_separation_of_variables`: 1D linear homogeneous PDEs with homogeneous boundary conditions: u = X(x)T(t), two ODEs **pending**
- `solve_pde_by_separation_of_variables_internal`: internal solver of `solve_pde_by_separation_of_variables` **pending**
- `classify_pde_heuristic`: inspects linearity, homogeneity, order, dimension and type to suggest methods **pending**
- `solve_pde_by_characteristics`: first-order PDE a u_x + b u_y = c: characteristic ODE dy/dx = b/a (constant coefficients only) **pending**
- `solve_pde_by_characteristics_internal`: internal solver of `solve_pde_by_characteristics` **pending**
- `solve_pde_by_greens_function`: Green-function integral solution **pending**
- `solve_pde_by_greens_function_internal`: internal solver of `solve_pde_by_greens_function` **pending**
- `solve_second_order_pde`: classify and dispatch second-order PDEs **pending**
- `solve_second_order_pde_internal`: internal solver of `solve_second_order_pde` **pending**
- `solve_wave_equation_1d_dalembert`: u = F(x-ct) + G(x+ct) **pending**
- `solve_wave_equation_1d_dalembert_internal`: internal solver of `solve_wave_equation_1d_dalembert` **pending**
- `solve_heat_equation_1d`: Fourier sine series u = sum A_n exp(-alpha n^2 pi^2 t/L^2) sin(n pi x/L) **pending**
- `solve_heat_equation_1d_internal`: internal solver of `solve_heat_equation_1d` **pending**
- `solve_laplace_equation_2d`: Solves the 2D Laplace equation `u_xx + u_yy = 0`. **pending**
- `solve_laplace_equation_2d_internal`: internal solver of `solve_laplace_equation_2d` **pending**
- `solve_wave_equation_3d`: Solves the 3D wave equation `u_tt = c²(u_xx + u_yy + u_zz)`. **pending**
- `solve_wave_equation_3d_internal`: internal solver of `solve_wave_equation_3d` **pending**
- `solve_heat_equation_3d`: Solves the 3D heat equation `u_t = α(u_xx + u_yy + u_zz)`. **pending**
- `solve_heat_equation_3d_internal`: internal solver of `solve_heat_equation_3d` **pending**
- `solve_laplace_equation_3d`: Solves the 3D Laplace equation `u_xx + u_yy + u_zz = 0`. **pending**
- `solve_laplace_equation_3d_internal`: internal solver of `solve_laplace_equation_3d` **pending**
- `solve_poisson_equation_2d`: Solves the 2D Poisson equation `u_xx + u_yy = f(x,y)`. **pending**
- `solve_poisson_equation_2d_internal`: internal solver of `solve_poisson_equation_2d` **pending**
- `solve_poisson_equation_3d`: Solves the 3D Poisson equation `u_xx + u_yy + u_zz = f(x,y,z)`. **pending**
- `solve_poisson_equation_3d_internal`: internal solver of `solve_poisson_equation_3d` **pending**
- `solve_helmholtz_equation`: Solves the Helmholtz equation `∇²u + k²u = 0`. **pending**
- `solve_helmholtz_equation_internal`: internal solver of `solve_helmholtz_equation` **pending**
- `solve_schrodinger_equation`: Solves the time-dependent Schrödinger equation `iℏ ∂ψ/∂t = -ℏ²/(2m) ∇²ψ + V(x)ψ`. **pending**
- `solve_schrodinger_equation_internal`: internal solver of `solve_schrodinger_equation` **pending**
- `solve_klein_gordon_equation`: Solves the Klein-Gordon equation `∂²φ/∂t² - c²∇²φ + m²c⁴/ℏ² φ = 0`. **pending**
- `solve_klein_gordon_equation_internal`: internal solver of `solve_klein_gordon_equation` **pending**
- `solve_burgers_equation`: Solves the 1D Burgers' equation `u_t + u*u_x = 0`. **pending**
- `solve_burgers_equation_internal`: internal solver of `solve_burgers_equation` **pending**
- `solve_with_fourier_transform`: transforms the spatial variable to reduce the PDE to an ODE **pending**
- `solve_with_fourier_transform_internal`: internal solver of `solve_with_fourier_transform` **pending**

## Integral transforms and convolution (28)

Legacy: Fourier = `defint(f*exp(-i w t), -oo, oo)`, Laplace = `defint(f*exp(-s t), 0, oo)`, Z = bilateral sum of `f*z^-n`; inverse Laplace = table lookup, then partial fractions, then a Bromwich contour fallback; inverse Fourier = `1/(2 pi)` times the integral; inverse Z = contour integral. Property helpers (shift, scaling, differentiation, convolution) were one-line formulas. Only the numeric FFT (`kernels::transforms`) and exact partial fractions (`rules::poly::apart`) exist.

### `src/symbolic/transforms.rs`

- `fourier_time_shift`: F{f(t-a)} = e^{-i w a} F(w) **pending**
- `fourier_frequency_shift`: F{e^{i a t} f} = F(w-a) **pending**
- `fourier_scaling`: F{f(a t)} = F(w/a)/|a| **pending**
- `fourier_differentiation`: F{f'} = i w F(w) **pending**
- `laplace_time_shift`: L{f(t-a)u(t-a)} = e^{-a s} F(s) **pending**
- `laplace_differentiation`: L{f'} = s F(s) - f(0) **pending**
- `laplace_frequency_shift`: L{e^{at} f} = F(s-a) **pending**
- `laplace_scaling`: L{f(a t)} = F(s/a)/a **pending**
- `laplace_integration`: L{int f} = F(s)/s **pending**
- `z_time_shift`: Z{x[n-k]} = z^-k X(z) **pending**
- `z_scaling`: Z{a^n x[n]} = X(z/a) **pending**
- `z_differentiation`: Z{n x[n]} = -z dX/dz **pending**
- `fourier_transform`: defint(f*exp(-i*w*t), t, -oo, oo) through the integrator **pending**
- `fourier_transform_internal`: defint(f*exp(-i*w*t), t, -oo, oo) through the integrator **pending**
- `inverse_fourier_transform`: 1/(2 pi) * defint(F*exp(i*w*t), w, -oo, oo) **pending**
- `inverse_fourier_transform_internal`: 1/(2 pi) * defint(F*exp(i*w*t), w, -oo, oo) **pending**
- `laplace_transform`: defint(f*exp(-s*t), t, 0, oo) through the integrator **pending**
- `laplace_transform_internal`: defint(f*exp(-s*t), t, 0, oo) through the integrator **pending**
- `inverse_laplace_transform`: table lookup (1/(s-a), w/(s^2+w^2), ...), else partial fractions term by term, else Bromwich contour integral **pending**
- `inverse_laplace_transform_internal`: table lookup (1/(s-a), w/(s^2+w^2), ...), else partial fractions term by term, else Bromwich contour integral **pending**
- `z_transform`: simplified bilateral sum of x[n]*z^-n **pending**
- `z_transform_internal`: simplified bilateral sum of x[n]*z^-n **pending**
- `inverse_z_transform`: 1/(2 pi i) contour integral of X(z) z^(n-1) by path integration **pending**
- `inverse_z_transform_internal`: 1/(2 pi i) contour integral of X(z) z^(n-1) by path integration **pending**
- `partial_fraction_decomposition`: rational function split over distinct and repeated denominator roots **partial**: `rules::poly::apart::apart` over Q (tests decompositions_recombine, repeated_irreducible_quadratic); not yet wired as an operator, no symbolic coefficients
- `partial_fraction_decomposition_internal`: rational function split over distinct and repeated denominator roots **partial**: `rules::poly::apart::apart` over Q (tests decompositions_recombine, repeated_irreducible_quadratic); not yet wired as an operator, no symbolic coefficients
- `convolution_fourier`: convolution theorem: F{f*g} = F{f} F{g} **pending**
- `convolution_laplace`: convolution theorem: L{f*g} = F(s) G(s) **pending**

## Topology (other agent) (24)

### `src/numerical/topology.rs`

- `find_connected_components`: Finds the connected components of a graph using Breadth-First Search (BFS). **pending**
- `vietoris_rips_complex`: Constructs a Vietoris-Rips simplicial complex from a set of points for a given radius. **pending**
- `betti_numbers_at_radius`: Computes the Betti numbers for a point cloud at a given radius. **pending**
- `compute_persistence`: Computes the persistent homology (persistence diagram) for a point cloud. **pending**
- `euclidean_distance`: Computes the Euclidean distance between two points. **pending**

### `src/symbolic/topology.rs`

- `new`: Creates a new `Simplex` from a slice of vertex indices. **pending**
- `dimension`: Returns the dimension of the simplex. **pending**
- `boundary`: Computes the boundary of the simplex (numerical version). **pending**
- `symbolic_boundary`: Computes the boundary of the simplex (symbolic version). **pending**
- `add_term`: Adds a simplex with a given coefficient to the chain. **pending**
- `add_simplex`: Adds a simplex and all its faces (sub-simplices) to the complex. **pending**
- `get_simplices_by_dim`: Returns a reference to the vector of simplices of a specific dimension. **pending**
- `get_boundary_matrix`: Constructs the k-th boundary matrix `∂_k` for the simplicial complex. **pending**
- `get_symbolic_boundary_matrix`: Constructs the k-th symbolic boundary matrix `∂_k` for the simplicial complex. **pending**
- `apply_boundary_operator`: Applies the k-th boundary operator `∂_k` to a k-chain. **pending**
- `apply_symbolic_boundary_operator`: Applies the k-th symbolic boundary operator `∂_k` to a symbolic k-chain. **pending**
- `compute_euler_characteristic`: Computes the Euler characteristic `χ` of the simplicial complex. **pending**
- `verify_boundary_property`: Verifies the fundamental property of boundary operators: **pending**
- `verify_coboundary_property`: Verifies the fundamental property of coboundary operators: **pending**
- `compute_homology_betti_number`: Computes the k-th Betti number, `β_k`, which is a topological invariant. **pending**
- `compute_cohomology_betti_number`: Computes the k-th cohomology Betti number, `β^k`, which is a topological invariant. **pending**
- `create_grid_complex`: Creates a 2D grid simplicial complex. **pending**
- `create_torus_complex`: Creates a 2D torus simplicial complex. **pending**
- `vietoris_rips_filtration`: Creates a Vietoris-Rips filtration from a set of points in a metric space. **pending**

## Group theory and Lie algebras (24)

### `src/symbolic/discrete_groups.rs`

- `cyclic_group`: Creates a cyclic group `C_n` of order `n`. **pending**
- `dihedral_group`: Creates a dihedral group `D_n` of order `2n`. **pending**
- `symmetric_group`: Creates a symmetric group `S_n` of order `n!`. **pending**
- `klein_four_group`: Creates the Klein four-group `V_4`. **pending**

### `src/symbolic/group_theory.rs`

- `new`: Creates a new group. **pending**
- `multiply`: Multiplies two group elements. **pending**
- `inverse`: Computes the inverse of a group element. **pending**
- `is_abelian`: Checks if the group is abelian (commutative). **pending**
- `element_order`: Computes the order of an element g (smallest k such that g^k = e). **pending**
- `conjugacy_classes`: Finds the conjugacy classes of the group. **pending**
- `center`: Finds the center of the group Z(G) = {z in G | zg = gz for all g in G}. **pending**
- `is_valid`: Checks if the representation is valid (homomorphism property). **pending**
- `character`: Computes the character of a representation. **pending**

### `src/symbolic/lie_groups_and_algebras.rs`

- `lie_bracket`: Computes the Lie bracket `[X, Y] = XY - YX` for matrix Lie algebras. **pending**
- `lie_bracket_internal`: internal solver of `lie_bracket` **pending**
- `exponential_map`: Computes the exponential map `e^X` for a Lie algebra element `X` using a Taylor series expansion. **pending**
- `adjoint_representation_group`: Computes the adjoint representation of a Lie group element `g` on a Lie algebra element `X`. **pending**
- `adjoint_representation_algebra`: Computes the adjoint representation of a Lie algebra element `X` on another Lie algebra element `Y`. **pending**
- `commutator_table`: Computes the commutator table for a Lie algebra. **pending**
- `check_jacobi_identity`: Checks if the basis of a Lie algebra satisfies the Jacobi identity. **pending**
- `so3_generators`: Returns the basis generators for the `so(3)` Lie algebra (infinitesimal rotations). **pending**
- `so3`: Creates the `so(3)` Lie algebra. **pending**
- `su2_generators`: Returns the basis generators for the `su(2)` Lie algebra. **pending**
- `su2`: Creates the `su(2)` Lie algebra. **pending**

## Computer graphics (other agent) (22)

### `src/symbolic/computer_graphics.rs`

- `translation_2d`: Generates a 3x3 2D translation matrix. **pending**
- `translation_3d`: Generates a 4x4 3D translation matrix. **pending**
- `rotation_2d`: Generates a 3x3 2D rotation matrix. **pending**
- `rotation_3d_x`: Generates a 4x4 3D rotation matrix around the X-axis. **pending**
- `rotation_3d_y`: Generates a 4x4 3D rotation matrix around the Y-axis. **pending**
- `rotation_3d_z`: Generates a 4x4 3D rotation matrix around the Z-axis. **pending**
- `scaling_2d`: Generates a 3x3 2D scaling matrix. **pending**
- `scaling_3d`: Generates a 4x4 3D scaling matrix. **pending**
- `perspective_projection`: Generates a 4x4 perspective projection matrix. **pending**
- `orthographic_projection`: Generates a 4x4 orthographic projection matrix. **pending**
- `look_at`: Generates a 4x4 "look at" view matrix. **pending**
- `evaluate`: Evaluates the Bezier curve at a given parameter `t`. **pending**
- `derivative`: Computes the derivative (tangent vector) of the Bezier curve at parameter `t`. **pending**
- `split`: Splits the Bezier curve at parameter `t` using De Casteljau's algorithm. **pending**
- `new`: Creates a new polygon from a list of vertex indices. **pending**
- `apply_transformation`: Applies a geometric transformation to the entire mesh. **pending**
- `compute_normals`: Computes the surface normal vectors for each polygon in the mesh. **pending**
- `triangulate`: Triangulates all polygons in the mesh into triangles. **pending**
- `shear_2d`: Generates a 3x3 2D shear matrix. **pending**
- `reflection_2d`: Generates a 3x3 2D reflection matrix across a line through the origin. **pending**
- `reflection_3d`: Generates a 4x4 3D reflection matrix across a plane through the origin. **pending**
- `rotation_axis_angle`: Generates a 4x4 3D rotation matrix around an arbitrary axis using Rodrigues' formula. **pending**

## Special functions (18)

### `src/symbolic/special.rs`

- `inverse_erfc`: Computes the inverse complementary error function, `erfc⁻¹(x)`. **pending**
- `bessel_k0`: Computes the modified Bessel function of the second kind, K₀(x). **pending**
- `bessel_k1`: Computes the modified Bessel function of the second kind, K₁(x). **pending**
- `ln_factorial`: Computes the logarithm of the factorial, `ln(n!)`. **pending**

### `src/symbolic/special_functions.rs`

- `polygamma`: Symbolic representation for the Polygamma function, `ψ⁽ⁿ⁾(z)`. **partial**: `kernels::special::polygamma_numerical` only; no `polygamma` operator
- `erfi`: Symbolic representation and smart constructor for the Imaginary Error Function, `erfi(z)`. **pending**
- `bessel_j`: Symbolic representation and smart constructor for the Bessel function of the first kind, `J_n(x)`. **partial**: `rules::special` `besselj`/`bessely`/`besseli` (test bessel_values; integer orders 0 and 1 only)
- `bessel_y`: Symbolic representation and smart constructor for the Bessel function of the second kind, `Y_n(x)`. **partial**: `rules::special` `besselj`/`bessely`/`besseli` (test bessel_values; integer orders 0 and 1 only)
- `bessel_i`: Symbolic representation for the Modified Bessel function of the first kind, `I_n(x)`. **partial**: `rules::special` `besselj`/`bessely`/`besseli` (test bessel_values; integer orders 0 and 1 only)
- `bessel_k`: Symbolic representation for the Modified Bessel function of the second kind, `K_n(x)`. **pending**
- `generalized_laguerre`: Symbolic representation for the Generalized Laguerre Polynomials, `L_n^α(x)`. **pending**
- `bessel_differential_equation`: Represents Bessel's differential equation: **pending**
- `legendre_differential_equation`: Represents Legendre's differential equation: **pending**
- `legendre_rodrigues_formula`: Represents Rodrigues' Formula for Legendre Polynomials: **pending**
- `laguerre_differential_equation`: Represents Laguerre's differential equation: **pending**
- `hermite_differential_equation`: Represents Hermite's differential equation: **pending**
- `hermite_rodrigues_formula`: Represents Rodrigues' Formula for Hermite Polynomials: **pending**
- `chebyshev_differential_equation`: Represents Chebyshev's differential equation: **pending**

## Calculus, integration and series (17)

### `src/numerical/convergence.rs`

- `aitken_acceleration`: Accelerates the convergence of a sequence using Aitken's delta-squared process. **pending**
- `find_sequence_limit`: Numerically finds the limit of a sequence by generating terms and applying acceleration. **pending**
- `richardson_extrapolation`: Performs Richardson extrapolation on a sequence of approximations. **partial**: private `richardson` inside `kernels::series::sum_to_infinity` (sums only, doubling checkpoints)
- `wynn_epsilon`: Applies Wynn's epsilon algorithm to accelerate the convergence of a sequence. **partial**: private `wynn_epsilon` inside `kernels::series::sum_to_infinity`; not callable on its own

### `src/numerical/series.rs`

- `evaluate_power_series`: Evaluates a power series at a point given its coefficients and center. **partial**: `kernels::polynomial::Polynomial::eval` for finite coefficient lists; no centre argument

### `src/symbolic/calculus.rs`

- `check_analytic`: Checks if a complex function `f(z)` is analytic by verifying the Cauchy-Riemann equations. **pending**
- `find_poles`: Finds the poles of a rational expression by solving for the roots of the denominator. **pending**
- `calculate_residue`: Calculates the residue of a complex function at a given pole. **partial**: `kernels::complex::residue` numeric closure kernel only
- `is_inside_contour`: This function is used in complex analysis, particularly with the Residue Theorem, to determine which poles of a function lie within a given integration path. **pending**
- `path_integrate`: Computes a path integral of a complex function over a given contour. **partial**: `rules::linalg` `line_integral` for real parametrised curves; `kernels::complex::contour_integral` numerically
- `improper_integral`: Calculates an improper integral from -infinity to +infinity using the residue theorem. **partial**: `defint` over `oo` limits, numeric only via `Quadrature` (test numeric_quadrature_without_a_closed_form); no residue-theorem symbolic result

### `src/symbolic/convergence.rs`

- `analyze_convergence`: Analyzes the convergence of a series given its general term `a_n`. **partial**: `rules::calculus` `converges` (ratio test only)

### `src/symbolic/integration.rs`

- `risch_norman_integrate`: Main entry point for Risch-Norman style integration. **partial**: `rules::calculus` `integral` staged heuristics; no Risch-Norman undetermined-coefficients ansatz
- `integrate_poly_exp`: Integrates the polynomial part of a transcendental function extension F(t). **partial**: `integral` by_parts stage handles polynomial*exp(ax) (test by_parts); no general exponential extension

### `src/symbolic/series.rs`

- `analyze_convergence`: Analyzes the convergence of a series using the Ratio Test. **partial**: `rules::calculus` `converges` (ratio test only)
- `asymptotic_expansion`: Computes the asymptotic expansion of an expression around a given point (e.g., infinity). **pending**
- `analytic_continuation`: Performs analytic continuation of a function represented by a power series. **pending**

## Polynomials, factorisation and algebraic geometry (17)

### `src/symbolic/cad.rs`

- `cad`: Computes the Cylindrical Algebraic Decomposition for a set of polynomials. **pending**

### `src/symbolic/cas_foundations.rs`

- `risch_integrate`:  **partial**: `rules::calculus` `integral` heuristic stages (no Risch algorithm)
- `cylindrical_algebraic_decomposition`:  **pending**
- `simplify_with_relations`: Simplifies an expression using a set of polynomial side-relations. **pending**
- `simplify_with_relations_internal`: internal solver of `simplify_with_relations` **pending**
- `normalize_with_relations`: Normalizes an expression to a canonical form using a set of polynomial side-relations. **pending**

### `src/symbolic/poly_factorization.rs`

- `factor_gf`: Factors a polynomial over a finite field. **partial**: private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API
- `poly_derivative_gf`: Computes the derivative of a polynomial over a prime field. **partial**: private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API
- `square_free_factorization_gf`: Performs square-free factorization of a polynomial over a prime field. **pending**
- `berlekamp_factorization`: Factors a square-free polynomial over a small prime field using Berlekamp's algorithm. **pending**
- `cantor_zassenhaus`: Factors a square-free polynomial over a large prime field using Cantor-Zassenhaus algorithm. **partial**: private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API
- `distinct_degree_factorization`: Performs Distinct-Degree Factorization (DDF) of a polynomial over a finite field. **partial**: private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API
- `poly_gcd_gf`: Computes the greatest common divisor (GCD) of two polynomials over a prime field. **partial**: private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API
- `poly_pow_mod`: Computes base^exp mod modulus for polynomials over a prime field. **partial**: private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API
- `poly_mul_scalar`: Helper to multiply a polynomial by a scalar `BigInt`. **partial**: private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API
- `poly_extended_gcd`: Polynomial Extended Euclidean Algorithm for `a(x)s(x) + b(x)t(x) = gcd(a(x), b(x))`. **partial**: private GF(p) helper inside `rules::poly::univariate`, reachable only through `factor` over Q; no public GF(p) API

### `src/symbolic/polynomial.rs`

- `leading_term`: Returns the term with the highest degree in the specified variable. **partial**: `poly::repr::Poly::terms` ordered by monomial (no dedicated accessor)

## Cryptography (other agent) (15)

### `src/symbolic/cryptography.rs`

- `is_infinity`: Returns true if this point is the point at infinity. **pending**
- `x`: Returns the x-coordinate if this is an affine point, None if infinity. **pending**
- `y`: Returns the y-coordinate if this is an affine point, None if infinity. **pending**
- `new`: Creates a new elliptic curve y^2 = x^3 + ax + b over GF(p). **pending**
- `is_on_curve`: Checks if a point is on the curve. **pending**
- `negate`: Negates a point on the curve (P -> -P). **pending**
- `double`: Doubles a point on the curve (2P). **pending**
- `add`: Adds two points on the curve. **pending**
- `scalar_mult`: Performs scalar multiplication (`k * P`) using the double-and-add algorithm. **pending**
- `generate_keypair`: Generates a new ECDH (Elliptic Curve Diffie-Hellman) key pair. **pending**
- `generate_shared_secret`: Generates a shared secret using one's own private key and the other party's public key. **pending**
- `point_compress`: Compresses a curve point to its x-coordinate and a sign bit. **pending**
- `point_decompress`: Decompresses a point from its x-coordinate and sign bit. **pending**
- `ecdsa_sign`: Signs a message hash using ECDSA. **pending**
- `ecdsa_verify`: Verifies an ECDSA signature. **pending**

## Integral equations and calculus of variations (14)

### `src/numerical/calculus_of_variations.rs`

- `evaluate_action`: Evaluates the action of a functional for a given path. **pending**
- `euler_lagrange`: d/dx(dL/dy') - dL/dy via diff **pending**

### `src/symbolic/calculus_of_variations.rs`

- `euler_lagrange`: d/dx(dL/dy') - dL/dy via diff **pending**
- `euler_lagrange_internal`: d/dx(dL/dy') - dL/dy via diff **pending**
- `solve_euler_lagrange`: feeds the Euler-Lagrange ODE to the ODE solver **pending**
- `solve_euler_lagrange_internal`: feeds the Euler-Lagrange ODE to the ODE solver **pending**
- `hamiltons_principle`: Euler-Lagrange equations for the action of a Lagrangian **pending**

### `src/symbolic/integral_equations.rs`

- `new`: Creates a new instance of a Fredholm integral equation of the second kind. **pending**
- `solve_neumann_series`: Solves the Fredholm integral equation using the method of successive approximations (Neumann Series). **pending**
- `solve_separable_kernel`: Solves a Fredholm integral equation of the second kind with a separable (or degenerate) kernel. **pending**
- `solve_successive_approximations`: Solves the Volterra integral equation using the method of successive approximations. **pending**
- `solve_by_differentiation`: Solves the Volterra equation by converting it into an Ordinary Differential Equation (ODE). **pending**
- `solve_airfoil_equation`: Solves the airfoil singular integral equation. **pending**
- `solve_airfoil_equation_internal`: internal solver of `solve_airfoil_equation` **pending**

## Output formats and plotting (13)

### `src/output/latex.rs`

- `to_latex`: Converts an expression to a LaTeX string. **pending**
- `to_latex_prec_with_parens`: Helper to add parentheses if needed. **pending**
- `to_greek`: Converts common Greek letter names to LaTeX. **pending**

### `src/output/plotting.rs`

- `plot_function_2d`: Plots a 2D function y = f(x) and saves it to a file. **pending**
- `plot_series_2d`: Plots multiple data series on a single 2D plot and saves it to a file. **pending**
- `plot_vector_field_2d`: Plots a 2D vector field and saves it to a file. **pending**
- `plot_surface_3d`: Plots a 3D surface z = f(x, y) and saves it to a file. **pending**
- `plot_surface_2d`: Plots a 2D array as a 3D surface plot and saves it to a file. **pending**
- `plot_parametric_curve_3d`: Plots a 3D parametric curve (x(t), y(t), z(t)) and saves it to a file. **pending**
- `plot_vector_field_3d`: Plots a 3D vector field and saves it to a file. **pending**
- `plot_3d_path_from_points`: Plots a 3D path from a series of points (x, y, z) and saves it to a file. **pending**
- `plot_heatmap_2d`: Plots a 2D heat map from a 2D array of data and saves it to a file. **pending**

### `src/output/typst.rs`

- `to_typst`: Converts an expression to a Typst string. **pending**

## Fractals and chaos (other agent) (12)

### `src/symbolic/fractal_geometry_and_chaos.rs`

- `new`: Creates a new Iterated Function System. **pending**
- `apply`: Applies the IFS to a point (symbolically) to generate the set of possible next points. **pending**
- `similarity_dimension`: Calculates the similarity dimension for a self-similar IFS. **pending**
- `new_mandelbrot_family`: Creates a new Mandelbrot/Julia system z -> z^2 + c. **pending**
- `iterate`: Iterates the system once: **pending**
- `orbit`: Computes the orbit of a point up to n iterations. **pending**
- `fixed_points`: Finds fixed points of the system: **pending**
- `stability_index`: Checks the stability of a fixed point z*. **pending**
- `find_fixed_points`: Calculates the fixed points of a 1D map f(x). **pending**
- `analyze_stability`: Analyzes the stability of a fixed point for a 1D map f(x). **pending**
- `lyapunov_exponent`: Calculates the symbolic Lyapunov exponent for a 1D chaotic map `x_{n+1} = f(x_n)`. **pending**
- `lorenz_system`: Returns the Lorenz System equations. **pending**

## Indefinite sums and products (11)

### `src/numerical/indefinite_sum.rs`

- `try_closed_form_sum`: Attempts to compute the indefinite sum of `f(var)` in closed form. **pending**
- `eval_antidiff`: Evaluates the symbolic anti-difference result `F(x_val)` numerically. **pending**
- `eval_normalized`: Evaluates the anti-difference with the Nörlund normalization: **pending**
- `new`: Creates a new Abel-Plana engine for the given expression. **pending**
- `eval`: Evaluates the Abel-Plana formula for the anti-difference at `x`. **pending**
- `eval_indefinite_product_numerical`: Computes the indefinite product F(x) = ∏_{k} f(k) numerically. **pending**
- `series_antidiff`: Computes Δ⁻¹ f(z) via a Taylor series expansion of f around a point `p`, applying the Hurwitz zeta formula termwise: **pending**
- `compute_taylor_coeffs_numerical`: Computes Taylor coefficients c_m = f^(m)(p)/m! numerically via forward differences. **pending**
- `eval_indefinite_sum_numerical`: Evaluates the indefinite sum F(x) - F(h) using a three-strategy cascade: **pending**
- `expr_contains_var`: Returns `true` if the expression contains the named variable. **pending**
- `extract_linear_coeff`: Attempts to extract the coefficient `a` if `expr = a * var` (or just `var`). **pending**

## Symbolic ODEs (8)

### `src/symbolic/ode.rs`

- `solve_ode_system`: Solves a system of coupled ordinary differential equations via E-Graph. **partial**: `rules::ode` `odeint` (numeric systems only; no symbolic system solver)
- `solve_ode_system_internal`: internal solver of `solve_ode_system` **partial**: `rules::ode` `odeint` (numeric systems only; no symbolic system solver)
- `solve_by_reduction_of_order`: Solves a second-order homogeneous linear ODE by reduction of order. **pending**
- `solve_by_reduction_of_order_internal`: internal solver of `solve_by_reduction_of_order` **pending**
- `solve_ode_by_series`: Solves an Ordinary Differential Equation using the power series method. **pending**
- `solve_ode_by_series_internal`: internal solver of `solve_ode_by_series` **pending**
- `solve_ode_by_fourier`: Solves a linear Ordinary Differential Equation using the Fourier Transform method via E-Graph. **pending**
- `solve_ode_by_fourier_internal`: internal solver of `solve_ode_by_fourier` **pending**

## Simulation scenarios and numeric solvers (test gaps) (7)

### `src/physics/physics_sim/geodesic_relativity.rs`

- `simulate_black_hole_orbits_scenario`: An example scenario that simulates several types of orbits around a black hole. **partial**: `sim::models::geodesic_relativity::simulate_black_hole_orbits_scenario` (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested)

### `src/physics/physics_sim/gpe_superfluidity.rs`

- `simulate_bose_einstein_vortex_scenario`: An example scenario that finds the ground state of a BEC, which may contain a vortex. **partial**: `sim::models::gpe_superfluidity::simulate_bose_einstein_vortex_scenario` (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested)

### `src/physics/physics_sim/ising_statistical.rs`

- `simulate_ising_phase_transition_scenario`: An example scenario that simulates the Ising model across a range of temperatures to observe the phase transition. **partial**: `sim::models::ising_statistical::simulate_ising_phase_transition_scenario` (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested)

### `src/physics/physics_sim/linear_elasticity.rs`

- `simulate_cantilever_beam_scenario`: An example scenario for a cantilever beam under a point load. **partial**: `sim::models::linear_elasticity::simulate_cantilever_beam_scenario` (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested)

### `src/physics/physics_sim/schrodinger_quantum.rs`

- `simulate_double_slit_scenario`: An example scenario simulating a wave packet hitting a double slit. **partial**: `sim::models::schrodinger_quantum::simulate_double_slit_scenario` (no test: scenario writes .csv/.npy files into the working directory; underlying run_* function is tested)

### `src/physics/physics_sm.rs`

- `ifft3d`: Performs a 3D IFFT. **partial**: `sim::physics_sm::ifft3d` (only covered by an #[ignore]d test that exposes a library bug)
- `solve_advection_diffusion_3d`: Solves the 3D advection-diffusion equation. **partial**: `sim::physics_sm::solve_advection_diffusion_3d` (only covered by an #[ignore]d test that exposes a library bug)

## Verification (proof) (7)

### `src/symbolic/proof.rs`

- `verify_equation_solution`: Verifies a solution to a single equation or a system of equations using numerical sampling. **pending**
- `verify_indefinite_integral`: Verifies an indefinite integral `F(x)` for an integrand `f(x)` by checking if `F'(x) == f(x)`. **pending**
- `verify_definite_integral`: Verifies a definite integral by comparing the symbolic result with numerical quadrature. **pending**
- `verify_ode_solution`: Verifies a solution to an ODE `G(x, y, y', y'', ...) = 0` by numerical sampling. **pending**
- `verify_matrix_inverse`: Verifies a matrix inverse `A⁻¹` by checking if `A * A⁻¹` is the identity matrix. **pending**
- `verify_derivative`: Verifies a symbolic derivative `f'(x)` by comparing it to a numerical differentiation. **pending**
- `verify_limit`: Verifies a symbolic limit `lim_{x->x0} f(x) = L`. **pending**

## Linear algebra, statistics and combinatorics (5)

### `src/nightly/matrix.rs`

- `with_backend`: Sets the backend for the matrix. **partial**: `kernels::matrix::Matrix::with_backend` (no direct test; used by the faer paths)

### `src/numerical/stats.rs`

- `skewness`: Computes the skewness of a slice of data. **partial**: `kernels::stats::skewness` (only covered by an #[ignore]d test that exposes a library bug)

### `src/symbolic/combinatorics.rs`

- `expand_binomial`: Expands an expression of the form `(a+b)^n` using the Binomial Theorem. **partial**: `expand((a+b)^n)` of the poly rules, literal n (symbolic n needs a sum operator)

### `src/symbolic/matrix.rs`

- `svd_decomposition`: Performs Singular Value Decomposition (SVD) of a matrix `A`. **partial**: `rules::linalg` `svd` (numeric, faer)

### `src/symbolic/stats_regression.rs`

- `polynomial_regression_symbolic`: Computes the symbolic coefficients for a polynomial regression `y = c0 + c1*x + ... **partial**: `polynomial_regression` (exact/float data only, symbolic data not supported)
