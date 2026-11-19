# Mathematical coverage and input policies

[中文](MATH_COVERAGE.zh.md) · [Embedding kernel](KERNEL.md)

Choose an algorithm by its mathematical domain, data representation and error contract. The live manifest lists the operators in a particular build. Rust algorithms and dynamic operators may expose different interfaces to the same implementation.

## Numerical inputs

The checked ABI reports malformed shapes and unsupported domains as an error or a documented IEEE result. This includes empty tensors, inconsistent extents and storage, nonfinite values and out-of-range control parameters. Low-level Rust accessors may require valid indices. Host callbacks have their own Rust panic contract.

| Area | Input and result policy |
|---|---|
| Dense linear algebra | Checked kernels distinguish nonfinite data, shape violations and singular factorizations. `Lu::solve_into` reuses a caller-provided output of the exact dimension and requires a finite right-hand side. Use `pinv`, `lstsq_svd`, `null_space` or `orth` with an appropriate singular-value threshold for rank-deficient problems. |
| Matrix functions | Symmetric `sqrtm` and `logm` use spectral methods. A symmetric positive semidefinite matrix can have a singular square root. The general real methods use Denman–Beavers or inverse scaling and squaring, with failures reported when the algorithm cannot produce the real principal result. |
| Statistics | Sample summaries and tests follow their per-operator missing and empty data policies. Use `percentile` with `p` in `[0,100]`. Survival-function paths retain small tail probabilities without subtracting a CDF rounded to one. Weighted empirical summaries require finite nonnegative weights and positive total weight. Zero-weight samples contribute no mass. |
| Real and complex special functions | Use documented values at poles and on branch cuts. Invalid numeric arguments can yield NaN. Complex `ln_gamma` is an analytic branch. `LambertW` takes a branch index. Complex Bessel J and I accept real order. Fresnel and the complex error-function family preserve scaled tails near floating-point limits. `bessel_i0e` removes the dominant exponential from real modified Bessel I₀ for large arguments. The Kaiser window uses this form to keep its normalization finite. |
| Exact integer algorithms | The 64-bit modular routines use exact wider intermediates. `is_prime` is deterministic for `u64`, and factorization uses Pollard–Brent. `checked_lcm` reports overflow. The discrete logarithm reduces common factors before the coprime baby-step giant-step search, whose table can require space proportional to the square root of the reduced modulus. The symbolic rule set supplies arbitrary-precision algorithms. |
| Quadrature | Match the method to the integrand. Tanh–sinh handles integrable endpoint singularities. Adaptive Gauss–Kronrod supplies an error estimate for smooth integrands. `levin_adaptive` targets oscillatory integrals with a supplied phase derivative. Simpson and Romberg enforce work bounds. |
| ODE, SDE and roots | Solvers expose method-specific evaluation or step bounds and typed failures. Stiff deterministic problems have implicit methods. SDE methods include drift-implicit options, with different admissible noise structures and a choice of Itô or Stratonovich calculus. Brent requires a sign bracket or an endpoint root. Read the exact solver signature for tolerances and returned diagnostics. |

## Structured and exact domains

- `kernels::dense` and `kernels::schur` provide real Schur, generalized Schur, Sylvester and Lyapunov paths. `kernels::lowrank` provides column-pivoted QR, randomized SVD, LOBPCG and Anderson acceleration. `lobpcg` solves ordinary symmetric problems. `lobpcg_generalized` accepts matrix-free stiffness and positive-definite mass operators, returns mass-orthonormal eigenvectors, and can skip search updates for converged pairs while retaining them in Rayleigh–Ritz.
- `kernels::sparse_exact` and the symbolic `sparse_symbolic` family cover elimination over rational, finite and symbolic fields, CRT determinants with a Hadamard bound, Dixon lifting with rational reconstruction and scalar Wiedemann methods. The modular minimal polynomial uses randomized projections. Rank over rational inputs uses rational elimination. Choose the dense exact algorithms when the task needs Smith or Hermite normal forms.
- `kernels::holonomic` and `holonomic_*` cover annihilators, closure under addition, multiplication, powers and differentiation, coefficient recurrences and exact zero tests at an ordinary expansion point. The supported coefficient path uses a single variable and rational Taylor data at zero. Definite sums use the separate telescoping rules. A generated closure annihilator satisfies the requested relation and can have higher order than a global minimal annihilator.
- Numeric finite fields include prime `GF(p)` and `GF(256)`. The symbolic `gf_*`, `gfx_*` and `gfp_factor` families add extension fields and symbolic factorization. The constructor and operations check their own field invariants.
- Numbered Bessel J and Y use recurrence in suitable parameter regions. Complex asymptotic expansions can lose relative precision close to zeros or with large order relative to the argument. Where residual or component-level accuracy matters, check it against independent values in the regime of interest.


## Signal processing

- `kernels::signal` implements orthonormal DCT-II and DST-II with their inverses using cached FFT plans. The legacy `signal::fft` and `signal::ifft` entry points and FFT convolution use the same plan and scratch pools as `kernels::fft`. Multidimensional real FFTs shrink the Hermitian axis before processing the other axes. Large independent FFT batches and STFT frames can run in parallel with deterministic ordering.
- `welch_psd` averages complete, optionally overlapping windowed frames into a one-sided power spectral density. `welch_cross_spectrum` also returns cross-spectral density using `conj(X) * Y` and magnitude-squared coherence. Coherence is zero in bins where either input has zero power. `inverse_stft` reconstructs a covered interval from full STFT spectra and rejects any sample with no nonzero window contribution. `upfirdn` applies explicit FIR taps without materializing the upsampled signal.
- The typed `num_signal_welch_psd` and `num_signal_welch_cross_spectrum` records carry separate dimensions for frequencies, density and coherence. Frequency bins use the sample rate units. Density is measured per unit frequency, and coherence is dimensionless.
- `analytic_signal`, `hilbert` and `analytic_envelope` use periodic boundaries. `zoom_fft` evaluates uniformly spaced frequencies in cycles per input sample via chirp convolution. The typed `num_signal_zoom_fft` operator takes and returns a vector of complex samples.
- `kernels::wavelets::wavelet_packet` decomposes both low and high branches into equal-width leaf bands, with the low branch preceding the high branch. `inverse_wavelet_packet` reconstructs a signal from that leaf matrix. `FrftPlan` reuses the eigenbasis of the discrete fractional Fourier transform for multiple orders on signals of the same length.

## Weighted graphs

`kernels::graph::Graph` stores directed arcs, including self-loops and parallel arcs. `Graph::from_edges` validates a batch and reserves each adjacency list before insertion. Model an undirected weighted graph with two equally weighted arcs per edge.

| Operation | Contract |
|---|---|
| Dijkstra | Single-source shortest paths with nonnegative weights. |
| Bellman–Ford | Single-source shortest paths with negative weights, rejecting negative cycles reachable from the source. |
| DAG shortest paths | Single-source shortest paths in linear time on directed acyclic graphs, including negative weights. |
| Floyd–Warshall | All-pairs shortest paths in row-major order. A negative matrix diagonal signals a negative cycle. |
| Topological sorting, strong and weak components | Edge weights are ignored. Strong component labels follow topological order. Weak component labels follow the smallest vertex in each component. The symbolic SCC operator shares the iterative strong-component traversal. |
| Prim's minimum spanning forest | Preserves all components and isolated vertices of an undirected input represented by symmetric arcs. |
| PageRank | Counts each outgoing arc separately and redistributes dangling mass uniformly. |
| Dinic maximum flow and minimum cut | Uses finite nonnegative capacities and distinct terminals. Returns the flow value and a Boolean mask for the source side of a minimum cut. |

The `num_native_graph_*` handles retain the physical dimension of edge weights. Bulk construction takes a two-column endpoint matrix and a rank-one weight tensor. Path records return `distances` in coherent edge units and `previous` as optional predecessor indices. Unreachable distances are positive infinity. Maximum flow returns its `value` in coherent capacity units and a dimensionless `source_side` mask. Native weighted graphs allow up to 256 vertices and 65,536 arcs. The maximum-flow call further limits `V² E` to 50 million. Consult the manifest for the callable names and typed schemas.

## Low-rank matrices

`kernels::lowrank::interpolative_decomposition` runs truncated column-pivoted QR to approximate a matrix using its original columns and a dimensionless coefficient matrix. The actual rank can fall below the requested maximum on rank-deficient input. `num_linear_interpolative_decomposition` accepts a native matrix handle and returns the selected column indices and the coefficients in original column order.

## Dimensions and grids

Symbolic heat, wave, Laplace and Poisson equations use one operator per family. The problem determines its spatial dimension, and compatibility aliases include `_1d`, `_2d` and `_3d`. Finite-difference families such as `heat_fd`, `wave_fd`, `poisson_fd` and `advect_fd` infer spatial rank from the input shape on the unit box. Boundary-aware variants accept Dirichlet or Neumann faces. Pure-Neumann Poisson checks source and outward-flux compatibility and returns a zero-mean result. The method-of-lines helper targets one-dimensional grids. Consult the explicit method for its supported boundary conditions.

`sim::dimensional` uses x-fastest flat fields with axes `[nx]`, `[nx, ny]` or `[nx, ny, nz]`. The shared transport paths support advection, diffusion and their combination. The Poisson paths provide relaxation and multigrid. Other solvers include wave, Euler, geometric fields and scenario-specific simulation drivers. Scenario dimensions, field layouts and method options are part of the live scenario catalogue. Verification compares against analytical solutions, dimensional reductions, symmetry, convergence order or conservation, using `tests/sim/dimension_coverage.rs` and `tests/kernels/dimension_coverage.rs`.

## Selecting a callable interface

Literal term operators under `rules::numeric` accept numbers, lists and method symbols. An unresolved input or failed precondition can leave the call symbolic. Use `Session::apply_diag` to inspect reported numerical rejections. The expression-level families cover robust statistics, Fourier and wavelet transforms, `coherence`, `wavelet_packet`, `inverse_wavelet_packet`, sparse iterative solves, exact integer convolution, tensor decompositions, sampled functions, polynomial roots and grid solvers. Where indices enter the interface they are zero-based. The native numerical manifest exposes the cosine and sine transforms, analytic signal, Hilbert transform, amplitude envelope and zoom FFT under `num_signal_*`.

The typed numerical ABI additionally exposes checked records, tensors and session-owned models. It supports reusable CHOL/LU and multigrid state, tensor-train algebra, complex transform data, checked homotopy certificates and ODE/SDE objects with solver-specific controls. Discover the exact callable names and signatures in `KernelManifest` or Aevia's `aevia rssn` output. Native methods that accept Rust callbacks or borrow host storage also remain available directly. [Kernel embedding contract](KERNEL.md) · [Numerical callbacks](NUMERICAL_CALLBACKS.md)

## Verification map

| Area | Representative tests |
|---|---|
| Dense linear algebra, matrix functions and null spaces | `tests/math_linear_algebra.rs`, `tests/kernels/dense.rs` |
| Complex functions and real tails | `tests/math_special_complex.rs`, `tests/math_special_edges.rs`, `kernels::fresnel` regressions |
| Primality, factors, modular roots and discrete logarithms | `tests/math_number_theory.rs` |
| Statistics, quadrature, interpolation and inverse problems | `tests/math_numeric_edges.rs`, `tests/math_robustness.rs`, `tests/kernels/robust.rs` |
| ABI domain, shape and work-budget behavior | `tests/abi_fuzz.rs`, `tests/analysis_native.rs`, `tests/source_policy_math.rs` |
| Integrals, limits, sums, series and symbolic solving | `tests/math_symbolic_classics.rs`, `tests/symbolic_battery.rs` |
| Scientific layouts and invariants | `tests/native_*.rs`, `tests/dimensioned_coverage.rs` |
| Signal analysis, FIR resampling, directed graph traversal | `tests/numerical_signal_graph.rs` |
| Real tensor FFTs, cosine and sine transforms, Hilbert, zoom FFT, coherence and wavelet packets | `tests/signal_transforms.rs` |
| Weighted paths, minimum cuts, graph isomorphism and native graph schemas | `tests/kernels/graph.rs`, `tests/native_discrete.rs` |
| Low-rank approximation and units of interpolation coefficients | `src/kernels/lowrank.rs`, `tests/native_linear.rs` |

The symbolic battery distinguishes exact-value checks from documented requests that remain unresolved in symbolic mode. Use the named test as the source of its expected behavior. Numerical answers also depend on precision, conditioning and the applicable error thresholds.
