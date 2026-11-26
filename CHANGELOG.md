# RSSN changelog · 变更日志

Release changes to the RSSN kernel, algorithms and embedding APIs. 本文记录 RSSN 内核、算法和嵌入 API 的版本变化。

## [3.0.0-beta.1]

### Engine and embedding

- The scientific kernel exposes a deterministic manifest, type and dimension contracts, structured results and typed numerical callbacks. Embedded hosts generate checked interfaces from the operator registry that executes the calls.
- Indexed rewrites and incremental matching reduce repeated search. Independent requests use private computation graphs with ordered result import. Optional JIT and device backends define their supported precision and fallback behavior.
- Workers prepare and evaluate symbolic batches on private graphs with bounded in-flight jobs. CPU dispatch reads the active Rayon pool size, and FFT routines share plans and scratch storage through `kernels::reuse`.
- Symbolic windows reuse extraction work across overlapping terms. Compiled CPU batches lease sized scratch buffers and skip parallel dispatch for small work. Concurrent FFT requests share a plan build for the same length and direction.
- Column-wise tapes read the selected source column for each gather. A prepared gather tape can take measured source columns with a common reference and return values in their original dimension.
- Prepared scalar programs can take one measured column per input, each with its own unit. Large batches run on the active CPU pool in fixed, ordered blocks and honor the embedding cancellation scope.
- Standard sessions share installed rules and their manifests until the global unit table changes. Per-session units, term graphs, native objects, and computation caches remain session-owned.
- The C ABI captures the shared standard manifest and reuses its JSON encoding. Each call still returns a separately owned string.
- The C ABI gains a host-callback surface: `rssn_callback_new` registers a C function with an explicit signature, `rssn_callback_new_for` copies the signature the operator declares at an argument position, `rssn_callback_signatures` lists what a position accepts, and `rssn_callback_free` retires a registration into a reusable slot. `rssn_value_call` now runs with the session's callback bank, so a `Value::Function` argument reaches the higher-order operators (`num_bfgs`, `num_brent`, the adaptive quadrature family, `num_bootstrap`, ...) from every binding instead of only from Rust. All 5411 operators are now reachable through the C ABI; the `backend_support.ffi` flag no longer excludes any of them. A callback reads its arguments as cells with the ordinary `rssn_value_*` readers, which now look through a `dyn`/`quantity` wrapper, and `rssn_value_new_bool`, `_new_complex`, `_new_record`, `_field_name` and `_tag` complete the readers and constructors it needs.
- Native object bindings derive runtime arity from the signature's argument list. Dense value preflight checks bounds without allocating a traversal buffer.
- Binary manifest, plugin and compute records use native bincode-next 3.2 traits. Signature and manifest bytes remain deterministic. Manifest verification reads the existing operator and unit tables in place and streams the canonical encoding into the hash. Compatibility tests compare the supported previous wire layouts with the native encoding. JSON and CBOR expose metadata to self-describing consumers.
- `signatures.overlay.json` now holds exactly the operators whose rule set declares no signature. The entries that repeated a declared overload are gone: they were invisible to the manifest (a declaration wins) and they broke the bijection `tests/signature_overlay.rs` enforces, which now fails on an entry for an unregistered operator and on an entry the manifest never reads.
- The manifest completes an operator's description from the rustdoc of the function behind it across generic type parameters and re-exported paths (`FdMode<F>::parallel`, `crate::sim::md_engine::FlexibleBarostat::energy`), so 551 more operators publish the description their implementation carries next to their family's contract sentence. `doc_path` names the function that supplied it, and an ambiguous path still resolves to nothing.
- `category` is never empty: an operator whose registration names none is grouped under the rule set that registered it. `tests/manifest_completeness.rs` states the invariants a host relies on — name, rule set, category, documentation, cost, keys, overload sources, and typed-or-explained — for every operator in the manifest.
- One shipped witness, `is_lie_algebra(sl_basis(4))`, no longer replays: the call answers `true` but takes 20.5 s against the miner's 20 s budget. It is dropped from `signatures.witnessed.json` rather than kept as evidence that cannot be reproduced, so `signature-mine --check` replays all 1047 witnesses of 727 operators. The overload it belonged to keeps its other witness.
- Four operators now state a result sharper than their symbolic contract — `arg(term) -> float`, `boltzmann_entropy(term) -> real`, `pde_order(term, term) -> int`, `todd_coxeter_index(term, term, term) -> int` — as *refinements*, not overloads: their parameters are the contract's own, so a second overload would make every call ambiguous. `infer` ignores refinements by design (`infer.rs`); a host with a checked result channel reads them, which is the only way a caller who passes handles learns that `arg` answers with a float rather than a term. Each entry carries the call that produced it, and `signature-mine --check` replays all 1051 witnesses of 731 operators.
- 132 operators whose handler builds one fixed `Value` kind on every success path now publish that result instead of only the symbolic contract: `graph_node_id(Obj<graph>, Int) -> Int`, `lie_rank(Symbol) -> Int`, `perm_group_is_abelian([[Int]]) -> Bool`, `gfp_divmod([Int], [Int], Int) -> [[Int]]`, `point_group_order(Symbol) -> Int`, `bdd_size([[Int]]) -> Int`, `crt([Int], [Int]) -> Int`, `pauli_x() -> [[Int]]`, `sc_homology_mod(Obj<simplicial_complex>, Int, Int) -> Int`. Every specific overload rests on two independent sources: the code path that builds the value, and a call the engine answered with that kind. A specific overload whose parameters would be indistinguishable from the symbolic contract is refused — that sharper result belongs in `refinements`, where it cannot make a call ambiguous — and an operator whose calls disagree with the code (it answers `false` on one path and a list on another) is left on the symbolic contract. The symbolic `(Term, …) -> Term` overload always stays last, so a call the operator does not accept still types as a term. Four of the 132 (`gf_inv`, `gfp_isirreducible`, `perm_order`, `perm_group_order`) state exactly the signature the Aevia bridge had curated by hand.

### Coverage: optimisation, statistics, control, meshes, data and precision

- Add `kernels::convex`: linear programs with two-sided bounds, convex quadratic programs by a dual active-set method, second-order cone programs by a Nesterov–Todd predictor–corrector interior point method, and semidefinite programs by ADMM with a Jacobi symmetric eigensolver. The linear solver is the existing tested two-phase simplex, so a bound shift and the upper-bound rows are the only new code on that path.
- Add `kernels::milp`: mixed-integer linear programming by best-bound branch and bound with pseudo-cost branching over that simplex, and a sequential quadratic program whose subproblem is solved by projected conjugate gradients on the working set. The subproblem never forms a dense factorisation, so a sparse problem costs `O(nnz * iterations)` per QP instead of the dense `O(n^3)` of `kernels::optim::sqp`, which stays as the reference for small dense problems.
- Add `kernels::glm` (Gaussian, Poisson, binomial, gamma and inverse-Gaussian families with seven links, IRLS with step halving, deviance, AIC and standard errors) and `kernels::timeseries` (autocorrelation, partial autocorrelation, Yule–Walker AR, Hannan–Rissanen ARIMA with a conditional-sum-of-squares refinement and forecasting, additive and multiplicative Holt–Winters, Kalman filtering and Rauch–Tung–Striebel smoothing).
- Add `kernels::survival` (Kaplan–Meier with Greenwood variance and a complementary-log-log band, Nelson–Aalen, the log-rank test, and Cox regression by Newton–Raphson with Breslow or Efron ties and Harrell's concordance) and `kernels::mcmc` (random-walk Metropolis–Hastings, Hamiltonian Monte Carlo with leapfrog integration, a Gibbs sweep, a deterministic PCG32 generator, split R-hat, effective sample size, autocorrelation and the Geweke score).
- Add `kernels::control`: state space with controllability, observability, pole placement checks and zero-order-hold discretisation by scaling and squaring, transfer functions with series, parallel, feedback, Bode, Nyquist and margin computation, continuous and discrete Riccati equations with LQR, LQG and Kalman gains, a condensed-QP MPC step, PID tuning from a first-order-plus-dead-time model, and ARX and eigensystem-realisation identification.
- Add `kernels::meshgen` (Bowyer–Watson Delaunay triangulation in 2D with Lawson flips, constrained triangulation by edge flipping, Ruppert-style quality refinement, an advancing front, and Delaunay tetrahedralisation with the insphere predicate) and `kernels::cad` (binary and ASCII STL, an ISO-10303-21 subset of STEP, welding, degenerate removal, consistent orientation, hole filling, manifold and closedness checks, and BSP-tree Boolean union, difference and intersection).
- Add `io::frame`: typed columns with nulls, a data frame with select, filter, sort, join and group-by aggregation, RFC 4180 CSV with type inference, and JSON in both array-of-objects and object-of-arrays forms. Add `io::parquet` (flat schema, PLAIN encoding, RLE/bit-packed definition levels, hand-written Thrift compact metadata, Snappy decompression), `io::hdf5` (superblocks 0–3, object headers 1 and 2, contiguous, compact and chunked datasets, deflate and shuffle, attributes, plus a writer for contiguous datasets) and `io::netcdf` (CDF-1 and CDF-2 with record dimensions, in and out). Every codec validates lengths and offsets before slicing and returns `Err` on a malformed file.
- Add `kernels::mp`: an arbitrary-precision tower. `Real` is `mant * 2^exp` with the mantissa of any width, every operation correctly rounded once from an exact `BigInt` result in five rounding modes; the elementary functions return balls whose radius bounds the series truncation and the rounding; `consts` computes `pi`, `e`, `ln 2`, `sqrt 2` and the golden ratio by integer series; `decimal` reads and writes decimal text exactly and renders a ball only when both endpoints agree on the digits. The `mp_*` operators carry their operands as decimal string literals, so precision travels with the value and no session state, global rounding mode or new Aevia value type is needed. `mp_to_float` is a single correctly rounded narrowing step.
- Expose the new families as table-driven operators (`qp_solve`, `socp_solve`, `sdp_solve`, `milp_solve`, `glm_fit`, `glm_predict`, `acf`, `pacf`, `difference`, `ar_fit`, `arima_fit`, `arima_forecast`, `holt_winters`, `holt_winters_mul`, `kalman_filter`, `care`, `dare`, `lqr`, `dlqr`, `lqg`, `kalman_gain`, `tf_state_space`, `ss_transfer`, `tf_poles`, `tf_bode`, `tf_nyquist`, `pid_tune`, `mpc_plan`, `arx_fit`, `kaplan_meier`, `nelson_aalen`, `log_rank`, `cox_fit`, `cox_predict`, `rhat`, `effective_sample_size`, `chain_autocorrelation`, `geweke`, `delaunay_2d`, `delaunay_3d`, `cdt_2d`, `refine_quality`, and the `mp_*` family). Each records why a call did not reduce.
- Extend the diagnostics channel beyond the table-driven families: `kernel::diag::classify` derives a reason for a heavy operator that recorded none, from the operator's arity, its typed overloads and its declared `DynamicReason`, so `Session::apply_diag` no longer answers `Unspecified` for a wrong argument count, a call no overload accepts, or a symbolic argument to a numeric-only operator.

### Performance: shared numerical infrastructure

- Add `kernels::par`: deterministic data-parallel helpers. Maps collect in index order and reductions fold fixed-size blocks in block order, so `RAYON_NUM_THREADS=1` and sixteen threads give the same bits; every helper has a sequential path below a grain size, and work items inherit the active cancellation scope.
- Add `kernels::reuse`: bounded scratch pools (per-thread slots over sharded shared storage) that lease empty buffers by minimum capacity and retain allocations subject to byte and count limits, plus plan caches for immutable built artefacts and a process-wide `global_pool()`. FFT plans are keyed by length and direction and concurrent requests for one key share a build.
- Add `execution::HostEnvironment`: the host's CPU count, SIMD width and L1/L2 tile sizes are probed once, and `should_split`, `should_split_work`, `choose_chunk_size` and `with_scoped_pool` decide when a kernel splits without oversubscribing an enclosing pool.
- Add batched modular inversion over a slice, slice-level modular add/sub/scale/axpy/dot, and an NTT plan cache with `poly_mul_ntt_mod`, differentially tested against schoolbook multiplication.
- Add u64 fast paths under the generic exact algorithms, each checked against its reference: modular row echelon for F4, a u64 S-polynomial merge for Gröbner bases, modular determinant and Dixon lifting, and a `gcd_u128` that falls back to the tested u64 gcd. `kernels::qlinalg` gains Bareiss elimination, modular determinants, characteristic polynomials, Hermite and Smith normal forms.
- Rewrite `apart` on a modular Z algorithm (rational reconstruction, f-adic expansion) with a pseudo-division path for symbolic parameters; add a parallel factor-block decomposition with a serial/parallel determinism test.
- Reuse scratch and matching state across a scheduler run: `MatchScratch` is shared by consecutive rewrites, the reachable set is rebuilt in place, tier sweeps reuse their bucket workspace, and incremental searches memoise the changed-node list per (trigger, epoch, depth).

### Consolidation: one modular-arithmetic implementation

- Add `kernels::modular`: `Montgomery` for an odd modulus, `Barrett` for any modulus of at least two, `mod_inv_u64`, `extended_gcd_u64` and the shared `add_mod_raw`/`sub_mod_raw`. `finite_field::{MontgomeryPrimeField, BarrettMod}`, `number_theory::MontU64`, `sparse_exact::{MontGfp, BarrettGfp}` and `qlinalg`'s private context are now aliases. The copies had drifted in their reduction contracts — one reduced its inputs before adding, one assumed `p < 2^63`, one handled the REDC carry and one did not — so a fix reached one copy and left the others.
- `sparse_exact::Gfp` carries a `Barrett` reducer instead of the bare modulus, so a sparse-elimination multiply is a multiply, a shift and at most two subtractions rather than a `u128` hardware division: 4.8–5.3x on the multiply (1.8 ns against 8.6–9.3 ns, release). `Barrett` implements `Domain` directly, and the Montgomery special cases in `Gfp::batch_inv` and `reduce_rational_matrix_mod_p` are gone.

### Consolidation: one implementation per primitive

- `is_prime` (six copies), `mod_inv_u64` (five), same-type `gcd` (four), `isqrt` (three) and Tonelli–Shanks (four) now delegate to the `number_theory` and `modular` implementations. Types that genuinely differ (`i64`, `BigInt`, `i128`) keep a wrapper but no longer their own algorithm.
- The merge caught a real bug: `number_exact::garner_step_with` used `mont.mul` where the canonical type spells `mul_mont`, so `crt([3, 5], [4, 9])` returned 35 instead of 23.

### Performance: foundational kernels

- Dense linear algebra: cache-blocked `matmul` that streams one panel of the right operand per output block, a parallel LU trailing update, column-wise Cholesky, a row-wise QR `R` update, and pooled reflector buffers. Bit-identical to the scalar oracle at every size and thread count; matmul 2.9–3.8x, QR 2.35–4.4x, LU and Cholesky up to 2.1–2.4x at the sizes tested.
- Sparse and interpolation: parallel `A^T x` per column, block-parallel Jacobi preconditioner setup and application, and parallel tensor-spline line reduction, each bit-identical across thread counts.
- Special functions and quadrature: bounded `(tag, order)` cache for Gauss–Legendre and Clenshaw–Curtis rules (105 us to 0.7 us at order 64, 18.5 ms to 1.9 us at 512), a cached Bernoulli table (489 us to 2.2 us), and element-wise parallel `*_batch` entry points (`erf` over 200k values: 981 ms to 77.6 ms).
- Polynomial arithmetic: Karatsuba threshold 33 to 40 and FFT threshold 129 to 48, both re-measured (mul 1.22x at 48, 1.84x at 64, 3.01x at 128); `eval_batch` advances four Horner chains at once (3x at eight points); Karatsuba leases its split buffers from `reuse::f64_scratch`.

### Accuracy: composition in double-double

- `Polynomial::compose` accumulates in double-double and rounds once. Composition is ill-conditioned when the intermediate terms `c_k Q(x)^k` dwarf the result, and the `f64` Horner and Brent-Kung strategies then rounded differently — up to 3.6e-16 relative apart on a non-contractive inner polynomial. Both now run over `Dd` and round to the same `f64` at every degree from 4 to 90.
- `mul_fft` documents its error (`O(log n)` epsilons against one or two for schoolbook). `mul_karatsuba` is the tighter-association option, `pow` is consistent with repeated `mul`, and a regression test bounds the FFT tier at 1e-11 relative against the schoolbook oracle.

### Numerical and symbolic algorithms

- Add the Risch–Bronstein integrator over a tower of monomial extensions (`rules::calculus::risch`). A tower is built from the exponentials and logarithms of the integrand (`theta_i' = a_i*theta_i + b_i` over the field below it, with `x` as the first level), and integration at each level runs Hermite reduction for repeated factors, the logarithmic part from the rational roots of `res_theta(D, A - z*D'*theta')`, and a bounded monomial ansatz for the polynomial part solved by one linear system over `Q`. The exponential place `theta = 0` is integrated as negative powers rather than as a residue, which is what `∫ dx/e^x = -e^{-x}` needs. Answers are checked exactly (`p' = h`) and then numerically like every other antiderivative; an integral that needs an algebraic residue field, or a tower deeper than three levels, is left to the other methods.
- Add stiff one-dimensional method-of-lines integration with time-dependent Dirichlet and Neumann boundaries, including measured callback programs for embedded hosts.
- Correct finite-temperature XXZ stochastic series expansion sampling. The corrected update is reproducible for a fixed seed and changes the trajectories produced by earlier builds.
- Add one-sided Welch spectral density, inverse STFT with overlap normalization, and polyphase FIR resampling. Add checked topological ordering and a nonrecursive strongly connected component API shared with symbolic graph operators.
- Add analytic signals, Hilbert transforms, zoom FFT, spectral coherence and wavelet packets. DCT and DST use FFT-based algorithms and fractional Fourier transforms can reuse their spectral basis.
- Add overlap-save convolution, autocorrelation, FIR filter design, frequency response, second-order-section filtering, polyphase resampling and stationary wavelet transforms with coefficient thresholding.
- Add batched quantiles and empirical CDF, Benjamini–Hochberg, Benjamini–Yekutieli and Holm adjustment, batched multivariate normal log densities, and atomic covariance batch updates with correlation output.
- Add signed shortest paths, weak components, maximum flow with its minimum cut, and column interpolative decomposition. The kernel manifest carries the numerical inputs, structured results and dimensional contracts.
- Add multi-source BFS, strong-component condensation, articulation points and bridges, and three sample time-series lag estimators. Large Welch workloads reduce fixed frame blocks in parallel. Dormand–Prince, conjugate gradients and restarted GMRES reuse inner-loop workspaces.
- Optimized dense multiplication uses row-major input and output views. Large products use the active Rayon pool.
- Euclidean norms, distance and vector directions scale intermediate squares at extreme ranges. Distance no longer materializes the difference vector, and Rust callers can normalize into their own buffers.
- Vector reflection writes one output buffer and preserves finite components when an intermediate correction overflows. Arnoldi scales its basis construction and reuses projection storage. Native dense eigensolvers write matrix-vector products into the iteration workspace.
- Bilinear interpolation uses the same scaled interval coordinates for single queries and prepared grids, including finite endpoints whose difference exceeds Float range.
- Prepared Newton–CG, polynomial antidifferences and geometric grid transfers retain reusable numerical state. Implicit Gauss–Legendre steps use a separate Jacobian at each collocation stage.
- The library adds logarithmic Bessel K for real orders and exact Chebyshev residues with full-width unsigned moduli. Sparse exact elimination reuses row storage and indexes pivot candidates on larger systems.
- Block sparse multiplication reuses row workspaces and accepts batches of dense right-hand sides. Repeated moderate Zoom FFT calls share the transformed chirp for the same grid.
- Drift-implicit SDE methods reuse Newton workspaces and provide a sampled path. PCHIP and Akima batch queries reuse segment coefficients. Online covariance accepts exact integer frequencies, and integer kernels provide a Möbius sieve and linear congruence solver.
- Sparse conjugate gradients rescale extreme right-hand sides. Fixed-step delay solvers reject time grids that cannot advance, and the wave stencil combines speed and step before squaring.

- Added weighted empirical quantiles and CDFs, discrete logarithms for composite moduli with non-invertible bases, and checked exact integer polynomial products. Improved LU factor reuse, eigensolver scaling, least-squares projections and parallel fractal correlation counts.
- The Rust LOBPCG API solves matrix-free generalized eigenproblems with a positive-definite mass operator and soft-locks converged vectors.
- Fixed-step ODE solvers reject invalid time grids and unrepresentable trajectory sizes before calling user functions. Event location reuses its state buffers.
- The single-pole Meijer G reduction computes the power and exponential on a common scale and preserves subnormal results. The commutative Platen solver reuses its correction workspace and rescales Brownian corrections that would overflow or underflow in intermediate products.
- Complex winding uses signed ray crossings with an exact comparison when product rounding obscures the sign. The numerical registry accepts measured closed polygons, and the argument-principle sampler evaluates each contour vertex once.
- Prepared ray–triangle intersections support measured coordinate references, finite or unbounded ray intervals, and batches with a shared direction. Near-edge tests use an exact sign check. The quadratic sieve validates its work limits and can return an ordered exact factor pair.
- Prepared Reed–Solomon encoders reuse feedback tables. Encoding, decoding and error counting share the 255-symbol codeword bound, and Rust callers can encode into an existing byte buffer.
- Exact factorial, binomial and ordered-selection counts use balanced integer products for large selections. Homotopy continuation checks degree sums and the Bézout path count before constructing start paths.
- Prepared inverse real FFTs reuse unpacking factors and report separate sample and bin counts. The measured interface preserves the spectrum dimension, including for subnormal samples.
- Linear regression borrows sample columns and preserves finite fits across extreme moment scales. Checked polynomial long division can recover a finite remainder after a product overflows. Adaptive Simpson quadrature reuses endpoint and midpoint samples across subdivisions.
- Matrix access and deserialization check both extents. Matrix transposition can reuse a correctly shaped destination. Nelder–Mead sorts vertex indices and reuses the storage for trial points.
- Persistent homology keeps the earliest filtration value for duplicate simplices, treats signed zeros as the same filtration level, and reuses the boundary reduction workspace.
- Double-double square root retains precision at extreme scales. The numerical registry provides double-double Euclidean norm and complex modulus, with a separate readout of the two components in their original physical units. Complex double-double division and reciprocal scale the intermediate products and preserve representable components.
- CSV and JSON array output use buffered writes. The CSV reader reuses its line storage and parses directly into the result array. Raster writers stream row-major pixels, and 16-bit PGM conversion uses a bounded reusable buffer.
- Binomial masses and beta densities handle support endpoints. A bounded binomial recurrence preserves low-count probabilities when the success probability is tiny. Binomial sampling uses discrete inversion with one uniform draw for nondegenerate parameters. Uniform density, CDF and sampling also support finite endpoints whose difference exceeds the floating-point range. Seeded binomial sample sequences change with the new algorithm.
- Improved extreme-tail probabilities for statistical tests and used thin SVD for ridge and OLS regression. Kendall correlation handles ties in O(n log n). Bootstrap procedures check sample and output sizes before allocating.
- Extended symbolic Mellin integrals, limits and infinite sums. Large-argument complex Faddeeva, error and Fresnel functions use scaled computations to avoid intermediate overflow.
- Poisson, heat and wave grid solvers accept per-face Dirichlet and Neumann conditions. Pure-Neumann Poisson checks the source–flux balance and returns a zero-mean solution. Structured kernel diagnostics identify rejected shapes, domains and convergence criteria.
- Three-dimensional shallow-water Euler and AMR runs accept a barotropic state in `(h, hu, hv, hw)`. The kernel manifest describes boundary conditions and tensor ranks for typed hosts.

- Fixed-step ODEs can return only the final state. Diagonal ETDRK4 shares a real-weight stepper across fixed and adaptive routes. Kerr–Newman trajectories can sample by stride and initialize either supported future-directed energy branch within the ergoregion.
- Finite-field homology computes all requested degrees in one pass over the boundary ranks. Binary elimination packs row bits, and odd-prime arithmetic covers the full unsigned word range. Permutation rotations accept signed arbitrary-precision exponents.
- Seeded Ising runs support online observables, and Potts sweeps reuse neighborhood weight tables. Complex Lambert W distinguishes both sides of its negative real cut, and inverse erfc retains accuracy for subnormal probabilities.
- The transverse magnetic FDTD solver reuses its field buffers and can return only the requested final electric field. The physical operator returns the field and elapsed time in their respective dimensions.
- The two-dimensional Schrödinger solver prepares splitting phases for repeated initial waves. Its measured terminal operator supports signed time increments and returns complex amplitudes with inverse-length dimensions.
- Cartesian finite-volume Euler and shallow-water solvers reuse their Runge–Kutta workspace and support the two-wave HLLE flux. Prepared triangular P1/P2 meshes reuse geometry and quadrature for independent diffusion–reaction problems. Extreme finite loads use scaled conjugate gradients.
- Prepared Cartesian Q1 Poisson solves reuse the sine-basis stiffness spectrum on uniform two- and three-dimensional grids. Prepared two-dimensional fast multipole evaluations reuse source geometry across strength vectors. Chain dynamics exposes measured accelerations in the original body order.
- AMR computes independent flux lines in parallel on large patches and validates refinement extents and final times. Prepared ideal-gas Riemann solutions reuse the wave pattern across queries, including vacuum and extreme common pressure and density scales.
- Checked explicit Euler, midpoint and Heun end at the requested time in either direction. Periodic one-dimensional DG advection can reuse a fixed mesh and velocity across time steps.
- Sparse MD neighbour searches scale with occupied cells. The one-dimensional method of lines avoids copying the full field at each right-hand-side evaluation. Truncated holonomic series can compute the first coefficients of large powers directly.
- Channel flow reuses the pressure system between steps. Closed fluid regions require compatible net flux and receive a pressure gauge. Rigid molecular dynamics precomputes its rotor composition before updating the state, and the cell-flow exponential uses closed forms for diagonal matrices and simple shear.
- Three-dimensional Laplace boundary-element solves with fixed boundary kinds can reuse the dense collocation operator and its factorization across boundary values.
- Boundary-element queries reuse already validated mesh and coefficient data and reject mixing interior and exterior solutions. Cut-cell callbacks retain the signs of subnormal distances. One-dimensional advection leaves its retained mesh intact if the measured output cannot be represented.
- The parallel CPU device splits large one-row matrix products over output columns while preserving the reduction order. Real algebraic field embeddings generate successive powers in the quotient basis and eliminate independent rows in large exact solves.

### Nested repetition, more encodings, compound data and nested groups

- `io::parquet` nests lists to any depth in both directions. The reader rebuilds `list<list<...>>` from the repetition and definition levels one level at a time, and the writer shreds a nested `Column::List` into the three-level schema the format prescribes. A list of lists is no longer rejected.
- `io::parquet` decodes `DELTA_BYTE_ARRAY` pages (prefix lengths plus suffixes), next to the delta, dictionary, byte-stream-split and RLE encodings it already handled.
- `io::parquet` writes version-2 data pages on request (`WriteOptions::data_page_version`), which keeps the levels outside the compressed payload, for flat and nested columns alike.
- `io::hdf5` reads compound datasets as one dataset per member (`/pair/a`, `/pair/b`): integer, floating-point and enumeration members become values, fixed-width strings become text, and members it cannot represent are skipped.
- `io::hdf5` reads dense attribute storage: the attribute info message points at a fractal heap and a version-2 B-tree, and each attribute message is decoded from the heap. Attributes written either way come back the same.
- `io::hdf5` writes nested groups: a dataset path such as `/a/b/temperature` creates every intermediate group, each with its own object header, local heap, group B-tree and symbol table nodes.
- Attribute messages holding variable-length strings are resolved through the global heap instead of failing the file.
- `tests/fixtures/` holds the binary fixtures the `io` tests read, together with the `generate.py` that reproduces them from `h5py`, `pyarrow` and `netCDF4`. No test is `#[ignore]`d for lack of a fixture.

### Repeated Parquet fields and a fuller HDF5 reader

- `io::parquet` reads and writes repeated fields. A `LIST`, a `MAP` and a bare repeated primitive all become `Column::List`, a flattened element column plus row offsets, so a null list, an empty list and a null element stay distinguishable. `MAP` is exposed as its `key` and `value` lists, one column each. The writer emits the three-level layout with correct definition and repetition levels for every combination of optional and required list and element.
- `io::frame` gained `Column::List`: `slice`, `new_list`, `list_bounds`, `list_row`, `list_values`. JSON reads and writes nested arrays as list columns, CSV writes a list cell as a JSON array, and `take`, `concat`, `sort`, `group` and the row-key encoding handle lists. `DType::List` is the new logical type.
- `io::hdf5` reads dense groups: fractal heap (`FRHP`), direct and indirect heap blocks including nested ones, the doubling table, and version-2 B-trees at any depth with the record-count and subtree-total field widths the format prescribes. A group of 5000 links opens.
- `io::hdf5` reads enumerations, variable-length strings through the global heap, and fixed-width strings, and writes fixed-width string datasets. `Dataset` carries `strings` next to `data`, with `Dataset::new_numeric` and `Dataset::new_strings` constructors. Attribute messages are written in the version-3 layout the format specifies (nine-byte header, unpadded name and terminator) and read back from files written by the reference library.
- `io::hdf5` skips a dataset whose datatype it cannot represent (compound, array, reference, chunked strings) instead of failing the whole file.

### Wider optimization and distribution surface

- `kernels::optimize` adds `find_root` over six bracketed scalar solvers (bisection, Brent, secant, TOMS 748, Newton, Halley, with central differences standing in for a missing derivative) and two more [`SolverKind`] entries: bounded L-BFGS (`Lbfgsb`) and bounded CMA-ES.
- `kernels::distributions_ext::StatrsDistribution` exposes the whole `statrs` univariate family — 29 distributions from Bernoulli to Weibull — behind one enum with `pdf`, `cdf`, `sf`, `mean`, `variance`, `entropy` and a `from_name` factory. `kernels::distributions::Distribution` keeps the hand-written kernels the rule layer evaluates.
- `kernels::special` adds `harmonic`, `inverse_beta_reg`, `multinomial` and `ln_binomial`.
- `signatures.overlay.json` pins an overload for every operator that does not declare one: the audit reports 5411 of 5411 typed, with no operator left on a derived-only signature.
- `io::netcdf` reads and writes `CDF` (64-bit data) alongside `CDF` and `CDF`. `numrecs` is a 32-bit count up to CDF-2 and 64-bit in CDF-5; dimension counts, name lengths, dimension identifiers, `vsize` and `begin` widen with the format. Files from all three formats open in the reference `netCDF4` library, and files it writes are read back here.

### Optimizer migration and dependency refresh

- `kernels::optimize` now runs on `basin` instead of `argmin`. `minimize` takes an `FnProblem` (cost closure, optional gradient closure, optional box bounds) and a `SolverKind`, and returns a plain `MinimizeResult`; `SolverKind` names 24 solvers across first-order (gradient descent, nonlinear CG, projected gradient, TNC), quasi-Newton (BFGS, L-BFGS), derivative-free (Nelder–Mead, BOBYQA, NEWUOA, LINCOA, COBYLA, DIRECT, MADS, SLSQP), global/stochastic (CMA-ES, differential evolution, global-best PSO, random search, Solis–Wets, steady-state GA, GBNM) and scalar (Brent, golden section) families. Box bounds are part of the problem, so a bounded solver reads them from `MinimizeOptions::bounds` and an unbounded one ignores them. `argmin` and its `argmin-math` shim are gone.
- Dependencies moved to their current releases: `ndarray` 0.17, `ndarray-npy` 0.11, `rand` 0.10 (one version instead of two), `gif` 0.14, `roxmltree` 0.21, `pollster` 1, `syn` 3 (the codegen tools), Cranelift 0.138 and `wgpu` 30 (the GPU backend was ported to the new instance descriptor, pipeline-layout and error-scope APIs). `num-bigint` stays on 0.4 because `num-rational` 0.4 — the newest release — pins it and `BigRational` is `Ratio<BigInt>`.
- `basin` is pulled in with default features off: the wrapper uses its plain `Vec<f64>` backend, so the optimizer adds no linear-algebra backend of its own.

### Hardening and interoperability

- `kernels::convex::socp` now stops on a scale-invariant merit (relative duality gap and residuals normalised by the data norms), keeps the best iterate instead of the last one, and ends a stalled run early instead of spending the remaining budget. The normal block is assembled into reused scratch, so an iteration no longer allocates once per variable. Large-data instances that could never reach an absolute `1e-10` gap now converge.
- `io::parquet` writes dictionary pages with `RLE_DICTIONARY` data pages when a column chunk repeats values (`WriteOptions::use_dictionary`), and gains a gzip codec in both directions, backed by a fixed-Huffman DEFLATE encoder and the existing inflate decoder. The reader also accepts Data Page v2 headers, including the split level/value compression and the RLE boolean encoding that version-2 pages use. Round-trip files were read back with Apache Arrow, and Arrow-written files (v1 and v2 pages, Snappy, gzip, dictionary and plain encodings) were read back here.
- `io::hdf5` writes object-header messages whose declared size includes the eight-byte padding the HDF5 library expects, and puts the sign-bit position in the float datatype message. Both the superblock-v0 and superblock-v2 writers now produce files that the HDF5 library itself opens, including numeric, fixed-length string and attribute data.
- `kernels::meshgen::delaunay` caches each triangle's circumcircle and decides membership with a floating-point filter, falling back to the exact predicate only within its error bound. The result is unchanged, and the per-insertion scan no longer pays for an exact expansion at every triangle.
- `io::parquet` reads and writes Zstandard and `LZ4_RAW` pages as well as Snappy and gzip, and exposes `DECIMAL` columns as `f64` divided by `10^scale` whatever their physical storage (`INT32`, `INT64`, `BYTE_ARRAY` or `FIXED_LEN_BYTE_ARRAY`). The LZ4 block codec is implemented here; Zstandard uses the reference `zstd` crate. Arrow-written files in every supported codec and both page versions are read back here, and Arrow reads the files written here.
- `io::hdf5` reads dense new-style groups (a fractal heap indexed by a version-2 B-tree), enumerations (as their underlying integers) and variable-length strings (resolved through the global heap). Files written by h5py 3.16 with `libver="latest"` are read back, including a group with a thousand links.
- `rules::calculus::risch` builds towers of up to five monomial extensions instead of three.

### Build and compatibility

- Source builds require Rust 1.100 or newer. The default Rust artifact is an rlib. Build C libraries explicitly with their generated headers. C archives include English and Chinese instructions for linking directly from the extracted directory. [Build guide](docs/BUILDING.md)
- C batch computation accepts uninitialized output buffers. The `ffi` build profile supports converting unwinding Rust panics to C error codes at guarded entry points.
- The manifest API version stays at 2. Inspect the manifest hash before sharing serialized contracts across independently built kernels. [Kernel embedding](docs/KERNEL.md)

### 内核与嵌入

- 原生对象绑定从同一参数列表生成签名与运行时参数个数。稠密值检查无需分配遍历缓冲区。
- 科学计算内核提供确定性的清单、类型与量纲契约、结构化结果以及类型化数值回调。宿主从实际执行调用的同一算子注册表生成受检查的接口。
- 重写索引和增量匹配减少重复搜索。独立请求使用各自的计算图，结果按输入顺序导入。可选 JIT 和设备后端明确支持的精度与回退行为。
- 工作线程在独立计算图中准备和执行符号批处理，并限制同时保留的任务数。CPU 调度读取当前 Rayon 池的大小，FFT 通过 `kernels::reuse` 共享计划并复用临时空间。
- 符号化简的树窗口复用相交子图的抽取结果。编译后的 CPU 批处理按容量借用临时缓冲区，小任务顺序执行。同尺寸、同方向的并发 FFT 请求共享计划构建。
- 列式计算带逐行读取被选中的来源列。预备的选列计算带可接收共用参考量的带量纲来源列，并保留结果的原始量纲。
- 标准会话共享已安装规则和清单，直至全局单位表变动。会话的自定义单位、项图、原生对象及计算缓存各自独立。
- `signatures.overlay.json` 现在只包含未声明签名的算子。重复声明签名的那批条目已移除：它们对清单不可见（声明优先），却破坏了 `tests/signature_overlay.rs` 的双向对应，该测试现在对未注册算子的条目和清单从不读取的条目都会失败。
- 清单可跨泛型参数与再导出路径（`FdMode<F>::parallel`、`crate::sim::md_engine::FlexibleBarostat::energy`）找到算子背后的函数，因此另有 551 个算子在所属族的契约句之外发布其实现所记载的描述。`doc_path` 记录提供文本的函数；有歧义的路径仍然不解析。
- `category` 不再为空：注册未指定类别的算子归入其规则集。`tests/manifest_completeness.rs` 为每个算子固化宿主依赖的不变量：名称、规则集、类别、说明、代价、键、重载来源，以及有类型或已解释。
- 已发布的见证 `is_lie_algebra(sl_basis(4))` 不再可重放：该调用返回 `true`，但耗时 20.5 秒，超出挖掘器 20 秒预算。它已从 `signatures.witnessed.json` 中移除，而不是作为无法复现的证据保留；`signature-mine --check` 现在可重放 727 个算子的 1047 条见证。其所属重载保留另一条见证。
- 四个算子现在声明比符号契约更精确的结果——`arg(term) -> float`、`boltzmann_entropy(term) -> real`、`pde_order(term, term) -> int`、`todd_coxeter_index(term, term, term) -> int`——以 *refinements* 而非重载的形式：它们的实参就是契约自身的实参，再加一个重载会让每次调用都产生歧义。`infer` 按设计忽略 refinements（`infer.rs`）；带结果检查通道的宿主会读取它们，这是传句柄的调用方得知 `arg` 返回浮点数而非项的唯一边界。每条都带有产生它的那次调用，`signature-mine --check` 可重放全部 1051 条见证、731 个算子。
- 132 个算子的处理函数在每条成功路径上都构造同一种 `Value`，现在直接发布该结果类型，而不只是符号契约：`graph_node_id(Obj<graph>, Int) -> Int`、`lie_rank(Symbol) -> Int`、`perm_group_is_abelian([[Int]]) -> Bool`、`gfp_divmod([Int], [Int], Int) -> [[Int]]`、`point_group_order(Symbol) -> Int`、`bdd_size([[Int]]) -> Int`、`crt([Int], [Int]) -> Int`、`pauli_x() -> [[Int]]`、`sc_homology_mod(Obj<simplicial_complex>, Int, Int) -> Int`。每个具体重载都基于两个独立来源：构造该值的代码路径，以及一次引擎以该类型作答的调用。若某具体重载的实参无法与符号契约区分则不予采用——那种更精确的结果属于 `refinements`，在那里不会让调用产生歧义；若算子的实际调用与代码判定不符（一条路径返回 `false`、另一条返回列表），则保留符号契约。符号重载 `(Term, …) -> Term` 始终排在最后，因此算子不接受的调用仍以项定型。132 个之中有 4 个（`gf_inv`、`gfp_isirreducible`、`perm_order`、`perm_group_order`）正是 Aevia 桥接层此前手工整理出的签名。
- C ABI 捕获共享的标准清单快照并复用 JSON 编码，每次调用仍返回独立所有权的字符串。
- 二进制清单、插件消息和计算记录采用 bincode-next 3.2 原生编码。签名和清单保持确定性编码，校验清单时直接读取现有的算子表和单位表，逐段将规范编码传给哈希器。兼容性测试将支持的旧线格式与原生编码进行比较。JSON 与 CBOR 继续提供自描述元数据。

### 数值与符号算法

- 横磁 FDTD 求解器复用场缓冲区，并可直接返回指定步数的最终电场。物理算子分别保留电场和经过时间的量纲。
- 二维薛定谔求解器为多个初始波函数复用分步相位。带量纲的末态算子支持带符号时间步，并返回逆长度量纲的复振幅。
- 笛卡尔有限体积欧拉和浅水求解器复用 Runge–Kutta 工作区，并支持双波 HLLE 通量。预备的 P1/P2 三角网格为独立扩散反应问题复用几何与求积数据。极端有限载荷使用缩放后的共轭梯度。
- 预备笛卡尔 Q1 Poisson 求解在二维和三维均匀网格上复用正弦基刚度谱。二维快速多极子法在多组源强度间复用源点几何。链式动力学按输入物体顺序返回带量纲的加速度。
- AMR 在较大的网格块上并行计算独立通量线，并校验细化范围和终止时间。预备理想气体 Riemann 求解在多次查询中复用波系，支持真空区及密度和压力的极端同比缩放。
- 受检的显式 Euler、中点法及 Heun 方法可沿任意时间方向到达请求终点。周期一维 DG 平流可在多个时间步之间复用固定网格与速度。
- 稀疏分子动力学近邻搜索的网格存储随占用单元数增加。一维线法在求值右端项时复用积分器状态。`holonomic_series` 可直接求高次幂的前若干系数。
- 通道流在各时间步间复用压力系统。封闭流体区域须满足净通量相容性，并固定压力基准。刚体分子动力学在更新状态前预计算转子组合，晶胞流的矩阵指数对对角矩阵和简单剪切使用闭式解。
- 边界类型固定时，三维 Laplace 边界元求解器可对不同边界值复用稠密配点算子与分解。
- 边界元求值复用已验证的网格与系数，并拒绝内域解和外域解混用。切割单元回调保留次正规距离的符号。一维平流无法恢复带量纲的输出时保留原网格状态。
- 并行 CPU 设备可按输出列划分大型单行矩阵乘积，并保持每个元素的累加顺序。实代数域嵌入在商代数的基上逐次生成原始元的幂，在大型精确求解中并行消去独立行。
- 新增多源 BFS、强连通分量缩点图、割点与桥，以及三种时间序列延迟统计量。大型 Welch 任务并行计算固定帧块并依序合并。Dormand–Prince、共轭梯度和重启 GMRES 复用迭代工作区。
- 优化后的稠密矩阵乘法直接使用行主序输入和输出。大型乘积使用当前 Rayon 线程池。
- 欧几里得范数、距离和向量方向计算会在极端范围缩放中间平方项。距离计算省去差向量的分配，Rust 调用方可将单位向量写入自备缓冲区。
- 向量反射直接写入输出缓冲区，并在中间校正项溢出时保留可表示的分量。Arnoldi 缩放基向量构造并复用投影存储。原生稠密特征求解器将矩阵向量乘积直接写入迭代工作区。
- 双线性插值在单点查询与预备网格中共用缩放后的区间坐标，包括端点有限但间距超过 Float 范围的情况。
- 预备的 Newton–CG、多项式原差分和几何网格传递保留可复用的计算状态。隐式 Gauss–Legendre 在各配点分别计算雅可比。
- 支持实阶 Bessel K 的对数求值，以及完整无符号模数范围内的精确 Chebyshev 余数。精确稀疏消元复用行存储，并为较大系统建立主元索引。

- 增加单边 Welch 功率谱密度、重叠归一化逆 STFT 和多相 FIR 重采样。增加受检查的拓扑排序，并让 Rust 和符号图算子共用非递归强连通分量实现。
- 增加解析信号、Hilbert 变换、zoom FFT、谱相干度和小波包。DCT、DST 采用 FFT 算法，分数阶 Fourier 变换支持复用谱基。
- 增加重叠保留卷积、自相关、FIR 滤波器设计、频率响应、二阶节滤波、多相重采样，以及支持系数阈值处理的平稳小波变换。
- 增加批量分位数和经验分布函数、Benjamini–Hochberg、Benjamini–Yekutieli 与 Holm 多重检验校正、多元正态批量对数密度，以及可输出相关矩阵的原子协方差批量更新。
- 增加含负边权的最短路径、弱连通分量、最大流及最小割、列插值分解。内核清单声明输入、结构化输出和量纲契约。
- 稀疏共轭梯度对极端尺度的右端项做缩放。固定步长时滞求解拒绝无法前进的时间网格，波动方程离散格式先合并波速和步长再计算平方。

- 增加加权经验分位数和累积分布、含不可逆底数的复合模数离散对数，以及带溢出检查的精确整数多项式乘法。改进 LU 分解复用、特征值计算的缩放、最小二乘投影和并行分形关联计数。
- Rust LOBPCG API 支持带正定质量算子的无矩阵广义特征值问题，并软锁定已收敛的向量。
- 固定步长 ODE 求解器在调用用户函数前检查时间网格和轨迹容量。事件定位复用状态缓冲区。
- Meijer G 的单极点公式在共同尺度上计算幂和指数，保留次正规结果。可交换 Platen 求解器复用校正工作区，并缩放会在中间乘积中上溢或下溢的布朗校正项。
- 复数绕数使用带符号的射线穿越计数，浮点乘积舍入无法区分符号时使用精确比较。数值注册表支持带量纲的闭合多边形，辩值原理的采样在每个轮廓顶点只求值一次。
- 预备的射线与三角形求交支持带量纲的坐标参考值、有限或无上界的射线区间，以及方向相同的批量射线。贴近边缘时检查精确符号。二次筛校验工作预算，并可返回有序的精确因子对。
- 预备的 Reed–Solomon 编码器复用反馈表。编码、解码和错误计数共用 255 个符号的码字上限，Rust 调用方可将码字写入自备字节缓冲区。
- 精确阶乘、二项式系数和有序选法数在大规模选择时使用平衡整数乘积。同伦延拓在构造起始路径之前检查次数和与 Bézout 路径数。
- 预备的逆实数 FFT 复用解包系数，分别报告样本数和频率分量数。带量纲的接口保留频谱量纲，支持次正规样本。
- 线性回归借用样本列，并在中间矩尺度极端时保留可表示的拟合结果。受检查的多项式长除法可在中间乘积溢出后恢复有限余式。自适应 Simpson 求积在各次细分之间复用端点和中点样本。
- 矩阵访问和反序列化检查两个轴的尺寸。矩阵转置可复用形状匹配的目标矩阵。Nelder–Mead 对顶点索引排序，并复用试探点的缓冲区。
- 持久同调对重复单形保留最早的过滤值，将带符号的零视为同一过滤层，并复用边界约化工作区。
- 双双精度平方根在极端尺度下保留精度。数值注册表提供双双精度的欧几里得模长和复数模长，支持按原量纲分别读取高、低分量。复双双精度除法与倒数缩放中间乘积，保留可表示的分量。
- CSV 和 JSON 数组采用缓冲写入。CSV 读取复用行缓冲区，并直接将解析值写入结果数组。光栅文件按行主序流式写入像素，16 位 PGM 转换复用有容量上限的缓冲区。
- 二项分布概率质量与 Beta 密度正确处理支撑集端点。成功概率极小时，有界递推可保持小计数的二项概率。非退化二项分布使用一次均匀随机数反演采样。即使有限端点之差超出浮点范围，均匀分布的密度、累积分布与采样仍可求值。新算法会改变给定种子的二项采样序列。
- 改进统计检验的极端尾部概率，岭回归与普通最小二乘采用薄 SVD。Kendall 相关系数以 O(n log n) 处理并列值。自助法在分配内存前检查样本与输出大小。
- 扩展符号 Mellin 积分、极限和无穷级数。复数 Faddeeva 函数、误差函数和菲涅耳积分在大参数下采用缩放计算，避免中间溢出。
- 泊松、热传导和波动网格求解器支持按面指定 Dirichlet 或 Neumann 条件。全 Neumann 泊松问题检查源项与通量平衡并返回零均值解。结构化内核诊断指出不符合要求的形状、定义域或收敛准则。
- 三维浅水 Euler 和自适应网格细化求解器接收 `(h, hu, hv, hw)` 形式的正压状态。内核清单声明宿主所需的边界条件和张量秩约束。

- 定步长 ODE 可以只返回末态。对角 ETDRK4 的定步与自适应路径共用实权重步进器。Kerr–Newman 轨道可按步距采样，并在能层内选择相容的未来指向能量分支。
- 有限域同调通过一次边界秩计算得到各次数结果。二元消元打包行位，奇素数运算覆盖完整无符号机器字范围。置换循环支持有符号任意精度指数。
- 带种子的 Ising 模拟支持在线观测，Potts 扫描复用邻域权重表。复数 Lambert W 区分负实割线两侧，反余误差函数保留次正规概率的尾部精度。

### 统一与性能

- 新增 `kernels::modular`，作为无除法模运算的唯一实现：`Montgomery`（奇模数）与 `Barrett`（任意不小于 2 的模数），并统一 `mod_inv_u64`、`extended_gcd_u64` 与共享的加减法。`finite_field::MontgomeryPrimeField`、`number_theory::MontU64`、`sparse_exact::MontGfp`、`finite_field::BarrettMod`、`sparse_exact::BarrettGfp` 以及 `qlinalg` 的私有上下文均改为别名；此前的副本在约减契约上已经分叉（有的先约减再相加，有的依赖 `p < 2^63`，有的处理 REDC 进位、有的不处理）。
- `sparse_exact::Gfp` 改为携带 `Barrett` 约减器，稀疏消元中的一次乘法由 128 位硬件除法变为一次乘法加移位：模乘本身提速 4.8–5.3 倍（1.8 ns 对 8.6–9.3 ns）。`Barrett` 直接实现 `Domain`，`batch_inv` 与 `reduce_rational_matrix_mod_p` 中的 Montgomery 特判随之删除。
- `is_prime`（六份）、`mod_inv_u64`（五份）、同类型 `gcd`（四份）、`isqrt`（三份）与 Tonelli–Shanks（四份）合并为单一实现。合并过程中发现并修复了 `number_exact` 的一处真实错误：`crt([3, 5], [4, 9])` 曾返回 35 而非 23。
- `Polynomial::compose` 改用双倍精度累加并只在最后舍入一次；两条求值路径（Horner 与 Brent–Kung）现在舍入到同一个 `f64`，原先相对偏差可达 3.6e-16。`mul_fft` 的误差量级（`O(log n)` 个 eps）写入文档，并有回归测试将其限制在 1e-11。
- 基础内核提速：分块矩阵乘法 2.9–3.8 倍、QR 2.35–4.4 倍、LU 与 Cholesky 最高 2.1–2.4 倍；求积规则按 `(tag, 阶数)` 缓存（64 阶 105 us 降至 0.7 us，512 阶 18.5 ms 降至 1.9 us）；多项式乘法在 128 次时 3.01 倍。

### 任意层嵌套、更多编码、复合类型与嵌套组

- `io::parquet` 的列表读写都支持任意层嵌套。读取端按层从重复层级与定义层级逐层还原 `list<list<...>>`，写出端把嵌套的 `Column::List` 拆成规范要求的三级结构；不再拒绝"列表的列表"。
- `io::parquet` 新增 `DELTA_BYTE_ARRAY` 页解码（前缀长度 + 后缀），与已有的 delta、字典、byte-stream-split、RLE 编码并列。
- `io::parquet` 可按需写出第二版数据页（`WriteOptions::data_page_version`），层级不再参与压缩，扁平列与嵌套列都支持。
- `io::hdf5` 把复合数据集按成员拆成多个数据集（`/pair/a`、`/pair/b`）：整数、浮点与枚举成员取值，定长字符串成员取文本，无法表示的成员跳过。
- `io::hdf5` 支持稠密属性存储：属性信息消息指向 fractal heap 与第二版 B 树，逐个解出属性消息。两种存储方式写出的属性读回结果一致。
- `io::hdf5` 支持写出嵌套组：`/a/b/temperature` 这样的路径会逐级创建中间组，每个组都有自己的对象头、局部堆、组 B 树与符号表节点。
- 属性消息中的变长字符串经全局堆解析，不再让整个文件失败。
- `tests/fixtures/` 存放 `io` 测试读取的二进制样例，并附带 `generate.py`，可用 `h5py`、`pyarrow` 与 `netCDF4` 重建。没有测试因为缺少样例而被 `#[ignore]`。

### Parquet 重复字段与更完整的 HDF5 读取

- `io::parquet` 支持读写重复字段。`LIST`、`MAP` 与裸重复原始类型都映射为 `Column::List`（扁平元素列加行偏移），因此空值列表、空列表与空值元素三者可区分。`MAP` 以 `key`、`value` 两列列表暴露。写出端按三级结构生成定义层级与重复层级，覆盖列表与元素的可选/必选四种组合。
- `io::frame` 新增 `Column::List`：`slice`、`new_list`、`list_bounds`、`list_row`、`list_values`。JSON 读写嵌套数组为列表列，CSV 把列表单元格写成 JSON 数组，`take`、`concat`、排序、分组与行键编码都处理列表。新增逻辑类型 `DType::List`。
- `io::hdf5` 支持稠密组：fractal heap（`FRHP`）、直接与间接堆块（含嵌套）、倍增表，以及任意深度的第二版 B 树（记录数字段与子树总数字段宽度按规范推导）。5000 个链接的组可正常打开。
- `io::hdf5` 支持枚举、经全局堆的变长字符串与定长字符串，并可写出定长字符串数据集。`Dataset` 在 `data` 之外新增 `strings`，并提供 `Dataset::new_numeric` 与 `Dataset::new_strings` 构造器。属性消息按规范的第三版布局写出（九字节头、名称与终止符不加填充），也能读回参考实现写出的文件。
- `io::hdf5` 对无法表示的数据类型（compound、array、reference、分块字符串）跳过该数据集，而不是让整个文件失败。

### 更宽的优化与分布接口

- `kernels::optimize` 新增 `find_root`，覆盖六种带区间的标量求根器（二分、Brent、割线、TOMS 748、Newton、Halley；缺导数时用中心差分代替），并新增两个 [`SolverKind`]：有界 L-BFGS（`Lbfgsb`）与有界 CMA-ES。
- `kernels::distributions_ext::StatrsDistribution` 用一个枚举暴露 `statrs` 的全部一元分布族（从 Bernoulli 到 Weibull 共 29 个），提供 `pdf`、`cdf`、`sf`、`mean`、`variance`、`entropy` 与 `from_name` 工厂。`kernels::distributions::Distribution` 保留规则层使用的手写内核。
- `kernels::special` 新增 `harmonic`、`inverse_beta_reg`、`multinomial` 与 `ln_binomial`。
- `signatures.overlay.json` 为**每个未自行声明签名**的算子固定了重载：审计报告 5411/5411 已定型，不再有算子停留在派生签名上。
- `io::netcdf` 在 `CDF`、`CDF` 之外支持读写 `CDF`（64 位数据）。`numrecs` 在 CDF-2 及以前是 32 位计数、CDF-5 为 64 位；维度计数、名称长度、维度标识、`vsize` 与 `begin` 随格式加宽。三种格式写出的文件都能被参考实现 `netCDF4` 打开，其写出的文件也能回读。

### 优化器迁移与依赖刷新

- `kernels::optimize` 由 `argmin` 迁移到 `basin`。`minimize` 接收 `FnProblem`（代价闭包、可选梯度闭包、可选盒约束）与 `SolverKind`，返回普通的 `MinimizeResult`；`SolverKind` 覆盖 24 个求解器：一阶（梯度下降、非线性共轭梯度、投影梯度、TNC）、拟牛顿（BFGS、L-BFGS）、无导数（Nelder–Mead、BOBYQA、NEWUOA、LINCOA、COBYLA、DIRECT、MADS、SLSQP）、全局/随机（CMA-ES、差分进化、全局最优 PSO、随机搜索、Solis–Wets、稳态遗传算法、GBNM）与标量（Brent、黄金分割）。盒约束属于问题本身：有界求解器从 `MinimizeOptions::bounds` 读取，无界求解器忽略。`argmin` 及其 `argmin-math` 垫片已移除。
- 依赖升级到当前版本：`ndarray` 0.17、`ndarray-npy` 0.11、`rand` 0.10（只保留一个版本）、`gif` 0.14、`roxmltree` 0.21、`pollster` 1、`syn` 3（代码生成工具）、Cranelift 0.138 与 `wgpu` 30（GPU 后端已适配新的 instance descriptor、pipeline layout 与 error scope API）。`num-bigint` 保留 0.4：最新的 `num-rational` 0.4 固定依赖它，而 `BigRational` 就是 `Ratio<BigInt>`。
- `basin` 以关闭默认特性的方式引入：封装层只使用其 `Vec<f64>` 后端，因此优化器不额外引入线性代数后端。

### 加固与互操作

- `kernels::convex::socp` 改用与量纲无关的收敛判据（相对对偶间隙，残差按数据范数归一），返回搜索过程中最优的迭代点而非最后一点，并在长期无进展时提前结束。法方程块复用暂存缓冲区，单次迭代不再按变量个数反复分配。数据量级较大的实例此前无法达到绝对 `1e-10` 间隙，现在可以正常收敛。
- `io::parquet` 在列块重复值较多时写出字典页与 `RLE_DICTIONARY` 数据页（`WriteOptions::use_dictionary`），并新增 gzip 编解码（自实现定长 Huffman DEFLATE 编码，复用已有 inflate 解码）。读取端支持 Data Page v2 头，包括级别与数值分段压缩以及 v2 页使用的 RLE 布尔编码。写出的文件已用 Apache Arrow 回读校验，Arrow 写出的 v1/v2 页、Snappy、gzip、字典与 PLAIN 编码文件也已在本实现中回读。
- `io::hdf5` 写出的对象头消息把八字节对齐填充计入声明长度，并在浮点数据类型消息中写入符号位位置。superblock v0 与 v2 两种写法的文件均可被 HDF5 库直接打开，包括数值、定长字符串与属性数据。
- `kernels::meshgen::delaunay` 缓存每个三角形的外接圆，用浮点过滤判定包含关系，仅在误差范围内回退到精确谓词。结果与纯精确判定一致，逐点插入不再为每个三角形支付一次精确展开。
- `io::parquet` 读写 Zstandard 与 `LZ4_RAW` 压缩页（原有 Snappy、gzip 之外），并把 `DECIMAL` 列按其物理存储（`INT32`、`INT64`、`BYTE_ARRAY`、`FIXED_LEN_BYTE_ARRAY`）统一暴露为除以 `10^scale` 的 `f64`。LZ4 块编解码在本模块实现，Zstandard 使用官方 `zstd` crate。Arrow 以各受支持编解码器、两种数据页版本写出的文件均可回读，本实现写出的文件也可由 Arrow 读取。
- `io::hdf5` 支持稠密新式组（fractal heap + 版本 2 B 树）、枚举（按底层整数暴露）与变长字符串（经全局堆解析）。h5py 3.16 以 `libver="latest"` 写出的文件（含一个上千链接的组）可回读。
- `rules::calculus::risch` 的塔深上限由 3 提升到 5 层。

### 构建与兼容性

- 源码构建需要 Rust 1.100 或更新版本。默认生成 Rust rlib。C 库按独立流程构建并携带生成的头文件。C 归档提供中英文说明，可从解压目录直接链接库和运行示例。[构建指南](docs/BUILDING.zh.md)
- C 批量计算允许传入未初始化的输出缓冲区。`ffi` 构建配置允许受保护的入口将可展开的 Rust panic 转为 C 错误码。
- 清单 API 版本保持为 2。不同构建的内核共享序列化契约前，应核对清单哈希。[内核嵌入](docs/KERNEL.md)

## [3.0.0-rc.1]

The release candidate establishes the graph rewriting engine, numerical ABI and standard operator registry on the 3.0 branch. The release notes above cover changes to this architecture through the beta series.

该候选版本在 3.0 分支建立图重写引擎、数值 ABI 与标准算子注册表。上方发布说明覆盖该架构在 beta 系列中的变更。
