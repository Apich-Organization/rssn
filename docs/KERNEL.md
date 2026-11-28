# The embedding kernel

[中文](KERNEL.zh.md) · [Dimensions and units](DIMENSIONS.md)

RSSN exposes two embedding modules. They share the graph, operator table and dimension arena.

| Module | Interface |
|---|---|
| `rssn::kernel` | Manifest, dynamic values, signatures, plugins, codecs and structured results |
| `rssn::dim` | Dimension terms, dimension arena, unit registry and physical constants |

## Kernel contracts

### 1. Encodings

`kernel::codec` is the only place that turns a value into bytes, and it always
uses `bincode-next`:

* `encode_binary` / `decode_binary` — native deterministic bincode-next 3.2, for
  hashing and for storage. The same value always produces the same bytes.
* `encode_cbor` / `decode_cbor` — deterministic CBOR, for interchange.
* `digest` / `short_digest` — BLAKE3 over the deterministic encoding.

`KernelManifest::to_json` and the C ABI's JSON accessors provide inspectable
metadata and results. The build embeds required `signatures.overlay.json`
and `signatures.witnessed.json`. Edit the overlay to declare a stronger
contract. Regenerate the witnessed table from replayed evidence. Missing files
fail the build. [Build guide](BUILDING.md)

Rust hosts with Serde types can also use `encode_deterministic` and
`decode_deterministic`. The native and Serde encoders share the deterministic
binary format for the supported kernel records. The native decoder expects one
complete value. Deterministic CBOR uses the self-describing Serde representation.

### 2. Operator identity

An operator is described once, by `OpDescriptor`, and that description becomes:

* a `Signature` (declared, from `signatures.overlay.json`, or inferred);
* an `OpKey` — `(rule_set, name, signature_hash)` — which is stable across
  builds and registration orders;
* an `OpInfo` in the manifest.

```rust
let manifest = session.manifest();
let sin = manifest.op("sin").expect("sin is a core operator");
assert!(!sin.overloads.is_empty());
let key = &sin.keys[0];                 // cache on (manifest_hash, key)
```

Every field of an `OpInfo` is filled for every operator, because a host that
generates an interface cannot special-case a hole. Two of them are filled by
default when a registration says nothing:

* `category` — the family a host groups operators by. A registration may name a
  finer one (`.category(..)`); otherwise it is the rule set's name, so
  `category` is never empty.
* `doc`, `doc_source` and `doc_path` — the description, where it came from, and
  the function whose rustdoc supplied it (next section).
* `since` — the release an operator appeared in. It is empty when the
  registration does not record one; a host must read an empty value as
  "unrecorded", not as "unreleased", and `stability` is what says whether an
  operator is settled.

`tests/manifest_completeness.rs` states these invariants and fails on an
operator that violates one.

### 3. Signature derivation

A host generating typed bindings cannot work from names. Every non-leaf
operator therefore has at least one signature, and it is *derived from what the
registration states*, never guessed:

| The operator declares | Therefore |
|---|---|
| a scalar evaluator (`eval`) | `n` reals in, a real out |
| `OpFlags::PREDICATE` | the result is `Bool` |
| a `Binder` at index `k` | argument `k` is the bound `Symbol` |
| `Arity::Variadic` | the last parameter repeats |
| `OpFlags::LEAF` | no call signature; a `dynamic_reason` instead |

and every non-leaf operator also has the **symbolic contract** — `n` terms in,
one term out. An operator with a scalar evaluator therefore has two overloads,
`(Real, …) -> Real` and `(Term, …) -> Term`, which is exactly what the engine
does: `sin(0.5)` folds, `sin(x)` stays a term. `SignatureSource` says where a
signature came from: `Declared` (the rule set wrote it), `Overlay` (a human
corrected it in `signatures.overlay.json`), `Derived` (the rules above),
`Witnessed` (derived, plus overloads the engine was seen to satisfy; next
section), or `None`.

The overlay is the type source for operators whose rule set declares no
signature, and only for those: a declared signature wins, so an overlay entry
for an operator that already declares one would be dead weight that silently
disagrees with the code the day the declaration changes. The file therefore
holds exactly the operators the manifest reports as `overlay`-sourced, and
`tests/signature_overlay.rs` fails when the two sets drift apart in either
direction — an entry for an operator that is not registered, or an entry the
manifest never reads.

A registration that wrote its own description keeps it; the rustdoc is appended
when that description is thin, so `doc_source` stays `declared` and
`doc_path` stays empty only when the operator's own text is what a host reads.

### 3a. Witnessed signatures

For operators such as `factorint` or `groebner`, a registration that supplies
only arity derives the broad symbolic contract `(Term, ..) -> Term`. Evidence
from successful graph computations can give a more specific result.

`kernel::sort` closes the gap with evidence the engine already produces. A
finished run leaves a request such as `factorint(360)` in the same equivalence
class as its answer `list(list(2, 3), list(3, 2), list(5, 1))`; that class is a
typing judgement. The module

1. defines a small lattice of **sorts** (`Int`, `Rat`, `Float`, `Sym`, `Expr`,
   `Bool`, `Prop`, `Obj<name>`, lists with a rectangularity flag, tuples) with a
   join, and classifies any e-class into one using only what the graph states
   (literals, list constructors, the `PREDICATE` flag, a non-scalar dimension
   rule);
2. with `RSSN_SORT_TRACE=<file>` set, logs one *fact* per operator application
   at the end of every top-level run: operator, argument sorts, result sort,
   and the call that showed it;
3. `signature-mine` treats the log as a source of *calls*, not of truths (a
   test may have loaded only some rule sets, and an application is never
   evidence about itself: `gcd(12, 18)` alone in its class says nothing about
   `gcd`). It runs every distinct call again on the full engine, keeps the
   ones that are answered, clusters them by the kinds of their arguments and
   joins each cluster (anti-unification).
4. Each candidate is then **probed**: the numbers in its witnesses' reduced
   arguments are shifted, shapes kept, and the call is run again. This is what
   stops `sin(0) = 0` from typing `sin` as `(Int) -> Int`, or one matrix of
   perfect squares from typing `matsqrt` as integer-valued; with fewer than
   three facts and a symbolic argument, exact number leaves of the result are
   widened to expressions as well.
5. A signature is kept when it says something about the *answer* and its
   parameter list is comparable with, or disjoint from, every other overload
   of the operator. One whose parameters equal an overload the registration
   already states is kept as a **refinement** of that overload's result
   (`holonomic_zero: (Term, Term) -> Bool` next to the symbolic contract),
   listed separately in the manifest and taking no part in resolution.

```text
RSSN_SORT_TRACE=sorts.tsv cargo test --lib       # the test suite is the corpus
cargo run --profile tools --bin signature-mine -- sorts.tsv
cargo run --profile tools --bin signature-mine -- --check   # replay the shipped table
```

The required `signatures.witnessed.json` table is generated and versioned. The
manifest records `signature_source: witnessed`, the source of each overload in
`overload_sources`, sharper results in `refinements`, and replayable calls in
`witnesses`.

What a witnessed signature promises is what a gradual type promises at a
dynamic boundary. It is evidence, not a totality proof: the parameter side may
be narrower than what the operator accepts (the symbolic contract stays as the
most general overload, so nothing becomes uncallable), and a host must validate
the result it decodes. Overload sets stay unambiguous by construction: a
witnessed signature is kept only when its parameter list is comparable with, or
disjoint from, every other overload (`sort::coherent`), and the overload
matcher breaks a tie between comparable overloads by pointwise dominance.

### 3b. Documentation provenance

`doc` is markdown, `doc_source` says where it came from, and `doc_path` names
the public function whose rustdoc completed it. A registration states the
function it exports in `rust_path`, which is the path a host would name the item
by; the harvest records the path the item is *defined* at. Three cases bridge
the two, and each requires a unique candidate — a description that might belong
to another function is worse than none:

1. the same path;
2. the same path with the generic arguments of its types removed
   (`FdMode<F>::parallel` against `FdMode::parallel`);
3. the same last two segments, or the same name in a submodule of the named
   module, when exactly one documented function matches — which is what makes a
   `rust_path` that names a public re-export
   (`crate::sim::md_engine::FlexibleBarostat::energy`) find the item defined in
   `md_cell`.

A registration that wrote its own description keeps it; the rustdoc is appended
when that description is thin, so `doc_source` stays `declared` and
`doc_path` stays empty only when the operator's own text is what a host reads.

### 4. Dynamic operator reasons

Every operator is either typed or explained. The explanation is a
`DynamicReason`:

| Reason | Meaning |
|---|---|
| `request` | reduced by rules; the result type is whatever the rule that fires produces |
| `variadic` | the arity is not fixed |
| `fully_polymorphic` | it accepts any type |
| `shape_dependent` | a container whose shape depends on the arguments |
| `opaque_payload` | it carries host data the kernel does not interpret |
| `not_yet_specified` | a declaration is missing; the audit treats this as a defect |

`cargo run --locked --bin signature-audit` reports the live operator counts by
signature source. Its exit status is nonzero if an operator has neither a
signature nor an accepted reason.

### 5. Determinism

A manifest built from the same operator table and the same rule sets is
byte-identical, in any process, in any registration order. The mechanism:

* every list in the manifest is sorted and deduplicated;
* no ids, addresses, timestamps or hash-map iteration orders appear;
* `manifest_hash` is BLAKE3 over the deterministic encoding of the manifest
  with the hash field zeroed.

A host keys caches on `(manifest_hash, OpKey, target)`.

## Calling an operator dynamically

```rust
use rssn::kernel::{CallCtx, OpError, Value};

let session = rssn::api::Session::new();

// By name, with values.
let value = session.call_value("sin", &[Value::Float(0.5)], &CallCtx::new())?;

// An answer that does not reduce comes back as a handle to a term.
let x = Value::Symbol("x".to_owned());
let derivative = session.call_value("diff", &[x.clone(), x], &CallCtx::new())?;
```

`CallCtx` carries the target (symbolic or numeric), the budget, the dimension
mode, assumptions, symbol bindings, substitutions, host callbacks, and an
optional trace. It is read-only for the operator being called; only the trace
is appended to.

Errors are structured, not strings:

```rust
match session.call_value("sin", &[], &CallCtx::new()) {
    Err(OpError::Arity { expected, got }) => { /* ... */ },
    Err(OpError::UnknownOp(key)) => { /* ... */ },
    _ => unreachable!(),
}
```

### Choosing between overloads

`call_value` resolves by arity and then by type, most specific first; a tie is
`OpError::Ambiguous` rather than a silent choice. A generated facade should
resolve once and then call through the key it got, so the static choice and the
executed operator cannot disagree:

```rust
let args = [Value::Float(0.5)];
let key = session.resolve_overload("sin", &args)?;   // OpKey
session.call_key(&key, &args, &CallCtx::new())?;     // exact dispatch
```

`Session::overloads(name)` lists every signature with its key. An operator with
no signature has no overloads and is reachable by name only; its key carries an
all-zero `signature_hash`.

The match lattice is `Exact` > `Widened` (an `Int` into a `Real`, a `Real` into
a `Complex`) > `Shaped` (same constructor, looser argument) > `Generic` (a type
variable or `Dyn`), and a call is scored by its *worst* parameter.

### Fitting and reusing a polynomial

A portable polynomial stores descending coefficients in the normalized coordinate
`x / coordinate_scale`. Fit the model once and pass it back to the kernel to
evaluate a vector of coordinates in a single call. Sample order is preserved.

```python
import rssn

with rssn.Session() as session:
    model = session.call_value(
        "num_polynomial_model_interpolate_lagrange", [0.0, 2.0, 4.0], [1.0, 5.0, 17.0], 2.0
    )
    values = session.call_value("num_polynomial_model_eval_batch", model, [3.0, 1.0])
    print([v.float() for v in values])  # [10.0, 2.0]
```

For measured data the positive scale and query coordinates share the units of
the fitted axis. The predictions carry the units of the ordinates. A model
built with `num_polynomial_model_from_roots` is monic in the normalized
coordinate. To chain two models, pass the outer model first to
`num_polynomial_model_compose`. Its coordinate dimension must match the
inner model’s ordinate dimension.

### Evaluating special functions in batches

Pass one rank-one Float tensor to `num_special_gamma_batch` and read back a
tensor in the same order. The kernel uses its parallel execution policy for
larger inputs. The log-gamma, digamma, error function, complementary error
function and Bessel variants share this interface.

```python
import rssn

with rssn.Session() as session:
    result = session.call_value("num_special_gamma_batch", [1.0, 2.0, 3.0, 4.0])
    print([item.float() for item in result])  # [1.0, 1.0, 2.0, 6.0]
```

Inputs and outputs are dimensionless. The value ABI reports a domain error
when the scalar kernel produces a nonfinite value, as it does at a gamma pole.

### Hyperbolic charts and Weierstrass periods

Use `num_uhp_to_disk` and `num_disk_to_uhp` to move complex points
between the open upper half-plane and the open Poincaré disk. The two
distance operators `num_uhp_distance` and `num_poincare_disk_distance`
measure the same geodesic length in curvature −1 coordinates. For
example, the half-plane points `i` and `4i` are separated by `ln(4)`.
Each input and the result are dimensionless. Boundary points and chart
images that round outside the open target chart return a domain error.

The Weierstrass derivative `num_weierstrass_p_prime_complex(z, omega1, omega2)`
takes two complex half-periods. The fundamental periods are
`2 * omega1` and `2 * omega2`. Evaluate the associated named fields
`g2` and `g3` with `num_weierstrass_invariants_complex(omega1, omega2)`.
These calls require a finite non-collinear pair. The derivative rejects
lattice poles. The calculation uses binary64, and a nonfinite result
from the underlying theta evaluation returns a domain error.

### Updating vectors and scanning batches

`num_vector_axpy(alpha, x, y)` computes `alpha*x + y`. Each scaled increment
must have the same dimension as its destination. `num_vector_axpby` scales
both vectors and requires the two products to share a dimension. Both calls
return an owned vector.

```python
import rssn

with rssn.Session() as session:
    updated = session.call_value("num_vector_axpy", 2.0, [1.0, 3.0], [5.0, 7.0])
    print([x.float() for x in updated])  # [7.0, 13.0]
```

The Hadamard product and outer product multiply the input dimensions.
`num_vector_hadamard_div` divides them. The fixed-block Euclidean norm and
inclusive prefix sum retain the input dimension, and `num_par_dot_fixed`
returns the product dimension over the common prefix. In the fixed-block
operators, the block size determines the floating-point grouping and
`grain` sets the threshold for parallel work. Keep the block size fixed when
comparing results across worker counts.

### Exact integer and modular arithmetic

Numerical vector operators under `num_finite_field_` use `Tensor<1, Int>`,
accept the full `u64` range and support any modulus of at least 2, including
composite moduli. Add, subtract, AXPY and dot product require equal vector
lengths. Strict batch inversion reports a domain error if any element is
noninvertible. `num_finite_field_batch_inverse_skip_zeros` leaves a zero in
each noninvertible position. Empty vectors are valid and an empty dot product
is zero.

`num_number_theory_gcd_u128`, `num_number_theory_mul_mod_u128` and
`num_number_theory_mod_pow_u128` accept the full unsigned 128-bit range.
Modular scalar operations require positive moduli. Inputs and outputs are
dimensionless exact integers. See the [binding contract](../bindings/api/CONTRACT.md)
for wide integer cells in each host language.

### Passing measured floats through the C ABI

Use `rssn_value_with_coherent_dimension` to attach a concrete dimension to an existing
Float or dense binary64 tensor with values in coherent units. For metres
pass `{"base":"Length"}` as UTF-8 JSON with its byte length. The result
has its own value cell in the same session. Rank, empty extents and the Float
representation survive the attachment. Supply a dimensioned zero where
an operation requires a measured argument.

The JSON uses the same dimension schema as the manifest. Convert named
units to coherent values before attaching the dimension. Use
`rssn_value_json` to inspect the resulting wire dimension.

### Measured tridiagonal systems

Use `num_dense_solve_spd_tridiagonal(diag, off, rhs)` for a symmetric
positive-definite tridiagonal system. Use
`num_dense_solve_tridiagonal(sub, diag, sup, rhs)` for a general
tridiagonal system with nonzero elimination pivots. Both take linear time
and storage in the system size. `diag` and `rhs` have length `n`. The
off-diagonals have length `max(n - 1, 0)`. The diagonals share one
dimension, and the solution has the right-hand-side dimension divided
by the matrix dimension. The general solver uses the Thomas algorithm
without pivoting. For a nonsingular system requiring a pivot, use
`num_dense_solve` with a full matrix.

### Measuring meshes

`num_meshgen_delaunay_2d_incremental(points)` accepts 3 to 4096 distinct
planar points and returns counter-clockwise triangles with zero-based vertex
indices. The coordinate tensor has one common dimension.
`num_meshgen_triangle_areas(points, triangles)` and
`num_meshgen_tetra_volumes(points, tetrahedra)` produce areas and volumes in
cell order, with squared and cubed coordinate units. Out-of-range vertex
indices give domain errors. Degenerate cells have zero measure.

`num_meshgen_min_angle_deg` reports the smallest angle in degrees.
`num_meshgen_boundary_faces` keeps faces incident to exactly one tetrahedron
and preserves their winding. Positively oriented tetrahedra give outward faces.

### Combinatorial counts and beta quantiles

`num_harmonic(0)` returns zero. `num_ln_binomial(n, k)` evaluates the
logarithm for `0 <= k <= n` and remains available when the coefficient
exceeds binary64. `num_multinomial(n, parts)` checks that the parts sum
exactly to `n` and reports overflow as a domain error. For arbitrary-precision
integer coefficients, use the symbolic `binomial` and `multinomial` operators.

`num_inverse_beta_reg(a, b, p)` accepts finite positive shape parameters and
a probability in `[0, 1]`. Endpoint and unit-shape cases use analytic
formulas. The general case brackets the solution in binary64, with precision
set by the forward probability evaluation. An unrepresentable evaluation
returns an error.

### Zero-phase filtering

`num_signal_filtfilt(b, a, samples)` applies a digital filter in both directions
with odd reflection and steady-state initialization. The coefficient vectors
`b` and `a` are dimensionless. The result has the original length and sample
units. Padding is limited to three times the filter order or one fewer than
the number of samples, whichever is smaller. A singleton uses the steady
response of the two passes. Empty coefficients, a zero leading denominator,
an undefined steady state, and non-finite values produce domain errors.

### Polynomial root bounds and the unit disk

`num_cauchy_root_bound`, `num_cauchy_sharp_root_bound`, and
`num_fujiwara_root_bound` accept finite complex coefficients in ascending
order. `num_enestrom_kakeya_bounds` accepts at least two strictly positive
real coefficients and returns a record with `lower` and `upper` radii.
The input polynomial uses a dimensionless coordinate. Its coefficients
may share ordinate units, which cancel from the bounds. A bound that
exceeds binary64 range produces a domain error through the value API.

For a normalized univalent function on the unit disk,
`num_koebe_growth_bounds(r)` and `num_koebe_distortion_bounds(r)` give
`lower` and `upper` bounds for the function and derivative magnitudes.
`num_koebe_rotation_bound(r)` gives the bound on the continuous derivative
argument in radians. The dimensionless radius must satisfy `0 <= r < 1`.

### Primitive roots and exact polynomial products

`num_number_theory_extended_gcd_u64(a, b)` returns the exact record
`{gcd, x, y}` with signed Bézout coefficients. `num_number_theory_mod_inv_u64`
returns the least nonnegative inverse for a coprime input and modulus at
least two. The Jacobi and Legendre operators distinguish an invalid
denominator from the valid value zero.

`num_number_theory_primitive_root` returns the least positive generator
for a cyclic multiplicative group of residues.
`num_finite_field_primitive_root_prime` requires a prime modulus and also
accepts two. `num_finite_field_poly_mul_ntt_mod(a, b, p)` returns the
coefficient convolution under an odd prime `p`. Use the same coefficient
order in both inputs. The next power of two of a nonempty result length
must divide `p - 1`. An empty input yields an empty tensor for a valid
prime. This value operator applies the standard native tensor size limit.

Use `num_finite_field_ntt_plan(modulus, length)` to validate an odd prime
and a power-of-two length that divides `modulus - 1`. The portable result
contains `modulus` and `length`. Pass it to
`num_finite_field_ntt_transform(values, plan, inverse)` with an equally
long tensor of unsigned 64-bit coefficients. The forward transform uses
the root derived from the smallest positive generator, and the inverse
divides by the length in the field. Both calls use ordinary residues.
Cached twiddles stay internal and a reconstructed plan is revalidated
on use. Length must fit the native value budget.

### Survival curves and Cox regression

`num_survival_kaplan_meier(time, event, confidence)` returns event times,
at-risk and event counts, probabilities, Greenwood variance and log-log
confidence bounds. The time tensor may carry a physical dimension.
`num_survival_median(curve)` returns the first time at which the curve
drops to one-half or below, or an empty option when it never does.

`num_survival_fit_cox(time, event, covariates, ties, max_iterations, tolerance)`
fits the partial likelihood on dimensionless covariate rows without an
intercept. Censored observations remain in their time's risk set.
Use `"breslow"` or `"efron"` for ties. The fitted record contains the
coefficients and convergence diagnostics. Pass it to
`num_survival_predict_cox(fit, covariates)` for relative hazards. Prediction
requires the same column order and count as the design.

### Rectangular CSR operations

`num_csr_transpose` preserves the dimension of stored values.
`num_csr_matmul_sparse` and `num_csr_matmul_dense` multiply the
dimensions of both inputs. A structural CSR record contains `rows`,
`cols`, `indptr`, `indices`, and `values`. Column indices are zero-based
and strictly increasing per row.

`num_csr_lsqr(matrix, rhs, damping, tolerance, max_iterations)` solves
rectangular least squares. The damping coefficient carries the same
dimension as the stored entries and the solution has RHS dimensions
divided by the matrix entry dimension. The return value contains the
solution, unaugmented relative residual, iteration count and convergence
flag. A least-squares optimum can have a nonzero residual.

### Sparse sequential quadratic programming

`num_sqp_solve(initial, shape, tolerance, trust_radius, max_iterations,
max_evaluations, evaluate, hessian)` solves a normalized nonlinear
program with typed callbacks. The evaluation packet holds the objective,
gradient, constraints and sparse Jacobian at one point. The Hessian callback
receives the point and both multiplier vectors and returns the upper
triangle of the Lagrangian Hessian. Sparse entry records contain zero-based
`rows`, `columns` and `values` tensors. Equality constraint rows precede
inequalities in the form `g(x) <= 0`.

Normalize decision variables, the objective and constraint functions before
building these callbacks. The solver calls the evaluation callback once
per consecutive point and shares the work budget across evaluation and
Hessian calls. It reports the final point, multipliers, KKT residual,
convergence and callback counts. The `num_sqp_oracle_*` operators let callers
inspect each typed oracle component before solving.
The callback interface accepts up to 128 coordinates and 128 constraints
with up to 16,384 triplets per sparse record.

### Exact decimal exchange

Use `num_mp_decimal_from_text(text, precision)` to round a decimal literal
directly to a dyadic with ties to even. The `from_text_in` variant takes a
rounding mode. Use `num_mp_decimal_to_text(value, digits)` to format an exact
dyadic and `to_text_in` for directed formatting. Digit arguments accept
1 to 10,000 significant digits. The decimal parse and enclosure routes accept
up to 100,000 bytes and a signed decimal scale within ±100,000 before removing
trailing zeros. A result outside the binary exponent range gives a domain error.

`num_mp_decimal_enclose` returns a dyadic ball that contains the exact input.
`num_mp_decimal_format_ball` returns text when both endpoints round to the
same result. The text-only `num_mp_decimal_round_text` accepts a decimal scale
within ±40,000 and rounds without a binary conversion.

### Exact binary multiprecision values

`num_mp_real_round(mantissa, exponent, precision, mode)` returns a
dimensionless record `{mantissa, exponent, precision}`. All fields are
exact integers. The value is `mantissa * 2^exponent`. A nonzero mantissa
has exactly `precision` bits in magnitude. Zero has exponent zero and
keeps the requested precision.

Choose 1 to 32,768 bits of precision. The constructor accepts input
mantissas up to 65,536 bits. Input and normalised exponents have magnitude
at most 262,144. Valid rounding modes are `nearest` (ties to even),
`zero`, `away`, `up` and `down`.

Use `num_mp_real_exp(value, precision)` for an exact dyadic input and
`num_mp_ball_exp(interval, precision)` for a closed dyadic interval. Both
return a dimensionless dyadic ball encompassing the result. The midpoint
has the requested bit precision, and the radius includes rounding and
series error. Exponents outside the supported range give a domain error.

Use `num_mp_real_zero` and `num_mp_real_one` for constants at a specified
precision. `num_mp_real_from_i64` rounds a signed 64-bit integer to nearest-even.
`num_mp_real_from_bigint` also accepts wider integers and an explicit rounding
mode. `num_mp_real_from_f64` rounds a finite binary64 input to the requested
precision, with both signed zeros represented by the same dyadic zero.
`num_mp_real_to_f64` rounds directly to binary64 with ties to even. It preserves
subnormals and signed underflow to zero and reports overflow as a domain error.
`num_mp_real_sign` returns -1, 0 or 1 from the exact mantissa.

`num_mp_real_compare(left, right)` returns -1, 0 or 1 for exact order
across precisions. `num_mp_real_add`, `num_mp_real_sub`, `num_mp_real_mul` and
`num_mp_real_div` accept `(left, right, precision, mode)` and round the
exact result once. Division rejects a zero divisor. Keep the record
fields as integers when passing a result back through the value API.

`num_mp_ball_add`, `num_mp_ball_sub`, `num_mp_ball_mul` and `num_mp_ball_div`
accept `(left, right, precision)` and enclose the sums, differences, products
and quotients of two closed intervals written as `{midpoint, radius}`. Both fields are exact
dyadic records as above, with `radius` nonnegative. The record denotes
`[midpoint - radius, midpoint + radius]`. The requested precision sets
the midpoint precision. The radius may use a different precision to
preserve the enclosure. The divisor interval must exclude zero, including
at its endpoints. Exponents outside the representable range give a domain
error.

A complex rectangle is a dimensionless record `{real, imag}` whose two
fields are closed dyadic intervals. Use `num_mp_complex_new` to combine
two intervals, `num_mp_complex_point` for exact coordinates, or
`num_mp_complex_from_real` to lift an interval with zero imaginary part.
Addition, subtraction, multiplication, division, reciprocal, squared
magnitude and signed integer power use the same precision range as real
interval arithmetic. Squared magnitude returns a dyadic interval.
Division, reciprocal and negative powers report a domain error when the
divisor enclosure includes the origin.

### Time-series statistics, reconstruction, and state-space models

`num_timeseries_integrate(differences, initial)` undoes the number of first
differences specified by the length of `initial`. The inputs share one
dimension and the reconstructed samples retain it.

`num_timeseries_rolling_moments(samples, window, ddof)` returns `mean` and
`variance` vectors for complete windows in start order. Means retain the
sample dimension and variances square it. Use `ddof=0` for population
variance or `ddof=1` for sample variance with a window of at least two.
A positive width larger than the series gives two empty vectors. Samples
must be finite and `ddof` must be smaller than the width.

`num_timeseries_kalman_filter(y, f, h, q, r, x0, p0, observation_scale)`
returns filtered and predicted means, the rank-three filtered covariance
history in time order, and log likelihood. The first prediction uses the
given prior. State components share one dimension, observations share
another, and the observation matrix maps between them. The positive
`observation_scale` uses observation units and defines the reference for
the dimensionless log likelihood. Covariance matrices use squared units
of the corresponding state or observations.

`num_timeseries_kalman_smooth(y, f, h, q, r, x0, p0)` returns one
smoothed state mean per observation in time order. The solve requires a
positive-definite predicted covariance at each backward step.

### Directed roots for exact dyadics and intervals

`num_mp_real_cbrt_in(value, precision, rounding)` rounds an exact real
cube root, including negative inputs. `num_mp_real_cbrt` and
`num_mp_real_hypot` give dyadic intervals enclosing the cube root and the
Euclidean norm of two exact dyadics. Use `num_mp_ball_sqrt` for the
square roots of a closed interval wholly in the nonnegative real numbers.

`num_mp_real_limits()` reports the supported bit precision and binary
exponent range. Dyadic and interval accessors preserve the exact integer
fields. Use `num_mp_real_rounding_mode(rounding)` to obtain the canonical
spelling of a supported rounding mode.

### Measured state-space simulation

`num_control_c2d_zoh(a, b, c, d, dt)` discretises a continuous
model under a zero-order hold and returns the portable record
`{a, b, c, d, dt}`. Specify a positive `dt` with time units. For
`n` states, `m` inputs and `p` outputs, pass matrices with shapes
`[n,n]`, `[n,m]`, `[p,n]` and `[p,m]`. A continuous `a`
has inverse-time units. Its discrete counterpart is dimensionless,
and the discrete `b` gains one factor of time.

`num_control_simulate_discrete(model, x0, inputs)` accepts the
record directly. `num_control_simulate_continuous(a, b, c, d, x0,
inputs, times)` also accepts a strictly increasing measured time
grid, with one input row per time. Both return one output row per input
row, starting with the output at the given initial state, including
the direct feedthrough `d`. Continuous simulation holds each input
until the following sample.

## Plugins

A plugin is a bundle of extensions registered through `Registry`:

```rust
use rssn::kernel::{Registry, SimplePlugin};

let plugin = SimplePlugin::new("acme/twice", "0.1.0", |reg: &mut Registry<'_>| {
    reg.register_op(
        rssn::graph::OpDescriptor::new("acme_twice", rssn::graph::Arity::Fixed(1))
            .signature(rssn::kernel::Signature::scalars(1, rssn::kernel::TypeExpr::Real)),
    )?;
    reg.register_rewrites(
        rssn::graph::Tier::Normalize,
        &["acme/twice-zero: acme_twice(0) => 0"],
    )?;
    Ok(())
});
session.register(Box::new(plugin))?;
```

Registration refuses an operator that declares neither a signature nor a
`DynamicReason`: a plugin cannot make the manifest incomplete. It also refuses
a name that clashes with an existing operator of a different shape; adding an
overload of the same name requires `.overload()` on the descriptor.

## Certified complex rectangle functions

`num_mp_complex_abs(value, precision)` returns a closed dyadic interval
for the modulus of a complex rectangle. `num_mp_complex_exp` returns a
rectangle enclosing the exponential when the imaginary midpoint is
within ±2^20 radians. `num_mp_complex_ln` returns a rectangle enclosing
the principal logarithm and rejects an input that touches the
non-positive real axis. Arithmetic outside the exponent or precision
limits returns an error.

The circular operators `num_mp_complex_sin`, `num_mp_complex_cos` and
`num_mp_complex_tan` and their hyperbolic counterparts
`num_mp_complex_sinh`, `num_mp_complex_cosh` and `num_mp_complex_tanh`
also take and return dyadic rectangles. The circular operators reduce the
real midpoint and the hyperbolic operators reduce the imaginary midpoint.
These angle reductions support ±2^20 radians. Tangents return a domain
error when the computed denominator rectangle includes zero. On the
imaginary axis for tangent and the real axis for hyperbolic tangent, the
bounded real hyperbolic tangent route handles large finite inputs.

`num_mp_complex_sqrt` encloses the principal root. On the negative real
axis the result has positive imaginary part. An input rectangle that
straddles the cut may enclose roots with both signs of the imaginary part.
`num_mp_complex_pow` takes a base and an exponent rectangle. A point
exponent that is an exact real integer uses the signed 64-bit integer
route, including zero. The general case computes the principal value via
the certified logarithm and exponential. Its base must avoid the
non-positive real axis, and the exponential angle range still applies.

The `num_mp_complex_{asin,acos,atan,asinh,acosh,atanh}` operators enclose
principal inverse functions of dimensionless dyadic rectangles. Asin and
acos include both endpoints for real intervals in [-1, 1], and asinh
includes the corresponding imaginary endpoints. Acosh includes the point
1 and rejects rectangles that meet the real axis below 1. Atan and atanh
reject rectangles that meet their respective imaginary and real cuts,
including their branch points. Operations return a domain error when an
intermediate enclosure reaches a singularity or exceeds the exponent range.

Use `num_mp_expression_enclose(text, precision)` when you need an exact
dyadic midpoint and outward radius from the bounded arithmetic parser.
The precision argument is a bit count. For example,
`sqrt(81) + hypot(3, 4) + cbrt(-8)` encloses the exact value 12.
Invalid syntax, unsafe operand intervals and out-of-range exponents
return domain errors.

`num_mp_ball_lipschitz(interval, mapped_midpoint, bound, precision)`
propagates the input radius through a supplied nonnegative Lipschitz
bound and encloses rounding of the mapped midpoint. The caller supplies
the exact image of the input midpoint and a bound that holds throughout
the interval.

The typed `num_mp_complex_erf` and `num_mp_complex_erfc` operators
enclose the image of a closed dimensionless dyadic rectangle. The Taylor
series uses a uniform geometric remainder bound. In the right sector
`Re(z) >= |Im(z)|`, the complementary function can use the first-omitted-term
bound for the asymptotic expansion to retain the small positive-real tail.
Very wide inputs or insufficient series budgets produce domain errors.
See [DLMF §7.6](https://dlmf.nist.gov/7.6) and
[DLMF §7.12](https://dlmf.nist.gov/7.12).

## Certified real trigonometry

`num_mp_real_sin`, `num_mp_real_cos` and `num_mp_real_tan` enclose the
trigonometric values of an exact dyadic angle in radians. The maximum
supported magnitude is 2^20 radians. Use `num_mp_real_sin_cos` for both
sine and cosine with one argument reduction. The matching `num_mp_ball_*`
operators take a closed angle interval whose midpoint is within the
supported range. Tangent returns a domain error when the computed cosine
enclosure includes zero.

`num_mp_complex_bessel_j(value, order, precision)` encloses principal Bessel J
over dyadic rectangles in both arguments. Exact signed integer orders of
magnitude at most 512 admit either side of the negative real axis. General
orders require a nonnegative real part throughout the order rectangle, at
most 512 requested bits, and an argument disjoint from the non-positive real
cut. At the exact origin, zero order gives one and orders with strictly
positive real part or nonzero exact integer order give zero. All other
arguments must satisfy `abs(Re(value)) + abs(Im(value)) <= 256` throughout
their rectangle. Use the returned radius to check how operand uncertainty
affects the result. The power-series tail and Stirling normalization use
outward bounds from [DLMF §10.2](https://dlmf.nist.gov/10.2#E2) and
[DLMF §5.11](https://dlmf.nist.gov/5.11.ii).

The typed gamma family uses the same certified Stirling expansion as
Bessel J. `num_mp_complex_gamma` and `num_mp_complex_digamma` accept
dyadic rectangles that avoid the non-positive integer poles.
`num_mp_complex_log_gamma` instead uses the analytic continuation real
on positive arguments and excludes the entire non-positive real axis.
All three accept at most 1024 requested bits. An exhausted recurrence
budget returns an error. Gamma of a supported positive integer uses the
exact factorial before outward rounding.

The typed zeta and polylogarithm operators return outward dyadic
enclosures. Use `num_mp_complex_riemann_zeta(order, precision)` with a
rectangle excluding the pole at one. The bounded Euler–Maclaurin path
supports real parts at least -32 and modulus at most 256. Integer orders
with closed formulas, including the negative even zeros, can use a
direct route. `num_mp_complex_polylog(value, order, precision)` requires
the whole argument rectangle strictly inside the unit disk and real order
at least -256. An exhausted series budget returns an error. Both accept
up to 512 requested bits. Their certified tails use the expansions in
[DLMF §25.2](https://dlmf.nist.gov/25.2#E9) and
[DLMF §25.12](https://dlmf.nist.gov/25.12#E10).

## Certified real hyperbolic functions

`num_mp_real_sinh`, `num_mp_real_cosh` and `num_mp_real_tanh` return
outward dyadic enclosures for exact real inputs. Their `num_mp_ball_*`
variants accept closed intervals. The cosine variant includes the minimum
of one if the interval crosses zero. The tangent preserves a finite
enclosure for large finite inputs as the function approaches ±1.

## Certified principal angles

`num_mp_real_atan` and `num_mp_ball_atan` enclose the inverse tangent of
an exact real value or a closed real interval. `num_mp_real_atan2(y, x, precision)`
encloses the principal angle in (-π, π]. The origin is outside its domain,
and the negative real axis has angle +π. The rectangle adapter
`num_mp_ball_atan2` rejects a box containing the origin or touching the
negative real axis unless the box lies wholly on that axis. The complex
rectangle argument adapter `num_mp_complex_arg` excludes the non-positive
real axis, including boxes wholly on the cut.

## Certified inverse sine and cosine

`num_mp_real_asin` and `num_mp_real_acos` enclose the principal inverse
sine and cosine of exact real dyadics in [-1, 1], including both
endpoints. The corresponding `num_mp_ball_*` variants enclose closed
intervals contained in that domain. They return domain errors when any
part of the input interval lies outside [-1, 1]. At one, inverse cosine
returns an exact zero-radius zero.

## Certified real powers

`num_mp_real_powi` and `num_mp_ball_powi` enclose powers with a signed
64-bit integer exponent. An exponent of zero gives one. Negative powers
reject base intervals containing zero. The rational variants
`num_mp_real_pow_ratio` and `num_mp_ball_pow_ratio` accept a signed
64-bit numerator and a positive 32-bit denominator. For a negative base
the reduced denominator must be odd.

`num_mp_real_pow` and `num_mp_ball_pow` take dyadic exponents. A
non-integer exponent on a nonzero real base requires a positive base. If
the interval exponent varies, the base interval must be strictly positive.

## Real logarithms and binary exponential

`num_mp_real_ln`, `num_mp_real_log2`, `num_mp_real_log10` and
`num_mp_real_exp2` return outward dyadic enclosures of exact real inputs.
The corresponding `num_mp_ball_*` operations accept closed dyadic
intervals. The logarithms require strictly positive inputs, including
both endpoints of an interval. The binary exponential preserves exact
powers of two when the input is an integer and the result is representable.

`num_mp_real_scale_decimal(value, exponent, precision)` scales a dyadic
by a signed integral power of ten and rounds the exact result once.
`num_mp_real_scale_decimal_in` accepts an explicit rounding direction.

## Multiprecision numerical constants

`num_mp_real_const_named(name, precision)` returns a certified nearest-even
dyadic for a mathematical constant. The names are `pi`, `tau`, `e`,
`ln2`, `sqrt2` and `phi`. Use `num_mp_real_const_enclose` to retain an
outward interval. `num_mp_real_const_pi_scaled(bits)` returns the exact
integer floor of π × 2^bits for 0 to 32,768 fractional bits. The name
catalogs report exactly which constants and expression functions the
multiprecision parser accepts.

## Physical constants

`dim::constants` carries the CODATA 2022 constants with an exact rational value,
a dimension and an SI unit, and distinguishes the ones that are **exact by
definition** (the seven SI defining constants and everything derived from them)
from the **measured** ones, which carry a standard uncertainty. `c` really is
the integer `299792458`; `hbar` is `h/(2π)` exactly, carried as a rational times
`π^-1`.

`KernelManifest.constants` lists them, so a host can generate
`Units.Constants` without a second table. Note that `G` is the gravitational
constant here and the gauss in the unit table — different namespaces, on
purpose.

## Loading a plugin from a shared library

`kernel::cdylib` defines a byte-level C ABI: one `bincode-next` buffer of
`Vec<Value>` in, one of `Result<Value, String>` out. A `#[repr(C)]` mirror of
`Value` would have to track every variant across two independently compiled
artifacts and would break silently the day one is added; the deterministic
encoding is already the stable description of a value.

It is behind the `plugin-cdylib` feature because a shared library runs
arbitrary code in the host process. `Features::plugin_kinds` reports whether the
build can load one, so a host that meets a `kind = "cdylib"` plugin in a build
without the capability refuses it with a precise diagnostic.
`tests/cdylib_plugin.rs` builds a plugin in C with `cc` and calls it end to end.

## Build capabilities

`KernelManifest.features` is what a host gates on before generating bindings:

| Field | Meaning |
|---|---|
| `dynamic_call`, `parallel_symbolic`, `higher_order` | the dynamic value ABI, the parallel symbolic phase, host callbacks |
| `jit`, `jit_branch`, `jit_index` | the optional JIT and what it lowers |
| `exact_units`, `affine_units`, `user_base_dimensions` | the exact unit layer |
| `structured_results` | `Solutions` / `OdeSolution` records cross the ABI |
| `sde` | the `sde` rule set ships, so the noise operators below are callable |
| `plugin_abi`, `plugin_kinds` | the plugin ABI, and how a plugin may reach this build |

The `sde` flag covers the SDE rule set. `sde_wiener(seed, m, steps, h)` returns
seeded Brownian increments of shape `steps × m` and dimension `sqrt(dim(h))`.
`sde_levy_areas(seed, h, dw, terms)` returns conditioned iterated integrals
with the step's dimension. `sde_scheme_info(name)` reports a named scheme's
strong order, stochastic calculus, commutativity, adaptivity and Lévy-area
requirements.

Separate native SDE bindings expose owned drift and diffusion tapes, retained
noise samples and actual integration methods. `num_analysis_sde_ensemble`
executes an ensemble with a shared evaluation budget. Rust drivers also accept
native closures. For host callback values and their validation, see the
[callback contract](CALLBACKS.md) and [numerical callback operators](NUMERICAL_CALLBACKS.md).

## Which values can enter the graph

`KernelManifest.value_kinds` lists every `Value` variant with
`in_term_graph: bool`. A variant with `in_term_graph = true` becomes a node of
the expression graph and can be passed to any operator; one with `false` reaches
only an operator with a dynamic evaluator, or comes back as a result.

`Int`, `Real`, `Float`, `Bool`, `String`, `Symbol`, `Unit`, `Quantity`,
`Tensor`, `Term` and `Dyn` enter the graph. `Complex` enters only when its
imaginary part is zero. `Sparse`, `Record`, `Variant`, `Function` and `Opaque`
do not, and the manifest says why in each `note`.

## Structured results

A result that is more than one number is a `Value::Record` with a declared
`FlattenSchema`. That is deliberate: it crosses the dynamic ABI unchanged, it
serialises with `bincode-next` for free, and the layout travels with the value.

| Type | Fields |
|---|---|
| `Solutions` | `variable`, `roots`, `multiplicities`, `conditions`, `complete` |
| `OdeSolutionView` | `t`, `y`, `flatten_order`, `events`, `nfev`, `rejected` |
| `SimResult` | `scenario`, `t`, `fields`, `schema` |

`Session::solve_structured` lifts the `solve` operator's `list(...)` into a
`Solutions`, and `OdeSolution::view` names an integrator's anonymous columns.

An empty `multiplicities` means the solver did not distinguish them; it does not
mean "all one" as a claim. `complete` is true only when the solver answered in
its exhaustive `list(...)` form; a single root found on the way is reported with
`complete = false`, because it is not a claim that the others do not exist.

`OdeSolutionView::dimensions` travels inside the record, so a host reading the
`Value` does not have to go back to Rust for it.

## Dimensions and units

`rssn::dim` owns the only unit table in the crate.

* `UnitDef` — name, aliases, `DimExpr` dimension, exact `Scale` (`rat * pi^pi`),
  optional affine offset, prefixability, display, category.
* `UnitRegistry` — `list`, `lookup`, `resolve`, `register`, `format`, `snapshot`.
  `StandardRegistry` is the SI table; `register` adds to a process-wide overlay.
* Registration is checked: a name or alias that already belongs to a different
  unit is refused, so a plugin cannot silently redefine the metre.
* Everything is exact. The scale of a foot is the rational `381/1250`, the scale
  of a degree is `pi/180`, the offset of `degC` is `27315/100`. A float would
  lose the distinction between `1 in` and `0.0254 m` after a few conversions.
* The table covers SI, the non-SI units accepted with it, the IAU astronomy
  units (`au`, `ly`, `pc`, `Rsun`, `Msun`, `Mearth`, `jup`), the CGS set
  (`dyn`, `erg`, `Gal`, `P`, `St`, `G`, `Oe`, `Mx`, `Ba`), and the common
  non-SI units of the Aevia design standard (`Da`/`u`, `kn`, `mph`, `acre`,
  and the dimensionless information units `bit`, `B`, `KiB`).
* `KernelManifest.units` carries `aliases`, `prefixable`, `category` and
  `since` as well as the name, dimension and scale, so a host can generate its
  unit table from the manifest alone.

A typed front end should declare dimensions with `Config::dimension_term`,
which takes a `DimExpr` and never parses text:

```rust
let config = Config::new().dimension_term("x", DimExpr::base("Currency"));
```

Unlike `Config::dimension(symbol, "Currency")` it cannot confuse a unit symbol
with a base dimension, and a base the arena has never seen is defined on the
way in.

`DimExpr` is the *type-level* mirror of a dimension: names instead of arena-local
ids, so it can go in a manifest and be read back by another process.
`DimExpr::to_term` interns unknown names on demand, which is how a host defines
a base dimension (`Currency`) the engine has never seen.

## Threads and handles

A `Session` is single-threaded (its store sits behind a `RefCell`, so it is not
`Sync`): use one session per thread. `Term<'s>` borrows a session and cannot
outlive it; `TermRef` is `Copy + Send + Sync` and owns nothing — it is a
`(session_id, node)` pair, refused by any session but the one that minted it.
That turns a use-after-free of a term handle into an error instead of a wrong
answer.

`fork` gives a session its own graph, which can be moved to another thread;
`merge` imports every node of another and translates dimension assumptions
between the two arenas, failing on a disagreement rather than picking a winner.
`Session::compute_all` is the supported way to compute several requests in
parallel inside one session.

`Session::with_native_object_scope` lets a Rust host run a temporary action in
the current session. The scope releases native objects created by the action
on return or unwind. Earlier objects are still available, and local handles
must stay inside the action.

## Measured CAD meshes

`num_cad_mesh_new` accepts measured Cartesian vertices and zero-based exact
triangle indices in separate three-column tensors. Mesh queries share a
cached edge-incidence summary. Signed volume has dimension `Length³` and
requires a closed consistently wound surface. Welding and repair accept a
length tolerance, return new meshes and report what changed.

`num_cad_mesh_boolean` combines two connected, outward-wound convex shells
with `union`, `difference`, or `intersection`. The result is a new measured
mesh. Empty intersections return an empty mesh. The operation returns an
error for unsupported geometry or exhausted work limits. In Aevia, use
`Rssn_Cad.mesh_boolean` with the same three operation names.

The STL operations accept and return in-memory ASCII text or byte tensors.
Pass the positive length represented by one coordinate unit in the file. A read model may repeat vertices across
triangles. Weld the mesh before examining topology. Release each handle
with `num_native_release`.

`num_cad_step_read_text` reads bounded planar STEP solids with straight
`LINE` edges and returns a tensor of mesh handles. Declared STEP length units
are converted to metres. Unitless coordinates are interpreted in metres.

## The C ABI

```c
char       *json = rssn_manifest_json();   /* rssn_string_free(json) */
uint64_t    hash = rssn_manifest_hash();
```

Both describe the same *standard* manifest snapshot, including the process-wide
units visible at the first manifest call. The JSON encoding is reused and each
call returns its own string. A host that installed plugins reads the richer
per-session manifest through the Rust API. The rest of the C interface is
summarised in [ARCHITECTURE.md](../ARCHITECTURE.md#c-interface).

## What the JIT can do

* Arithmetic, `sqrt`, `abs`, `floor`, `ceil`, comparisons, `select` and integer
  powers are native instructions.
* `Inst::Gather` is the load/index instruction: `inputs[base + offset]`, with
  the index clamped before the address is computed, so an out-of-range index
  reads NaN and never memory outside the buffer.
* An operator can declare how it compiles with `OpDescriptor::jit(Lowering)`:
  `Intrinsic`, `Extern1`, `Extern2` or `Gather`. Without a declaration, the
  backend decides: a scalar evaluator is called through a shim, everything else
  runs in the interpreter.
* The cache key is `mix(tape_structural_hash, manifest_hash ^ target_key())`,
  so two builds for different micro-architectures, or with different operator
  semantics, never share machine code.

## Prepared modular arithmetic

`num_modular_montgomery_new(modulus)` prepares an odd unsigned 64-bit modulus
greater than one. `num_modular_barrett_new(modulus)` also supports even moduli.
Both constructors reuse the same Rust reducer throughout the lifetime of a
session-owned handle. Use the `_new_prime` constructors when you need to check
primality. Barrett supports the prime 2.

The scalar operators accept exact unsigned words and return integer residues.
An inverse or quotient returns `None` when the denominator is not coprime to
the modulus. A bad input, wrong reducer type or stale handle returns an error.
Montgomery encoded-word operations require inputs below that plan's modulus.
Convert with `num_modular_montgomery_to_mont` or its 128-bit counterpart and
convert the answer back with `num_modular_montgomery_from_mont`. Pass the same
plan that encoded your residues. Raw Montgomery reduction accepts a u128 below
`modulus * 2^64`. Barrett reduction supports the full u128 range.

An immutable reducer retains a few machine words in the session's bounded
native-object storage. Keep its handle for repeated calculations and release it
with `num_native_release` when finished. Language bindings can find these methods
in the live manifest and use the same native-object ABI.

## Adding an operator

For numerical functions supported by the composable codecs, add a Rust source
path, ABI name, namespace and dimension policy to `codec-policy.tsv`. Run
`bazel run //:refresh_contracts` and use the emitted
overload through the same manifest and value-call APIs.
[Binding generation and supported contract shapes](../tools/rssn-api-gen/README.md)

For operators registered directly, including symbolic rules and callback or
native object adapters, follow the registration steps below.

Typed native object adapters derive the signature and arity from the same
argument list. At runtime the adapter checks arity before decoding. Values
are validated in depth-first order, with limits on depth, element count,
and tensor shape before the native method runs.

1. Register it with `Installer::op` or `Registry::register_op`.
2. Declare a `Signature` (`.signature(..)`), or a `DynamicReason`
   (`.dynamic_reason(..)`).
3. Add `.doc(..)`, `.category(..)` and `.since(..)` — the manifest is the
   documentation a host reads.
4. Run `bazel test //:contracts_test` and the operator tests.

## Prepared transfer functions

`num_control_transfer_new` prepares a normalized continuous-time SISO model
from finite numerator and denominator coefficient vectors in descending powers
of the Laplace variable. Use `_new_proper` when you need a state-space step
response. The native handle can be reused for frequency evaluation, zeros,
stability, series, parallel and negative feedback connections. Pole and zero
queries reuse results within that handle's session.

Choose a common reference time and express both coefficients and samples in
that scale before calling the ABI. Coefficients, frequency, time and response
are dimensionless normalized values at this boundary. An absent margin on the
internal log-frequency grid is `None`, while a pole or numerical failure
produces a checked error. Magnitude margin is in dB and phase margin is in
degrees. Release a handle with `num_native_release` when finished.

## Normalized PID and identification

`num_control_pid_response` uses a proper prepared transfer model and a
normalized time grid. `num_control_simulate_pid` accepts sampled error, a
positive interval and optional ordered output limits. Both use parallel gains
`kp`, `ki`, `kd` and an optional positive `filter_time`. Scope the gains and
samples to the same reference time and error scales. Output saturation freezes
the integrator when its next contribution would push it farther past the bound.

`num_control_era` fits a discrete single-input model from normalized input
samples and a sample-by-output tensor. The first input sample must be nonzero.
Call `num_control_era_matrices` to inspect the fitted A, B and C matrices.
`num_control_impulse_response` starts at a zero first sample and returns one
column per output. The model has no direct feedthrough. Release it with
`num_native_release`.

## Prepared discrete MPC

`num_control_mpc_prepare` takes normalized discrete A and B matrices, a
positive sample time, and a record of horizons, diagonal weights and optional
component bounds. Input weights must be positive. The prepared handle stores
the prediction maps, Hessian and constraint rows, which are reused across
`num_control_mpc_closed_loop` calls. Each sample shifts the constraints to the
current state. An infeasible or non-finite solve returns an error.

The native route limits matrix size, horizons and rollout work. Each
trajectory includes the initial state. Release the plan with
`num_native_release` when finished.

### Finite numerical transforms

The `num_transform_*_checked` exports report a value, an absolute error
estimate, `converged`, and the actual evaluation count. The Mellin, inverse
Mellin, Hilbert, Hankel and Stieltjes callbacks use dimensionless transformed
coordinates and retain the callback amplitude dimension. For scales with
physical units, normalize the coordinates and apply the analytic scale factor
to the result. Adaptive error estimates cover only the finite interval.
Laplace contours and Gaver–Stehfest compare two resolutions, so their error
estimate is empirical. Expose analytic callbacks on the entire chosen contour
and choose shifts to the right of the relevant singularities.

Mellin–Barnes uses a pole-free vertical segment. Meijer G additionally requires
the vertical line to separate the two numerator pole families. Both finite
estimates exclude the omitted contour tails.
The Borel polynomial preserves the common coefficient dimension and expects
a dimensionless series argument.

## Verified complex disks

Complex disks use `[real, imaginary, radius]` on the last tensor axis. The
special-function operators `num_complex_disk_erf`, `num_complex_disk_gamma`,
`num_complex_disk_log_gamma` and `num_complex_disk_bessel_j` take normalized
coordinates and return an optional enclosure of the entire input disk.
An enclosure that reaches a pole, crosses a branch cut, or exceeds the
supported finite range can produce `None`.

Root certification uses `num_complex_disk_certify_root` for a scalar jet
and `num_complex_disk_certify_system` for a normalized box of up to 64
disks. Supply enclosure callbacks for the full disk or box, including the
derivative or Jacobian. The scalar test needs two callback invocations,
and the system test needs three. The optional result contains a certified
root enclosure when the Krawczyk test succeeds.

`num_complex_disk_certify_rouche` proves the strict inequality
`|g(z)| < |f(z)|` along a circle. Both analytic callbacks enclose their
images over each arc-covering disk. Use 4 to 4,096 arcs and a budget of at
least two callback invocations per arc. `false` means the test did not
establish the inequality.
