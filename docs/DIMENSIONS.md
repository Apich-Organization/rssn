# Dimensions as a first-class typing layer

This document is the implementation plan for integrating dimensional analysis
into the whole engine, not just the `units` rule set. It is the contract for
the core work and for the operator-coverage workstream that follows it.

## 1. Goal

Every scalar term has a dimension, and every operator states how the
dimensions of its arguments determine the dimension of its value. A request
`compute(f(x), ...)` therefore reports, besides the answer term and numeric
value, the dimension of the answer; in strict mode an ill-typed expression
(`x + y` with incompatible dimensions, `exp(3 m)`, `ln(1 m)`, a dimensionful
exponent) is rejected with a hard error.

The dimension is an inferred attribute of an e-class, not part of the node
payload: terms are hash-consed and the same term may be dimensionless in one
session and a length in another. Dimensions are computed by an abstract
interpretation over the term DAG, cached per e-class, with dimension
variables for unknown symbols and unification for dimension constraints.

## 2. Decisions

These were fixed before implementation:

| Decision | Choice |
|---|---|
| Default mode | `Advisory` (report and diagnose, never change the answer) |
| Symbolic exponents | **Required**: `x^n`, `x^(a+b)`, `x^(1/2)`, `prod(f(k), k=1..n)` |
| Angles | Dimensionless (`rad = 1`, `deg = pi/180`), with an `Angle` kind for UX |
| Affine temperatures | `Point`/`Vector` boundary checks in the unit layer; core sees `Temperature` |
| Undeclared symbols | Unknown dimension variable (reported as unknown, not assumed dimensionless) |
| `Strict` conflicts | Hard `ComputeError::Dimension` |
| User-defined base dimensions | Supported (`Config::define_base_dimension`) |
| Bypass | Supported (`Config::bypass_dimension_checks`, `DimRule::Any`, dimension `any`) |
| Operator coverage | Every operator must carry a semantic rule; a scanning test enforces it |

## 3. Algebra of dimensions

### 3.1 Base dimensions

Base dimensions are interned names. The seven SI bases are registered in every
arena: `Length`, `Mass`, `Time`, `Current`, `Temperature`, `Amount`,
`Luminosity`. Additional bases (for example `Currency`, `Bit`, `Count`) can be
registered per session/graph.

### 3.2 Exponents

Exponents are rational functions over **atoms**:

* a rational number;
* a symbol that occurs as an exponent (`n`);
* an opaque term that cannot be reduced to arithmetic (for example
  `sin(a)`), interned by its printed form.

The rational-function representation is needed for `x^(a+b)`, `x^(a*b)`,
`x^(1/n)`, `(x^a)^b = x^(a*b)`; it also makes `x^a * x^b = x^(a+b)` and
`x^a / x^a = 1` automatic. Equality of two exponent expressions is decided by
cross-multiplication of numerator and denominator (no polynomial GCD is
required). Atoms are interned in the graph's dimension arena, so dimensions
are printable and storable.

A concrete dimension is a product of `base^exponent`; an inferred dimension is
a product of `base^exponent` and `dimvar^exponent`. The exponents of
dimension variables are normally rational constants; symbolic exponents on
variables are allowed in the representation but block solving (see 3.4).

### 3.3 Special values

* `Any` — polymorphic scalar: literal `0`, `oo`, the empty list (no scalar
  leaves, so its element dimension is unconstrained), and operators
  registered with a bypass rule. Unifies with everything.
* `NonMetric` — booleans, predicates, dimension labels (`dimension(q)`),
  unresolved containers. Compatible only with `NonMetric` (or `Any`).
* `Dimensionless` — the empty product `1`.
* Dimension variables — one per undeclared symbol / unknown function value.

### 3.4 Unification

`unify(a, b)` computes `q = a / b` and solves:

1. For each dimension variable with a constant rational exponent `r != 0`,
   bind the variable to `rest^(-1/r)` (abelian-group unification). The
   substitution is applied transitively; occurs check.
2. Variables with symbolic exponents cannot be solved; if `q` is not
   structurally `1`, report a symbolic constraint.
3. After removing variables, every base exponent must be zero. A non-zero
   constant exponent is a mismatch; a non-zero symbolic exponent is a symbolic
   constraint (`L^(1-n) = 1`, i.e. `n = 1`) and is reported as such.

This is exactly the theory of dimensions as a `Q`-vector space, extended by a
field of scalar exponent expressions. The symbolic part is conservative: it
never invents a solution that is not implied.

## 4. Units

The unit layer supplies dimensions plus exact scales and affine offsets.

* `quantity(v, u)` is an inert term. Its dimension is the dimension of the
  unit expression `u`; its numeric value in coherent SI units is `v * scale`.
* Unit expressions are parsed by the `units` rule set. The dimension core
  reads them through a resolver installed on the `quantity` operator, so
  `graph` does not depend on `rules`.
* An affine unit (`degC`, `degF`) is a point of a one-dimensional affine
  space. The unit resolver reports the `affine` flag, and the dimension core
  sees `Temperature` for a point and for a difference alike. The `units`
  kernel refuses to scale, power or add two points (`quantity(1, degC) * 2`,
  `quantity(1, degC) + quantity(2, degC)` stay unreduced) instead of quietly
  converting them through kelvin; inference reports `DimError::Affine` for
  the same shapes, so Strict mode rejects them. A point plus a difference
  written in `K` is still a temperature (`quantity(20, degC) + quantity(2,
  K)`). A `degC` inside a product or a power makes the unit unresolvable and
  is reported as `UnknownUnit`. The point/vector distinction is therefore a
  boundary check, not a full affine type (see §13).
* Dimensionless scaled units (`%`, `ppm`, `rad`, `deg`) keep their scale and an
  `Angle`/`Ratio` kind; `sin(quantity(30, deg))` converts to `sin(pi/6)`.
* Derived named units (`N`, `Pa`, `J`, ...) are a printing preference: the
  dimension is the same product of bases.
* User-defined base dimensions are registered per graph before parsing the
  request. Text parsing uses units first and bases second; `base:Name`
  forces a base dimension and resolves a custom base whose name is also a
  unit symbol (`base:T` for a custom base `T`, while bare `T` stays the
  tesla).

## 5. Operator signatures

A `DimRule` is attached to an operator as an attribute (like `Partials`), so
rules that re-register an operator do not fight over the descriptor:

```rust
pub enum DimRule {
    Unspecified,                 // coverage gap: forbidden by the scanner
    NonMetric,                   // booleans, labels
    Container,                   // a tuple/decomposition: non-metric, but
                                 // compatible with its expanded literal
    Any,                         // explicit bypass
    Same(usize),                 // result = dim of argument i
    Sum,                         // all arguments equal; result = common
    Product,                     // sum of exponents
    Quotient,                    // difference of exponents
    Power,                       // exponent argument must be dimensionless
    Root(u32),                   // exponent divided by n
    Dimensionless,               // all arguments dimensionless; result 1
    Inverse,                     // 1 / dim(argument 0)
    Square,                      // dim(argument 0)^2
    Diff { f: usize, x: usize }, // dim(f) / dim(x)
    Integral { f: usize, x: usize }, // dim(f) * dim(x)
    SumOver { f: usize, index: usize },
    ProductOver { f: usize, index: usize }, // dim(f) ^ count
    Family {                     // a distribution: scale/rate source plus
        result: FamilySource,    // dimensionless shape parameters
        same: &'static [usize],
        dimensionless: &'static [usize],
    },
    MatrixFunction {             // matfun(f, x, A): homogeneous degree k
        function: usize,         // gives dim(f) * entry^k
        variable: usize,
        matrix: usize,
    },
    Quantity,                    // quantity(v, u): the unit layer supplies it
    Custom(fn(&DimCtx<'_>, &[DimTerm]) -> Result<DimTerm, DimError>),
}
```

`matexp` and the contour/Cauchy operators use `Custom` rules
(`dim::matexp_dimension`, `dim::contour_dimension`,
`dim::cauchy_derivative_dimension`) because their rules are predicates on a
combination of arguments (`dim(t) * dim(A) = 1`, `dim(z) = 1`) rather than a
data-expressible formula; the manifest therefore records them as
domain-specific rules.

Operators that take a binder (`defint`, `sum`, `product`, `limit`, ...) use the
binding rules above: the bound symbol's dimension is looked up in the symbol
assumptions; the bounds must agree with the index dimension; `defint` requires
both limits to have the dimension of the integration variable.

A curated table in `graph::dim` supplies rules for the core operator names
(`add`, `mul`, `pow`, `abs`, `exp`, `ln`, `sin`, `sqrt`, `diff`, `integral`,
`sum`, `product`, `list`, `eq`, ...). Rule sets may override them with
`Installer::dim`. Operators without a rule are treated as unknown dimensions
(never wrong), and are reported by the coverage scanner.

### 5.1 Special cases

* `pow(x, y)`: `y` is unified with dimensionless; the result is `x^y` with `y`
  converted to an exponent expression. If `x` is `Any`, the result is `Any`.
* `product(f, k, a, b)`: `dim(f)^(b-a+1)`; the count is an exponent expression.
* `matexp(A, t)`: only the product `dim(t) * dim(A)` must be dimensionless;
  `A` may be a rate and `t` a time (`exp(A t)`), and an unknown product is
  accepted as unknown rather than guessed. The result is `1`, and a product
  with a concrete dimension is a `Mismatch`. `matfun` keeps the stricter
  single-total-degree rule above.
* Complex-analysis operators operate on a **dimensionless complex plane**:
  every operator that consumes a complex variable (`contour_integral`,
  `cauchy_integral`, `cauchy_derivative`,
  `contour_integral_residue_theorem`, `residue`, `poles`, `zeros_of`,
  `singularity`, `branch_cut`, `check_analytic`, `count_zeros_poles`,
  `winding_number`, `continue_along`, `radius_of_convergence`,
  `mobius_apply`, `is_inside_contour`) requires it to be dimensionless and
  reports `NonDimensionless` otherwise (blocking in Strict mode). A contour
  integral `∮ f dz` still carries the integrand's dimension, because `dz`
  contributes `dim(z) = 1`; `cauchy_derivative` keeps its old
  `dim(f) / dim(z)` form when `z` is still an undeclared symbol, which
  reduces to `dim(f)` once `z` is known to be dimensionless. The branch
  wrappers (`log_branch`, `sqrt_branch`, `root_branch`, `power_branch`, the
  inverse-trig/hyperbolic branches) are `Dimensionless`: every argument,
  including the branch index, is a pure number and the value is one. This
  keeps `power_branch(z, w, k)` consistent with its defining rewrite
  `exp(w log_branch(z, k))`, which a raw `Same(0)` rule contradicted. The
  two count operators no longer unify their non-scalar contour argument with
  `1`, which used to report a spurious conflict on every request, and
  `poles`, `zeros_of`, `singularity` and `branch_cut` are `Container` values:
  their list answers unify with the request instead of conflicting with it.
* `sqrt`/`cbrt`: `Root(2)`/`Root(3)`, so `sqrt(m^2) = m`.
* `abs`, `floor`, `ceil`, `round`, `re`, `im`, `conj`: `Same(0)`.
* `sign`: dimensionless.
* `atan2(y, x)`: arguments must have equal dimensions; result dimensionless.
* `exp`, `ln`, `log`, trigonometric, hyperbolic and inverse functions:
  argument dimensionless, result dimensionless (inverse trigonometric results
  are dimensionless angles).
* `list`: unifies the scalar leaves; heterogeneous lists yield `Any` with a
  diagnostic rather than a hard error (matrices are homogeneous lists of
  lists). An empty list is `Any` (its element dimension is unconstrained).
* `big_o`: `Any` (an order term follows the series it belongs to).
* Meijer G: `Custom`. If every `b_j` has the same real part, the result is
  `z^b`; otherwise `z` must be dimensionless. The rule reports a
  non-dimensionless argument when the exponents are incompatible.
* Physical constants: the coverage workstream attaches their dimensions
  (`c`, `G`, `hbar`, `k_B`, `N_A`, `e`, `m_e`, `eps_0`, `mu_0`, ...).
* Unknown functions `f(x)`: result is a fresh variable `?f(...)`; a
  `Config::dimension_of("f", "m/s^2")` declaration fixes it.

## 6. Inference

`DimInfer` is an abstract interpretation over the DAG:

```
class_dim(class):
    if memo has class: return it
    if visiting class: return Any          # cycle through an equality
    visiting.insert(class)
    result = Any
    for member in graph.members(class):
        d = derive(member)                 # children via class_dim
        result = unify(result, d)          # conflict -> diagnostic
    visiting.remove(class); memo[class] = result
    return result
```

`derive(node)`:

* literal: `Any` for zero, dimensionless otherwise;
* symbol: assumption, else a variable keyed by the symbol; undeclared symbols
  stay unknown and are reported in Advisory mode;
* other operators: the `DimRule` attribute (or the core fallback), applied to
  the dimensions of the children.

The class-level unification is what makes the analysis useful for soundness:
all members of an e-class must have the same dimension. A conflict means a
rewrite changed the dimension unless the dimension is unknown/`Any`.

Analysis is run on the input and the final answer of every `compute` in
Advisory/Strict mode. Intermediate tiers are not analysed by default. The
`graph::soundness` hook is a *numeric* refutation of rewrite rules, not a
dimensional one. Dimensional agreement of rewrites is checked in two layers:
`tests/dimensions.rs::standard_rewrites_keep_a_consistent_dimension` pins a
representative set (`x + x`/`2*x`, `x^2/x`, `expand`/`factor`, trigonometric
identities, unit conversions), and
`tests/dimensions.rs::per_rewrite_random_dimension_scan` instantiates every
declarative rewrite of `rules::standard()` under a deterministic fixed-seed set
of dimension assignments (three structured and three pseudo-random per rule)
and requires both sides to agree when both are well-typed. The scan counts its
skips separately and asserts its own non-vacuity, so a rule that stops being
exercised fails the test rather than silently passing.

Errors are collected, not thrown. In `Advisory` they become
`Answer::diagnostics`; in `Strict` the errors on the request or answer become
`ComputeError::Dimension`, except the non-conflicts: `Unannotated` (a coverage
gap), `MixedList` (a structural note) and an `UnknownUnit` inside a term (an
unknown, not a conflict). Those stay warnings even in Strict; an unknown unit
in requested dimension text is a declaration error and is blocking.

## 7. Public API

```rust
Config::new()                       // DimensionMode::Advisory
    .dimensions(DimensionMode::Strict)
    .dimension("x", "m")            // unit text, also "Length/Time", "any"
    .dimension_of("f", "m/s^2")     // unknown function value
    .define_base_dimension("Currency")
    .bypass_dimension_checks();     // explicit opt-out of hard errors

Answer {
    term, value, error, reduced, report,
    dimension: Option<Dimension>,   // inferred dimension of `term`
    diagnostics: Vec<DimError>,     // Advisory findings
}

Term::dimension() -> Dimension       // session-level assumptions
Term::dimension_report() -> DimOutcome

DimError::code() -> &'static str     // stable kind: "mismatch", ...
DimError::is_blocking() -> bool      // a warning rather than a conflict
```

`Dimension` is a printable snapshot (owned strings), independent of the graph:
`Length*Mass/Time^2`, `Length^n`, `1`, `any`, `boolean`, or a partially known
dimension with unresolved variables shown as `?x`.

`Answer::value` remains the coherent-SI numeric value; for a dimensionful
quantity it is the value after applying the unit scale. `Term::eval` keeps its
current meaning for dimensionless terms.

## 8. Coverage protocol

Requirement: **every operator registered by the standard rule sets has a
semantic dimension rule** (not `Unspecified`).

* `DimRule` is stored as an operator attribute; `graph::dim::default_rule`
  covers the core names and `Installer::dim` overrides/adds.
* The unit test `graph::dim::tests::every_standard_operator_has_a_dimension_rule`
  installs `rules::standard()`, iterates every operator and asserts that no
  operator is missing a rule or has `Unspecified`. Three companion tests keep
  the scanner honest: `the_coverage_scanner_detects_a_missing_rule` proves it
  reports an operator registered without a rule,
  `only_documented_operators_bypass_dimensions` rejects an accidental `Any`
  (the allowlist is `oo`, `infinity`, `big_o`, `hessian_definiteness`), and
  `no_operator_has_an_accidental_rule` pins the core table.
* The per-operator work is delegated to the coverage workstream, rule set by
  rule set: add `i.dim(op, rule)` in the install function (or extend the
  central table for families such as Bessel functions), add positive/negative
  tests, and remove the operator from the scanner's report. The scanner is the
  definition of done: it must be unignored and green when the coverage PRs
  land.
* `Unspecified` falls back to an unknown dimension, so an uncovered operator
  is never silently wrong; it is reported as a coverage gap.

Interim state: the scanner test is present and marked
`#[ignore = "operator coverage workstream"]` while the core lands. Removing
the `#[ignore]` is the last step of the coverage workstream; it is a checklist
item, not an architecture phase.

## 9. Numeric values

Dimensions do not change numeric evaluation of dimensionless terms. For
`quantity`, the dimension layer uses the unit resolver; `Answer::value` is the
coherent-SI value and `Answer::dimension` says what it is. Integration with
`Term::eval_in` and simulation inputs is part of the coverage workstream.

## 10. Files

* `src/graph/dim.rs` — arena, `Dim`, `DimTerm`, `Dimension`, `DimRule`,
  `DimError` (including the stable `code()`), `DimensionMode`, unit/dimension
  text parsing, core rule table, `Graph` integration.
* `src/graph/dim_expr.rs` — rational functions over atoms (symbolic
  exponents), polynomial arithmetic, display.
* `src/graph/dim_infer.rs` — `DimInfer`, class memo, unification, diagnostics.
* `src/dim/mod.rs`, `src/dim/units.rs`, `src/dim/constants.rs` — the public
  `rssn::dim` facade over the graph kernel: re-exports, the one unit table
  (`UnitDef`, `UnitRegistry`, `StandardRegistry`, exact `Scale`s and affine
  offsets) and the physical constants with their dimensions.
* `src/graph/store.rs` — `dims: RefCell<DimArena>`, symbol dimension
  assumptions, `Graph::dimension`, `Graph::assume_dimension`,
  `Graph::define_base_dimension`.
* `src/graph/rule.rs` — `Installer::dim`; `OpTable` iteration for the scanner.
* `src/rules/units.rs` — install the `quantity` unit resolver and dimension
  rules; dimension text parsing uses it.
* `src/api/mod.rs` — `Config` fields, `Answer` fields, `ComputeError::Dimension`,
  `Term::dimension`, application of assumptions in `prepare`, analysis in
  `run_job`/`finish`.
* `src/kernels/dim.rs` (optional) — numeric helpers for affine units if the
  unit layer needs them.

## 11. Testing

1. Algebra: exponent field arithmetic, symbolic exponents, normalization,
   display, unification (solvable/unsolvable/occurs check).
2. Core signatures: every core operator, positive/negative; the scaling-law
   property test (scale a symbol by `λ^d`, require the value to scale by
   `λ^dim`).
3. Units: exact conversions, affine temperatures, angles, custom bases.
4. Strict/Advisory/Off: same terms for Advisory/Off, hard errors only in
   Strict; bypass works.
5. Rewrite soundness: a representative set of standard rewrites (`x + x`,
   `x^2/x`, `expand`/`factor`, trigonometric identities, unit conversions)
   pinned in `tests/dimensions.rs`, plus
   `per_rewrite_random_dimension_scan`, a deterministic fixed-seed scan over
   every declarative rewrite of `rules::standard()` (three structured and
   three pseudo-random dimension assignments per rule; coverage counters and
   an empty refutation list are asserted).
6. Coverage scanner: zero `Unspecified` (definition of done).
7. API/FFI: `Answer.dimension`, `Term::dimension`, diagnostics, serialization.
8. Machine-readable diagnostics: `DimError::code` pins every variant; the
   `matexp` product rule and the dimensionless complex plane of the whole
   complex-analysis family (including branch wrappers) have positive and
   negative tests. The tests are wired into the build gate in §14.

## 12. Checklist

- [x] `graph::dim` core (algebra, arena, unification, rules, errors).
- [x] `Graph` integration, assumptions, base dimensions, unit resolver.
- [x] `Installer::dim`, core rule table.
- [x] `Config`/`Answer`/`Term`/`ComputeError` integration.
- [x] Tests 1–8 for the core (`graph::dim` unit tests and `tests/dimensions*.rs`).
- [x] Operator coverage workstream: every standard operator has a rule,
      scanner green (test unignored).
- [x] Docs (`ARCHITECTURE.md`, README, `rssn.h`/`rssn.hpp` status code).
- [x] Serialization of `Dimension` and the C API dimension accessor
      (`rssn_dimension`, `rssn_compute_json`); diagnostics are text (§13.8).
- [x] Machine-readable `DimError::code()` (stable snake_case kinds), pinned by
      tests; the codes are not yet surfaced through the C API (§13.8).
- [x] `matexp(A, t)` only requires the product `dim(t) * dim(A)` to be
      dimensionless, so a rate and a time are accepted.
- [x] Every complex-analysis operator that consumes a complex variable
      requires it to be dimensionless; branch wrappers are dimensionless
      (which fixed the `power_branch` inconsistency the scan found).
- [x] Per-rewrite random dimension scan over the whole standard rule set,
      deterministic and self-checking (coverage + empty refutation list).

## 13. Known limitations

These are genuine limits of the current layer, stated so that they are not
mistaken for supported behaviour. Each is pinned by a test (a diagnostic, a
hard error, or a test of the current behaviour). The list is ordered by
priority: typing gaps that can make a dimensionful value look dimensionally
correct first, then representation and binding ergonomics.

Closed in the current audit: diagnostics have a stable machine-readable kind
(`DimError::code`, §7); `matexp(A, t)` accepts a dimensionful `A`/`t` whose
product is dimensionless; every complex-analysis operator that consumes a
complex variable requires it to be dimensionless; branch wrappers are typed as
the dimensionless-plane functions they are (the scan found and this closed the
`power_branch` inconsistency); and the highest-priority item, the missing
per-rewrite random dimensional scan, is now a deterministic test
(`tests/dimensions.rs::per_rewrite_random_dimension_scan`, §6).

1. **A structured result is a typed container, not a tuple type (priority:
   medium).** `lu`, `qr`, `svd`, `eigenvects`, `jordan`, `jordan_blocks`,
   `spectral` and `ldl` are typed `DimRule::Container` (`non-metric`): the
   operator knows the parts have different dimensions, and an
   Advisory/Strict computation reports no `MixedList` for a structured
   request. The dimensions of the parts are documented and tested
   (`L`/`P`/`Q`/eigenvectors dimensionless, `U`/`R`/`S`/eigenvalues the entry
   dimension, projectors dimensionless, `d` the pivots), but a `Term` built
   directly from the expanded literal (for example `Term::dimension` on the
   stored answer outside `compute`) still analyses it as nested `list`s and
   can report `MixedList`; there is no re-usable `(A, B, C)` product type.
   `vandermonde` and `companion` stay `NonMetric` because their entries mix
   powers of the input.
2. **Affine points are boundary-checked, not fully typed (priority: medium).**
   The kernel and the dimension rule reject scaling, powers and point+point,
   and report `DimError::Affine` (blocking in Strict). But a point and a
   difference are both `Temperature` in the type layer: `quantity(t, K)` is
   read as a difference when added to a `degC` point, `convert` still shifts
   through kelvin, and `point - point` is left unreduced rather than typed as
   a difference. A full point/vector distinction would need a separate affine
   dimension kind.
3. **`matfun(f, x, A)` is homogeneous only for a single total degree
   (priority: low).** If `f` is a polynomial in `x` whose terms all have the
   same degree `k` and whose coefficients do not mention `x`, `f(A)` has
   entries `dim(f) * entry^k` (`x^2` gives `A^2`; `x^3 + x` is not
   homogeneous). A non-homogeneous `f` (including `exp`, `sin` and other
   non-polynomial functions) requires a dimensionless matrix and is reported
   with `NonDimensionless` otherwise. Symbolic exponents (`x^n`) are not
   analysed.
4. **A refused term-level unit is `any` plus a non-blocking diagnostic
   (priority: low).** `quantity(1, parsec)` or `quantity(1, degC*m)` has no
   resolvable unit; the resolver returns `None`, the dimension is `any`, the
   `units` kernel leaves the request unreduced, and `DimError::UnknownUnit`
   reports the unit text and the term. It is an *unknown*, not a conflict, so
   Strict mode answers with the diagnostic rather than failing. A unit in
   requested dimension text (`Config::dimension("x", "parsec")`) is a
   declaration error: `UnknownUnit` with no context is blocking.
5. **Strict mode checks the input as written (priority: low).** A term can be
   ill-typed although it simplifies to a well-typed one: `ln(exp(x))` with
   `x` a length is rejected, because `exp(x)` requires a dimensionless
   argument, even though the term equals `x`. This is inherent in typing the
   input and the answer independently; Advisory mode still reports the same
   finding without failing the computation.
6. **A class containing the literal zero is polymorphic (priority: low).**
   This is what makes kernels that simplify intermediate expressions to zero
   safe (a zero-equivalent expression such as `x - x` has a length dimension
   syntactically, but the class is `any`). It also means Strict mode does not
   catch an ill-typed subterm that is multiplied by zero.
7. **Unit text wins over base-dimension text (priority: low).**
   `dimension_from_text` first asks the unit resolver, so `"T"` means tesla,
   `"C"` coulomb, `"F"` farad. A custom base registered under one of those
   names is only reachable with the explicit `base:Name` syntax; bare text is
   never silently reinterpreted. `base:` also works for the SI bases
   (`base:Length`), with the registered name matched exactly (the aliases are
   not applied). A trailing `/` is tolerated (`"Length/"` is `Length`).
8. **Diagnostics are not yet surfaced as codes through the C API (priority:
   low).** `DimError::code()` gives every finding a stable snake_case kind
   (`mismatch`, `non_dimensionless`, `symbolic_constraint`, `unknown_unit`,
   `affine`, `mixed_list`, `unannotated`, `parse`), documented and pinned by
   tests, so a Rust caller no longer has to parse the message. The C API
   still exposes `Answer.dimension` and the diagnostics as strings in
   `rssn_compute_json`; a binding that needs to branch on the kind has to
   parse the text until a later round serializes the code.

## 14. Landing status

The core landed before the six algorithm branches were merged; the coverage
workstream ran on the merged tree. Test commands for the core:

```
cargo test --lib graph::dim
cargo test --test dimensions
cargo test --test dimensions_core
cargo test --test dimensions_special
cargo test --test dimensions_ode_pde
cargo test --test dimensions_physics
```

`graph::dim::tests::every_standard_operator_has_a_dimension_rule` installs
`rules::standard()` and asserts that no operator has `Unspecified`; it is no
longer ignored. The four coverage workstreams (core, special/math, ODE/PDE and
physics/misc) annotated the remaining operators; the complex-analysis and
verify operators were completed centrally in the core signature table.

The limitations audit closed the machine-readable codes, the `matexp` product
rule and the complex-variable/plane checks. They are pinned by
`dim::tests::dim_error_codes_pin_every_variant`,
`dim::tests::matexp_requires_only_a_dimensionless_product`,
`dim::tests::branch_wrappers_are_dimensionless_plane_functions`,
`dim::tests::contour_operators_require_a_dimensionless_complex_variable`,
`tests/dimensions.rs::diagnostics_carry_stable_codes`,
`tests/dimensions_special.rs::matrix_functions_carry_the_degree_of_an_entry`,
`tests/dimensions_special.rs::complex_analysis_operators_check_their_complex_variable`
and
`tests/dimensions_special.rs::contour_operators_require_a_dimensionless_complex_variable`.

The rewrite-soundness audit added
`tests/dimensions.rs::per_rewrite_random_dimension_scan`: it instantiates every
declarative rewrite of `rules::standard()` (261 rules in the current tree)
under three structured (all dimensionless, all length, all time) and three
pseudo-random dimension assignments per rule, generated from a fixed per-rule
seed. A run currently reports 261/261 rewrites exercised, 882 dimensioned
instantiations checked, 683 skipped because the *pattern* is dimensionally
ill-typed under the assignment (outside the rewrite's domain), 1 skipped
because the *right* side is domain-restricted, 0 guarded out, and an empty
refutation list; the test asserts the counters, so a rule that stops being
exercised fails rather than passing vacuously. The scan found one genuine
inconsistency: branch wrappers were typed `Same(0)`, which made
`power_branch(z, w, k)` disagree with its defining rewrite
`exp(w log_branch(z, k))`. They are now `Dimensionless` — every argument,
including the branch index, is a pure number — and every complex-analysis
operator that consumes a complex variable checks that it is dimensionless.
The same pass fixed the spurious `non-metric vs 1` conflict on
`poles`/`zeros_of`/`singularity`/`branch_cut` (now `Container`, so their list
answers unify) and the spurious `NonDimensionless` on the non-scalar contour
argument of `count_zeros_poles`/`winding_number`.
