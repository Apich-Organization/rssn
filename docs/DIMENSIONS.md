# Dimensions and units in RSSN

[中文](DIMENSIONS.zh.md) · [The embedding kernel](KERNEL.md)

RSSN infers the dimensions of symbolic computations and carries dimensioned quantities through the kernel ABI. This document describes the symbolic rules, the quantity representation and the evidence a host can inspect.

## 1. Modes

`DimensionMode::Advisory` is the default. `compute` analyses the request and its answer, returns the inferred answer dimension and includes dimension findings in `Answer::diagnostics`. `DimensionMode::Strict` rejects a blocking finding with `ComputeError::Dimension`. `Config::bypass_dimension_checks` keeps diagnostics and allows the computation to continue after a conflict. Use `DimensionMode::Off` to skip the analysis.

Strict mode checks the expression supplied by the caller. If `x` is a length, `ln(exp(x))` is rejected because `exp` requires a dimensionless argument, even though the composition could be simplified to `x`.

## 2. Bases and exponents

The default bases are `Length`, `Mass`, `Time`, `Current`, `Temperature`, `Amount` and `Luminosity`. Additional names can be registered for a graph. A dimension is a product of bases raised to exponents. Exponents support exact rationals, variables and rational functions over symbolic atoms. The name of a symbolic exponent and an uninterpreted expression such as `sin(a)` can each serve as an atom. Equality of rational expressions is checked by cross-multiplication.

Dimension variables record undeclared symbols and function results. For unification, the engine divides the two dimensions and isolates each variable with a nonzero rational exponent. It substitutes the resulting product transitively. Residual concrete exponents give a mismatch, and residual symbolic exponent equations give a constraint. A variable with a symbolic exponent stays unresolved when the equation does not determine it.

`1` denotes the dimensionless empty product. `NonMetric` denotes predicates, labels and containers without a single scalar dimension. `Any` is the polymorphic dimension used for zero, infinity, empty lists and explicit bypass rules. A term class containing zero may use `Any` even if a different member has a concrete dimension.

## 3. Units and conversion

A unit carries a dimension, scale and optional affine offset. Coherent SI values are `value × scale + offset`. Exact rational or π-valued scales preserve computations such as `sin(quantity(30, deg)) = sin(π/6)`. Named derived units such as `N` and `Pa` are presentations of products of the same bases.

Angles, ratios and information units have dimension `1`. The unit metadata still identifies angles and ratios. `bit`, `B` and `KiB` have scale factors `1`, `8` and `8192` relative to bits. Unit names are resolved before base names in a dimension declaration. `base:T` names a registered base called `T`, while bare `T` is tesla. `base:Length` uses the same explicit form for an SI base.

The `units` rules provide the resolver used by the symbolic graph. A host may supply a full `UnitValue` without registering it. Materialization retains its dimensions, exact scale, offset and display metadata. A registered name is reused only when the entire supplied definition matches. Unit products can combine registered and anonymous atoms. Custom bases survive coherent conversion and rational powers. Materialization rejects invalid scale denominators, incompatible affine metadata, symbolic dimensions on unit definitions and π-valued offsets.

## 4. Affine quantities

`degC` and `degF` identify points. A point may be converted through its coherent origin or shifted by a difference in the same dimension. For example, `quantity(20, degC) + quantity(2, K)` adds a temperature difference to a point. The core dimension of both quantities is `Temperature`. The unit rules reject point scaling, point exponentiation and point addition. A point inside a compound unit fails unit resolution. Point subtraction remains a symbolic request. Code that needs an explicit point-to-point difference can first convert both points to coherent Kelvin values and subtract the readings.

## 5. Operator rules

Each standard operator registers a `DimRule`. A rule may preserve an argument dimension, multiply or divide dimensions, constrain equal dimensions, take a root, impose dimensionlessness or delegate to a domain-specific Rust function. Rules that take an index or integration variable also use its declared dimension. Operator registration and the manifest carry these rules. The current Rust variants are defined in [`src/graph/dim.rs`](../src/graph/dim.rs).

| Expression | Dimension rule |
|---|---|
| `pow(x, y)` | `y` is dimensionless. The result carries the dimension of `x` raised to `y`. |
| `diff(f, x)` · `integral(f, x)` | Divide or multiply the dimension of `f` by the dimension of `x`. Definite bounds use the dimension of `x`. |
| `product(f, k, a, b)` | Raise the dimension of `f` to `b - a + 1`. Bounds and index agree. |
| `atan2(y, x)` | Match dimensions of both arguments. The angle is dimensionless. |
| `sde_wiener(seed, m, steps, h)` | Use the square root of `h`'s dimension. |
| `sde_levy_areas(seed, h, dw, terms)` | Use `h`'s dimension. The iterated integral has the dimension of a pair of increments. |
| `matexp(A, t)` | Require `dim(A) × dim(t) = 1`. This accepts a rate matrix times a duration. |
| `matfun(f, x, A)` | A homogeneous polynomial of degree `k` contributes `entry^k` times the coefficient dimension. Other functions require dimensionless entries. |
| Complex-plane operators | The complex variable and branch indexes are dimensionless. A contour integral keeps its integrand dimension. |

## 6. Inference

The analysis walks the term DAG and memoizes dimension evidence per equivalence class. Constants and user assumptions seed the analysis. Each node derives a candidate from its children's dimensions and its operator rule. Class members are unified. During a cycle, the in-progress class provides `Any` until its traversal completes. A computation analyses the request and final answer. Registered assertions such as the dimension of `f(x)` help when the structure alone leaves a variable unresolved.

The engine reports unresolved dimensions, symbolic constraints and conflicts separately. `Unannotated`, `MixedList` and a term-level `UnknownUnit` are advisory findings even in Strict mode. Unknown unit text in a dimension declaration is a blocking declaration error. Downstream hosts that require a concrete result should check the outcome and diagnostics before extracting a quantity.

## 7. Embedding API

```rust
let config = Config::new()
    .dimensions(DimensionMode::Strict)
    .dimension("x", "m")
    .dimension_of("f", "m/s^2")
    .define_base_dimension("Currency");
```

`Config::dimension_term` accepts a typed `DimExpr`. `Term::dimension` returns the rendered dimension under the current assumptions. `Term::dimension_report` also exposes its analysis outcome. `Answer::dimension` stores the final dimension, and `Answer::value` is the numeric value in coherent SI. A `Dimension` can outlive the graph. Descriptions look like `Length*Mass/Time^2`, `Length^n`, `1`, `any`, `boolean` or `?x` for an unresolved variable. `DimError::code()` and `DimError::is_blocking()` let a host classify findings. The C API exposes the result and per-thread diagnostic accessors. [C and embedding contracts](KERNEL.md)

## 8. Coverage contract

`Installer::dim` attaches a dimension rule to a registered operator. The standard-operator scan rejects missing rules and undocumented bypasses. A deterministic rewrite suite exercises standard declarative rules against structured and pseudo-random dimension assignments. It counts skipped assignments, requires exercised rules and compares the dimensions when both sides are well-typed. Per-domain tests pin matrix, complex-plane, stochastic and structured-value behavior. [Architecture and verification](../ARCHITECTURE.md)

## 9. Dynamic values

The kernel preserves explicit dimensions in `Value::Dyn`. An existing quantity inside a dynamic value must agree with its dimension witness. Nonempty homogeneous numeric tensors retain their element dimension on graph leaves and on readback. Mixed quantities remain structured symbolic terms. Hosts should retain typed tensors when an empty tensor carries trailing extents or a nontrivial element dimension. The nested-list term representation cannot reconstruct metadata without leaves.

Core arithmetic and comparisons offer symbolic signatures as well as literal numeric signatures. A symbolic call can accept quantities and term handles. The selected signature describes the admissible argument kinds, and dimension validation determines whether the result is a concrete quantity.

## 10. Systems of units

`KernelManifest::unit_systems` exposes four systems. Each fixes the value of selected physical quantities to one and chooses generators for the remaining dimensions.

| System | Unity quantities | Remaining generators |
|---|---|---|
| `natural` | `c`, `hbar`, `k_B`, `epsilon_0` in the Heaviside–Lorentz convention | Energy, Amount, Luminosity |
| `planck` | `c`, `hbar`, `G`, `k_B`, `4 pi epsilon_0` | Amount, Luminosity |
| `atomic` | `hbar`, `m_e`, `e`, `4 pi epsilon_0` | Temperature, Amount, Luminosity |
| `geometrized` | `c`, `G` | Length, Current, Temperature, Amount, Luminosity |

`UnitSystem::reduce` solves the decomposition over rational numbers. The unity dimensions and the generators must span the seven SI bases independently. The remaining dimension is still an SI dimension. Divide the coherent value by the corresponding powers of the unity quantities to obtain the reduced value. `Reduction::exact_scale` gives a symbolic scale for integer powers. `factor_f64`, `to_system` and `from_system` also handle fractional powers. For example, an electron mass in the natural system has energy dimension and a value of approximately `0.51099895 MeV`.

## 11. Low-level numerical kernels

Algorithms that operate on ordinary `f64` storage use the caller's coherent values. Typed numerical bindings supply dimensional evidence around those kernels. See [numerical callbacks](NUMERICAL_CALLBACKS.md) for dimensioned inputs, tolerances and results.

## 12. Run a quantity example

Run `cargo run --locked --example units_and_dimensions` from the repository root. The example shows unit conversion and the diagnostic difference between Advisory and Strict modes. [Example source](../examples/units_and_dimensions.rs)

## 13. Structured results and domain boundaries

Factorizations such as `lu`, `qr` and `svd` return structures with different dimensions in different fields. In `svd`, singular values carry the matrix-entry dimension and singular vectors are dimensionless. A structured operator uses `DimRule::Container`; independently constructing a heterogeneous nested list can report `MixedList` because its entries do not share a scalar dimension. A host that needs field-level dimensions should use the typed kernel contract for the factorization.

A term with an unknown unit remains symbolic and reports `UnknownUnit`. Provide a registered unit or pass a complete `UnitValue` for an anonymous unit. For matrix functions with an unresolved symbolic polynomial degree, use dimensionless entries or first substitute a degree that determines a homogeneous polynomial. For validation before multiplication by a polymorphic zero, analyse the subexpression that must satisfy the dimension contract on its own.
