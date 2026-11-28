# RSSN source binding generation

The source generator builds without the RSSN library, even if the checked-in
fragments are missing or invalid. From the RSSN repository root:

```sh
bazel test //:api_generator_test //:generated_sources_check_test
bazel run //:refresh_contracts
```

Use `bazel test //:contracts_test` to check the generated fragments together
with the linked operator manifest. The underlying generator also accepts
`--bindings-report` or `--structured-report` to print a single fragment.
An optional final argument specifies an alternate RSSN source root.
Policies describe reviewed semantic contracts; supported Rust function and
record shapes come from syntax and are checked by generated Rust bindings.
The Aevia-side `rssn-api-audit` separately compares the live manifest with the
public Rust inventory; generation itself does not depend on that manifest.

`contracts.tsv` is the reviewed source for every legacy numerical dimensional
policy. It names argument and nested record bindings, equality constraints,
and scalar or record result formulas. The generated contract is interpreted
by both manifest signature construction and runtime dimension checking; adding
an enum case without a corresponding policy fails compilation. Structured
result policies can be referenced by name to avoid repeating field formulas.

`codec-policy.tsv` registers concrete kernel functions through composable codecs.
Rows contain only source path, ABI name, namespace and reviewed dimension policy.
`--codecs-report`, `--codecs-write` and `--codecs-check` derive the full call,
borrowed-argument conversion, arity and result type from the public Rust function.
Generated code shares the same codecs for signatures and values; Rust checks each
kernel call. Every policy declares its physical relationships; dimensionless
inputs are required only when mathematically appropriate. Optional numeric
constraints apply to Some payloads and skip None; mixed result records retain
separate field units. `extent_pow(D,arg,axis,multiplier)` declares a result power
determined from a validated input tensor extent (determinants and densities).
Its symbolic manifest result uses the checked dynamic route, and runtime
decoration resolves the exact integer exponent.
Tuple results are records (`item0`, `item1`, ...); optional results are tagged
`None`/`Some`; fixed-width point vectors use rank-two tensors. Structured result
codecs still come from the separately reviewed `structured-policy.tsv`.

The generator builds without linking RSSN. Callback algorithms have dedicated
serial adapters with explicit evaluation budgets in `numerical_callbacks`.
They preserve host errors and never pretend a generic Rust closure is a plain
value codec. Public source APIs and dynamic manifest coverage are audited
separately in Aevia's `rssn-api-audit` tool.

Reviewed fixed SI laws can use `base("Length")` and the other six canonical
base names from `src/dim/units.rs`. Unknown names and computed strings fail
generation; physical constants never become implicit dimensionless inputs.
Concrete public `sim` scalar functions undergo the same source checks as
`kernels` and `backend` functions. Unsafe, generic, conditional, private and
mutable-reference functions require a separate adapter.

Sampled functional-analysis data use `{item0: coordinates, item1: ordinates}`.
This maps `Vec<(f64,f64)>` to independently dimensioned vectors, and bases to
matching matrices. Norms include the coordinate integration measure, projection
preserves the first function's ordinate units, and L2 normalization returns
inverse square-root coordinate units. Axis grids must be finite and increasing.

`num_stats_percentile` explicitly copies input before calling the native mutable
slice API; read-only robust quantile also supports all Hyndman-Fan definitions.
Scope decisions: `classical::damped_harmonic_oscillator`
is not added: phase is not consistently represented in its critical/overdamped
branches. Stateful simulation models, save-to-file scenarios, backend internals
and generic callbacks are not implied to be covered by the scalar policy set.
