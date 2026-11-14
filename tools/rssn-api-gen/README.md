# RSSN source binding generation

The source generator builds without the RSSN library, even if the checked-in
fragments are missing or invalid. From the RSSN repository root:

```sh
cargo run -p rssn-api-gen -- --codecs-check
cargo run -p rssn-api-gen -- --bindings-check
cargo run -p rssn-api-gen -- --structured-check
cargo run -p rssn-api-gen -- --contracts-check
cargo run -p rssn-api-gen -- --codecs-write
cargo run -p rssn-api-gen -- --bindings-write
cargo run -p rssn-api-gen -- --structured-write
cargo run -p rssn-api-gen -- --contracts-write
cargo test -p rssn-api-model -p rssn-api-gen
```

Use `--bindings-report` or `--structured-report` to review the generated
output. An optional final argument specifies an alternate RSSN source root.
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
kernel call. The new raw-number policy deliberately requires dimensionless inputs.
Existing quantity-aware contracts continue to apply to their original exports.
Tuple results are records (`item0`, `item1`, ...); optional results are tagged
`None`/`Some`; fixed-width point vectors use rank-two tensors. Structured result
codecs still come from the separately reviewed `structured-policy.tsv`.

The generator builds without linking RSSN. Callback algorithms have dedicated
serial adapters with explicit evaluation budgets in `numerical_callbacks`.
They preserve host errors and never pretend a generic Rust closure is a plain
value codec. Public source APIs and dynamic manifest coverage are audited
separately in Aevia's `rssn-api-audit` tool.
