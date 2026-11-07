# RSSN source binding generation

The source generator builds without the RSSN library, even if the checked-in
fragments are missing or invalid. From the RSSN repository root:

```sh
cargo run -p rssn-api-gen -- --bindings-check
cargo run -p rssn-api-gen -- --structured-check
cargo run -p rssn-api-gen -- --contracts-check
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
