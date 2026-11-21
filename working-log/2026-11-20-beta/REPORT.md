# Beta verification report

The final gate ran on Windows against the shared source archive SHA-256 `1b390e200c8d3de37c2f33e9f5c4914acc1bd315a7f06b9621c1ed5dea473597`. Both hosts verified their staged sources before and after testing. The RSSN catalog included 1,148 files. The Aevia catalog included 960 files and one intentional deletion. The [source catalog](evidence/tested-snapshot/source-catalog.json) records the tested paths and hashes.

| Windows gate | Passed | Failed | Ignored |
| --- | ---: | ---: | ---: |
| RSSN workspace, all features | 6,115 | 0 | 3 |
| Aevia language workspace, all features | 3,216 | 0 | 5 |
| Aevia platform workspace, all features | 282 | 0 | 0 |
| Aevia registry documentation, no default features | 32 | 0 | 0 |
| **Total** | **9,645** | **0** | **8** |

Each gate exited 0. The [RSSN terminal receipt](evidence/windows/rssn/terminal-receipt.json) and [Aevia terminal receipt](evidence/windows/aevia/terminal-receipt.json) include log digests. The [target table](evidence/windows/targets.tsv) maps 494 direct integration roots to source hashes and log lines across 548 target result records. Service fixtures were verified separately on Debian. Some Windows service tests return early when their fixtures are unavailable. [Earlier fixture-backed Debian run](evidence/debian/platform-with-fixtures-earlier.log)

The first Windows platform run exposed an index replacement race. The old executable reproduced it eight times. The client now closes the staging handle before publication and retries transient Windows replacement errors using the same staging file. The existing destination stays readable during retry. The isolated Windows package run passed 70 tests, and all 24 concurrent-writer repeats and six held-reader repeats passed. Both regression tests passed again in the final full platform suite. [Baseline](evidence/regressions/windows-cache-baseline.json) · [Repair verification](evidence/regressions/windows-cache-repair.json)

## Architecture and coverage

- RSSN detects available CPU workers and honors the active Rayon pool. Coarse independent work uses the shared scheduling policy. Bounded scratch leases and reusable plans reduce allocations. Current changes include column partitioning for wide one-row CPU matrix products and independent-row elimination for exact algebraic fields.
- Aevia uses its worker budget and the same RSSN runtime for native calls. Homogeneous parallel results use contiguous buffers. Pinned native calls borrow large contiguous Float inputs to avoid building an intermediate neutral tensor, while preserving the common witness and signature checks.
- Aevia generates bindings and references from the RSSN manifest. This candidate regenerated 88 bilingual files covering 5,353 operators. All 44 English outputs matched the checked-in files. The 25 runnable manual examples passed separately through AST and IR evaluation. The local audit checked 149 Markdown documents, 1,030 references, and 129 Cargo selectors. [Generated reference hashes](evidence/docs/bilingual-reference-hashes.json) · [Manual execution log](evidence/docs/manuals-both-engines.log) · [Local link audit](evidence/docs/markdown-links.log)
- Debian built the current default-feature sources with configured Clippy and passed 64 RSSN and 83 cross-repository Aevia integration tests. Aevia also passed 125 type-checker units, 272 type-checker integration tests, and 75 resolver units. [RSSN receipt](evidence/debian/rssn-default-receipt.json) · [Aevia integration log](evidence/debian/aevia-runtime-integrations.log)

| Project | Module groups | Targeted evidence | Partial evidence |
| --- | ---: | ---: | ---: |
| RSSN | 345 | 141 | 204 |
| Aevia | 116 | 82 | 34 |
| **Total** | **461** | **223** | **238** |

The [coverage ledger](coverage/modules.tsv) records every group. Partial rows may have source review, a smoke test, or execution of selected algorithms. Per-operator boundary and performance work remains in the [next steps](HANDOFF.md#do-next). The session allowed three simultaneous subagents, so the work used successive module batches.

## Performance

These are paired Debian measurements of specific inputs.

| Workload | Baseline | Candidate |
| --- | ---: | ---: |
| Aevia type check of a pre-parsed 4,096-element tensor literal | 1,246.991 ms | 48.428 ms |
| Aevia pinned measured norm, 100,000 Float values | 32.044 ms | 22.016 ms |
| Aevia Fortran NPY conversion, 512 × 512 × 4, two workers | 16.145 ms | 6.753 ms |
| RSSN CPU product, 1 × 64 × 65,537, two workers | 2.373 ms | 1.571 ms |

[Type checking](evidence/performance/aevia-type-aliases.tsv) · [Native arguments](evidence/performance/aevia-native-arguments.tsv) · [NPY](evidence/performance/aevia-npy-fortran.log) · [CPU matrix product](evidence/performance/rssn-cpu-wide-matmul.tsv)

## Evidence after formatting

Nightly rustfmt changed 10 RSSN and five Aevia Rust files. Cleanup also removed one trailing blank line from an Aevia editor fixture. Every other committed source hash matches the tested snapshot. The full gate predates this cleanup. [Rustfmt patch](evidence/after-format/formatting.patch) · [All cleanup deltas](evidence/after-format/committed-source-delta.tsv) · [Format checks](evidence/after-format/format-check.log) · [Cargo metadata checks](evidence/after-format/metadata-check.log)

From the `evidence` directory, run `sha256sum -c SHA256SUMS` to check the archived receipts. Raw Windows test logs keep their original line endings. The coverage ledgers retain historical artifact names from the original archive so their scope can be traced. The evidence folder contains the final receipts and the supporting runs cited above.
