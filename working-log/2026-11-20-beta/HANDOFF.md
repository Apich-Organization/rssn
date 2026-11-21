# RSSN and Aevia beta handoff

The beta round closed on November 20, 2026 at 20:30 Eastern. Both Windows hosts finished the full test plan on the same frozen source. The cleanup then ran nightly formatting and committed the changes in module order. Start with the source revisions below and use the [verification report](REPORT.md) for the scope of each result.

| Checkout | Source revision |
| --- | --- |
| RSSN, `v0.3` | `9952d45b9a9c6b3c6c128fed0f2d3cf7d56cbae6` |
| Aevia, `v0.1` | `174dc51bccbca4ecaa2ead55bc73686f30e80dea` |

RSSN's subsequent working-log commit contains only this bundle. Keep the two checkouts adjacent because Aevia resolves RSSN as a local path dependency. The [commit ledger](commits.tsv) records the module boundaries. The Windows results apply to the frozen candidate, before 15 Rust formatting changes and removal of one extra blank line in an editor fixture. The [committed source delta](evidence/after-format/committed-source-delta.tsv) lists every difference. Cargo nightly format checks, locked offline metadata checks, and the manual structure check passed after cleanup.

## Do next

1. **Verify the release artifacts.** Run the existing release dry-run for the paired source revisions. Inspect the optimized/LTO C archives, confirm exports from the packaged headers, and build and run the C and C++ clients against the extracted static and dynamic libraries. Check the archive hashes and signing or attestation steps before publishing. The earlier Windows C/C++ smoke covered development-profile artifacts. [RSSN build and release procedure](../../docs/BUILDING.md)
2. **Deepen the weakest RSSN modules first.** Start with measured grid and continuum ABI, stateful particles and relativity, then exact polynomial and quantifier elimination. For each leaf, identify missing named operator calls, pair analytic or exact controls with rejection and rollback cases, and record manifests and performance evidence alongside any extension. Use `native_execution`, `native_continuum`, `native_grids`, `native_legacy_physics`, `native_particles`, `real_qe`, `f4_f5`, and `symbolic_resultant` as starting Cargo test selectors. [Coverage ledger](coverage/modules.tsv)
3. **Improve Aevia's dynamic RSSN bindings.** Replace dynamic signatures with exact type, shape, and dimension contracts where the source declaration supports them. Review Rust-only interfaces for a portable representation before exposing them. Keep the RSSN manifest authoritative. From Aevia, run `cargo run --locked -p rssn-api-audit -- --verify`, then use the built CLI with `tools/check_stdlib_docs.py --aevia` and `tools/check_manuals.py --aevia` to validate generated references and executable examples.
4. **Measure before adding parallelism.** For another parallel candidate, compare serial and 1/2/4-worker behavior with identical data and conversion code, preserve accumulation order and the first error, and measure elapsed time plus allocation traffic. Use the active RSSN pool and existing scratch or plan leases. Deep Aevia tensor checking still visits child shapes repeatedly during unification. The WGPU cold-validation investigation also needs an isolated cross-backend comparison before any dependency change. [Measured baselines](REPORT.md#performance)
5. **Keep service fixtures explicit.** The final Windows platform run covers the cache regression. PostgreSQL and GPG fixtures were unavailable on those hosts. Rerun affected platform service paths with the required fixtures on Debian when changing them. [Earlier fixture-backed Debian run](evidence/debian/platform-with-fixtures-earlier.log)

Use scoped tests during development and schedule the full release gate once the next candidate is frozen. Reuse the warmed Cargo targets and schedule one build at a time per target directory. Before removing compiled artifacts to recover space, identify idle legacy caches and keep the newest verified artifacts and their receipts.
