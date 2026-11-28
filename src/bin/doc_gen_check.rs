//! Verifies that every operator in the kernel manifest has non-empty
//! documentation, and reports coverage grouped by `rule_set`.
//!
//! ```text
//! cargo run --bin doc_gen_check                  # human-readable report (fails if any op lacks doc)
//! cargo run --bin doc_gen_check -- --summary     # per-rule-set summary only
//! cargo run --bin doc_gen_check -- --json        # structured JSON report
//! cargo run --bin doc_gen_check -- --allow-missing # exit 0 even when docs are missing
//! ```
//!
//! The exit status is `0` when every operator has non-empty documentation (or
//! when `--allow-missing` is passed), `1` when any operator lacks
//! documentation under strict checking, and `2` when the check itself could
//! not run.

use std::collections::BTreeMap;
use std::process::ExitCode;

use rssn::graph::Engine;
use rssn::graph::Graph;
use rssn::kernel::KernelManifest;
use serde::Serialize;

/// Entry point for `doc_gen_check`.
fn main() -> ExitCode {
    let mut json = false;
    let mut summary = false;
    let mut allow_missing = false;

    for argument in std::env::args().skip(1) {
        match argument.as_str() {
            | "--json" => json = true,
            | "--summary" => summary = true,
            | "--allow-missing" => allow_missing = true,
            | "--help" | "-h" => {
                print_help();
                return ExitCode::SUCCESS;
            },
            | other => {
                eprintln!("doc-gen-check: unknown argument `{other}`");
                print_help();
                return ExitCode::from(2);
            },
        }
    }

    let sets = rssn::rules::standard();
    let mut graph = Graph::new();
    if Engine::install(&mut graph, &sets).is_err() {
        eprintln!("doc-gen-check: the standard rule sets do not install");
        return ExitCode::from(2);
    }
    let manifest = KernelManifest::build(&graph, &sets, &[]);

    let mut total_covered: usize = 0;
    let mut total_missing: usize = 0;
    let mut per_set: BTreeMap<&str, RuleSetStats> = BTreeMap::new();

    for op in &manifest.ops {
        let entry = per_set.entry(op.rule_set.as_str()).or_default();
        if op.doc.trim().is_empty() {
            total_missing = total_missing.saturating_add(1);
            entry.missing = entry.missing.saturating_add(1);
            entry.missing_ops.push(op.name.clone());
        } else {
            total_covered = total_covered.saturating_add(1);
            entry.covered = entry.covered.saturating_add(1);
        }
    }

    let total_ops = manifest.ops.len();

    if json {
        let rule_sets = per_set
            .iter()
            .map(|(name, stats)| {
                let set_total = stats.covered.saturating_add(stats.missing);
                JsonRuleSetReport {
                    rule_set: (*name).to_string(),
                    covered: stats.covered,
                    missing: stats.missing,
                    coverage_percent: percent(stats.covered, set_total),
                    missing_ops: stats.missing_ops.clone(),
                }
            })
            .collect();

        let report = JsonReport {
            rssn_version: manifest.rssn_version.clone(),
            api_version: manifest.api_version.to_string(),
            total_ops,
            covered: total_covered,
            missing: total_missing,
            coverage_percent: percent(total_covered, total_ops),
            rule_sets,
        };

        match serde_json::to_string_pretty(&report) {
            | Ok(rendered) => println!("{rendered}"),
            | Err(error) => {
                eprintln!("doc-gen-check: failed to serialize JSON report: {error}");
                return ExitCode::from(2);
            },
        }
    } else {
        println!(
            "rssn {} · api {}",
            manifest.rssn_version, manifest.api_version
        );
        println!(
            "{total_ops} operators in {} rule set(s): {total_covered} documented ({:.1}%), {total_missing} missing",
            per_set.len(),
            percent(total_covered, total_ops),
        );

        println!("\nper rule set:");
        for (name, stats) in &per_set {
            let set_total = stats.covered.saturating_add(stats.missing);
            println!(
                "  {name:<24} {:>5} documented {:>5} missing ({:.1}%)",
                stats.covered,
                stats.missing,
                percent(stats.covered, set_total),
            );
        }

        if !summary && total_missing > 0 {
            println!("\nmissing documentation:");
            for (name, stats) in &per_set {
                if stats.missing_ops.is_empty() {
                    continue;
                }
                println!("  [{name}]");
                for op_name in &stats.missing_ops {
                    println!("    {op_name}");
                }
            }
        }
    }

    if total_missing == 0 {
        ExitCode::SUCCESS
    } else if allow_missing {
        eprintln!(
            "doc-gen-check: {total_missing} operator(s) lack documentation (allowed via --allow-missing)"
        );
        ExitCode::SUCCESS
    } else {
        eprintln!("doc-gen-check: {total_missing} operator(s) lack documentation");
        ExitCode::FAILURE
    }
}

/// Prints CLI usage information.
fn print_help() {
    println!(
        "doc-gen-check — verify that every manifest operator has documentation\n\
         \n\
         USAGE:\n    doc_gen_check [--json] [--summary] [--allow-missing]\n\
         \n\
         OPTIONS:\n\
         \x20   --json           print a structured JSON coverage report\n\
         \x20   --summary        print only overall and per-rule-set counts\n\
         \x20   --allow-missing  exit 0 even when operators lack documentation\n\
         \x20   --help           this message"
    );
}

/// `part` as a percentage of `whole`, or `0.0` when `whole` is zero.
fn percent(
    part: usize,
    whole: usize,
) -> f64 {
    if whole == 0 {
        return 0.0;
    }
    let part = u32::try_from(part).unwrap_or(u32::MAX);
    let whole = u32::try_from(whole).unwrap_or(u32::MAX);
    f64::from(part) * 100.0 / f64::from(whole)
}

/// Per-rule-set documentation coverage counters.
#[derive(Default)]
struct RuleSetStats {
    /// Operators with non-empty `doc`.
    covered: usize,
    /// Operators with empty or whitespace-only `doc`.
    missing: usize,
    /// Names of operators missing `doc`.
    missing_ops: Vec<String>,
}

/// Structured JSON output for `--json`.
#[derive(Serialize)]
struct JsonReport {
    /// `rssn` crate version from the manifest.
    rssn_version: String,
    /// Kernel API version from the manifest.
    api_version: String,
    /// Total number of operators in `manifest.ops`.
    total_ops: usize,
    /// Number of operators with non-empty `doc`.
    covered: usize,
    /// Number of operators with empty or whitespace-only `doc`.
    missing: usize,
    /// Overall documentation coverage percentage.
    coverage_percent: f64,
    /// Coverage breakdown grouped by `rule_set`.
    rule_sets: Vec<JsonRuleSetReport>,
}

/// Per-rule-set entry in [`JsonReport`].
#[derive(Serialize)]
struct JsonRuleSetReport {
    /// Name of the rule set.
    rule_set: String,
    /// Number of documented operators in this rule set.
    covered: usize,
    /// Number of undocumented operators in this rule set.
    missing: usize,
    /// Coverage percentage for this rule set.
    coverage_percent: f64,
    /// Operator names in this rule set that lack documentation.
    missing_ops: Vec<String>,
}
