//! Reports how much of the operator table is typed, and fails if any
//! operator is neither typed nor explained.
//!
//! ```text
//! cargo run --bin signature-audit            # human-readable report
//! cargo run --bin signature-audit -- --json  # the whole manifest as JSON
//! cargo run --bin signature-audit -- --quiet # only failures, exit code only
//! ```
//!
//! # Why this is a binary and not only a test
//!
//! The rule the audit enforces is a property of a *build*, and the build a
//! downstream host cares about is the one it linked, possibly with plugins it
//! registered. Running it by hand against that build is the only way to know
//! whether the manifest a host will read is complete.
//!
//! The exit status is `0` when every operator is typed or explained, `1` when
//! some are not, and `2` when the audit itself could not run. A CI job can
//! therefore use it as a gate without parsing the output.

use std::collections::BTreeMap;
use std::process::ExitCode;

use rssn::graph::Engine;
use rssn::graph::Graph;
use rssn::kernel::KernelManifest;
use rssn::kernel::SignatureSource;

fn main() -> ExitCode {
    let mut json = false;
    let mut quiet = false;
    let mut verbose = false;
    for argument in std::env::args().skip(1) {
        match argument.as_str() {
            | "--json" => json = true,
            | "--quiet" | "-q" => quiet = true,
            | "--verbose" | "-v" => verbose = true,
            | "--help" | "-h" => {
                print_help();
                return ExitCode::SUCCESS;
            },
            | other => {
                eprintln!("signature-audit: unknown argument `{other}`");
                print_help();
                return ExitCode::from(2);
            },
        }
    }

    let Some(manifest) = build() else {
        eprintln!("signature-audit: the standard rule sets do not install");
        return ExitCode::from(2);
    };

    if json {
        println!("{}", manifest.to_json());
        return ExitCode::SUCCESS;
    }

    let audit = Audit::of(&manifest);
    if !quiet {
        audit.report(verbose);
    }
    if audit.unexplained.is_empty() {
        ExitCode::SUCCESS
    } else {
        // Printed even under `--quiet`: a gate that fails silently is useless.
        eprintln!(
            "signature-audit: {} operator(s) are neither typed nor explained",
            audit.unexplained.len()
        );
        for name in &audit.unexplained {
            eprintln!("  {name}");
        }
        eprintln!(
            "signature-audit: declare a `Signature`, or a `DynamicReason` saying why not"
        );
        ExitCode::FAILURE
    }
}

fn print_help() {
    println!(
        "signature-audit — check that every operator is typed or explained\n\
         \n\
         USAGE:\n    signature-audit [--json] [--quiet] [--verbose]\n\
         \n\
         OPTIONS:\n\
         \x20   --json      print the whole manifest as JSON and stop\n\
         \x20   --quiet     print nothing unless the audit fails\n\
         \x20   --verbose   also list the operators that were inferred\n\
         \x20   --help      this message"
    );
}

/// Installs the standard rule sets and builds the manifest.
fn build() -> Option<KernelManifest> {
    let sets = rssn::rules::standard();
    let mut graph = Graph::new();
    Engine::install(&mut graph, &sets).ok()?;
    Some(KernelManifest::build(&graph, &sets, &[]))
}

/// What the audit found.
struct Audit<'m> {
    manifest: &'m KernelManifest,
    declared: usize,
    inferred: usize,
    explained: usize,
    /// Names of operators that are neither typed nor explained.
    unexplained: Vec<String>,
    /// Operators whose signature was derived rather than declared.
    inferred_names: Vec<String>,
    /// Per rule set: (typed, explained, unexplained).
    per_set: BTreeMap<&'m str, (usize, usize, usize)>,
    /// Every distinct `dynamic_reason`, with how many operators use it.
    reasons: BTreeMap<&'m str, usize>,
}

impl<'m> Audit<'m> {
    fn of(manifest: &'m KernelManifest) -> Self {
        let mut audit = Self {
            manifest,
            declared: 0,
            inferred: 0,
            explained: 0,
            unexplained: Vec::new(),
            inferred_names: Vec::new(),
            per_set: BTreeMap::new(),
            reasons: BTreeMap::new(),
        };
        for op in &manifest.ops {
            let entry = audit.per_set.entry(op.rule_set.as_str()).or_default();
            if !op.overloads.is_empty() {
                match op.signature_source {
                    | SignatureSource::Declared | SignatureSource::Overlay => audit.declared += 1,
                    | SignatureSource::Derived => {
                        audit.inferred += 1;
                        audit.inferred_names.push(op.name.clone());
                    },
                    // An operator with overloads always has a source; a
                    // `None` here would be a manifest bug, and counting it as
                    // typed is the harmless reading.
                    | SignatureSource::None => audit.declared += 1,
                }
                entry.0 += 1;
            } else if let Some(reason) = &op.dynamic_reason {
                audit.explained += 1;
                entry.1 += 1;
                *audit.reasons.entry(reason.as_str()).or_default() += 1;
            } else {
                audit.unexplained.push(op.name.clone());
                entry.2 += 1;
            }
        }
        audit
    }

    fn report(
        &self,
        verbose: bool,
    ) {
        let total = self.manifest.ops.len();
        let typed = self.declared + self.inferred;
        // The identity is the 32-byte digest; `hash_u64` is the first eight
        // bytes read little-endian, which is *not* the first eight bytes in
        // hex order. Printing both, labelled, stops a host pinning the wrong
        // one.
        println!(
            "rssn {} · api {}",
            self.manifest.rssn_version, self.manifest.api_version
        );
        println!(
            "manifest hash {} (u64 {:016x})",
            self.manifest.hash_hex(),
            self.manifest.hash_u64()
        );
        println!(
            "{total} operators in {} rule set(s): {typed} typed ({:.1}%), {} explained, {} unexplained",
            self.manifest.rule_sets.len(),
            percent(typed, total),
            self.explained,
            self.unexplained.len()
        );
        println!(
            "  declared: {}   derived: {}   overlay: {}",
            self.declared,
            self.inferred,
            self.manifest
                .ops
                .iter()
                .filter(|op| op.signature_source == SignatureSource::Overlay)
                .count()
        );
        println!(
            "  units: {}   base dimensions: {}   constants: {} ({} exact)   plugins: {}",
            self.manifest.units.len(),
            self.manifest.base_dimensions.len(),
            self.manifest.constants.len(),
            self.manifest.constants.iter().filter(|c| c.exact).count(),
            self.manifest.plugins.len()
        );
        println!(
            "  plugin kinds: {}   value kinds: {} ({})",
            self.manifest.features.plugin_kinds.join(", "),
            self.manifest.value_kinds.len(),
            self.manifest
                .value_kinds
                .iter()
                .filter(|kind| kind.in_term_graph)
                .count()
        );

        println!("\nper rule set:");
        for (name, (typed, explained, unexplained)) in &self.per_set {
            println!(
                "  {name:<24} {typed:>5} typed {explained:>5} explained {unexplained:>4} unexplained"
            );
        }

        if !self.reasons.is_empty() {
            println!("\nwhy operators are untyped:");
            for (reason, count) in &self.reasons {
                println!("  {count:>5}  {reason}");
            }
        }

        if verbose && !self.inferred_names.is_empty() {
            println!("\nderived (no declaration, no overlay):");
            for name in &self.inferred_names {
                println!("  {name}");
            }
        }
    }
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
