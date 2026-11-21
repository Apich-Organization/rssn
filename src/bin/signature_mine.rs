//! Mines `signatures.witnessed.json` from sort traces.
//!
//! ```text
//! # 1. record: run any workload with tracing on (the test suite is the corpus)
//! RSSN_SORT_TRACE=/tmp/sorts.tsv cargo test --lib
//! # 2. mine: cluster, replay every witness, write the table
//! cargo run --profile tools --bin signature-mine -- /tmp/sorts.tsv
//! cargo run --profile tools --bin signature-mine -- --check   # replay the shipped table
//! ```
//!
//! The first step costs nothing to maintain: every operator that has a test
//! is exercised with real arguments, and the engine logs the sort of each
//! request next to the sort of its answer. The second step keeps only what a
//! fresh engine reproduces, so the table can be regenerated at any time and a
//! stale entry cannot survive. See `rssn::kernel::sort` for the design.
//!
//! Exit status: `0` on success, `1` when `--check` finds a witness that no
//! longer replays, `2` when the tool could not run.

use std::collections::BTreeMap;
use std::process::ExitCode;
use std::time::Duration;

use rssn::graph::Arity;
use rssn::graph::Budget;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Saturate;
use rssn::kernel::Signature;
use rssn::kernel::sort;
use rssn::kernel::sort::Candidate;
use rssn::kernel::sort::Fact;
use rssn::kernel::sort::Replay;
use rssn::kernel::sort::Sort;

const OUTPUT: &str = "signatures.witnessed.json";

fn main() -> ExitCode {
    let mut check = false;
    let mut output = OUTPUT.to_owned();
    let mut inputs = Vec::new();
    let mut args = std::env::args().skip(1);
    while let Some(argument) = args.next() {
        match argument.as_str() {
            | "--check" => check = true,
            | "--output" | "-o" => {
                let Some(path) = args.next() else {
                    eprintln!("signature-mine: --output needs a path");
                    return ExitCode::from(2);
                };
                output = path;
            },
            | "--help" | "-h" => {
                println!(
                    "usage: signature-mine [--output FILE] TRACE...\n       signature-mine --check"
                );
                return ExitCode::SUCCESS;
            },
            | other => inputs.push(other.to_owned()),
        }
    }

    let sets = rssn::rules::standard();
    let mut template = Graph::new();
    let engine = match Engine::install(&mut template, &sets) {
        | Ok(engine) => engine,
        | Err(error) => {
            eprintln!("signature-mine: the standard rule sets do not install: {error}");
            return ExitCode::from(2);
        },
    };
    // A witness is a call a test already answered; a generous but bounded
    // budget keeps one pathological replay from stalling the run.
    let budget = Budget {
        time_limit: Some(Duration::from_secs(20)),
        ..Budget::default()
    };
    // The budget is checked between engine iterations; a single kernel can
    // still take arbitrarily long on arguments nobody tested (probing shifts
    // the numbers of a witness). Each run therefore happens on its own thread
    // with a deadline. A run that misses it is reported and treated as
    // unanswered; its thread is abandoned, which is acceptable for an offline
    // tool and impossible to do better without cancellation in the engine.
    let engine = std::sync::Arc::new(engine);
    let run = move |graph: &mut Graph, root: NodeId| {
        let mut owned = graph.clone();
        let (sender, receiver) = std::sync::mpsc::channel();
        let (engine, budget) = (std::sync::Arc::clone(&engine), budget.clone());
        std::thread::spawn(move || {
            engine.run(&mut owned, &[root], &Env::symbolic(), &Saturate, &budget);
            let _ = sender.send(owned);
        });
        match receiver.recv_timeout(Duration::from_secs(45)) {
            | Ok(done) => *graph = done,
            | Err(_) => {
                eprintln!(
                    "signature-mine: no answer within 45 s, skipped: {}",
                    graph.display(root)
                );
            },
        }
    };

    if check {
        return check_shipped(&template, &run);
    }
    if inputs.is_empty() {
        eprintln!("signature-mine: give at least one trace file (see --help)");
        return ExitCode::from(2);
    }

    let mut facts: Vec<Fact> = Vec::new();
    let mut malformed = 0_usize;
    for path in &inputs {
        let text = match std::fs::read_to_string(path) {
            | Ok(text) => text,
            | Err(error) => {
                eprintln!("signature-mine: cannot read {path}: {error}");
                return ExitCode::from(2);
            },
        };
        for line in text.lines() {
            match Fact::parse(line) {
                | Some(fact) => facts.push(fact),
                | None => malformed += 1,
            }
        }
    }
    // A fact about an operator this build does not have is not evidence.
    facts.retain(|fact| template.ops().lookup(&fact.op).is_some());
    facts.sort_by(|a, b| (&a.op, &a.params, &a.witness).cmp(&(&b.op, &b.params, &b.witness)));
    facts.dedup();

    let existing = |name: &str, arity: usize, variadic: bool| -> Vec<Signature> {
        template
            .ops()
            .lookup_all(name)
            .iter()
            .map(|&id| template.ops().get(id))
            .filter(|desc| {
                match desc.arity {
                    | Arity::Fixed(n) => !variadic && n as usize == arity,
                    | Arity::Variadic => variadic,
                }
            })
            .flat_map(|desc| {
                // Only what the registration itself states: the shipped
                // witnessed overloads are being replaced, not extended.
                if !desc.overloads.is_empty() {
                    desc.overloads.clone()
                } else {
                    rssn::kernel::infer::derive(desc).unwrap_or_default()
                }
            })
            .collect()
    };
    let binder = |name: &str| -> Option<usize> {
        template
            .ops()
            .lookup(name)
            .and_then(|id| template.ops().get(id).binder)
            .map(|b| usize::from(b.var))
    };
    // Operators with a declared or overlay signature are already typed by
    // hand; mining them would only add noise next to a deliberate choice.
    facts.retain(|fact| {
        template.ops().lookup_all(&fact.op).iter().all(|&id| {
            let desc = template.ops().get(id);
            desc.overloads.is_empty()
                && rssn::kernel::infer::signatures(desc).1 != rssn::kernel::SignatureSource::Overlay
        })
    });

    let (table, stats) = sort::mine(&template, &facts, &existing, &binder, &run);
    let json = match serde_json::to_string_pretty(&table) {
        | Ok(json) => json,
        | Err(error) => {
            eprintln!("signature-mine: cannot encode the table: {error}");
            return ExitCode::from(2);
        },
    };
    if let Err(error) = std::fs::write(&output, format!("{json}\n")) {
        eprintln!("signature-mine: cannot write {output}: {error}");
        return ExitCode::from(2);
    }
    println!(
        "signature-mine: {} facts ({} malformed lines), {} candidates",
        stats.facts, malformed, stats.candidates
    );
    println!(
        "  kept {} signatures for {} operators ({} refine the result of a registered overload)",
        stats.kept, stats.ops, stats.refinements
    );
    println!(
        "  {} calls were not answered by the full engine, {} candidates would be ambiguous or \
         redundant, {} results were widened by probing",
        stats.unconfirmed, stats.incoherent, stats.widened
    );
    println!("  wrote {output}");
    ExitCode::SUCCESS
}

/// Replays every shipped witness; reports the ones that no longer hold.
fn check_shipped(
    template: &Graph,
    run: &dyn Fn(&mut Graph, NodeId),
) -> ExitCode {
    let table = sort::shipped();
    let mut failures: BTreeMap<String, String> = BTreeMap::new();
    let mut replayed = 0_usize;
    for (name, entries) in &table.ops {
        for entry in entries {
            let Some(result) = Sort::parse(&entry.result) else {
                failures.insert(name.clone(), format!("unreadable sort `{}`", entry.result));
                continue;
            };
            let params: Option<Vec<Sort>> = entry.params.iter().map(|p| Sort::parse(p)).collect();
            let Some(params) = params else {
                failures.insert(name.clone(), "unreadable parameter sort".to_owned());
                continue;
            };
            let candidate = Candidate {
                op: name.clone(),
                variadic: entry.signature.variadic.is_some(),
                params,
                result,
                facts: Vec::new(),
            };
            for witness in &entry.witnesses {
                let fact = Fact {
                    op: name.clone(),
                    variadic: candidate.variadic,
                    params: candidate.params.clone(),
                    result: candidate.result.clone(),
                    witness: witness.input.clone(),
                };
                replayed += 1;
                match sort::replay(template, &candidate, &fact, run) {
                    | Replay::Confirmed(_) => {},
                    | Replay::Refuted(why) | Replay::Unusable(why) => {
                        failures.insert(name.clone(), format!("{}: {why}", witness.input));
                    },
                }
            }
        }
    }
    if failures.is_empty() {
        println!(
            "signature-mine: {replayed} witnesses of {} operators replay",
            table.ops.len()
        );
        return ExitCode::SUCCESS;
    }
    for (name, why) in &failures {
        eprintln!("signature-mine: {name}: {why}");
    }
    eprintln!(
        "signature-mine: {} of {} operators have a stale witness; re-run the miner",
        failures.len(),
        table.ops.len()
    );
    ExitCode::from(1)
}
