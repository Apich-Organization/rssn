//! Opt-in per-rule cost accounting.
//!
//! Setting `RSSN_PROFILE_RULES=1` makes the scheduler and the windows record,
//! for every rewrite searched and every kernel run, how long it took and how
//! often it produced something. [`report`] renders the table, most expensive
//! first. With the variable unset every hook is one relaxed load of a cached
//! flag and nothing is recorded.
//!
//! The table is process-wide and keyed by the rule's address, so it is only
//! meaningful within one engine program; it exists to tell a performance
//! investigation which rules to look at, not to be shipped in output.

use super::op::OpTable;
use super::pattern::Pat;
use std::fmt::Write as _;
use std::sync::Mutex;
use std::sync::OnceLock;
use std::time::Duration;
use std::time::Instant;

/// Whether `RSSN_PROFILE_RULES` is set (read once).
#[must_use]
pub fn enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("RSSN_PROFILE_RULES").is_some_and(|v| v != "0"))
}

#[derive(Default)]
struct Entry {
    label: String,
    calls: u64,
    hits: u64,
    nanos: u128,
}

type Table = std::collections::HashMap<(&'static str, usize), Entry>;

fn table() -> &'static Mutex<Table> {
    static TABLE: OnceLock<Mutex<Table>> = OnceLock::new();
    TABLE.get_or_init(|| Mutex::new(Table::new()))
}

/// A started measurement; `None` when profiling is off.
#[must_use]
pub fn start() -> Option<Instant> {
    enabled().then(Instant::now)
}

/// Adds one measurement for the rule identified by `(kind, key)`; `label`
/// is called only the first time the rule is seen.
pub fn record(
    started: Option<Instant>,
    kind: &'static str,
    key: usize,
    hits: usize,
    label: impl FnOnce() -> String,
) {
    let Some(started) = started else {
        return;
    };
    let elapsed: Duration = started.elapsed();
    let Ok(mut table) = table().lock() else {
        return;
    };
    let entry = table.entry((kind, key)).or_insert_with(|| {
        Entry {
            label: label(),
            ..Entry::default()
        }
    });
    entry.calls = entry.calls.saturating_add(1);
    entry.hits = entry.hits.saturating_add(hits as u64);
    entry.nanos = entry.nanos.saturating_add(elapsed.as_nanos());
}

/// A readable rendering of a pattern, with operator names.
#[must_use]
pub fn describe(
    pat: &Pat,
    ops: &OpTable,
) -> String {
    let mut out = String::new();
    write_pat(pat, ops, &mut out);
    out
}

fn write_pat(
    pat: &Pat,
    ops: &OpTable,
    out: &mut String,
) {
    match pat {
        | Pat::Var(v) => {
            let _ = write!(out, "?{v}");
        },
        | Pat::Num(n) => {
            let _ = write!(out, "{n}");
        },
        | Pat::Exact(node) => {
            let _ = write!(out, "#{}", node.raw());
        },
        | Pat::Node(op, children) => {
            out.push_str(&ops.get(*op).name);
            out.push('(');
            for (i, child) in children.iter().enumerate() {
                if i > 0 {
                    out.push_str(", ");
                }
                write_pat(child, ops, out);
            }
            out.push(')');
        },
    }
}

/// The table so far, most expensive first, and resets it.
#[must_use]
pub fn report() -> String {
    let Ok(mut table) = table().lock() else {
        return String::new();
    };
    let mut rows: Vec<((&'static str, usize), Entry)> = table.drain().collect();
    rows.sort_by(|a, b| b.1.nanos.cmp(&a.1.nanos));
    let total: u128 = rows.iter().map(|(_, e)| e.nanos).sum();
    let mut out = format!(
        "{:>10} {:>6} {:>8} {:>8}  rule\n",
        "ms", "%", "calls", "hits"
    );
    for ((kind, _), entry) in rows.iter().take(60) {
        let ms = entry.nanos as f64 / 1e6;
        let share = if total == 0 {
            0.0
        } else {
            entry.nanos as f64 * 100.0 / total as f64
        };
        let _ = writeln!(
            out,
            "{ms:>10.3} {share:>5.1}% {:>8} {:>8}  [{kind}] {}",
            entry.calls, entry.hits, entry.label
        );
    }
    out
}
