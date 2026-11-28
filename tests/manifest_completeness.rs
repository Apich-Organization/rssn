//! Every operator must be *completely* described by the manifest.
//!
//! `signature-audit` answers "is every operator typed or explained", and
//! `doc-gen-check` answers "does every operator have documentation". This test
//! answers the remaining question a generated host asks before it trusts the
//! manifest: is every field of every operator present and internally
//! consistent, or did some registration quietly leave a hole?
//!
//! A hole is not a cosmetic problem. `category` is how a host groups the
//! operator table, `doc` is what it shows a user, `keys` is what it caches on,
//! and `overload_sources` is how it tells a declared contract from one the
//! engine was merely seen to satisfy. An operator missing one of them is an
//! operator a generated interface has to special-case.
//!
//! The check is deliberately structural. It does not re-derive the signatures
//! (the audit and the overlay test do that); it states the invariants that hold
//! for *every* operator whatever its rule set, so that a new family cannot be
//! added without meeting them.

use std::collections::BTreeMap;
use std::sync::Arc;

use rssn::kernel::KernelManifest;

fn standard() -> Arc<KernelManifest> {
    rssn::Session::new().manifest()
}

#[test]
fn every_operator_names_itself_its_family_and_its_category() {
    let manifest = standard();
    let mut problems = Vec::new();
    for op in &manifest.ops {
        if op.name.trim().is_empty() {
            problems.push(String::from("an operator has an empty name"));
        }
        if op.rule_set.trim().is_empty() {
            problems.push(format!("{}: empty rule_set", op.name));
        }
        if op.category.trim().is_empty() {
            problems.push(format!("{}: empty category", op.name));
        }
        if op.doc.trim().is_empty() {
            problems.push(format!("{}: empty doc", op.name));
        }
        if op.cost == 0 {
            problems.push(format!("{}: cost is zero", op.name));
        }
    }
    assert!(
        problems.is_empty(),
        "{} incomplete operator(s):\n  {}",
        problems.len(),
        problems.join("\n  ")
    );
}

#[test]
fn every_operator_is_typed_or_explained() {
    let manifest = standard();
    let mut problems = Vec::new();
    for op in &manifest.ops {
        if op.overloads.is_empty() {
            // The only acceptable reason for having no signature.
            if op.dynamic_reason.is_none() {
                problems.push(format!("{}: no overload and no dynamic_reason", op.name));
            }
        } else if op.dynamic_reason.is_some() {
            // A typed operator has nothing to explain.
            problems.push(format!("{}: has overloads and a dynamic_reason", op.name));
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

#[test]
fn the_key_overload_and_source_lists_agree() {
    let manifest = standard();
    let mut problems = Vec::new();
    for op in &manifest.ops {
        if op.keys.len() != op.overloads.len() {
            problems.push(format!(
                "{}: {} key(s) for {} overload(s)",
                op.name,
                op.keys.len(),
                op.overloads.len()
            ));
        }
        if op.overload_sources.len() != op.overloads.len() {
            problems.push(format!(
                "{}: {} source(s) for {} overload(s)",
                op.name,
                op.overload_sources.len(),
                op.overloads.len()
            ));
        }
        let distinct: BTreeMap<_, _> = op.keys.iter().map(|key| (key, ())).collect();
        if distinct.len() != op.keys.len() {
            problems.push(format!("{}: duplicate keys", op.name));
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}

#[test]
fn the_manifest_is_sorted_and_unique_by_name() {
    let manifest = standard();
    for pair in manifest.ops.windows(2) {
        let (left, right) = (&pair[0], &pair[1]);
        assert!(
            left.name < right.name,
            "operators out of order: {} then {}",
            left.name,
            right.name
        );
    }
}

#[test]
fn every_documented_path_belongs_to_the_operator_it_documents() {
    // `doc_path` is the function whose rustdoc completed `doc`. It is set only
    // when the harvest found the function, so it must always name something,
    // and it must agree with the `rust_path` the registration stated up to the
    // generic arguments and the re-exporting module (both are resolved by
    // `kernel::rustdoc`).
    let manifest = standard();
    let mut problems = Vec::new();
    for op in &manifest.ops {
        let Some(doc_path) = &op.doc_path else {
            continue;
        };
        if doc_path.trim().is_empty() {
            problems.push(format!("{}: empty doc_path", op.name));
            continue;
        }
        if let Some(rust_path) = &op.rust_path {
            // The last segment is the function's own name, which carries no
            // generic arguments; a resolved path may differ in the module it
            // names (a re-export) but never in the function.
            let stated = rust_path.rsplit("::").next().unwrap_or(rust_path);
            let found = doc_path.rsplit("::").next().unwrap_or(doc_path);
            if stated != found {
                problems.push(format!(
                    "{}: doc_path {doc_path} does not name the same function as rust_path {rust_path}",
                    op.name
                ));
            }
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
