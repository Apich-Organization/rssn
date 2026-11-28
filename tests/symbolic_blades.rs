//! Exact products and signature boundaries for symbolic Clifford algebra.

use std::collections::BTreeMap;

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::geometric_algebra::geometric_algebra;

type Signature = (u32, u32, u32);

fn reduce(expression: &str) -> Option<String> {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[geometric_algebra()]).unwrap();
    let root = graph.parse(expression).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{expression}");
    Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .map(|value| graph.display(value).to_string())
}

fn mv(
    (positive, negative, null): Signature,
    terms: &[(u32, i64)],
) -> String {
    let terms = terms
        .iter()
        .map(|&(basis, c)| format!("list({basis}, {c})"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("mv(list({positive}, {negative}, {null}), list({terms}))")
}

/// Reference via sorted basis words and explicit contractions.
fn reference(
    (positive, negative, null): Signature,
    left: u32,
    right: u32,
) -> (u32, i64) {
    let dimension = positive + negative + null;
    let word = |mask| (0..dimension).filter(move |&k| mask & (1_u32 << k) != 0);
    let mut word = word(left).chain(word(right)).collect::<Vec<_>>();
    let mut coefficient = 1;
    for right in 1..word.len() {
        let mut i = right;
        while i > 0 && word[i - 1] > word[i] {
            word.swap(i - 1, i);
            coefficient = -coefficient;
            i -= 1;
        }
    }
    let mut result = 0;
    while let Some(index) = word.pop() {
        if word.last() == Some(&index) {
            word.pop();
            if index >= positive + negative {
                return (0, 0);
            }
            if index >= positive {
                coefficient = -coefficient;
            }
        } else {
            result |= 1 << index;
        }
    }
    (result, coefficient)
}

#[test]
fn symbolic_blade_signature_overflow_is_rejected_before_products() {
    for bad in [
        "mv_basis(list(4294967295,2,0),1)",
        "mv_scalar(list(4294967295,1,0),3)",
        "mv_pseudoscalar(list(4294967295,0,2))",
        "mv_scalar(list(0,0,13),5)",
        "mv_scalar(list(0,4294967296,0),5)",
        "mv_blade(list(4,3,5),4096,3)",
        "mv_basis(list(4,3,5),13)",
    ] {
        assert_eq!(reduce(bad), None, "{bad}");
    }
    assert_eq!(
        reduce("mv_basis(list(4,3,5),12)"),
        Some(mv((4, 3, 5), &[(2048, 1)]))
    );
    assert_eq!(
        reduce("ga_gp(mv_scalar(list(0,0,0),2),mv_scalar(list(0,0,0),3))"),
        Some(mv((0, 0, 0), &[(0, 6)]))
    );
}

#[test]
fn symbolic_blade_products_follow_the_metric_and_orientation() {
    for (sig, left, right) in [
        ((2, 1, 1), 9, 5),
        ((2, 1, 1), 5, 9),
        ((2, 1, 1), 6, 4),
        ((2, 1, 1), 9, 8),
        ((4, 3, 5), 2048, 2048),
        ((7, 5, 0), 4095, 4095),
        ((0, 0, 12), 3, 2048),
        ((0, 12, 0), 1025, 2048),
    ] {
        let (basis, c) = reference(sig, left, right);
        let lhs = mv(sig, &[(left, 2)]);
        let rhs = mv(sig, &[(right, 3)]);
        let expected = if c == 0 {
            mv(sig, &[])
        } else {
            mv(sig, &[(basis, 6 * c)])
        };
        assert_eq!(reduce(&format!("ga_gp({lhs},{rhs})")), Some(expected));
    }
}

#[test]
fn symbolic_blade_dense_mixed_product_matches_exact_word_contractions() {
    let sig = (2, 2, 1);
    let a = [(0, 1), (3, -2), (6, 4), (10, 6), (11, 7), (17, -9), (28, 6)];
    let b = [(1, 1), (4, -11), (9, 2), (13, -3), (21, 5), (30, 3)];
    let mut expected = BTreeMap::new();
    for &(left, lc) in &a {
        for &(right, rc) in &b {
            let (basis, sign) = reference(sig, left, right);
            *expected.entry(basis).or_insert(0) += lc * rc * sign;
        }
    }
    let expected = expected
        .into_iter()
        .filter(|&(_, c)| c != 0)
        .collect::<Vec<_>>();
    let lhs = mv(sig, &a);
    let rhs = mv(sig, &b);
    assert_eq!(
        reduce(&format!("ga_gp({lhs},{rhs})")),
        Some(mv(sig, &expected))
    );
}

#[test]
fn symbolic_blade_products_keep_symbolic_coefficients_in_either_orientation() {
    let sig = "list(2,1,1)";
    let x = format!("mv_blade({sig},9,u)");
    let y = format!("mv_blade({sig},5,v)");
    assert_eq!(
        reduce(&format!("ga_gp({x},{y})")).as_deref(),
        Some("mv(list(2, 1, 1), list(list(12, u*v)))")
    );
    assert_eq!(
        reduce(&format!("ga_gp({y},{x})")).as_deref(),
        Some("mv(list(2, 1, 1), list(list(12, -u*v)))")
    );
}
