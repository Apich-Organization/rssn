//! Exercises the exact symbolic inversion count and permutation sign rules.

use num_bigint::BigInt;
use rssn::graph::Budget;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Number;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::rules::combinatorics::combinatorics;

fn run(
    graph: &mut Graph,
    engine: &Engine,
    root: NodeId,
) -> Option<BigInt> {
    engine.run(
        graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
    match graph.number_of(root)? {
        | Number::Int(n) => Some(n.clone()),
        | _ => panic!("permutation statistics must return an exact integer"),
    }
}

fn inversion_reference(perm: &[usize]) -> u64 {
    let mut count = 0;
    for j in 0..perm.len() {
        for i in 0..j {
            if perm[i] > perm[j] {
                count += 1;
            }
        }
    }
    count
}

fn check(
    perm: &[usize],
    expected: u64,
) {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[combinatorics()]).unwrap();
    let items = perm
        .iter()
        .map(|&r| graph.int(i64::try_from(r + 1).unwrap()))
        .collect::<Vec<_>>();
    let list = graph.node(core::LIST, &items);
    for (name, result) in [
        ("inversions", BigInt::from(expected)),
        (
            "permutation_sign",
            BigInt::from(if expected.is_multiple_of(2) {
                1
            } else {
                -1
            }),
        ),
    ] {
        let op = graph.ops().lookup(name).unwrap();
        let root = graph.node(op, &[list]);
        assert_eq!(
            run(&mut graph, &engine, root),
            Some(result),
            "{name}; n={}",
            perm.len()
        );
    }
}

#[test]
fn combinatorics_inversions_matches_quadratic_cases_across_rank_boundaries() {
    let mut seed = 5_934_678_125_u64;
    for n in [0, 1, 2, 7, 63, 64, 65, 127, 128, 129, 255, 256, 257] {
        let mut perm = (0..n).collect::<Vec<_>>();
        for p in (1..n).rev() {
            seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let j = (seed >> 32) as usize % (p + 1);
            perm.swap(p, j);
        }
        check(&perm, inversion_reference(&perm));
    }
}

#[test]
fn combinatorics_inversions_respects_order_and_parity_for_wide_runs() {
    let n = 4097;
    let mut perm = (0..n).collect::<Vec<_>>();
    check(&perm, 0);
    perm.swap(2047, 2048);
    check(&perm, 1);
    perm.sort_unstable_by(|a, b| b.cmp(a));
    check(&perm, 8_390_656);
    perm.swap(2047, 2048);
    check(&perm, 8_390_655);
}

#[test]
fn combinatorics_inversions_keeps_invalid_rank_lists_symbolic() {
    for op in ["inversions", "permutation_sign"] {
        for items in [
            "list(2, 2)",
            "list(0, 1)",
            "list(1, 3)",
            "list(-1, 2)",
            "list(1, 3/2)",
            "list(1, x)",
            "3",
        ] {
            let mut graph = Graph::new();
            let engine = Engine::install(&mut graph, &[combinatorics()]).unwrap();
            let root = graph.parse(&format!("{op}({items})")).unwrap();
            assert!(run(&mut graph, &engine, root).is_none(), "{op}({items})");
        }
    }
}
