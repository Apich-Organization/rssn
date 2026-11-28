//! Explicit truth-table and large exact model-count references.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::logic::logic;

fn model_count(formula: &str) -> Option<String> {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[logic()]).unwrap();
    let request = graph.parse(&format!("count_models({formula})")).unwrap();
    engine.run(
        &mut graph,
        &[request],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{formula}");
    Extractor::new(&graph, &[request], &ClosedForm)
        .build(&mut graph, request)
        .map(|count| graph.display(count).to_string())
}

#[test]
fn symbolic_bdd_count_matches_an_independent_truth_table() {
    let expected = (0_u8..16)
        .filter(|&assignment| {
            let [a, b, c, d] = [1, 2, 4, 8].map(|bit| assignment & bit != 0);
            (a || b) && (b != c) && (d == (!a || c))
        })
        .count();
    assert_eq!(
        model_count("and(or(a,b),xor(b,c),iff(d,implies(a,c)))"),
        Some(expected.to_string())
    );
    assert_eq!(
        model_count("and(iff(w,implies(z,x)),xor(y,x),or(y,z))"),
        Some(expected.to_string())
    );
}

#[test]
fn symbolic_bdd_count_retains_exact_counts_through_sixty_four_atoms() {
    let atoms = (0..64).map(|i| format!("v{i:02}")).collect::<Vec<_>>();
    let joined = atoms.join(",");
    assert_eq!(
        model_count(&format!("or({joined})")),
        Some(((1_u128 << 64) - 1).to_string())
    );
    assert_eq!(
        model_count(&format!("xor({joined})")),
        Some((1_u128 << 63).to_string())
    );
    assert_eq!(model_count(&format!("and({joined})")).as_deref(), Some("1"));
    assert_eq!(model_count(&format!("or({joined},v64)")), None);
}

#[test]
fn symbolic_bdd_count_keeps_empty_and_single_atom_conventions() {
    assert_eq!(model_count("true").as_deref(), Some("1"));
    assert_eq!(model_count("false").as_deref(), Some("0"));
    assert_eq!(model_count("a").as_deref(), Some("1"));
    assert_eq!(model_count("not(a)").as_deref(), Some("1"));
}
