//! Exact finite-difference checks for the public Gosper operator.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Number;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::rules;

fn closed(
    graph: &mut Graph,
    engine: &Engine,
    root: NodeId,
) -> Option<NodeId> {
    engine.run(
        graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
    Extractor::new(graph, &[root], &ClosedForm).build(graph, root)
}

#[test]
fn gosper_matches_exact_polynomial_and_factorial_differences() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    let index = graph.sym("k");
    for term in ["k^7", "k*factorial(k)"] {
        let query = graph.parse(&format!("antidifference({term},k)")).unwrap();
        let answer =
            closed(&mut graph, &engine, query).unwrap_or_else(|| panic!("unresolved {term}"));
        for point in [1_i64, 4, 7] {
            let current = graph.int(point);
            let next = graph.int(point + 1);
            let previous = graph.substitute(answer, index, current);
            let following = graph.substitute(answer, index, next);
            let minus_one = graph.int(-1);
            let negative = graph.node(core::MUL, &[minus_one, previous]);
            let diff = graph.node(core::ADD, &[following, negative]);
            let value = closed(&mut graph, &engine, diff).unwrap();
            let expected = if term == "k^7" {
                point.pow(7)
            } else {
                point * (1..=point).product::<i64>()
            };
            assert_eq!(
                graph.number_of(value),
                Some(&Number::from(expected)),
                "{term} at {point}: {}",
                graph.display(value)
            );
        }
    }
}

#[test]
fn gosper_leaves_an_unpaired_rational_pole_unresolved() {
    // The two simple poles are in distinct integer shift orbits. A rational
    // difference has zero sum of the simple residues along each orbit.
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    let root = graph.parse("antidifference(1/(k^2+k+1),k)").unwrap();
    let operator = graph.op(root);
    let result = closed(&mut graph, &engine, root);
    assert!(result.is_none_or(|value| graph.op(value) == operator));
}
