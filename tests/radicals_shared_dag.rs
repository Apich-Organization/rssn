//! Regression tests for applying radical requests to hash-consed expressions.

use rssn::graph::Arity;
use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::OpDescriptor;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::rules::radicals::radicals;

fn operation(
    graph: &mut Graph,
    engine: &Engine,
    name: &str,
    term: NodeId,
) -> NodeId {
    let op = graph.ops().lookup(name).unwrap();
    let root = graph.node(op, &[term]);
    engine.run(
        graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
    Extractor::new(graph, &[root], &ClosedForm)
        .build(graph, root)
        .unwrap()
}

fn assert_value(
    graph: &Graph,
    root: NodeId,
    expected: f64,
) {
    let value = graph.eval(root, &Env::default()).unwrap();
    assert!((value - expected).abs() < 1e-12, "{value} vs {expected}");
}

#[test]
fn radicals_shared_dag_preserves_shared_children_for_both_requests() {
    for (request, expression, expected) in [
        ("denest_sqrt", "(3+2*2^(1/2))^(1/2)", 1.0 + 2.0_f64.sqrt()),
        (
            "simplify_radicals",
            "1/((3+2*2^(1/2))^(1/2))",
            2.0_f64.sqrt() - 1.0,
        ),
    ] {
        let mut graph = Graph::new();
        let engine = Engine::install(&mut graph, &[rssn::rules::arith(), radicals()]).unwrap();
        let pair = graph
            .ops_mut()
            .register(OpDescriptor::new("same_pair", Arity::Fixed(2)))
            .unwrap();
        let mut input = graph.parse(expression).unwrap();
        for _ in 0..20 {
            input = graph.node(pair, &[input, input]);
        }
        let original_size = graph.len();
        let mut found = operation(&mut graph, &engine, request, input);
        for _ in 0..20 {
            assert_eq!(graph.op(found), pair);
            let &[left, right] = graph.children(found) else {
                panic!("expected binary pair")
            };
            assert_eq!(left, right);
            found = left;
        }
        assert_value(&graph, found, expected);
        assert!(
            graph.len() < original_size + 100,
            "the result should retain the shared graph"
        );
    }
}

#[test]
fn radicals_shared_dag_handles_a_deep_single_path() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rssn::rules::arith(), radicals()]).unwrap();
    let marker = graph
        .ops_mut()
        .register(OpDescriptor::new("marker", Arity::Fixed(1)))
        .unwrap();
    let mut input = graph.parse("(5-2*6^(1/2))^(1/2)").unwrap();
    for _ in 0..2048 {
        input = graph.node(marker, &[input]);
    }
    let mut found = operation(&mut graph, &engine, "denest_sqrt", input);
    for _ in 0..2048 {
        assert_eq!(graph.op(found), marker);
        let &[child] = graph.children(found) else {
            panic!("expected unary marker")
        };
        found = child;
    }
    assert_value(&graph, found, 3.0_f64.sqrt() - 2.0_f64.sqrt());
}

#[test]
fn radicals_shared_dag_keeps_the_order_and_values_of_unshared_arguments() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rssn::rules::arith(), radicals()]).unwrap();
    let triple = graph
        .ops_mut()
        .register(OpDescriptor::new("triple", Arity::Fixed(3)))
        .unwrap();
    let a = graph.parse("(3+2*2^(1/2))^(1/2)").unwrap();
    let b = graph.parse("(2+5*3^(1/2))^(1/2)").unwrap();
    let c = graph.parse("1/(1+2^(1/2))").unwrap();
    let root = graph.node(triple, &[a, b, c]);
    let result = operation(&mut graph, &engine, "simplify_radicals", root);
    assert_eq!(graph.op(result), triple);
    let &[first, second, third] = graph.children(result) else {
        panic!("expected triple")
    };
    assert_value(&graph, first, 1.0 + 2.0_f64.sqrt());
    assert_value(&graph, second, (2.0 + 5.0 * 3.0_f64.sqrt()).sqrt());
    assert_value(&graph, third, 2.0_f64.sqrt() - 1.0);
    assert!(graph.op(third) != core::POW || graph.children(third).len() != 2);
}
