//! Symbolic unit evaluation on shared terms and at the recursion boundary.

use rssn::graph::Arity;
use rssn::graph::Budget;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::OpDescriptor;
use rssn::graph::Saturate;
use rssn::rules::units::units;

fn run(
    graph: &mut Graph,
    engine: &Engine,
    root: NodeId,
) {
    engine.run(
        graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
}

#[test]
fn dimension_evaluates_a_deep_shared_scalar_expression() {
    let mut graph = Graph::new();
    let opaque = graph
        .ops_mut()
        .register(OpDescriptor::new("opaque_pair", Arity::Fixed(2)))
        .unwrap();
    let engine = Engine::install(&mut graph, &[units()]).unwrap();
    let quantity = graph.ops().lookup("quantity").unwrap();
    let dimension = graph.ops().lookup("dimension").unwrap();
    let start = graph.sym("x");
    let unit = graph.sym("m");
    let expected = graph.sym("Length");
    let value = (0..24).fold(start, |node, _| graph.node(opaque, &[node, node]));
    let quantity = graph.node(quantity, &[value, unit]);
    let request = graph.node(dimension, &[quantity]);
    run(&mut graph, &engine, request);
    assert_eq!(graph.find(request), graph.find(expected));
}

#[test]
fn scalar_depth_limit_is_checked_after_reusing_a_shallower_subexpression() {
    let mut graph = Graph::new();
    let unary = graph
        .ops_mut()
        .register(OpDescriptor::new("opaque_unary", Arity::Fixed(1)))
        .unwrap();
    let pair = graph
        .ops_mut()
        .register(OpDescriptor::new("opaque_pair", Arity::Fixed(2)))
        .unwrap();
    let engine = Engine::install(&mut graph, &[units()]).unwrap();
    let dimension = graph.ops().lookup("dimension").unwrap();
    let quantity = graph.ops().lookup("quantity").unwrap();
    let x = graph.sym("x");
    let meter = graph.sym("m");
    let expected = graph.sym("Length");
    let shared = graph.node(unary, &[x]);
    let deep = (0..46).fold(shared, |node, _| graph.node(unary, &[node]));
    let accepted = graph.node(pair, &[shared, deep]);
    let accepted = graph.node(quantity, &[accepted, meter]);
    let accepted = graph.node(dimension, &[accepted]);
    run(&mut graph, &engine, accepted);
    assert_eq!(graph.find(accepted), graph.find(expected));

    let too_deep = graph.node(unary, &[deep]);
    let rejected = graph.node(pair, &[shared, too_deep]);
    let rejected = graph.node(quantity, &[rejected, meter]);
    let rejected = graph.node(dimension, &[rejected]);
    run(&mut graph, &engine, rejected);
    assert_ne!(graph.find(rejected), graph.find(expected));
}
