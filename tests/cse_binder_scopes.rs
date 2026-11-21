//! Regression tests for preserving lexical scopes during CSE.

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
use rssn::rules::cse::cse;

fn run(
    graph: &mut Graph,
    engine: &Engine,
    op: &str,
    inputs: &[NodeId],
) -> NodeId {
    let list = graph.node(core::LIST, inputs);
    let op = graph.ops().lookup(op).unwrap();
    let root = graph.node(op, &[list]);
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

fn parts(
    graph: &Graph,
    result: NodeId,
) -> (Vec<NodeId>, Vec<NodeId>) {
    let [bindings, outputs] = graph.children(result) else {
        panic!("expected CSE lists")
    };
    (
        graph.children(*bindings).to_vec(),
        graph.children(*outputs).to_vec(),
    )
}

#[test]
fn cse_counts_large_expanded_trees_exactly_and_visits_repeated_outputs() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    let pair = graph
        .ops_mut()
        .register(OpDescriptor::new("pair", Arity::Fixed(2)))
        .unwrap();
    let mut root = graph.sym("x");
    let depth = 73_usize;
    for _ in 0..depth {
        root = graph.node(pair, &[root, root]);
    }
    let cost = run(&mut graph, &engine, "cse_cost", &[root, root]);
    let expanded = (num_bigint::BigUint::from(1_u8) << (depth + 1)) - 2_u8;
    assert_eq!(
        graph.display(cost),
        format!("list({expanded}, {})", depth + 1)
    );
    let plans = run(&mut graph, &engine, "cse", &[root, root]);
    let (definitions, results) = parts(&graph, plans);
    assert_eq!(definitions.len(), depth - 1);
    assert_eq!(results.len(), 2);
}

#[test]
fn cse_cost_treats_a_temporary_reference_inside_an_unshared_node_as_a_leaf() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    for name in ["inner", "outer", "second"] {
        graph
            .ops_mut()
            .register(OpDescriptor::new(name, Arity::Fixed(1)))
            .unwrap();
    }
    graph
        .ops_mut()
        .register(OpDescriptor::new("common", Arity::Fixed(2)))
        .unwrap();
    let a = graph.parse("outer(inner(common(x,y)))").unwrap();
    let b = graph.parse("second(common(x,y))").unwrap();
    let cost = run(&mut graph, &engine, "cse_cost", &[a, b]);
    assert_eq!(graph.display(cost), "list(5, 4)");
}

#[test]
fn exact_cost_keeps_nested_tree_occurrences_inside_a_lexical_scope() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    let pair = graph
        .ops_mut()
        .register(OpDescriptor::new("pair", Arity::Fixed(2)))
        .unwrap();
    let local = graph
        .ops_mut()
        .register(OpDescriptor::new("local", Arity::Fixed(2)).binder(0, 0b10))
        .unwrap();
    let x = graph.sym("x");
    let mut body = x;
    for _ in 0..81 {
        body = graph.node(pair, &[body, body]);
    }
    let root = graph.node(local, &[x, body]);
    let cost = run(&mut graph, &engine, "cse_cost", &[root]);
    let exact = num_bigint::BigUint::from(1_u8) << 81_usize;
    assert_eq!(graph.display(cost), format!("list({exact}, {exact})"));
}

#[test]
fn cse_binder_scopes_does_not_expose_a_bound_symbol_or_hide_body_cost() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    graph
        .ops_mut()
        .register(OpDescriptor::new("local", Arity::Fixed(2)).binder(0, 0b10))
        .unwrap();
    let first = graph.parse("local(x,(x+y)^2)").unwrap();
    let second = graph.parse("local(x,(x+y)^3)").unwrap();
    let x = graph.interner().find_symbol("x").unwrap();
    let result = run(&mut graph, &engine, "cse", &[first, second]);
    let (bindings, outputs) = parts(&graph, result);
    assert!(bindings.is_empty());
    assert_eq!(outputs, [first, second]);
    assert!(!graph.free_symbols(graph.find(result)).contains(&x));
    let cost = run(&mut graph, &engine, "cse_cost", &[first, second]);
    // Each unchanged binder has its own operation, power and sum.
    assert_eq!(graph.display(cost), "list(6, 6)");
}

#[test]
fn cse_binder_scopes_still_shares_unscoped_occurrences_and_preserves_body_identity() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    graph
        .ops_mut()
        .register(OpDescriptor::new("local", Arity::Fixed(3)).binder(0, 0b010))
        .unwrap();
    let first = graph.parse("local(t1,t1*(y+z)^2,(y+z)^2)").unwrap();
    let second = graph.parse("local(t1,t1*(y+z)^3,(y+z)^3)").unwrap();
    let result = run(&mut graph, &engine, "cse", &[first, second]);
    let (bindings, outputs) = parts(&graph, result);
    let [definition] = bindings.as_slice() else {
        panic!("expected a shared unscoped sum")
    };
    let &[name, expression] = graph.children(*definition) else {
        panic!("expected a definition")
    };
    assert_eq!(graph.display(name), "t2");
    assert_eq!(graph.display(expression), "y + z");
    assert_eq!(outputs.len(), 2);
    for (output, original) in outputs.into_iter().zip([first, second]) {
        let original_children = graph.children(original).to_vec();
        let rewritten = graph.children(output);
        assert_eq!(&rewritten[..2], &original_children[..2]);
        assert_ne!(rewritten[2], original_children[2]);
        let reconstructed = graph.substitute(output, name, expression);
        assert_eq!(graph.display(reconstructed), graph.display(original));
    }
}

#[test]
fn cse_binder_scopes_keeps_nested_binders_local() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    graph
        .ops_mut()
        .register(OpDescriptor::new("local", Arity::Fixed(2)).binder(0, 0b10))
        .unwrap();
    let first = graph.parse("local(x,local(y,(x+y+z)^2))").unwrap();
    let second = graph.parse("local(x,local(y,(x+y+z)^3))").unwrap();
    let result = run(&mut graph, &engine, "cse", &[first, second]);
    let (bindings, outputs) = parts(&graph, result);
    assert!(bindings.is_empty());
    assert_eq!(outputs, [first, second]);
    assert_eq!(
        graph.free_symbols(graph.find(result)),
        [graph.interner().find_symbol("z").unwrap()],
    );
}

#[test]
fn cse_binder_scopes_respects_the_highest_scope_bit_and_a_later_binder_slot() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    let op = graph
        .ops_mut()
        .register(OpDescriptor::new("wide_local", Arity::Fixed(34)).binder(32, 1 << 31))
        .unwrap();
    let x = graph.sym("x");
    let zero = graph.int(0);
    let mut inputs = [zero; 34];
    inputs[31] = graph.parse("(x+y)^2").unwrap();
    inputs[32] = x;
    let first = graph.node(op, &inputs);
    inputs[31] = graph.parse("(x+y)^3").unwrap();
    let second = graph.node(op, &inputs);
    let result = run(&mut graph, &engine, "cse", &[first, second]);
    let (bindings, outputs) = parts(&graph, result);
    assert!(bindings.is_empty());
    assert_eq!(outputs, [first, second]);
}
