//! Shared symbolic terms in transform and functional-application rules.

use rssn::graph::Arity;
use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::OpDescriptor;
use rssn::graph::OpFlags;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::rules::functional::functional;
use rssn::rules::transforms::transforms;

fn shared_coefficient(
    graph: &mut Graph,
    depth: usize,
    pending: bool,
) -> NodeId {
    let pair = graph
        .ops_mut()
        .register(OpDescriptor::new("pair", Arity::Fixed(2)))
        .unwrap();
    let mut left = graph.sym("a");
    let mut right = graph.sym("b");
    if pending {
        let request = graph
            .ops_mut()
            .register(OpDescriptor::new("waiting", Arity::Fixed(1)).flags(OpFlags::HEAVY))
            .unwrap();
        right = graph.node(request, &[right]);
    }
    for _ in 0..depth {
        (left, right) = (
            graph.node(pair, &[left, right]),
            graph.node(pair, &[right, left]),
        );
    }
    graph.node(pair, &[left, right])
}

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
fn symbolic_transform_keeps_a_shared_coefficient_and_exact_rational_factor() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[transforms()]).unwrap();
    let coefficient = shared_coefficient(&mut graph, 20, false);
    let decay = graph.parse("exp(-t)").unwrap();
    let f = graph.node(core::MUL, &[coefficient, decay]);
    let t = graph.sym("t");
    let s = graph.sym("s");
    let laplace = graph.ops().lookup("laplace").unwrap();
    let root = graph.node(laplace, &[f, t, s]);
    run(&mut graph, &engine, root);
    let answer = Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .unwrap();
    let inverse = graph.parse("1/(s+1)").unwrap();
    let expected = graph.node(core::MUL, &[coefficient, inverse]);
    run(&mut graph, &engine, expected);
    assert!(graph.same(answer, expected));
    assert!(
        graph
            .children(answer)
            .iter()
            .any(|&factor| graph.same(coefficient, factor))
    );
}

#[test]
fn symbolic_transform_preserves_the_pinned_derivative_theorem() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[transforms()]).unwrap();
    let root = graph.parse("laplace(diff(y(t),t),t,s)").unwrap();
    run(&mut graph, &engine, root);
    let answer = graph.pinned(graph.find(root)).unwrap();
    assert_eq!(graph.display(answer), "s*laplace(y(t), t, s) - y(0)");
    assert_eq!(
        Extractor::new(&graph, &[root], &ClosedForm).build(&mut graph, root),
        Some(answer)
    );
}

#[test]
fn symbolic_functional_multiplier_distinguishes_concrete_and_pending_subterms() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[functional()]).unwrap();
    let coefficient = shared_coefficient(&mut graph, 12, false);
    let three = graph.int(3);
    let apply = graph.ops().lookup("op_apply").unwrap();
    let root = graph.node(apply, &[coefficient, three]);
    run(&mut graph, &engine, root);
    let answer = Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .unwrap();
    let expected = graph.node(core::MUL, &[coefficient, three]);
    run(&mut graph, &engine, expected);
    assert!(graph.same(answer, expected));

    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[functional()]).unwrap();
    let coefficient = shared_coefficient(&mut graph, 12, true);
    let three = graph.int(3);
    let apply = graph.ops().lookup("op_apply").unwrap();
    let root = graph.node(apply, &[coefficient, three]);
    run(&mut graph, &engine, root);
    assert!(
        Extractor::new(&graph, &[root], &ClosedForm)
            .build(&mut graph, root)
            .is_none()
    );
}
