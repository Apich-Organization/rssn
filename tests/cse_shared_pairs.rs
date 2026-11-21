//! Regression tests for symbolic common-subexpression elimination.

use std::collections::HashSet;

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::rules::cse::cse;

fn run(
    graph: &mut Graph,
    engine: &Engine,
    op: &str,
    argument: NodeId,
) -> NodeId {
    let request = graph.ops().lookup(op).unwrap();
    let root = graph.node(request, &[argument]);
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

fn cse_parts(
    graph: &Graph,
    result: NodeId,
) -> (Vec<NodeId>, Vec<NodeId>) {
    let [definitions, outputs] = graph.children(result) else {
        panic!("expected CSE result")
    };
    (
        graph.children(*definitions).to_vec(),
        graph.children(*outputs).to_vec(),
    )
}

fn equality(
    graph: &Graph,
    equation: NodeId,
) -> (NodeId, NodeId) {
    assert_eq!(graph.op(equation), core::EQ);
    let [lhs, rhs] = graph.children(equation) else {
        panic!("expected equation")
    };
    (*lhs, *rhs)
}

fn verify_reconstruction(
    graph: &mut Graph,
    definitions: &[NodeId],
    results: &[NodeId],
    originals: &[NodeId],
) {
    let all_temporaries: HashSet<_> = definitions
        .iter()
        .map(|&definition| graph.as_symbol(equality(graph, definition).0).unwrap())
        .collect();
    assert_eq!(all_temporaries.len(), definitions.len());
    for (a, b, c) in [(1.5, -0.25, 2.0), (-0.5, 0.25, -1.0)] {
        let mut env = Env::default();
        for (name, value) in [("x", a), ("y", b), ("t3", c)] {
            env.bind(graph.interner_mut().symbol(name), value);
        }
        let original_env = env.clone();
        let mut available = HashSet::new();
        for &definition in definitions {
            let (lhs, rhs) = equality(graph, definition);
            for symbol in graph.free_symbols(graph.find(rhs)) {
                assert!(!all_temporaries.contains(symbol) || available.contains(symbol));
            }
            let symbol = graph.as_symbol(lhs).unwrap();
            env.bind(symbol, graph.eval(rhs, &env).unwrap());
            available.insert(symbol);
        }
        assert_eq!(results.len(), originals.len());
        for (&rewritten, &original) in results.iter().zip(originals) {
            let (after_lhs, after_rhs) = equality(graph, rewritten);
            let (before_lhs, before_rhs) = equality(graph, original);
            assert_eq!(graph.display(after_lhs), graph.display(before_lhs));
            assert_eq!(
                graph.eval(after_rhs, &env).unwrap().to_bits(),
                graph.eval(before_rhs, &original_env).unwrap().to_bits(),
            );
        }
    }
}

#[test]
fn cse_shared_pairs_skips_output_and_input_symbols_with_ordered_definitions() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    let source = graph
        .parse("list(t1=(1+(x+y)^2)^2+(1+(x+y)^2)^3,t2^2=(x+y)^4+(1+(x+y)^2)^4+t3)")
        .unwrap();
    let originals = graph.children(source).to_vec();
    let result = run(&mut graph, &engine, "cse", source);
    let (definitions, results) = cse_parts(&graph, result);
    assert!(definitions.len() >= 2);
    for (i, &definition) in definitions.iter().enumerate() {
        let (lhs, _) = equality(&graph, definition);
        assert_eq!(graph.display(lhs), format!("t{}", i + 4));
    }
    verify_reconstruction(&mut graph, &definitions, &results, &originals);
}

fn unique_operands(
    graph: &mut Graph,
    stem: &str,
    count: usize,
) -> Vec<NodeId> {
    (0..count)
        .map(|i| graph.sym(&format!("{stem}{i}")))
        .collect()
}

#[test]
fn cse_shared_pairs_wide_disjoint_groups_preserve_cost_and_expression_order() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    let a = unique_operands(&mut graph, "a", 192);
    let b = unique_operands(&mut graph, "b", 192);
    let roots = [
        graph.node(core::ADD, &a),
        graph.node(core::ADD, &b),
        graph.node(core::MUL, &a),
    ];
    let source = graph.node(core::LIST, &roots);
    let cost = run(&mut graph, &engine, "cse_cost", source);
    assert_eq!(graph.display(cost), "list(573, 573)");
    let result = run(&mut graph, &engine, "cse", source);
    let (definitions, results) = cse_parts(&graph, result);
    assert!(definitions.is_empty());
    assert_eq!(
        results
            .iter()
            .map(|&node| graph.display(node))
            .collect::<Vec<_>>(),
        roots.map(|node| graph.display(node))
    );
}

#[test]
fn cse_shared_pairs_preserves_repeated_operands_and_existing_binary_pair() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[cse()]).unwrap();
    let x = graph.sym("x");
    let y = graph.sym("y");
    let p = graph.node(core::ADD, &[x, y]);
    let sensor = graph.sym("sensor");
    let observed = graph.node(core::MUL, &[p, sensor]);
    let mut a = unique_operands(&mut graph, "a", 35);
    let mut b = unique_operands(&mut graph, "b", 35);
    a.extend([x, x, y]);
    b.extend([x, y, y]);
    let items = [
        graph.node(core::ADD, &a),
        graph.node(core::ADD, &b),
        observed,
    ];
    let source = graph.node(core::LIST, &items);
    let result = run(&mut graph, &engine, "cse", source);
    let (definitions, results) = cse_parts(&graph, result);
    assert_eq!(definitions.len(), 1);
    let (name, rhs) = equality(&graph, definitions[0]);
    assert_eq!(graph.display(rhs), "x + y");
    let mut env = Env::default();
    for terms in [&a, &b] {
        for &v in terms {
            env.bind(graph.as_symbol(v).unwrap(), 0.125);
        }
    }
    for (node, value) in [(x, 1.5), (y, -0.25), (sensor, 2.0)] {
        env.bind(graph.as_symbol(node).unwrap(), value);
    }
    env.bind(
        graph.as_symbol(name).unwrap(),
        graph.eval(rhs, &env).unwrap(),
    );
    for (new, old) in results.into_iter().zip(items) {
        assert_eq!(
            graph.eval(new, &env).unwrap().to_bits(),
            graph.eval(old, &env).unwrap().to_bits()
        );
    }
}
