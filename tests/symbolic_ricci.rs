//! Analytic checks for symbolic Ricci contraction and curvature orientation.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Saturate;
use rssn::rules::geometry::geometry;

fn reduce(
    graph: &mut Graph,
    engine: &Engine,
    text: &str,
) -> Option<NodeId> {
    let root = graph.parse(text).unwrap();
    engine.run(
        graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{text}");
    Extractor::new(graph, &[root], &ClosedForm).build(graph, root)
}

fn environment(
    graph: &mut Graph,
    bindings: &[(&str, f64)],
) -> Env {
    let mut env = Env::default();
    for &(name, value) in bindings {
        env.bind(graph.interner_mut().symbol(name), value);
    }
    env
}

fn check_value(
    graph: &Graph,
    env: &Env,
    term: NodeId,
    expected: f64,
) {
    let value = graph.eval(term, env).unwrap();
    assert!((value - expected).abs() < 1e-12, "{value} vs {expected}");
}

fn check_tensor<const N: usize>(
    graph: &Graph,
    env: &Env,
    term: NodeId,
    expected: [[f64; N]; N],
) {
    let rows = graph.children(term);
    assert_eq!(rows.len(), N);
    for (&row, wanted) in rows.iter().zip(expected) {
        let entries = graph.children(row);
        assert_eq!(entries.len(), N);
        for (&entry, value) in entries.iter().zip(wanted) {
            check_value(graph, env, entry, value);
        }
    }
}

#[test]
fn symbolic_ricci_of_a_sphere_times_a_line_matches_positive_analytic_curvature() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[geometry()]).unwrap();
    let metric = "list(list(4,0,0),list(0,4*sin(t)^2,0),list(0,0,1))";
    let vars = "list(t,p,z)";
    let ricci = reduce(&mut graph, &engine, &format!("ricci({metric},{vars})")).unwrap();
    let scalar = reduce(
        &mut graph,
        &engine,
        &format!("ricci_scalar({metric},{vars})"),
    )
    .unwrap();
    let einstein = reduce(&mut graph, &engine, &format!("einstein({metric},{vars})")).unwrap();
    for angle in [0.3_f64, 0.7, 1.1] {
        let env = environment(&mut graph, &[("t", angle), ("p", 0.2), ("z", 1.0)]);
        check_tensor(
            &graph,
            &env,
            ricci,
            [[1.0, 0.0, 0.0], [0.0, angle.sin().powi(2), 0.0], [0.0; 3]],
        );
        check_value(&graph, &env, scalar, 0.5);
        check_tensor(
            &graph,
            &env,
            einstein,
            [[0.0; 3], [0.0; 3], [0.0, 0.0, -0.25]],
        );
    }
}

fn rank_four(
    graph: &Graph,
    node: NodeId,
    indices: [usize; 4],
) -> NodeId {
    indices
        .into_iter()
        .fold(node, |parent, index| graph.children(parent)[index])
}

#[test]
fn symbolic_ricci_contracts_both_orientations_of_the_public_riemann_tensor() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[geometry()]).unwrap();
    // For e^(2u) I in three dimensions and u = x + y, Ric is
    // du ⊗ du - |du|² I. Its off-diagonal entries exercise both orientations.
    let exp = "exp(2*(x+y))";
    let metric = format!("list(list({exp},0,0),list(0,{exp},0),list(0,0,{exp}))");
    let vars = "list(x,y,z)";
    let ricci = reduce(&mut graph, &engine, &format!("ricci({metric},{vars})")).unwrap();
    let riemann = reduce(&mut graph, &engine, &format!("riemann({metric},{vars})")).unwrap();
    let expected = [[-1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 0.0, -2.0]];
    let env = environment(&mut graph, &[("x", 0.125), ("y", -0.25), ("z", 0.0)]);
    check_tensor(&graph, &env, ricci, expected);
    for (s, row) in expected.iter().enumerate() {
        for (v, &wanted) in row.iter().enumerate() {
            let mut contracted = 0.0;
            for k in 0..3 {
                let direct = rank_four(&graph, riemann, [k, s, k, v]);
                let reverse = rank_four(&graph, riemann, [k, s, v, k]);
                let value = graph.eval(direct, &env).unwrap();
                check_value(&graph, &env, reverse, -value);
                contracted += value;
            }
            assert!((contracted - wanted).abs() < 1e-12);
        }
    }
}

#[test]
fn symbolic_ricci_of_a_non_diagonal_metric_transforms_as_a_covariant_tensor() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[geometry()]).unwrap();
    // Pull e^(2x) I back by the shear (x,y,z) -> (x,x+y,z).
    let exp = "exp(2*x)";
    let metric = format!("list(list(2*{exp},{exp},0),list({exp},{exp},0),list(0,0,{exp}))");
    let vars = "list(x,y,z)";
    let result = reduce(&mut graph, &engine, &format!("ricci({metric},{vars})")).unwrap();
    // The Ricci tensor in the unsheared coordinates is diag(0,-1,-1).
    let expected = [[-1.0, -1.0, 0.0], [-1.0, -1.0, 0.0], [0.0, 0.0, -1.0]];
    for point in [0.125, -0.25] {
        let env = environment(&mut graph, &[("x", point), ("y", 0.2), ("z", 1.0)]);
        check_tensor(&graph, &env, result, expected);
    }
}

#[test]
fn symbolic_ricci_preserves_singular_and_incompatible_shape_refusals() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[geometry()]).unwrap();
    for source in [
        "ricci(list(list(1,0),list(0,0)),list(x,y))",
        "ricci(list(list(1,0),list(0,1)),list(x,y,z))",
        "ricci(list(list(1,0,0),list(0,1,0)),list(x,y))",
    ] {
        assert!(reduce(&mut graph, &engine, source).is_none(), "{source}");
    }
}
