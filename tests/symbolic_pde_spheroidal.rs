//! Regular radial values used by the symbolic spheroidal PDE solver.

use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Graph;

fn graph() -> Graph {
    let mut graph = Graph::new();
    Engine::install(&mut graph, &[rssn::rules::pde()]).unwrap();
    graph
}

fn radial(
    graph: &mut Graph,
    args: [f64; 5],
) -> f64 {
    let op = graph.ops().lookup("spheroidal_radial").unwrap();
    let nodes = args.map(|value| graph.float(value));
    let term = graph.node(op, &nodes);
    graph.eval(term, &Env::numeric(0.0)).unwrap_or(f64::NAN)
}

#[test]
fn laplace_radial_modes_agree_with_analytic_profiles() {
    let mut graph = graph();
    for x in [1.0, 1.0001, 1.7, 4.9, 1.13] {
        let linear = radial(&mut graph, [1.0, 0.0, 1.0, 0.0, x]);
        assert!((linear - x).abs() < 2e-7 * x, "ξ={x}: {linear} vs {x}");
        let order_one = radial(&mut graph, [1.0, 1.0, 1.0, 0.0, x]);
        let expected = ((x - 1.0) * (x + 1.0)).sqrt();
        assert!((order_one - expected).abs() < 2e-8 * expected.max(1.0));
    }
    for x in [0.0, 0.0001, 0.6, 5.0, 0.2] {
        let linear = radial(&mut graph, [-1.0, 0.0, 1.0, 0.0, x]);
        assert!((linear - x).abs() < 2e-8 * x.max(1.0));
        // The sign convention also accepts finite zero as oblate.
        assert_eq!(
            radial(&mut graph, [0.0, 0.0, 1.0, 0.0, x]).to_bits(),
            linear.to_bits()
        );
    }
}

#[test]
fn radial_queries_reject_bad_controls_and_keep_a_finite_prefix() {
    let mut graph = graph();
    for args in [
        [f64::NAN, 0.0, 1.0, 0.0, 0.4],
        [1.0, 0.0, 1.0, 0.0, 0.99],
        [-1.0, 0.0, 1.0, 0.0, -0.1],
        [-1.0, 0.0, 1.0, f64::INFINITY, 0.1],
        [1.0, 1.5, 2.0, 0.0, 1.2],
        [1.0, 2.0, 1.0, 0.0, 1.2],
        [1.0, 1e308, 1e308, 0.0, 1.2],
        [1.0, 0.0, 1.0, 0.0, 1e308],
    ] {
        assert!(radial(&mut graph, args).is_nan());
    }
    let early = radial(&mut graph, [1.0, 0.0, 3.0, -1000.0, 1.0001]);
    assert!(early.is_finite());
    assert!(radial(&mut graph, [1.0, 0.0, 3.0, -1000.0, 40.0]).is_nan());
    assert_eq!(
        radial(&mut graph, [1.0, 0.0, 3.0, -1000.0, 1.0001]).to_bits(),
        early.to_bits()
    );
    assert_eq!(radial(&mut graph, [1.0, 0.0, 3.0, -1000.0, 1.0]), 1.0);
}
