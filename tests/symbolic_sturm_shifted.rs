//! Shifted Euler–Cauchy evolution through the public symbolic PDE solver.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Evaluated;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Saturate;

fn solve(source: &str) -> (Graph, Engine, Option<NodeId>) {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rssn::rules::standard()).unwrap();
    let request = graph.parse(source).unwrap();
    engine.run(
        &mut graph,
        &[request],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    let solved = Extractor::new(&graph, &[request], &ClosedForm).build(&mut graph, request);
    (graph, engine, solved)
}

fn value_at(
    graph: &mut Graph,
    engine: &Engine,
    answer: NodeId,
    x: f64,
    t: f64,
) -> f64 {
    let equation = graph.display(answer);
    assert!(!equation.contains("sl_mode"), "{equation}");
    let &[_, rhs] = graph.children(answer) else {
        panic!("expected an explicit equation: {equation}");
    };
    let mut variables = Env::numeric(1e-10);
    for (name, value) in [("x", x), ("t", t)] {
        let symbol = graph.sym(name);
        variables.bind(graph.symbol_of(symbol).unwrap(), value);
    }
    if let Some(value) = graph.eval(rhs, &variables) {
        return value;
    }
    engine.run(graph, &[rhs], &variables, &Evaluated, &Budget::default());
    graph
        .approx(graph.find(rhs))
        .unwrap_or_else(|| panic!("unable to evaluate {equation}"))
        .mid
}

#[test]
fn shifted_log_charts_make_both_dirichlet_branches_explicit() {
    for (quadratic, bounds, profile, side) in [
        (
            "x^2 - 6*x + 9",
            "u(4,t)=0, u(3+exp(pi),t)=0",
            "sin(ln(x-3))",
            1.0,
        ),
        (
            "(x-3)^2",
            "u(3-exp(pi),t)=0, u(2,t)=0",
            "sin(ln(3-x))",
            -1.0,
        ),
    ] {
        let request = format!(
            "pdsolve(diff(u(x,t),t) = ({quadratic})*diff(diff(u(x,t),x),x) + (x-3)*diff(u(x,t),x), u(x,t), list({bounds}, u(x,0)={profile}))"
        );
        let (mut graph, engine, solution) = solve(&request);
        let solution = solution.expect("Euler–Cauchy solution");
        let equation = graph.display(solution);
        assert!(!equation.contains("defint"), "{equation}");
        for (distance, t) in [(1.0_f64, 0.5_f64), (2.35, 0.23), (6.0, 0.75)] {
            let value = value_at(&mut graph, &engine, solution, 3.0 + side * distance, t);
            assert!((value - (-t).exp() * distance.ln().sin()).abs() < 1e-10);
        }
    }
}

#[test]
fn shifted_negative_robin_boundary_uses_the_signed_distance() {
    let source = "pdsolve(diff(u(x,t),t) = (x-3)^2*diff(diff(u(x,t),x),x) + (x-3)*diff(u(x,t),x), u(x,t), list(at(diff(u(x,t),x),x,2)-u(2,t)=0, u(x,0)=1/(3-x)))";
    let (mut graph, engine, solution) = solve(source);
    let solution = solution.expect("shifted Robin solution");
    let (x, t) = (0.65_f64, 0.4_f64);
    let value = value_at(&mut graph, &engine, solution, x, t);
    assert!((value - t.exp() / (3.0 - x)).abs() < 1e-7);
}
