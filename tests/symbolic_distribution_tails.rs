//! Tail probabilities of continuous distributions in the expression graph.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::NodeId;
use rssn::graph::Saturate;

fn reduced(
    graph: &mut Graph,
    engine: &Engine,
    expression: &str,
    env: &Env,
) -> NodeId {
    let root = graph.parse(expression).unwrap();
    engine.run(graph, &[root], env, &Saturate, &Budget::default());
    Extractor::new(graph, &[root], &ClosedForm)
        .build(graph, root)
        .unwrap_or(root)
}

#[test]
fn concrete_cdf_preserves_tiny_but_representable_lower_tails() {
    let cases = [
        ("cdf(normal(0.0, 1.0), -10.0)", 7.619_853_024_160_526e-24),
        (
            "cdf(lognormal(0.0, 1.0), 0.00004539992976248485)",
            7.619_853_024_160_526e-24,
        ),
        ("cdf(exponential(2.0), 1e-20)", 2e-20),
        ("cdf(weibull(3.0, 2.0), 2e-8)", 1e-24),
        ("cdf(cauchy(0.0, 1.0), -1e20)", 3.183_098_861_837_907e-21),
        (
            "cdf(normal(1.5e308, 1.5e308), -1.5e308)",
            0.022_750_131_948_179_21,
        ),
    ];
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rssn::rules::stats()]).unwrap();
    let env = Env::numeric(0.0);
    for (expression, expected) in cases {
        let result = reduced(&mut graph, &engine, expression, &env);
        let observed = graph.number_of(result).unwrap().to_f64();
        assert!(
            (observed / expected - 1.0).abs() < 3e-13,
            "{expression}: {observed:e}"
        );
    }
}

#[test]
fn cdf_keeps_exact_values_and_symbolic_formulas() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rssn::rules::stats()]).unwrap();
    let symbolic = Env::symbolic();
    let center = reduced(&mut graph, &engine, "cdf(normal(0, 1), 0)", &symbolic);
    assert_eq!(graph.display(center).to_string(), "1/2");
    let result = reduced(&mut graph, &engine, "cdf(normal(0, 1), x)", &symbolic);
    assert!(graph.display(result).to_string().contains("erfc"));
    let mut numeric = Env::numeric(0.0);
    numeric.bind(graph.interner_mut().symbol("x"), -10.0);
    let tail = graph.eval(result, &numeric).unwrap();
    assert!((tail / 7.619_853_024_160_526e-24 - 1.0).abs() < 3e-13);

    for (expression, expected) in [
        ("cdf(cauchy(0, 1), -1)", "1/4"),
        ("cdf(cauchy(0, 1), 0)", "1/2"),
        ("cdf(cauchy(0, 1), 1)", "3/4"),
        ("cdf(cauchy(2, 3/2), 1/2)", "1/4"),
        ("cdf(cauchy(2, 3/2), 7/2)", "3/4"),
    ] {
        let result = reduced(&mut graph, &engine, expression, &symbolic);
        assert_eq!(graph.display(result).to_string(), expected, "{expression}");
    }

    let cauchy = reduced(&mut graph, &engine, "cdf(cauchy(0, 1), x)", &symbolic);
    numeric.bind(graph.interner_mut().symbol("x"), -1e20);
    let tail = graph.eval(cauchy, &numeric).unwrap();
    assert!((tail / 3.183_098_861_837_907e-21 - 1.0).abs() < 3e-13);
}

#[test]
fn cdf_preserves_support_and_declines_invalid_parameters() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rssn::rules::stats()]).unwrap();
    let env = Env::numeric(0.0);
    for expression in [
        "cdf(lognormal(0.0, 1.0), -1.0)",
        "cdf(weibull(3.0, 2.0), -1.0)",
    ] {
        let result = reduced(&mut graph, &engine, expression, &env);
        assert_eq!(graph.number_of(result).unwrap().to_f64(), 0.0);
    }
    let request = reduced(&mut graph, &engine, "cdf(normal(0.0, -1.0), -5.0)", &env);
    assert_eq!(graph.ops().get(graph.op(request)).name.as_ref(), "cdf");
}
