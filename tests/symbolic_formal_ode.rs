//! Numerical evaluation of formal solutions through the symbolic ODE operators.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::kernels::integrate::gauss_kronrod_any;
use rssn::rules;

fn euler_solution(
    operator: &str,
    terms: usize,
    branch: usize,
) -> Vec<f64> {
    // Differentiate x²y' + (1 + x)y = 1 to give the homogeneous equation.
    let equation = "x^2*diff(diff(y(x),x),x)+(1+3*x)*diff(y(x),x)+y(x)";
    let source = format!("{operator}({equation},y(x),0,{terms},1/4,{branch})");
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    let root = graph.parse(&source).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{source}");
    let value = Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .unwrap_or_else(|| panic!("unresolved {source}"));
    assert_eq!(graph.op(value), core::LIST, "{source}");
    graph
        .children(value)
        .iter()
        .map(|&child| graph.number_of(child).unwrap().to_f64())
        .collect()
}

#[test]
fn euler_integral_agrees_with_asymptotic_and_borel_evaluation() {
    // Integrate the finite geometric expansion of 1/(1 + t/4) against
    // exp(-t). The kth term is (-1/4)^k k!, and the remainder at k = 5
    // has the sign of that term and magnitude strictly below 5! / 4^5.
    let reference = gauss_kronrod_any(
        |t| (-t).exp() / (1.0 + t / 4.0),
        0.0,
        f64::INFINITY,
        1e-13,
        1000,
    )
    .value;
    let asymptotic = euler_solution("ode_asymptotic", 12, 0);
    let [re, im, error] = asymptotic[..] else {
        panic!("{asymptotic:?}");
    };
    assert!((re - 7.0 / 8.0).abs() < 1e-14 && im.abs() < 1e-14);
    assert!((error - 15.0 / 128.0).abs() < 1e-14, "{asymptotic:?}");
    assert!(0.0 < re - reference && re - reference < error);

    let borel = euler_solution("ode_borel", 12, 0);
    let [re, im, error, stokes] = borel[..] else {
        panic!("{borel:?}");
    };
    assert!((re - reference).abs() < 5e-13 && im.abs() < 1e-14);
    assert!(error.is_finite() && error >= 0.0);
    assert_eq!(stokes.to_bits(), 0.0_f64.to_bits());
}

#[test]
fn exhausted_asymptotic_series_reports_last_supplied_term() {
    let short = euler_solution("ode_asymptotic", 5, 0);
    let [re, im, error] = short[..] else {
        panic!("{short:?}");
    };
    assert!((re - 7.0 / 8.0).abs() < 1e-14 && im.abs() < 1e-14);
    assert!((error - 3.0 / 32.0).abs() < 1e-14, "{short:?}");

    // The other solution is exp(1/x) / x and has a terminating series.
    let terminating = euler_solution("ode_asymptotic", 12, 1);
    let [re, im, error] = terminating[..] else {
        panic!("{terminating:?}");
    };
    assert!((re - 4.0 * 4.0_f64.exp()).abs() < 1e-11 && im.abs() < 1e-14);
    assert_eq!(error.to_bits(), 0.0_f64.to_bits());
}
