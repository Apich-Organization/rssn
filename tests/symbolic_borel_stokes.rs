//! Public Borel–Padé median and lateral jump for a simple Stokes pole.

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

#[test]
fn borel_stokes_median_agrees_with_principal_value_and_residue() {
    // Differentiating -x²y' + (1 - x)y = 1 gives this equation.
    // Its regular formal branch has coefficients k! and a Borel pole on the ray.
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    let root = graph
        .parse("ode_borel(x^2*diff(diff(y(x),x),x)+(3*x-1)*diff(y(x),x)+y(x),y(x),0,60,1/20,0)")
        .unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
    let result = Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .unwrap();
    assert_eq!(graph.op(result), core::LIST);
    let &[re, im, error, stokes] = graph.children(result) else {
        panic!("{}", graph.display(result));
    };
    let number = |node| graph.number_of(node).unwrap().to_f64();
    assert_eq!(graph.number_of(stokes).unwrap().to_i64(), Some(1));

    // In the principal value integral exp(-t)/(1 - t/a), pair points at
    // t = a ± u for 0 < u < a to remove the pole before quadrature.
    let a = 20.0_f64;
    let paired = gauss_kronrod_any(
        |u| a * (u - a).exp() * -(-2.0 * u).exp_m1() / u,
        0.0,
        a,
        1e-13,
        1000,
    );
    let tail = gauss_kronrod_any(|u| a * (-a - u).exp() / u, a, f64::INFINITY, 1e-13, 1000);
    let median = paired.value - tail.value;
    let half_jump = std::f64::consts::PI * a * (-a).exp();
    assert!(
        (number(re) - median).abs() < 2e-8,
        "{}",
        graph.display(result)
    );
    assert!(number(im).abs() < 1e-13, "{}", graph.display(result));
    assert!(
        (number(error) - half_jump).abs() < 2e-8,
        "{}",
        graph.display(result)
    );
}
