//! Frobenius metadata and numerical selection of conjugate embeddings.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Number;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::rules;

#[test]
fn bessel_half_order_preserves_both_embeddings_and_exponents() {
    // x²y'' + xy' + (x² - 1/4)y = 0 has branches exp(±ix) / sqrt(x),
    // both represented over the same field Q(i).
    let equation = "x^2*diff(diff(y(x),x),x)+x*diff(y(x),x)+(x^2-1/4)*y(x)";
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    let roots = [
        graph
            .parse(&format!("ode_frobenius({equation},y(x),oo,6)"))
            .unwrap(),
        graph
            .parse(&format!("ode_asymptotic({equation},y(x),oo,6,1/2,0)"))
            .unwrap(),
        graph
            .parse(&format!("ode_asymptotic({equation},y(x),oo,6,1/2,1)"))
            .unwrap(),
    ];
    engine.run(
        &mut graph,
        &roots,
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
    let [report, positive, negative] = roots.map(|root| {
        Extractor::new(&graph, &[root], &ClosedForm)
            .build(&mut graph, root)
            .unwrap()
    });

    assert_eq!(graph.op(report), core::LIST);
    let &[solution, exponents, recurrence] = graph.children(report) else {
        panic!("{}", graph.display(report));
    };
    assert_eq!(graph.op(solution), core::EQ);
    assert_eq!(graph.op(recurrence), core::LIST);
    assert!(graph.children(recurrence).is_empty());
    let &[first, second] = graph.children(exponents) else {
        panic!("{}", graph.display(exponents));
    };
    let half = Number::fraction(1, 2).unwrap();
    let exponent_parts = [first, second].map(|pair| {
        assert_eq!(graph.op(pair), core::LIST);
        let &[exponent, rho] = graph.children(pair) else {
            panic!("{}", graph.display(pair));
        };
        assert_eq!(graph.number_of(rho), Some(&half));
        graph.find(exponent)
    });
    assert_ne!(exponent_parts[0], exponent_parts[1]);

    let x = 0.5_f64;
    for (value, sign) in [(positive, 1.0), (negative, -1.0)] {
        assert_eq!(graph.op(value), core::LIST);
        let &[re, im, error] = graph.children(value) else {
            panic!("{}", graph.display(value));
        };
        let number = |node| graph.number_of(node).unwrap().to_f64();
        assert!((number(re) - x.cos() / x.sqrt()).abs() < 1e-13);
        assert!((number(im) - sign * x.sin() / x.sqrt()).abs() < 1e-13);
        assert_eq!(number(error).to_bits(), 0.0_f64.to_bits());
    }
}
