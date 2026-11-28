//! Eigenring factorization through the public symbolic ODE entry point.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules;

#[test]
fn eigenring_recovers_two_independent_airy_equations() {
    // A common left multiple of D² - x and D² - 2x, checked by differentiating each equation twice.
    let equation = concat!(
        "dsolve(2*x^4*y(x)+(2-3*x^3)*diff(diff(y(x),x),x)",
        "-2*x*diff(diff(diff(y(x),x),x),x)",
        "+x^2*diff(diff(diff(diff(y(x),x),x),x),x),y(x))"
    );
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    let root = graph.parse(equation).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
    let answer = Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .unwrap();
    let answer = graph.display(answer);
    assert!(answer.contains("airy"), "{answer}");
    for coefficient in ["C1", "C2", "C3", "C4"] {
        assert!(
            answer.contains(coefficient),
            "missing {coefficient} in {answer}"
        );
    }
}
