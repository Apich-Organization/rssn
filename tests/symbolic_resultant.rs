//! Exact scalar and parameterized Sylvester resultants through the public graph.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Number;
use rssn::graph::Saturate;
use rssn::rules;
use rssn::rules::poly::repr::Gens;
use rssn::rules::poly::repr::Limits;
use rssn::rules::poly::repr::from_term;

#[test]
fn resultant_uses_ascending_coefficients_and_preserves_symbolic_coefficients() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    for (query, expected) in [
        ("resultant(6*x^2+13*x-5,4*x+3,x)", -182),
        ("resultant(6*x^2+13*x-5,3*x-1,x)", 0),
        ("discriminant(6*x^2+13*x-5,x)", 289),
    ] {
        let term = graph.parse(query).unwrap();
        engine.run(
            &mut graph,
            &[term],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        let root = Extractor::new(&graph, &[term], &ClosedForm)
            .build(&mut graph, term)
            .unwrap();
        assert_eq!(
            graph.number_of(root).and_then(Number::to_i64),
            Some(expected),
            "{query}"
        );
    }

    let query = "resultant(x^9+y,x^8+1,x)";
    let term = graph.parse(query).unwrap();
    let expected = graph.parse("1+y^8").unwrap();
    engine.run(
        &mut graph,
        &[term],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    let root = Extractor::new(&graph, &[term], &ClosedForm)
        .build(&mut graph, term)
        .unwrap();
    let mut generators = Gens::default();
    let actual = from_term(&mut graph, &mut generators, root, Limits::default()).unwrap();
    let expected = from_term(&mut graph, &mut generators, expected, Limits::default()).unwrap();
    assert_eq!(actual, expected, "{query}");
}

#[test]
fn resultant_rejects_coefficients_that_depend_on_the_eliminated_variable() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    for query in ["resultant(x+sin(x),x-1,x)", "discriminant(x^2+sin(x),x)"] {
        let term = graph.parse(query).unwrap();
        let request = graph.op(term);
        engine.run(
            &mut graph,
            &[term],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        if let Some(root) = Extractor::new(&graph, &[term], &ClosedForm).build(&mut graph, term) {
            assert_eq!(graph.op(root), request, "{query}: {}", graph.display(root));
        }
    }
}
