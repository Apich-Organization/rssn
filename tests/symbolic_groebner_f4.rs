//! Public exact F4 bases when the leading coefficient changes with the field.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::rules;
use rssn::rules::poly::repr::Gens;
use rssn::rules::poly::repr::Limits;
use rssn::rules::poly::repr::from_term;

#[test]
fn f4_elimination_matches_exact_rational_and_prime_field_references() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    let x = graph.sym("x");
    let y = graph.sym("y");
    for (query, expected) in [
        (
            "groebner_f4(list(2*x+y-1,x^2+y^2-1),list(x,y))",
            "list(y^2-2/5*y-3/5,x+1/2*y-1/2)",
        ),
        (
            "groebner_f4(list(2*x+y-1,x^2+y^2-1),list(x,y),5)",
            "list(y+4,x)",
        ),
    ] {
        // Eliminating x gives 5y² - 2y - 3, whose quadratic term vanishes in GF(5).
        let root = graph.parse(query).unwrap();
        let reference = graph.parse(expected).unwrap();
        engine.run(
            &mut graph,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(graph.conflicts().is_empty(), "{query}");
        let basis = Extractor::new(&graph, &[root], &ClosedForm)
            .build(&mut graph, root)
            .unwrap();
        assert_eq!(graph.op(basis), core::LIST, "{query}");
        let actual = graph.children(basis).to_vec();
        let reference = graph.children(reference).to_vec();
        assert_eq!(actual.len(), reference.len(), "{query}");
        let mut gens = Gens::default();
        gens.index(&graph, x);
        gens.index(&graph, y);
        for (term, expected) in actual.into_iter().zip(reference) {
            let actual = from_term(&mut graph, &mut gens, term, Limits::default()).unwrap();
            let expected = from_term(&mut graph, &mut gens, expected, Limits::default()).unwrap();
            assert_eq!(actual, expected, "{query}");
        }
    }
}
