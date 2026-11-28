//! Public Taylor expansion of inverse hyperbolic functions.

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
fn inverse_hyperbolic_taylor_operators_return_exact_polynomials() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &rules::standard()).unwrap();
    let x = graph.sym("x");
    for (source, odd_coefficients) in [
        (
            "taylor(asinh(x),x,0,9)",
            [(1, 1), (-1, 6), (3, 40), (-5, 112), (35, 1152)],
        ),
        (
            "taylor(atanh(x),x,0,9)",
            [(1, 1), (1, 3), (1, 5), (1, 7), (1, 9)],
        ),
    ] {
        let root = graph.parse(source).unwrap();
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
            .unwrap_or_else(|| panic!("unresolved {source}"));
        let mut generators = Gens::default();
        let variable = generators.index(&graph, x);
        let polynomial = from_term(&mut graph, &mut generators, answer, Limits::default()).unwrap();
        assert_eq!(generators.len(), 1, "{source}: {}", graph.display(answer));
        let actual = polynomial.univariate_in(variable).unwrap();
        let expected: Vec<_> = odd_coefficients
            .into_iter()
            .flat_map(|(n, d)| [Number::from(0), Number::fraction(n, d).unwrap()])
            .collect();
        assert_eq!(actual, expected, "{source}");
    }
}
