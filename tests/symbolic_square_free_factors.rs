//! Square-free factors through the symbolic graph and the typed kernel.

use rssn::Session;
use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Number;
use rssn::graph::Saturate;
use rssn::graph::op::core;
use rssn::kernel::CallCtx;
use rssn::kernel::TypeExpr;
use rssn::kernel::Value;
use rssn::kernel::infer::term;
use rssn::rules;
use rssn::rules::poly::repr::Gens;
use rssn::rules::poly::repr::Limits;
use rssn::rules::poly::repr::Poly;
use rssn::rules::poly::repr::from_term;

#[test]
fn square_free_factors_group_equal_multiplicities_and_keep_the_leading_coefficient() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rules::poly()]).unwrap();
    let x = graph.sym("x");
    let mut generators = Gens::default();
    generators.index(&graph, x);
    let to_poly = |graph: &mut Graph, generators: &mut Gens, node| {
        from_term(graph, generators, node, Limits::default()).unwrap()
    };
    for (expression, unit, expected) in [
        ("-6*x^2*(x^2+1)^3", (-6, 1), &[("x", 2), ("x^2+1", 3)][..]),
        (
            "(2*x+1)^4*(3*x-1)^4/5",
            (1296, 5),
            &[("x^2+x/6-1/6", 4)][..],
        ),
        ("7/3", (7, 3), &[][..]),
        ("x^24", (1, 1), &[("x", 24)][..]),
    ] {
        let query = format!("square_free_factors({expression},x)");
        let source = graph.parse(&query).unwrap();
        let original = graph.parse(expression).unwrap();
        engine.run(
            &mut graph,
            &[source],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(graph.conflicts().is_empty(), "{query}");
        let root = Extractor::new(&graph, &[source], &ClosedForm)
            .build(&mut graph, source)
            .unwrap();
        assert_eq!(graph.op(root), core::LIST, "{query}");
        let &[actual_unit, factors] = graph.children(root) else {
            panic!("{query}")
        };
        let actual_unit = graph.number_of(actual_unit).unwrap();
        assert_eq!(
            Some(actual_unit.clone()),
            Number::fraction(unit.0, unit.1),
            "{query}"
        );
        let mut reconstructed = Poly::constant(actual_unit.clone());
        assert_eq!(graph.op(factors), core::LIST, "{query}");
        let rows = graph.children(factors).to_vec();
        assert_eq!(rows.len(), expected.len(), "{query}");
        for (row, &(factor, multiplicity)) in rows.into_iter().zip(expected) {
            assert_eq!(graph.op(row), core::LIST, "{query}");
            let &[actual_factor, actual_multiplicity] = graph.children(row) else {
                panic!("{query}")
            };
            assert_eq!(
                graph
                    .number_of(actual_multiplicity)
                    .and_then(Number::to_i64),
                Some(i64::from(multiplicity)),
                "{query}"
            );
            let factor = graph.parse(factor).unwrap();
            let actual = to_poly(&mut graph, &mut generators, actual_factor);
            assert_eq!(
                actual,
                to_poly(&mut graph, &mut generators, factor),
                "{query}"
            );
            reconstructed = reconstructed
                .mul(&actual.pow(multiplicity, 1024).unwrap(), 1024)
                .unwrap();
        }
        assert_eq!(
            reconstructed,
            to_poly(&mut graph, &mut generators, original),
            "{query}"
        );
    }
}

#[test]
fn square_free_factors_leaves_unsupported_requests_unevaluated() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rules::poly()]).unwrap();
    for query in [
        "square_free_factors(0,x)",
        "square_free_factors(1/(x+1),x)",
        "square_free_factors((x-a)^2,x)",
        "square_free_factors(x^2,x+1)",
    ] {
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

#[test]
fn square_free_factors_declares_and_dispatches_the_symbol_overload() {
    let session = Session::with_rules(&[rules::poly()]).unwrap();
    let overloads = session.overloads("square_free_factors");
    let [(key, signature)] = overloads.as_slice() else {
        panic!("{overloads:?}")
    };
    assert_eq!(signature.params, vec![term(), TypeExpr::Symbol]);
    assert_eq!(signature.ret, term());
    let args =
        ["-6*x^2*(x^2+1)^3", "x"].map(|input| session.parse(input).unwrap().to_value().unwrap());
    assert_eq!(
        session
            .resolve_overload("square_free_factors", &args)
            .unwrap(),
        *key
    );
    let Value::Term(handle) = session.call_key(key, &args, &CallCtx::new()).unwrap() else {
        panic!("expected term")
    };
    let answer = session.term_ref(handle).unwrap().to_string();
    assert_eq!(answer, "list(-6, list(list(x, 2), list(x^2 + 1, 3)))");
}
