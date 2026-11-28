//! Exact integer coefficients of symbolic cells, boundaries, and integrals.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::geometry::geometry;

const SEGMENT: &str = "cell(list(t),list(t),list(list(0,1)))";
const BIG: &str = "340282366920938463463374607431768211457";

fn reduce(expression: &str) -> Option<String> {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[geometry()]).unwrap();
    let root = graph.parse(expression).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{expression}");
    Extractor::new(&graph, &[root], &ClosedForm)
        .build(&mut graph, root)
        .map(|value| graph.display(value).to_string())
}

#[test]
fn symbolic_boundary_negates_the_minimum_machine_integer_exactly() {
    let chain = format!("chain(list(list(-9223372036854775808,{SEGMENT})))");
    assert_eq!(
        reduce(&format!("boundary({chain})")).as_deref(),
        Some(
            "chain(list(list(9223372036854775808, cell(list(0), list(), list())), list(-9223372036854775808, cell(list(1), list(), list()))))"
        )
    );
    assert_eq!(
        reduce(&format!(
            "integrate_form(form(list(list(x,list()))),list(x),boundary({chain}))"
        ))
        .as_deref(),
        Some("-9223372036854775808")
    );
}

#[test]
fn symbolic_chain_integrates_an_arbitrary_integer_weight_exactly() {
    let chain = format!("chain(list(list({BIG},{SEGMENT})))");
    let antiderivative = "form(list(list(x,list())))";
    assert_eq!(
        reduce(&format!(
            "integrate_form({antiderivative},list(x),boundary({chain}))"
        )),
        Some(BIG.to_string())
    );
    let one_form = "form(list(list(x,list(1))))";
    assert_eq!(
        reduce(&format!("integrate_form({one_form},list(x),{chain})")),
        Some(format!("{BIG}/2"))
    );
}

#[test]
fn symbolic_boundary_keeps_each_face_and_its_original_orientation_order() {
    let chain = format!("chain(list(list(0,{SEGMENT}),list(3,{SEGMENT}),list(-1,{SEGMENT})))");
    assert_eq!(
        reduce(&format!("boundary({chain})")).as_deref(),
        Some(
            "chain(list(list(-3, cell(list(0), list(), list())), list(3, cell(list(1), list(), list())), list(1, cell(list(0), list(), list())), list(-1, cell(list(1), list(), list()))))"
        )
    );
    assert_eq!(
        reduce(&format!("boundary(chain(list(list(0,{SEGMENT}))))")).as_deref(),
        Some("chain(list())")
    );
}

#[test]
fn symbolic_chain_rejects_non_integer_multiplicities() {
    for weight in ["2.0", "3/2", "y"] {
        let chain = format!("chain(list(list({weight},{SEGMENT})))");
        assert_eq!(reduce(&format!("boundary({chain})")), None);
        assert_eq!(
            reduce(&format!(
                "integrate_form(form(list(list(1,list(1)))),list(x),{chain})"
            )),
            None
        );
    }
}
