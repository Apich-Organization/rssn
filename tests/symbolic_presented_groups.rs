//! Deterministic Todd–Coxeter action and a verified finite group presentation.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::discrete::discrete;

#[test]
fn presented_symmetric_group_has_a_three_coset_action() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[discrete()]).unwrap();
    for (source, expected) in [
        (
            "todd_coxeter(list(a,b),list(list(a^2),list(b^3),list(a,b,a,b)),list(list(a)))",
            "list(3, list(list(1, 3, 2), list(2, 3, 1)))",
        ),
        (
            "todd_coxeter(list(a,b),list(list(1,1),list(2,2,2),list(1,2,1,2)),list(list(1)))",
            "list(3, list(list(1, 3, 2), list(2, 3, 1)))",
        ),
        (
            "todd_coxeter_index(list(a,b),list(list(a^2),list(b^3),list(a,b,a,b)),list(list(a)))",
            "3",
        ),
        (
            "group_order(fp_group(list(a,b),list(list(a^2),list(b^3),list(a,b,a,b))))",
            "6",
        ),
        (
            "group_is_valid(fp_group(list(a,b),list(list(a^2),list(b^3),list(a,b,a,b))))",
            "true",
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
        assert!(graph.conflicts().is_empty(), "{source}");
        let value = Extractor::new(&graph, &[root], &ClosedForm)
            .build(&mut graph, root)
            .unwrap_or_else(|| panic!("unresolved {source}"));
        assert_eq!(graph.display(value), expected, "{source}");
    }
}

#[test]
fn presented_order_remains_available_above_the_cayley_table_limit() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[discrete()]).unwrap();
    let order = graph
        .parse("fp_group_order(list(a),list(list(a^721)))")
        .unwrap();
    let table = graph.parse("fp_group(list(a),list(list(a^721)))").unwrap();
    let roots = [order, table];
    engine.run(
        &mut graph,
        &roots,
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty());
    let number = Extractor::new(&graph, &[order], &ClosedForm)
        .build(&mut graph, order)
        .unwrap();
    assert_eq!(graph.number_of(number).unwrap().to_i64(), Some(721));
    let request = graph.op(table);
    let rendered = Extractor::new(&graph, &[table], &ClosedForm).build(&mut graph, table);
    assert!(rendered.is_none_or(|term| graph.op(term) == request));
}
