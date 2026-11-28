//! Exact real-root ordering when two distinct rational roots round to one float.

use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Number;
use rssn::graph::Saturate;
use rssn::rules;

#[test]
fn exact_root_of_uses_index_first_and_orders_nearby_positive_roots() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[rules::poly()]).unwrap();
    // (x - 1) * (2^60 * x - (2^60 + 1)). Both roots round to 1.0.
    for (index, numerator, denominator) in [
        (0, 1, 1),
        (1, 1_152_921_504_606_846_977, 1_152_921_504_606_846_976),
    ] {
        let query = format!(
            "root_of({index},1152921504606846977,-2305843009213693953,1152921504606846976)"
        );
        let term = graph.parse(&query).unwrap();
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
        let actual = graph.number_of(root).and_then(Number::to_rational);
        let expected =
            Number::fraction(numerator, denominator).and_then(|value| value.to_rational());
        assert_eq!(actual, expected, "{query}: {}", graph.display(root));
    }
}
