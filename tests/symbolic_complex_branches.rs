//! Periodic roots through the public symbolic complex rule set.

use std::collections::HashMap;

use num_complex::Complex64;
use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::complex::complex;

#[test]
fn symbolic_root_reduces_exact_indices_before_numeric_evaluation() {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[complex()]).unwrap();
    let z = graph.sym("z");
    let symbol = graph.symbol_of(z).unwrap();
    let index = format!("1{}1", "0".repeat(399));
    let degree = format!("1{}", "0".repeat(400));
    let half_degree = format!("5{}", "0".repeat(399));
    let cubic = 3.0_f64.sqrt();
    for (source, radicand, want) in [
        (
            "root_branch(z,3,10^400+1)".to_string(),
            8.0,
            Complex64::new(-1.0, -cubic),
        ),
        (
            format!("root_branch(z,3,{index})"),
            8.0,
            Complex64::new(-1.0, -cubic),
        ),
        (
            format!("root_branch(z,3,-{index})"),
            8.0,
            Complex64::new(-1.0, cubic),
        ),
        (
            format!("sqrt_branch(z,-{index})"),
            -4.0,
            Complex64::new(0.0, -2.0),
        ),
        (
            format!("root_branch(z,2,-{index})"),
            -4.0,
            Complex64::new(0.0, -2.0),
        ),
        (
            format!("root_branch(z,{degree},{half_degree})"),
            8.0,
            Complex64::new(-1.0, 0.0),
        ),
    ] {
        let root = graph.parse(&source).unwrap();
        engine.run(
            &mut graph,
            &[root],
            &Env::symbolic(),
            &Saturate,
            &Budget::default(),
        );
        assert!(graph.conflicts().is_empty(), "{source}");
        let result = Extractor::new(&graph, &[root], &ClosedForm)
            .build(&mut graph, root)
            .unwrap();
        let values = HashMap::from([(symbol, Complex64::new(radicand, 0.0))]);
        let actual = graph.eval_complex(result, &values).unwrap();
        assert!(
            (actual - want).norm() < 1e-13,
            "{source}: {actual} vs {want}"
        );
    }
}
