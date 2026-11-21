//! Tests the symbolic holonomic series rule through the real graph engine.

use num_bigint::BigInt;
use rssn::graph::Budget;
use rssn::graph::ClosedForm;
use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Extractor;
use rssn::graph::Graph;
use rssn::graph::Saturate;
use rssn::rules::holonomic::holonomic;

fn reduce(source: &str) -> Option<String> {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[holonomic()]).unwrap();
    let root = graph.parse(source).unwrap();
    engine.run(
        &mut graph,
        &[root],
        &Env::symbolic(),
        &Saturate,
        &Budget::default(),
    );
    assert!(graph.conflicts().is_empty(), "{source}");
    let answer = Extractor::new(&graph, &[root], &ClosedForm).build(&mut graph, root)?;
    Some(graph.display(answer))
}

fn choose(
    n: u64,
    k: u64,
) -> BigInt {
    (0..k).fold(BigInt::from(1u32), |c, j| c * (n - j) / (j + 1))
}

fn polynomial(
    terms: u64,
    mut coefficient: impl FnMut(u64) -> BigInt,
) -> String {
    (0..terms)
        .map(|i| format!("{}*x^{i}", coefficient(i)))
        .collect::<Vec<_>>()
        .join(" + ")
}

fn check(
    source: &str,
    expected: &str,
) {
    let observed = reduce(source).unwrap_or_else(|| panic!("unreduced {source}"));
    assert_eq!(
        observed,
        reduce(expected).expect("reference is a plain expression"),
        "{source}"
    );
}

#[test]
fn holonomic_series_rules_keeps_large_polynomial_powers_truncated() {
    for (exponent, terms) in [
        (8, 8),
        (64, 8),
        (256, 8),
        (4096, 6),
        (u64::from(u32::MAX), 4),
        (1 << 32, 4),
    ] {
        let src = format!("holonomic_series((1+x)^{exponent}, x, {terms})");
        check(&src, &polynomial(terms, |k| choose(exponent, k)));
    }
    check(
        "holonomic_series(((1+x)^17)^23, x, 6)",
        &polynomial(6, |k| choose(391, k)),
    );
    check("holonomic_series(x^4096, x, 8)", "0");
    check("holonomic_series(x^8, x, 8)", "0");
}

#[test]
fn holonomic_series_rules_keeps_rational_bases_and_signed_power_semantics() {
    check(
        "holonomic_series((1/(1-x))^64, x, 6)",
        &polynomial(6, |k| choose(63 + k, k)),
    );
    check(
        "holonomic_series((1-x)^(-4), x, 6)",
        &polynomial(6, |k| choose(3 + k, k)),
    );
    check("holonomic_series(x^0, x, 1)", "1");
    check("holonomic_series((1+x)^16, x, 0)", "0");
    check("holonomic_series(x, x, 4000)", "x");
}

#[test]
fn holonomic_series_rules_leaves_invalid_limits_and_nonanalytic_bases_unreduced() {
    for source in [
        "holonomic_series((1+x)^4096, x, 4001)",
        "holonomic_series((1+x)^4096, x, -1)",
        "holonomic_series((1+x)^4096, x, 7/2)",
        "holonomic_series((1+x)^4096, x, k)",
        "holonomic_series((1+x)^4096, 2, 8)",
        "holonomic_series(((1-x)/x)^100, x, 6)",
        "holonomic_series(((1-x)/x)^100, x, 0)",
        "holonomic_series(exp(x+1)^100, x, 6)",
        "holonomic_series(exp(x^4097), x, 6)",
    ] {
        assert!(reduce(source).is_none(), "{source}");
    }
}
