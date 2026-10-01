//! Helpers shared by the rule sets' unit tests.

// Each domain uses a different subset.
#![allow(dead_code)]

use crate::graph::Budget;
use crate::graph::ClosedForm;
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Evaluated;
use crate::graph::Extractor;
use crate::graph::Facts;
use crate::graph::Graph;
use crate::graph::RuleSet;
use crate::graph::Saturate;
use crate::graph::SizeCost;

/// Runs `sets` on `src` symbolically and returns the best term as text,
/// together with whether it is free of heavy operators.
pub(crate) fn reduce_with(
    sets: &[RuleSet],
    src: &str,
    assume: &[(&str, Facts)],
) -> (String, bool) {
    let mut g = Graph::new();
    let engine = Engine::install(&mut g, sets).unwrap_or_else(|e| panic!("{e}"));
    for (name, facts) in assume {
        let symbol = g.interner_mut().symbol(name);
        g.assume(symbol, *facts);
    }
    let root = g.parse(src).unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
    engine.run(&mut g, &[root], &Env::symbolic(), &Saturate, &Budget::default());
    assert!(g.conflicts().is_empty(), "unsound derivation for `{src}`: {:?}", g.conflicts());
    // The invariant check is quadratic; keep it to graphs where that is
    // affordable.
    if g.len() < 3_000 {
        assert_eq!(g.validate(), Ok(()));
    }
    match Extractor::new(&g, &[root], &ClosedForm).build(&mut g, root) {
        | Some(node) => (g.display(node), true),
        | None => {
            let node = Extractor::new(&g, &[root], &SizeCost).build(&mut g, root).unwrap_or(root);
            (g.display(node), false)
        },
    }
}

/// The closed form `sets` reduce `src` to. Panics if heavy operators remain.
pub(crate) fn simplify(
    sets: &[RuleSet],
    src: &str,
) -> String {
    let (text, reduced) = reduce_with(sets, src, &[]);
    assert!(reduced, "`{src}` was not fully reduced: {text}");
    text
}

/// The numeric value `sets` give `src` under `bindings`, with its error
/// estimate. NaN if no value was found.
pub(crate) fn numeric(
    sets: &[RuleSet],
    src: &str,
    bindings: &[(&str, f64)],
    tolerance: f64,
) -> (f64, f64) {
    let mut g = Graph::new();
    let engine = Engine::install(&mut g, sets).unwrap_or_else(|e| panic!("{e}"));
    let root = g.parse(src).unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
    let mut env = Env::numeric(tolerance);
    for (name, value) in bindings {
        env.bind(g.interner_mut().symbol(name), *value);
    }
    engine.run(&mut g, &[root], &env, &Evaluated, &Budget::default());
    assert!(g.conflicts().is_empty(), "unsound derivation for `{src}`: {:?}", g.conflicts());
    g.approx(g.find(root)).map_or((f64::NAN, f64::NAN), |b| (b.mid, b.rad))
}

/// Evaluates the term `src` directly (no rules) under `bindings`.
pub(crate) fn eval(
    sets: &[RuleSet],
    src: &str,
    bindings: &[(&str, f64)],
) -> f64 {
    let mut g = Graph::new();
    assert!(Engine::install(&mut g, sets).is_ok());
    let root = g.parse(src).unwrap_or_else(|e| panic!("cannot parse `{src}`: {e}"));
    let mut env = Env::numeric(0.0);
    for (name, value) in bindings {
        env.bind(g.interner_mut().symbol(name), *value);
    }
    g.eval(root, &env).unwrap_or(f64::NAN)
}
