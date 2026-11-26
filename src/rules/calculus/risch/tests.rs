//! Tests of the Risch–Bronstein tower integrator.
//!
//! Each case is run through the tower entry point directly, so a failure
//! points at the algorithm rather than at the heuristic search around it,
//! and the answer is checked by differentiating it back.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::panic,
    clippy::pedantic
)]

use super::Differentiate;
use super::Functions;
use super::integrate as tower;
use crate::graph::Cx;
use crate::graph::Engine;
use crate::graph::Env;
use crate::graph::Graph;
use crate::rules::calculus;
use crate::rules::testing::reduce_with;

/// Runs the tower integrator alone and checks the answer numerically.
fn check_tower(src: &str) {
    check_tower_at(src, &[0.37, 0.91, 1.63, 2.3]);
}

/// As [`check_tower`], sampling at caller-chosen points (a deep logarithm
/// tower is only real above `e^e`, `e^(e^e)`, ...).
fn check_tower_at(
    src: &str,
    points: &[f64],
) {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[calculus()]).unwrap();
    let integrand = graph.parse(src).unwrap();
    let x = graph.parse("x").unwrap();
    let symbol = graph.symbol_of(x).unwrap();
    let integral_op = graph.ops().lookup("integral").unwrap();
    let functions = *graph.ops().attr::<Functions>(integral_op).unwrap();
    let env = Env::symbolic();
    let answer = {
        let mut cx = Cx {
            graph: &mut graph,
            env: &env,
            engine: &engine,
        };
        tower(&mut cx, functions, integrand, x)
            .unwrap_or_else(|| panic!("the tower found no antiderivative of {src}"))
    };
    let derivative = Differentiate { diff: functions.diff }.derive(&mut graph, answer, symbol, x);
    let mut checked = 0;
    for &point in points {
        let mut values = Env::numeric(0.0);
        values.bind(symbol, point);
        let (Some(want), Some(got)) = (
            graph.eval(integrand, &values),
            graph.eval(derivative, &values),
        ) else {
            continue;
        };
        if want.is_finite() && got.is_finite() {
            assert!(
                (want - got).abs() <= 1e-6 * (1.0 + want.abs().max(got.abs())),
                "{src}: the derivative is {got}, not {want}"
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "{src}: no sample point was usable");
}

/// The tower declines an integrand that is not elementary in the tower.
fn declines(src: &str) {
    let mut graph = Graph::new();
    let engine = Engine::install(&mut graph, &[calculus()]).unwrap();
    let integrand = graph.parse(src).unwrap();
    let x = graph.parse("x").unwrap();
    let integral_op = graph.ops().lookup("integral").unwrap();
    let functions = *graph.ops().attr::<Functions>(integral_op).unwrap();
    let env = Env::symbolic();
    let mut cx = Cx {
        graph: &mut graph,
        env: &env,
        engine: &engine,
    };
    assert!(
        tower(&mut cx, functions, integrand, x).is_none(),
        "{src} was decided although it should not be"
    );
}

#[test]
fn rational_functions() {
    check_tower("x^3");
    check_tower("1/x");
    check_tower("1/(x + 1)");
    check_tower("(2*x + 1)/(x^2 + x + 1)");
    check_tower("2*x/(x^2 + 1)");
    check_tower("x/(x^2 + 1)^2");
    // Irrational residues: the logarithmic part needs an algebraic
    // extension, so the module declines instead of inventing a form.
    declines("1/(x^2 + 1)");
}

#[test]
fn exponential_extension() {
    check_tower("exp(x)");
    check_tower("1/exp(x)");
    check_tower("exp(-2*x)");
    check_tower("x^2*exp(x)");
    check_tower("1/(exp(x) + 1)");
    check_tower("1/(exp(x) - 1)");
    check_tower("exp(x)/(exp(x) + 1)^2");
    // `exp(x)/x` is the exponential integral: not elementary.
    declines("exp(x)/x");
}

#[test]
fn logarithmic_extension() {
    check_tower("ln(x)");
    check_tower("ln(x)/x");
    check_tower("ln(x)^3");
    check_tower("x*ln(x)^2");
    check_tower("1/(x*ln(x))");
    // `1/ln(x)` is the logarithmic integral: not elementary.
    declines("1/ln(x)");
}

#[test]
fn nested_towers() {
    // ln(ln(x)): the tower is Q(x) < Q(x)(log x) < Q(x)(log x)(log log x).
    check_tower("ln(ln(x))/(x*ln(x))");
    check_tower("1/(x*ln(x)*ln(ln(x)))");
    // exp(x)*ln(exp(x) + 1) is elementary in a three-level tower.
    check_tower("exp(x)*ln(exp(x) + 1)");
    // Four nested logarithms: the tower is Q(x) < log x < log log x <
    // log log log x < log log log log x, and the antiderivative needs the
    // next level again. The integrand is only real above e^(e^e).
    check_tower_at(
        "1/(x*ln(x)*ln(ln(x))*ln(ln(ln(x)))*ln(ln(ln(ln(x)))))",
        &[20.0, 25.0, 30.0],
    );
}

#[test]
fn the_engine_answers_nested_integrals() {
    // The heuristics around the tower cannot reach these; the answer must
    // come back reduced.
    for src in [
        "integral(ln(ln(x))/(x*ln(x)), x)",
        "integral(1/(x*ln(x)*ln(ln(x))), x)",
        "integral(exp(x)*ln(exp(x) + 1), x)",
    ] {
        let (answer, reduced) = reduce_with(&[calculus()], src, &[]);
        assert!(reduced, "{src} was not integrated: {answer}");
    }
}
