//! Dimension rules of the `ode` and `pde` rule sets through the public
//! session API.
//!
//! The checks are made in `Strict` mode, where a conflicting dimension is a
//! hard error, and repeated in `Advisory` mode to show that the same terms are
//! accepted there. Equations, lists of symmetries and classification reports
//! are non-metric values: their display is `non-metric`.

use rssn::graph::DimensionMode;
use rssn::prelude::*;

/// Computes `src` with explicit dimension assumptions and returns the printed
/// dimension together with the blocking (real conflict) diagnostics.
fn compute(
    src: &str,
    mode: DimensionMode,
    assumptions: &[(&str, &str)],
) -> (String, Vec<String>) {
    let session = Session::new();
    let term = session
        .parse(src)
        .unwrap_or_else(|error| panic!("cannot parse `{src}`: {error}"));
    let mut config = Config::new().dimensions(mode);
    for &(name, text) in assumptions {
        config = config.dimension(name, text);
    }
    let answer = session
        .compute(term, &config)
        .unwrap_or_else(|error| panic!("`{src}` ({mode:?}) failed: {error}"));
    let dimension = answer
        .dimension
        .as_ref()
        .map(ToString::to_string)
        .unwrap_or_default();
    let blocking = answer
        .diagnostics
        .iter()
        .filter(|error| error.is_blocking())
        .map(ToString::to_string)
        .collect();
    (dimension, blocking)
}

#[test]
fn dsolve_returns_a_non_metric_equation() {
    // First order: `y' = a y`; second order: `y'' = -y`. Both answers are the
    // equation object `y(x) = ...`, so their dimension is non-metric.
    for src in [
        "dsolve(diff(y(x), x) = a*y(x), y(x))",
        "dsolve(diff(diff(y(x), x), x) = -y(x), y(x))",
    ] {
        for mode in [DimensionMode::Strict, DimensionMode::Advisory] {
            let (dimension, blocking) = compute(src, mode, &[("a", "1"), ("x", "1")]);
            assert_eq!(dimension, "non-metric", "{src} ({mode:?})");
            assert!(blocking.is_empty(), "{src} ({mode:?}): {blocking:?}");
        }
    }
}

#[test]
fn an_equation_object_stays_non_metric() {
    for mode in [DimensionMode::Strict, DimensionMode::Advisory] {
        let (dimension, blocking) =
            compute("diff(y(x), x) = a*y(x)", mode, &[("a", "1"), ("x", "1")]);
        assert_eq!(dimension, "non-metric", "({mode:?})");
        assert!(
            blocking.is_empty(),
            "equation object ({mode:?}): {blocking:?}"
        );
    }
}

#[test]
fn pde_solutions_are_equations() {
    // Heat, wave and Laplace, all in dimensionless coordinates so that the
    // transport coefficients (set to one) are consistent.
    let heat = "pdsolve(diff(u(x, t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(1, t) = 0, u(x, 0) = sin(pi*x)))";
    let wave = "pdsolve(diff(diff(u(x, t), t), t) = diff(diff(u(x, t), x), x), u(x, t), list(u(0, t) = 0, u(pi, t) = 0, u(x, 0) = sin(x), at(diff(u(x, t), t), t, 0) = 0))";
    let laplace = "pdsolve(diff(diff(u(x, y), x), x) + diff(diff(u(x, y), y), y) = 0, u(x, y), list(u(0, y) = 0, u(1, y) = 0, u(x, 0) = 0, u(x, 1) = sin(pi*x)))";
    let assumptions = [("x", "1"), ("y", "1"), ("t", "1"), ("u", "1")];
    for (name, src) in [("heat", heat), ("wave", wave), ("laplace", laplace)] {
        for mode in [DimensionMode::Strict, DimensionMode::Advisory] {
            let (dimension, blocking) = compute(src, mode, &assumptions);
            assert_eq!(dimension, "non-metric", "{name} ({mode:?})");
            assert!(blocking.is_empty(), "{name} ({mode:?}): {blocking:?}");
        }
    }
}

#[test]
fn eigenfunction_and_greens_function_dimensions() {
    // The first Dirichlet eigenvalue of `-y''` on [0, 1] is `pi^2`: a pure
    // number, as the numeric Sturm-Liouville operator returns.
    let (dimension, blocking) = compute(
        "sl_eigen(1, 0, 1, 0, 1, 0, 1, 4, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1)",
        DimensionMode::Strict,
        &[],
    );
    assert_eq!(dimension, "1");
    assert!(blocking.is_empty(), "sl_eigen: {blocking:?}");

    // The Dirichlet Green's function of the ball of radius 2 m is
    // `1/(4 pi |x - s|) - 1/(4 pi rho)`, hence `1/Length`.
    let (dimension, blocking) = compute(
        "pde_greens_ball(quantity(2, m), 3)",
        DimensionMode::Strict,
        &[
            ("x", "m"),
            ("y", "m"),
            ("z", "m"),
            ("s_x", "m"),
            ("s_y", "m"),
            ("s_z", "m"),
        ],
    );
    assert_eq!(dimension, "1/Length");
    assert!(blocking.is_empty(), "pde_greens_ball: {blocking:?}");
}
