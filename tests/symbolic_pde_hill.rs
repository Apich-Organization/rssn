//! Periodic Fourier-Galerkin eigensystems exposed by the PDE rule set.

use rssn::graph::Engine;
use rssn::graph::Env;
use rssn::graph::Graph;
use std::f64::consts::TAU;

fn graph() -> Graph {
    let mut graph = Graph::new();
    Engine::install(&mut graph, &[rssn::rules::pde()]).unwrap();
    graph
}

fn evaluate(
    graph: &mut Graph,
    name: &str,
    fixed: &[f64],
    potentials: &[f64],
) -> f64 {
    let op = graph.ops().lookup(name).unwrap();
    let args: Vec<_> = fixed
        .iter()
        .chain(potentials)
        .map(|&value| graph.float(value))
        .collect();
    let term = graph.node(op, &args);
    graph.eval(term, &Env::numeric(0.0)).unwrap_or(f64::NAN)
}

#[test]
fn hill_constant_potential_preserves_extreme_periods_and_exact_multiplicity() {
    let mut graph = graph();
    let potential = [2.0; 41];
    let huge = 1e308;
    let lowest = evaluate(&mut graph, "hill_eigen", &[1.0, huge, 16.0], &potential);
    assert_eq!(lowest, 2.0);
    let constant = evaluate(
        &mut graph,
        "hill_mode",
        &[1.0, huge * 0.75, huge, 16.0],
        &potential,
    );
    assert!((constant * huge.sqrt() - 1.0).abs() < 1e-14);

    for (index, expected) in [(1.0, 2.0), (2.0, 3.0), (3.0, 3.0), (4.0, 6.0), (5.0, 6.0)] {
        assert!(
            (evaluate(&mut graph, "hill_eigen", &[index, TAU, 16.0], &potential) - expected).abs()
                < 1e-14
        );
    }

    // At the Nyquist frequency the eight sampled sines vanish, and the
    // sampled cosine has twice the usual mean square in the real basis.
    let nyquist = [2.0; 8];
    assert_eq!(
        evaluate(&mut graph, "hill_eigen", &[8.0, TAU, 4.0], &nyquist),
        16.0
    );
    assert_eq!(
        evaluate(&mut graph, "hill_eigen", &[9.0, TAU, 4.0], &nyquist),
        20.0
    );
}

#[test]
fn hill_handles_nontrivial_fourier_modes_and_odd_sample_counts() {
    let mut graph = graph();
    let count = 257;
    let potential: Vec<f64> = (0..count)
        .map(|i| 2.0 * (2.0 * TAU * i as f64 / count as f64).cos())
        .collect();
    let expected = [-0.455_138_604_1, -0.110_248_816_0, 1.859_108_072_5];
    for (index, expected) in expected.into_iter().enumerate() {
        let value = evaluate(
            &mut graph,
            "hill_eigen",
            &[(index + 1) as f64, TAU, 24.0],
            &potential,
        );
        assert!((value - expected).abs() < 1e-7);
    }
    let mut norm = 0.0;
    let mut overlap = 0.0;
    for i in 0..97 {
        let x = TAU * i as f64 / 97.0;
        let mode = evaluate(&mut graph, "hill_mode", &[1.0, x, TAU, 24.0], &potential);
        let other = evaluate(&mut graph, "hill_mode", &[2.0, x, TAU, 24.0], &potential);
        let repeated = evaluate(
            &mut graph,
            "hill_mode",
            &[1.0, x - TAU, TAU, 24.0],
            &potential,
        );
        assert!((mode - repeated).abs() < 2e-14);
        norm += mode * mode;
        overlap += mode * other;
    }
    assert!((norm * TAU / 97.0 - 1.0).abs() < 1e-12);
    assert!((overlap * TAU / 97.0).abs() < 1e-12);
}

#[test]
fn hill_rejects_invalid_periods_harmonics_and_samples() {
    let mut graph = graph();
    let potential = [1.0; 8];
    for fixed in [
        [1.0, 0.0, 4.0],
        [1.0, f64::INFINITY, 4.0],
        [1.0, TAU, 2.75],
        [1.0, TAU, 65.0],
        [10.0, TAU, 4.0],
    ] {
        assert!(evaluate(&mut graph, "hill_eigen", &fixed, &potential).is_nan());
    }
    assert!(evaluate(&mut graph, "hill_eigen", &[1.0, TAU, 4.0], &potential[..7]).is_nan());
    let mut invalid = potential;
    invalid[7] = f64::NAN;
    assert!(evaluate(&mut graph, "hill_eigen", &[1.0, TAU, 4.0], &invalid).is_nan());
    assert!(
        evaluate(
            &mut graph,
            "hill_mode",
            &[1.0, f64::NAN, TAU, 4.0],
            &potential
        )
        .is_nan()
    );
}
