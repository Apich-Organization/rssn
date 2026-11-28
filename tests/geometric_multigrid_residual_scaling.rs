//! Scale-safe geometric residuals with exact dyadic and dense stencil references.

use num_rational::BigRational;
use num_traits::ToPrimitive;
use rssn::kernels::multigrid::GeometricMg;
use rssn::kernels::multigrid::MgOptions;

fn row_weights(
    n: usize,
    row: usize,
) -> impl Iterator<Item = i64> {
    let diagonal = (4 * (n + 1).pow(2) + 3) as i64;
    let off_diagonal = -((n + 1).pow(2) as i64);
    (0..n * n).map(move |col| {
        match (row % n).abs_diff(col % n) + (row / n).abs_diff(col / n) {
            | 0 => diagonal,
            | 1 => off_diagonal,
            | _ => 0,
        }
    })
}

fn exact_relative(
    n: usize,
    rhs: &[f64],
    solution: &[f64],
) -> f64 {
    let exact_u: Vec<_> = solution
        .iter()
        .map(|&v| BigRational::from_float(v).unwrap())
        .collect();
    let mut residual = BigRational::from_integer(0.into());
    let mut right = residual.clone();
    for (row, &rhs) in rhs.iter().enumerate() {
        let product = row_weights(n, row)
            .zip(&exact_u)
            .filter(|(coefficient, _)| *coefficient != 0)
            .map(|(coefficient, value)| value * BigRational::from_integer(coefficient.into()))
            .sum::<BigRational>();
        let rhs = BigRational::from_float(rhs).unwrap();
        let difference = &rhs - product;
        residual += &difference * &difference;
        right += &rhs * &rhs;
    }
    (residual / right).to_f64().unwrap().sqrt()
}

#[test]
fn solve_and_fmg_have_scale_invariant_residual_histories() {
    let n = 7;
    let grid = GeometricMg::new(2, n, 3.0).unwrap();
    let solution: Vec<_> = (0..grid.size())
        .map(|i| ((37 * i % 41) as f64 - 20.0) / 32.0)
        .collect();
    let rhs: Vec<f64> = (0..grid.size())
        .map(|row| {
            row_weights(n, row)
                .zip(&solution)
                .map(|(coefficient, &value)| coefficient as f64 * value)
                .sum()
        })
        .collect();
    let options = MgOptions {
        tol: 2e-11,
        max_cycles: 60,
        ..MgOptions::default()
    };
    let expected_solve = grid.solve(&rhs, None, &options);
    let expected_fmg = grid.fmg(&rhs, 1, &options);
    for exp in [-600, 600] {
        let scale = 2.0_f64.powi(exp);
        let scaled_rhs: Vec<_> = rhs.iter().map(|&v| v * scale).collect();
        let solve = grid.solve(&scaled_rhs, None, &options);
        let fmg = grid.fmg(&scaled_rhs, 1, &options);
        assert!((solve.history[0] - 1.0).abs() < 4e-16);
        assert!(solve.history.len() > 1);
        for (result, expected) in [(solve, &expected_solve), (fmg, &expected_fmg)] {
            assert!(result.converged, "{exp}: {:?}", result.history);
            assert_eq!(result.history.len(), expected.history.len());
            for (&actual, &expected) in result.history.iter().zip(&expected.history) {
                assert!(actual.is_finite());
                assert!((actual - expected).abs() < 2e-15);
            }
            for (&actual, &expected) in result.x.iter().zip(&solution) {
                assert!((actual / scale - expected).abs() < 5e-10);
            }
            let relative = exact_relative(n, &scaled_rhs, &result.x);
            assert!((result.history.last().unwrap() - relative).abs() < 2e-15);
        }
    }
}

#[test]
fn convergence_at_and_on_either_side_of_a_dyadic_tolerance() {
    let grid = GeometricMg::new(1, 1, 0.0).unwrap();
    let tol = 2.0_f64.powi(-18);
    let delta = 2.0_f64.powi(-51);
    let options = MgOptions {
        tol,
        max_cycles: 0,
        ..MgOptions::default()
    };
    for exp in [-600, 0, 600] {
        let rhs = 2.0_f64.powi(exp);
        for (ratio, converged) in [(tol - delta, true), (tol, true), (tol + delta, false)] {
            // The 1D one-point stencil is exactly [8]. Every operation
            // in this fixture is dyadic, including the residual near tol.
            let x = (rhs / 8.0) * (1.0 - ratio);
            let result = grid.solve(&[rhs], Some(&[x]), &options);
            assert_eq!(result.history, [ratio]);
            assert_eq!(result.converged, converged);
        }
    }
}

#[test]
fn fmg_stopping_matches_the_exact_residual_near_tolerance() {
    let n = 3;
    let grid = GeometricMg::new(2, n, 3.0).unwrap();
    let rhs = [1.0, 5.0, 3.0, -4.0, 2.0, -1.0, 1.0, 4.0, 6.0];
    let options = MgOptions {
        tol: 0.0,
        max_cycles: 0,
        ..MgOptions::default()
    };
    for exp in [-600, 0, 600] {
        let scale = 2.0_f64.powi(exp);
        let scaled = rhs.map(|value| value * scale);
        let initial = grid.fmg(&scaled, 1, &options);
        let exact = exact_relative(n, &scaled, &initial.x);
        assert!(exact > 1e-4);
        assert!((initial.history[0] - exact).abs() < 1e-15);
        for (tol, converged) in [(exact - 1e-13, false), (exact + 1e-13, true)] {
            let result = grid.fmg(&scaled, 1, &MgOptions { tol, ..options });
            assert_eq!(result.history, initial.history);
            assert_eq!(result.converged, converged);
        }
    }
}

#[test]
fn absolute_norm_and_finite_relative_ratio_survive_extreme_components() {
    let grid = GeometricMg::new(1, 3, 0.0).unwrap();
    for exp in [-600, 600] {
        let scale = 2.0_f64.powi(exp);
        let residual = grid.residual_norm(&[0.0; 3], &[3.0 * scale, -4.0 * scale, 0.0]);
        assert_eq!(residual, 5.0 * scale);
    }
    assert_eq!(
        grid.residual_norm(&[0.0; 3], &[f64::MAX, f64::MAX, 0.0]),
        f64::INFINITY
    );

    let grid = GeometricMg::new(1, 1, 0.0).unwrap();
    let rhs = f64::MAX;
    let x = -rhs / 8.0;
    assert_eq!(grid.residual_norm(&[x], &[rhs]), f64::INFINITY);
    let measured = grid.solve(
        &[rhs],
        Some(&[x]),
        &MgOptions {
            max_cycles: 0,
            ..MgOptions::default()
        },
    );
    assert_eq!(measured.history, [2.0]);
    assert!(!measured.converged);
}

#[test]
fn distinguishes_exact_zero_subnormal_and_undefined_residuals() {
    let grid = GeometricMg::new(1, 1, 0.0).unwrap();
    let options = MgOptions {
        tol: 0.0,
        max_cycles: 0,
        ..MgOptions::default()
    };
    let smallest = f64::from_bits(1);
    let exact_zero = grid.solve(&[0.0], None, &options);
    assert_eq!(exact_zero.history, [0.0]);
    assert!(exact_zero.converged);
    let homogeneous = grid.solve(&[0.0], Some(&[smallest]), &options);
    assert_eq!(homogeneous.history, [2.0_f64.powi(-49)]);
    assert_eq!(grid.residual_norm(&[smallest], &[0.0]), 8.0 * smallest);
    assert!(!homogeneous.converged);
    let subnormal = grid.solve(&[smallest], None, &options);
    assert_eq!(subnormal.history, [1.0]);
    assert!(!subnormal.converged);

    let lax = MgOptions {
        tol: f64::INFINITY,
        ..options
    };
    for (rhs, x) in [(0.0, f64::INFINITY), (f64::NAN, 0.0), (1.0, f64::NAN)] {
        let result = grid.solve(&[rhs], Some(&[x]), &lax);
        assert!(result.history[0].is_nan());
        assert!(!result.converged);
    }
    let overflow = grid.solve(&[smallest], Some(&[1.0]), &lax);
    assert_eq!(overflow.history, [f64::INFINITY]);
    assert!(!overflow.converged);

    // The large components of this initial iterate already solve their rows.
    // The only residual, at the opposite boundary, is the smallest positive
    // Float and would vanish if the relative residual were rounded to zero.
    let grid = GeometricMg::new(1, 7, 0.0).unwrap();
    let large = 2.0_f64.powi(1000);
    let rhs = [128.0 * large, -64.0 * large, 0.0, 0.0, 0.0, 0.0, smallest];
    let initial = [large, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let underflow = grid.solve(&rhs, Some(&initial), &options);
    assert_eq!(underflow.history, [smallest]);
    assert!(!underflow.converged);
}
