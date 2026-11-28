//! X-fastest layout and boundary regressions for the fixed square-grid V-cycle.

use rayon::ThreadPoolBuilder;
use rssn::sim::physics_mtm::solve_poisson_2d_multigrid;
use std::f64::consts::PI;

#[test]
fn validates_shape_before_allocating_or_iterating() {
    assert!(solve_poisson_2d_multigrid(1, &[0.0], 1).is_err());
    assert!(solve_poisson_2d_multigrid(4, &[0.0; 16], 1).is_err());
    assert!(solve_poisson_2d_multigrid(3, &[0.0; 8], 1).is_err());
    assert!(solve_poisson_2d_multigrid(3, &[0.0; 10], 0).is_err());

    let oversized_dyadic_side = (1usize << (usize::BITS - 1)) + 1;
    assert!(solve_poisson_2d_multigrid(oversized_dyadic_side, &[], 0).is_err());
}

#[test]
fn boundary_only_and_zero_cycles_stay_zero() {
    for cycles in [0, 1, usize::MAX] {
        let u = solve_poisson_2d_multigrid(2, &[f64::NAN; 4], cycles).unwrap();
        assert!(u.iter().all(|value| value.to_bits() == 0));
    }
    let u = solve_poisson_2d_multigrid(9, &[f64::NAN; 81], 0).unwrap();
    assert!(u.iter().all(|value| value.to_bits() == 0));
}

#[test]
fn coarsest_interior_matches_the_single_discrete_equation() {
    let mut rhs = [f64::NAN; 9];
    rhs[4] = -24.0;
    let u = solve_poisson_2d_multigrid(3, &rhs, 1).unwrap();

    // The one unknown obeys 4 * u / h^2 = f, with h = 1/2.
    for (index, value) in u.into_iter().enumerate() {
        let expected: f64 = if index == 4 { -1.5 } else { 0.0 };
        assert_eq!(value.to_bits(), expected.to_bits());
    }
}

#[test]
fn asymmetric_discrete_eigenmode_has_the_expected_layout_and_converges() {
    let n = 65;
    let h = 1.0 / (n - 1) as f64;
    let mut exact = vec![0.0; n * n];
    let mut rhs = vec![f64::NAN; n * n];
    // For x-fastest indexing, this mode oscillates twice along x, three times
    // along y. Its five-point eigenvalue uses discrete angles, so the reference
    // has no spatial discretization error.
    let lambda = 4.0 * ((PI * h).sin().powi(2) + (1.5 * PI * h).sin().powi(2)) / h.powi(2);
    for y in 1..n - 1 {
        for x in 1..n - 1 {
            let k = y * n + x;
            exact[k] = (2.0 * PI * h * x as f64).sin() * (3.0 * PI * h * y as f64).sin();
            rhs[k] = lambda * exact[k];
        }
    }
    let check = |cycles| {
        let u = solve_poisson_2d_multigrid(n, &rhs, cycles).unwrap();
        assert!(u.iter().all(|value| value.is_finite()));
        u.iter()
            .zip(&exact)
            .map(|(u, target)| (u - target).abs())
            .fold(0.0, f64::max)
    };
    let one = check(1);
    let four = check(4);
    let ten = check(10);
    assert!(one > four && four > ten, "{one} {four} {ten}");
    assert!(ten < 1e-8, "{ten}");
}

#[test]
fn active_pool_size_does_not_change_any_node() {
    let n = 513;
    let mut rhs: Vec<_> = (0..n * n)
        .map(|i| ((i % 123) as f64 - 60.0) * 0.125)
        .collect();
    rhs[0] = f64::INFINITY;
    let solve = |threads| {
        ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| solve_poisson_2d_multigrid(n, &rhs, 3).unwrap())
    };
    let serial = solve(1);
    let parallel = solve(2);
    assert!(serial.iter().all(|value| value.is_finite()));
    assert!(
        serial
            .iter()
            .zip(parallel)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}
