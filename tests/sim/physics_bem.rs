//! Boundary element method for the 2D Laplace equation (ported from
//! `physics_bem_test.rs`).

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_bem::*;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 32,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

/// Unit square, counter-clockwise from the origin, with `n` nodes per side.
/// Left side u = 100, right side u = 0, top and bottom insulated.
fn conduction_problem(n: usize) -> (Vec<(f64, f64)>, Vec<BoundaryCondition<f64>>) {
    let mut points = Vec::new();
    let mut bcs = Vec::new();
    for i in 0..n {
        points.push((i as f64 / n as f64, 0.0));
        bcs.push(BoundaryCondition::Flux(0.0));
    }
    for i in 0..n {
        points.push((1.0, i as f64 / n as f64));
        bcs.push(BoundaryCondition::Potential(0.0));
    }
    for i in 0..n {
        points.push((1.0 - i as f64 / n as f64, 1.0));
        bcs.push(BoundaryCondition::Flux(0.0));
    }
    for i in 0..n {
        points.push((0.0, 1.0 - i as f64 / n as f64));
        bcs.push(BoundaryCondition::Potential(100.0));
    }
    (points, bcs)
}

fn elements_of(points: &[(f64, f64)]) -> Vec<Element2D> {
    let n = points.len();
    (0..n)
        .map(|i| {
            Element2D::new(
                Vector2D::new(points[i].0, points[i].1),
                Vector2D::new(points[(i + 1) % n].0, points[(i + 1) % n].1),
            )
        })
        .collect()
}

#[test]
fn vector_helpers_and_element_geometry() {
    let v = Vector2D::new(3.0, 4.0);
    assert_eq!(v.norm(), 5.0);
    let w = (v + Vector2D::new(1.0, 1.0)) - Vector2D::new(0.0, 2.0);
    assert_eq!((w.x, w.y), (4.0, 3.0));
    assert_eq!((v * 2.0).x, 6.0);
    let e = Element2D::new(Vector2D::new(0.0, 0.0), Vector2D::new(2.0, 0.0));
    assert_eq!(e.length, 2.0);
    assert_eq!((e.midpoint.x, e.midpoint.y), (1.0, 0.0));
    // Counter-clockwise traversal => the normal points to the right of the direction of travel, i.e. outward.
    assert_eq!((e.normal.x, e.normal.y), (0.0, -1.0));
    let d = Element2D::new(Vector2D::new(1.0, 1.0), Vector2D::new(1.0, 3.0));
    assert_eq!((d.normal.x, d.normal.y), (1.0, 0.0));
    assert!((Vector3D::new(2.0, 3.0, 6.0).norm() - 7.0).abs() < 1e-12);
}

#[test]
fn rectangle_conduction_matches_the_linear_profile() {
    let n = 10;
    let (points, bcs) = conduction_problem(n);
    let (u, q) = solve_laplace_bem_2d(&points, &bcs).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!((u.len(), q.len()), (4 * n, 4 * n));
    // Prescribed values are returned untouched.
    assert_approx_eq!(u[n + n / 2], 0.0, 1e-12);
    assert_approx_eq!(u[3 * n + n / 2], 100.0, 1e-12);
    assert_approx_eq!(q[n / 2], 0.0, 1e-12);
    // Exact solution u = 100 (1 - x): the insulated sides read 100 (1 - x_mid) at their elements.
    for i in 0..n {
        let x_mid = (i as f64 + 0.5) / n as f64;
        assert!((u[i] - 100.0 * (1.0 - x_mid)).abs() < 5.0, "bottom node {i}: {} vs {}", u[i], 100.0 * (1.0 - x_mid));
    }
    // The flux out of the cold side is u_x = -100 (pointing outward: q = du/dn = -du/dx = 100 at x = 1 ... sign per convention).
    let cold: f64 = q[n..2 * n].iter().sum::<f64>() / n as f64;
    assert!(cold.abs() > 50.0 && cold.abs() < 150.0, "mean flux on the cold side {cold}");
    // Global balance: net flux through the boundary vanishes (harmonic function).
    let net: f64 = q.iter().enumerate().map(|(i, v)| v * elements_of(&points)[i].length).sum();
    assert!(net.abs() < 5.0, "net flux {net}");
}

#[test]
fn evaluated_interior_potential_follows_the_linear_profile() {
    let (points, bcs) = conduction_problem(10);
    let (u, q) = solve_laplace_bem_2d(&points, &bcs).unwrap_or_else(|e| panic!("{e}"));
    let elements = elements_of(&points);
    for (x, want) in [(0.25, 75.0), (0.5, 50.0), (0.75, 25.0)] {
        let pot = evaluate_potential_2d((x, 0.5), &elements, &u, &q);
        assert!((pot - want).abs() < 6.0, "u({x}, 0.5) = {pot}, expected {want}");
    }
}

#[test]
fn cylinder_scenario_reproduces_the_harmonic_function_x() {
    // u = x on the unit circle has the exact interior solution u = x and boundary flux du/dn = x.
    let (u, q) = simulate_2d_cylinder_scenario().unwrap_or_else(|e| panic!("{e}"));
    assert_eq!((u.len(), q.len()), (40, 40));
    for i in 0..40 {
        let a = 2.0 * std::f64::consts::PI * i as f64 / 40.0;
        assert!((u[i] - a.cos()).abs() < 1e-12);
        // The flux is defined on the element between node i and i + 1: its midpoint angle is a + pi/40.
        let mid = a + std::f64::consts::PI / 40.0;
        assert!((q[i] - mid.cos()).abs() < 0.1, "q[{i}] = {} vs {}", q[i], mid.cos());
    }
    // The flux integrates to zero over the closed boundary.
    assert!(q.iter().sum::<f64>().abs() < 1e-6);
}

#[test]
fn mismatched_inputs_and_singular_systems_are_errors() {
    assert!(solve_laplace_bem_2d(&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], &[BoundaryCondition::Potential(1.0)]).is_err());
    // Pure Neumann data leaves the potential undetermined up to a constant.
    let square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let neumann = [BoundaryCondition::Flux(0.0); 4];
    assert!(solve_laplace_bem_2d(&square, &neumann).is_err());
}

#[test]
fn three_dimensional_bem_is_a_placeholder() {
    assert_eq!(solve_laplace_bem_3d(), Ok(()));
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_constant_potential_has_zero_flux(u_val in 0.0..1000.0f64) {
        let points = vec![(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
        let bcs = vec![BoundaryCondition::Potential(u_val); 4];
        let (u, q) = solve_laplace_bem_2d(&points, &bcs).map_err(TestCaseError::fail)?;
        for i in 0..4 {
            prop_assert!((u[i] - u_val).abs() < 1e-10);
            prop_assert!(q[i].abs() < 1e-5);
        }
    }

    #[test]
    fn prop_linear_boundary_data_is_reproduced_inside(a in -5.0..5.0f64, b in -5.0..5.0f64) {
        // u = a x + b y is harmonic; prescribe it on a circle and test the interior value at the centre offset.
        let n = 48;
        let (mut points, mut bcs) = (Vec::new(), Vec::new());
        for i in 0..n {
            let t = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            let (x, y) = (t.cos(), t.sin());
            points.push((x, y));
            bcs.push(BoundaryCondition::Potential(a * x + b * y));
        }
        let (u, q) = solve_laplace_bem_2d(&points, &bcs).map_err(TestCaseError::fail)?;
        let elements = elements_of(&points);
        let pot = evaluate_potential_2d((0.2, -0.3), &elements, &u, &q);
        let want = a * 0.2 - b * 0.3;
        prop_assert!((pot - want).abs() < 0.08 * (a.abs() + b.abs() + 1.0), "{pot} vs {want}");
    }
}
