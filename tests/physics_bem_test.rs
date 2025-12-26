//! Unit and property-based tests for the physics BEM module.

use assert_approx_eq::assert_approx_eq;
use proptest::prelude::*;
use rssn::physics::physics_bem::*;

#[test]

fn test_solve_laplace_bem_2d_rectangle()
{

    // Steady state heat conduction in a rectangle
    // u=100 at x=0, u=0 at x=1, q=0 at y=0 and y=1
    let mut points = Vec::new();

    let mut bcs = Vec::new();

    let n_per_side = 10;

    // Bottom: y=0, x from 0 to 1
    for i in 0 .. n_per_side {

        points.push((
            (i as f64)
                / n_per_side as f64,
            0.0,
        ));

        bcs.push(
            BoundaryCondition::Flux(
                0.0,
            ),
        );
    }

    // Right: x=1, y from 0 to 1
    for i in 0 .. n_per_side {

        points.push((
            1.0,
            (i as f64)
                / n_per_side as f64,
        ));

        bcs.push(BoundaryCondition::Potential(0.0));
    }

    // Top: y=1, x from 1 to 0
    for i in 0 .. n_per_side {

        points.push((
            1.0 - (i as f64)
                / n_per_side as f64,
            1.0,
        ));

        bcs.push(
            BoundaryCondition::Flux(
                0.0,
            ),
        );
    }

    // Left: x=0, y from 1 to 0
    for i in 0 .. n_per_side {

        points.push((
            0.0,
            1.0 - (i as f64)
                / n_per_side as f64,
        ));

        bcs.push(BoundaryCondition::Potential(100.0));
    }

    let (u, _q) = solve_laplace_bem_2d(
        &points,
        &bcs,
    )
    .unwrap();

    // Check boundary values
    assert_approx_eq!(
        u[n_per_side + n_per_side / 2],
        0.0,
        1e-1
    );

    assert_approx_eq!(
        u[3 * n_per_side
            + n_per_side / 2],
        100.0,
        1e-1
    );

    // Midpoints on insulated sides should be around 50
    assert!(
        u[n_per_side / 2] > 40.0
            && u[n_per_side / 2] < 60.0
    );

    assert!(
        u[2 * n_per_side
            + n_per_side / 2]
            > 40.0
            && u[2 * n_per_side
                + n_per_side / 2]
                < 60.0
    );
}

#[test]

fn test_evaluate_potential_2d() {

    let mut points = Vec::new();

    let mut bcs = Vec::new();

    let n_per_side = 10;

    for i in 0 .. n_per_side {

        points.push((
            (i as f64)
                / n_per_side as f64,
            0.0,
        ));

        bcs.push(
            BoundaryCondition::Flux(
                0.0,
            ),
        );
    }

    for i in 0 .. n_per_side {

        points.push((
            1.0,
            (i as f64)
                / n_per_side as f64,
        ));

        bcs.push(BoundaryCondition::Potential(0.0));
    }

    for i in 0 .. n_per_side {

        points.push((
            1.0 - (i as f64)
                / n_per_side as f64,
            1.0,
        ));

        bcs.push(
            BoundaryCondition::Flux(
                0.0,
            ),
        );
    }

    for i in 0 .. n_per_side {

        points.push((
            0.0,
            1.0 - (i as f64)
                / n_per_side as f64,
        ));

        bcs.push(BoundaryCondition::Potential(100.0));
    }

    let (u, q) = solve_laplace_bem_2d(
        &points,
        &bcs,
    )
    .unwrap();

    let n = points.len();

    let elements : Vec<_> = (0 .. n)
        .map(|i| {

            Element2D::new(
                Vector2D::new(
                    points[i].0,
                    points[i].1,
                ),
                Vector2D::new(
                    points[(i + 1) % n]
                        .0,
                    points[(i + 1) % n]
                        .1,
                ),
            )
        })
        .collect();

    // Internal point at (0.5, 0.5) should be around 50
    let pot = evaluate_potential_2d(
        (0.5, 0.5),
        &elements,
        &u,
        &q,
    );

    assert!(pot > 45.0 && pot < 55.0);
}

proptest! {
    #[test]
    fn prop_solve_laplace_bem_2d_symmetry(u_val in 0.0..1000.0f64) {
        let points = vec![
            (0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)
        ];
        let bcs = vec![
            BoundaryCondition::Potential(u_val),
            BoundaryCondition::Potential(u_val),
            BoundaryCondition::Potential(u_val),
            BoundaryCondition::Potential(u_val),
        ];

        let (u, q) = solve_laplace_bem_2d(&points, &bcs).unwrap();

        for i in 0..4 {
            prop_assert!((u[i] - u_val).abs() < 1e-10);
            prop_assert!(q[i].abs() < 1e-5); // Flux should be zero for constant potential
        }
    }
}
