//! Mixed boundary behavior on the minimum and ordinary 1D method-of-lines grids.

use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::pde::Boundary;
use rssn::kernels::pde::Grid1D;
use rssn::kernels::pde::PdeError;
use rssn::kernels::pde::StiffMethod;
use rssn::kernels::pde::TimedBoundary1D;
use rssn::kernels::pde::method_of_lines_1d;
use rssn::kernels::pde::method_of_lines_1d_timed;
use std::cell::RefCell;
use std::sync::Arc;
use std::sync::Mutex;

#[test]
fn every_boundary_pair_preserves_an_affine_steady_state() {
    let options = OdeOptions {
        rtol: 1e-9,
        atol: 1e-11,
        ..OdeOptions::default()
    };
    for n in [3, 9, 17] {
        let grid = Grid1D { x0: -0.5, x1: 0.5, n };
        let affine = |x| 1.5 + 2.0 * x;
        for left_fixed in [false, true] {
            for right_fixed in [false, true] {
                for method in [StiffMethod::Radau5, StiffMethod::Bdf(3)] {
                    let calls = Arc::new(Mutex::new(Vec::new()));
                    let boundary = |fixed, end, value| {
                        let calls = Arc::clone(&calls);
                        let sample = move |x, y| {
                            calls.lock().unwrap().push((end, x, y));
                            value
                        };
                        if fixed {
                            Boundary::dirichlet(sample)
                        } else {
                            Boundary::neumann(sample)
                        }
                    };
                    let bc = [
                        boundary(
                            left_fixed,
                            0,
                            if left_fixed {
                                affine(-0.5)
                            } else {
                                -2.0
                            },
                        ),
                        boundary(
                            right_fixed,
                            1,
                            if right_fixed {
                                affine(0.5)
                            } else {
                                2.0
                            },
                        ),
                    ];
                    let sampled = RefCell::new(Vec::new());
                    let initial = |x| {
                        sampled.borrow_mut().push(x);
                        affine(x)
                    };
                    let f =
                        |_, _, _, gradient, laplacian| 0.25 * laplacian - 0.125 * (gradient - 2.0);
                    let out = method_of_lines_1d(&f, &grid, &bc, &initial, 0.1, method, &options)
                        .unwrap();
                    assert_eq!(out.field.x, grid.nodes());
                    assert!(out.field.max_error(affine) <= 1e-11, "{method:?}/{n}");
                    assert!(out.steps > 0 && out.nfev > 0);
                    assert_eq!(
                        *sampled.borrow(),
                        &out.field.x[usize::from(left_fixed)..n - usize::from(right_fixed)]
                    );
                    let both_calls = calls.lock().unwrap();
                    let expected = if !left_fixed && right_fixed {
                        [(1, 0.5, 0.0), (0, -0.5, 0.0)]
                    } else {
                        [(0, -0.5, 0.0), (1, 0.5, 0.0)]
                    };
                    assert_eq!(both_calls.as_slice(), expected);
                }
            }
        }
    }
}

#[test]
fn invalid_time_or_grid_precedes_sampling() {
    let boundary = Boundary::dirichlet(|_, _| panic!("invalid request sampled boundary"));
    let bc = [boundary.clone(), boundary];
    let initial = |_| panic!("invalid request sampled initial data");
    let f = |_, _, _, _, _| panic!("invalid request sampled RHS");
    let valid = Grid1D {
        x0: 0.0,
        x1: 1.0,
        n: 3,
    };
    let invalid = Grid1D { n: 2, ..valid };
    for (grid, end) in [(valid, 0.0), (valid, f64::NAN), (invalid, 0.1)] {
        assert!(matches!(
            method_of_lines_1d(
                &f,
                &grid,
                &bc,
                &initial,
                end,
                StiffMethod::Radau5,
                &OdeOptions::default()
            ),
            Err(PdeError::InvalidInput(_))
        ));
    }
}

#[test]
fn time_dependent_dirichlet_and_neumann_data_reproduce_a_manufactured_field() {
    let grid = Grid1D {
        x0: 0.0,
        x1: 1.0,
        n: 9,
    };
    let options = OdeOptions {
        rtol: 1e-9,
        atol: 1e-11,
        ..OdeOptions::default()
    };
    let initial = |x| x;
    let end = 0.35;
    let oracle = |time, x, _, gradient, laplacian| {
        x + 2.0 + 0.125 * (gradient - (1.0 + time)) + 0.0625 * laplacian
    };
    let left_value = |time| 2.0 * time;
    let left_gradient = |time| -1.0 - time;
    let right_value = |time| 1.0 + 3.0 * time;
    let right_gradient = |time| 1.0 + time;
    let expected = |x| (1.0 + end) * x + 2.0 * end;
    for left in [
        TimedBoundary1D::Dirichlet(&left_value),
        TimedBoundary1D::Neumann(&left_gradient),
    ] {
        for right in [
            TimedBoundary1D::Dirichlet(&right_value),
            TimedBoundary1D::Neumann(&right_gradient),
        ] {
            for method in [StiffMethod::Radau5, StiffMethod::Bdf(3)] {
                let result = method_of_lines_1d_timed(
                    &oracle,
                    &grid,
                    &[left, right],
                    &initial,
                    end,
                    method,
                    &options,
                )
                .unwrap();
                assert_eq!(result.field.t, end);
                assert!(result.field.max_error(expected) < 2e-8, "{method:?}");
            }
        }
    }
}

#[test]
fn time_dependent_boundaries_reject_nonfinite_stages_even_if_the_rhs_ignores_them() {
    let grid = Grid1D {
        x0: 0.0,
        x1: 1.0,
        n: 3,
    };
    let boundaries = [
        TimedBoundary1D::Dirichlet(&|_| 0.0),
        TimedBoundary1D::Neumann(&|_| f64::NAN),
    ];
    let rhs =
        |_, _, _, _, _| panic!("nonfinite endpoints must be rejected before the PDE callback");
    for method in [StiffMethod::Radau5, StiffMethod::Bdf(3)] {
        assert!(matches!(
            method_of_lines_1d_timed(
                &rhs,
                &grid,
                &boundaries,
                &|_| 0.0,
                0.25,
                method,
                &OdeOptions::default()
            ),
            Err(PdeError::NonFinite)
        ));
    }
}
