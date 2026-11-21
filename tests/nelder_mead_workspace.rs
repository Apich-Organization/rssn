//! Candidate ordering and accuracy of the Nelder–Mead workspace.

use rssn::kernels::optim::nelder_mead;
use std::cell::RefCell;

#[test]
fn expansion_evaluates_reflection_before_expansion() {
    let visits = RefCell::new(Vec::new());
    let result = nelder_mead(
        |point| {
            let x = point[0];
            visits.borrow_mut().push(x);
            (x - 4.0).powi(2)
        },
        &[1.0],
        1.0,
        0.0,
        1,
    );
    assert_eq!(visits.into_inner(), [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(result.x, [4.0]);
    assert_eq!(result.evaluations, 4);
    assert_eq!(result.iterations, 1);
    assert!(!result.converged);
}

#[test]
fn an_accepted_contraction_keeps_the_previous_best_on_a_tie() {
    let visits = RefCell::new(Vec::new());
    let result = nelder_mead(
        |point| {
            let x = point[0];
            visits.borrow_mut().push(x);
            if x.to_bits() == 1.0_f64.to_bits() {
                2.0
            } else if x.to_bits() == 3.0_f64.to_bits() {
                1.5
            } else {
                1.0
            }
        },
        &[1.0],
        1.0,
        0.0,
        1,
    );
    assert_eq!(visits.into_inner(), [1.0, 2.0, 3.0, 2.5]);
    assert_eq!(result.x, [2.0]);
    assert_eq!(result.evaluations, 4);
    assert!(!result.converged);
}

#[test]
fn tied_vertices_keep_the_same_shrink_anchor_on_every_iteration() {
    let visits = RefCell::new(Vec::new());
    let result = nelder_mead(
        |point| {
            visits.borrow_mut().push(point[0]);
            7.0
        },
        &[1.0],
        1.0,
        0.0,
        2,
    );
    assert_eq!(
        visits.into_inner(),
        [1.0, 2.0, 0.0, 1.5, 1.5, 0.5, 1.25, 1.25],
    );
    assert_eq!(result.x, [1.0]);
    assert_eq!(result.evaluations, 8);
    assert_eq!(result.iterations, 2);
    assert!(!result.converged);
}

#[test]
fn empty_and_zero_budget_starts_keep_the_initial_evaluation_contract() {
    let visits = RefCell::new(Vec::new());
    let objective = |point: &[f64]| {
        visits.borrow_mut().push(point.to_vec());
        point.iter().map(|x| x * x).sum()
    };
    let empty = nelder_mead(objective, &[], 0.5, 0.0, 20);
    assert_eq!(visits.take(), [Vec::<f64>::new()]);
    assert!(empty.x.is_empty());
    assert_eq!(empty.evaluations, 1);
    assert_eq!(empty.iterations, 0);
    assert!(empty.converged);

    let initial = [-2.0, 0.0, 4.0];
    let result = nelder_mead(objective, &initial, 0.5, 0.0, 0);
    assert_eq!(
        visits.into_inner(),
        [
            initial,
            [-3.0, 0.0, 4.0],
            [-2.0, 0.25, 4.0],
            [-2.0, 0.0, 6.0],
        ],
    );
    assert_eq!(result.x, initial);
    assert_eq!(result.evaluations, 4);
    assert_eq!(result.iterations, 0);
    assert!(!result.converged);
}

#[test]
fn a_weighted_quadratic_recovers_its_planted_minimum() {
    let minimum = [1.75, -3.0, 0.5, 2.25];
    let weights = [0.25, 1.0, 4.0, 16.0];
    let result = nelder_mead(
        |point| {
            point
                .iter()
                .zip(minimum)
                .zip(weights)
                .map(|((x, center), weight)| weight * (x - center).powi(2))
                .sum()
        },
        &[0.0; 4],
        1.0,
        1e-13,
        2_000,
    );
    assert!(result.converged, "{result:?}");
    assert!(result.fx < 1e-12, "{result:?}");
    for (coordinate, target) in result.x.iter().zip(minimum) {
        assert!((coordinate - target).abs() < 1e-6, "{result:?}");
    }
}
