//! Tests for the optimisation kernels.

use rssn::kernels::dense::Mat;
use rssn::kernels::optim::{
    LpStatus, Relation, augmented_lagrangian, bfgs, differential_evolution, lbfgs, nelder_mead,
    simplex, simulated_annealing, trust_region,
};

fn rosen(x: &[f64]) -> f64 {
    100.0 * (x[1] - x[0] * x[0]).powi(2) + (1.0 - x[0]).powi(2)
}

fn rosen_grad(x: &[f64]) -> Vec<f64> {
    vec![
        -400.0 * x[0] * (x[1] - x[0] * x[0]) - 2.0 * (1.0 - x[0]),
        200.0 * (x[1] - x[0] * x[0]),
    ]
}

fn rosen_hess(x: &[f64]) -> Mat {
    Mat {
        rows: 2,
        cols: 2,
        data: vec![
            1200.0 * x[0] * x[0] - 400.0 * x[1] + 2.0,
            -400.0 * x[0],
            -400.0 * x[0],
            200.0,
        ],
    }
}

#[test]
fn nelder_mead_rosenbrock() {
    let r = nelder_mead(rosen, &[-1.2, 1.0], 0.1, 1e-14, 5000);
    assert!(r.converged);
    assert!((r.x[0] - 1.0).abs() < 1e-4 && (r.x[1] - 1.0).abs() < 1e-4, "{r:?}");
}

#[test]
fn bfgs_with_and_without_gradient() {
    let r = bfgs(rosen, Some(rosen_grad), &[-1.2, 1.0], 1e-8, 200);
    assert!(r.converged, "{r:?}");
    assert!((r.x[0] - 1.0).abs() < 1e-6 && (r.x[1] - 1.0).abs() < 1e-6);
    let r = bfgs(rosen, None::<fn(&[f64]) -> Vec<f64>>, &[-1.2, 1.0], 1e-5, 300);
    assert!((r.x[0] - 1.0).abs() < 1e-3, "{r:?}");
}

#[test]
fn lbfgs_high_dimension() {
    let n = 50;
    let f = |x: &[f64]| x.iter().enumerate().map(|(i, v)| (i as f64 + 1.0) * (v - 1.0).powi(2)).sum::<f64>();
    let g = |x: &[f64]| x.iter().enumerate().map(|(i, v)| 2.0 * (i as f64 + 1.0) * (v - 1.0)).collect();
    let r = lbfgs(f, Some(g), &vec![0.0; n], 8, 1e-8, 500);
    assert!(r.converged, "{r:?}");
    assert!(r.x.iter().all(|v| (v - 1.0).abs() < 1e-6));
    let r2 = lbfgs(rosen, Some(rosen_grad), &[-1.2, 1.0], 10, 1e-8, 500);
    assert!((r2.x[0] - 1.0).abs() < 1e-5);
}

#[test]
fn trust_region_dogleg() {
    let r = trust_region(rosen, rosen_grad, Some(rosen_hess), &[-1.2, 1.0], 1e-8, 200);
    assert!(r.converged, "{r:?}");
    assert!((r.x[0] - 1.0).abs() < 1e-6);
    let r = trust_region(rosen, rosen_grad, None::<fn(&[f64]) -> Mat>, &[-1.2, 1.0], 1e-6, 500);
    assert!((r.x[0] - 1.0).abs() < 1e-4, "{r:?}");
}

#[test]
fn global_methods_are_deterministic() {
    let rastrigin = |x: &[f64]| {
        10.0 * x.len() as f64
            + x.iter().map(|v| v * v - 10.0 * (2.0 * std::f64::consts::PI * v).cos()).sum::<f64>()
    };
    let bounds = [(-5.12, 5.12); 3];
    let a = differential_evolution(rastrigin, &bounds, 40, 400, 7);
    let b = differential_evolution(rastrigin, &bounds, 40, 400, 7);
    assert_eq!(a, b);
    assert!(a.fx < 1e-6, "{a:?}");
    let s = simulated_annealing(|x| (x[0] - 3.0).powi(2) + (x[1] + 1.0).powi(2), &[0.0, 0.0], 1.0, 5.0, 0.999, 20000, 1);
    assert!((s.x[0] - 3.0).abs() < 0.05 && (s.x[1] + 1.0).abs() < 0.05, "{s:?}");
    let s2 = simulated_annealing(|x| (x[0] - 3.0).powi(2) + (x[1] + 1.0).powi(2), &[0.0, 0.0], 1.0, 5.0, 0.999, 20000, 1);
    assert_eq!(s, s2);
}

#[test]
fn augmented_lagrangian_constraints() {
    // min x^2 + y^2 s.t. x + y = 1 -> (0.5, 0.5)
    let r = augmented_lagrangian(
        |x| x[0] * x[0] + x[1] * x[1],
        |x| vec![x[0] + x[1] - 1.0],
        |_| vec![],
        &[0.0, 0.0],
        1e-8,
        30,
    );
    assert!(r.converged, "{r:?}");
    assert!((r.x[0] - 0.5).abs() < 1e-5 && (r.x[1] - 0.5).abs() < 1e-5);
    // min (x-2)^2 s.t. x <= 1 -> x = 1
    let r = augmented_lagrangian(|x| (x[0] - 2.0).powi(2), |_| vec![], |x| vec![x[0] - 1.0], &[0.0], 1e-8, 30);
    assert!((r.x[0] - 1.0).abs() < 1e-5, "{r:?}");
}

#[test]
fn simplex_linear_programs() {
    // max 3x + 2y s.t. x + y <= 4, x + 3y <= 6, x <= 3  -> (3, 1), value 11
    let lp = simplex(
        &[-3.0, -2.0],
        &[vec![1.0, 1.0], vec![1.0, 3.0], vec![1.0, 0.0]],
        &[Relation::Le, Relation::Le, Relation::Le],
        &[4.0, 6.0, 3.0],
    );
    match lp {
        LpStatus::Optimal { x, value } => {
            assert!((x[0] - 3.0).abs() < 1e-9 && (x[1] - 1.0).abs() < 1e-9);
            assert!((value + 11.0).abs() < 1e-9);
        }
        other => panic!("{other:?}"),
    }
    // min x + y s.t. x + 2y >= 4, 3x + y >= 6 -> (1.6, 1.2) value 2.8
    let lp = simplex(
        &[1.0, 1.0],
        &[vec![1.0, 2.0], vec![3.0, 1.0]],
        &[Relation::Ge, Relation::Ge],
        &[4.0, 6.0],
    );
    match lp {
        LpStatus::Optimal { x, value } => {
            assert!((x[0] - 1.6).abs() < 1e-9 && (x[1] - 1.2).abs() < 1e-9);
            assert!((value - 2.8).abs() < 1e-9);
        }
        other => panic!("{other:?}"),
    }
    // equality + infeasible + unbounded
    let lp = simplex(&[1.0, 2.0], &[vec![1.0, 1.0]], &[Relation::Eq], &[5.0]);
    assert!(matches!(lp, LpStatus::Optimal { value, .. } if (value - 5.0).abs() < 1e-9));
    let lp = simplex(&[1.0], &[vec![1.0], vec![1.0]], &[Relation::Le, Relation::Ge], &[1.0, 2.0]);
    assert_eq!(lp, LpStatus::Infeasible);
    let lp = simplex(&[-1.0], &[vec![1.0]], &[Relation::Ge], &[1.0]);
    assert_eq!(lp, LpStatus::Unbounded);
}
