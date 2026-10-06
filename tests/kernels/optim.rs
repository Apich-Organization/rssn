//! Tests for the optimisation kernels.

use rssn::kernels::dense::Mat;
use rssn::kernels::optim::{
    LpStatus, Pricing, Relation, augmented_lagrangian, bfgs, differential_evolution, lbfgs, nelder_mead,
    simplex, simplex_with, simulated_annealing, trust_region,
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

fn lp_value(r: &rssn::kernels::optim::LpResult) -> f64 {
    match &r.status {
        LpStatus::Optimal { value, .. } => *value,
        other => panic!("{other:?}"),
    }
}

#[test]
fn simplex_beale_cycling_example() {
    // Beale's classic example on which Dantzig's rule with naive tie-breaking cycles
    let c = [-0.75, 20.0, -0.5, 6.0];
    let a = [
        vec![0.25, -8.0, -1.0, 9.0],
        vec![0.5, -12.0, -0.5, 3.0],
        vec![0.0, 0.0, 1.0, 0.0],
    ];
    let rel = [Relation::Le; 3];
    let b = [0.0, 0.0, 1.0];
    for pr in [Pricing::Dantzig, Pricing::SteepestEdge, Pricing::Bland] {
        let r = simplex_with(&c, &a, &rel, &b, pr);
        assert!((lp_value(&r) + 1.25).abs() < 1e-9, "{pr:?}");
        assert!(r.iterations < 200, "{pr:?}: {}", r.iterations);
    }
    assert!(matches!(simplex(&c, &a, &rel, &b), LpStatus::Optimal { .. }));
}

#[test]
fn simplex_klee_minty_and_pricing_rules_agree() {
    // Klee-Minty cube in dimension 6: max sum 2^(n-j) x_j, optimum 5^n at x_n = 5^n
    let n = 6;
    let c: Vec<f64> = (0..n).map(|j| -(2.0_f64.powi((n - 1 - j) as i32))).collect();
    let mut a = Vec::new();
    let mut b = Vec::new();
    for i in 0..n {
        let mut row = vec![0.0; n];
        for j in 0..i {
            row[j] = 2.0_f64.powi((i - j + 1) as i32);
        }
        row[i] = 1.0;
        a.push(row);
        b.push(5.0_f64.powi(i as i32 + 1));
    }
    let rel = vec![Relation::Le; n];
    for pr in [Pricing::Dantzig, Pricing::SteepestEdge, Pricing::Bland] {
        let r = simplex_with(&c, &a, &rel, &b, pr);
        assert!((lp_value(&r) + 5.0_f64.powi(n as i32)).abs() < 1e-6, "{pr:?}");
    }
}

#[test]
fn simplex_pricing_reduces_pivots_on_random_lps() {
    // random dense feasible LPs: max c.x, A x <= b with positive data
    let mut s = 4242u64;
    let mut rnd = || {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((s >> 33) as f64) / f64::from(1u32 << 31)
    };
    let (mut it_bland, mut it_dantzig, mut it_steep) = (0, 0, 0);
    for _ in 0..8 {
        let (m, n) = (30, 40);
        let c: Vec<f64> = (0..n).map(|_| -(0.5 + rnd())).collect();
        let a: Vec<Vec<f64>> = (0..m).map(|_| (0..n).map(|_| 0.1 + rnd()).collect()).collect();
        let b: Vec<f64> = (0..m).map(|_| 5.0 + 10.0 * rnd()).collect();
        let rel = vec![Relation::Le; m];
        let rb = simplex_with(&c, &a, &rel, &b, Pricing::Bland);
        let rd = simplex_with(&c, &a, &rel, &b, Pricing::Dantzig);
        let rs = simplex_with(&c, &a, &rel, &b, Pricing::SteepestEdge);
        let v = lp_value(&rb);
        assert!((lp_value(&rd) - v).abs() < 1e-7 * v.abs().max(1.0));
        assert!((lp_value(&rs) - v).abs() < 1e-7 * v.abs().max(1.0));
        assert_eq!(rb.bland_pivots, rb.iterations);
        it_bland += rb.iterations;
        it_dantzig += rd.iterations;
        it_steep += rs.iterations;
    }
    assert!(it_dantzig < it_bland, "dantzig {it_dantzig} bland {it_bland}");
    assert!(it_steep < it_bland, "steepest {it_steep} bland {it_bland}");
}

#[test]
fn simplex_two_phase_with_pricing() {
    // equality-constrained transportation-like problem solved by every rule
    let c = [4.0, 6.0, 5.0, 3.0, 2.0, 7.0];
    let a = [
        vec![1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        vec![0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
    ];
    let rel = [Relation::Eq; 5];
    let b = [30.0, 20.0, 15.0, 20.0, 15.0];
    let mut vals = Vec::new();
    for pr in [Pricing::Dantzig, Pricing::SteepestEdge, Pricing::Bland] {
        let r = simplex_with(&c, &a, &rel, &b, pr);
        vals.push(lp_value(&r));
    }
    assert!((vals[0] - vals[1]).abs() < 1e-9 && (vals[1] - vals[2]).abs() < 1e-9);
    // the optimum is no worse than this feasible plan (x = 10, 5, 15, 5, 15, 0)
    assert!(vals[0] <= 4.0 * 10.0 + 6.0 * 5.0 + 5.0 * 15.0 + 3.0 * 5.0 + 2.0 * 15.0 + 7.0 * 0.0 + 1e-9);
}
