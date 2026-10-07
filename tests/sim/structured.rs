//! Tests of the dimension-generic structured-grid layer.

use rssn::sim::boundary::{BcValue, BoundaryCondition as Bc, Face, FaceBcs};
use rssn::sim::structured::{
    AdaptiveOptions, Grid, HeatProblem, LinearSolver, NonlinearDiffusion, WaveProblem, heat_adaptive, heat_fixed,
    nonlinear_diffusion, nonlinear_diffusion_adaptive, solve_helmholtz, wave_leapfrog,
};
use std::f64::consts::PI;

fn exact(dim: usize, x: [f64; 3]) -> f64 {
    let e = (x[0] + if dim > 1 { x[1] } else { 0.0 }).exp();
    if dim == 3 { e * (2.0 * PI * x[2]).cos() } else { e }
}

/// `-Δ u` of [`exact`].
fn minus_lap(dim: usize, x: [f64; 3]) -> f64 {
    let c = match dim {
        1 => 1.0,
        2 => 2.0,
        _ => 2.0 - 4.0 * PI * PI,
    };
    -c * exact(dim, x)
}

fn mixed_bcs(dim: usize) -> FaceBcs {
    let u = move |x: [f64; 3]| exact(dim, x);
    if dim == 1 {
        return FaceBcs::dirichlet_zero()
            .with(Face::XMin, Bc::Neumann(BcValue::space(move |x| -u(x))))
            .with(Face::XMax, Bc::Robin { alpha: 2.0, beta: 1.0, g: BcValue::space(move |x| 3.0 * u(x)) });
    }
    let mut b = FaceBcs::dirichlet_zero()
        .with(Face::XMin, Bc::Dirichlet(BcValue::space(u)))
        .with(Face::XMax, Bc::Dirichlet(BcValue::space(u)))
        .with(Face::YMin, Bc::Neumann(BcValue::space(move |x| -u(x))))
        .with(Face::YMax, Bc::Robin { alpha: 2.0, beta: 1.0, g: BcValue::space(move |x| 3.0 * u(x)) });
    if dim == 3 {
        b = b.with_axis(2, &Bc::Periodic);
    }
    b
}

fn poisson_error(dim: usize, n: usize) -> f64 {
    let bcs = mixed_bcs(dim);
    let mut ns = vec![n + 1; dim];
    if dim == 3 {
        ns[2] = n;
    }
    let g = Grid::new(&ns, &vec![0.0; dim], &vec![1.0; dim], &bcs).unwrap();
    let u = solve_helmholtz(&g, &bcs, 0.0, |x| minus_lap(dim, x), LinearSolver::Auto).unwrap();
    g.max_error(&u, |x| exact(dim, x))
}

#[test]
fn poisson_mixed_bcs_second_order_all_dimensions() {
    for (dim, ns) in [(1, [16, 32, 64]), (2, [8, 16, 32]), (3, [4, 8, 16])] {
        let e: Vec<f64> = ns.iter().map(|&n| poisson_error(dim, n)).collect();
        for w in e.windows(2) {
            let order = (w[0] / w[1]).log2();
            assert!(order > 1.8, "dim {dim}: order {order}, errors {e:?}");
        }
    }
}

#[test]
fn helmholtz_indefinite_and_iterative_solvers_agree() {
    let bcs = FaceBcs::dirichlet_zero();
    let g = Grid::unit(2, 25, &bcs).unwrap();
    // -Δu - 4 u = f with u = sin(pi x) sin(pi y)
    let sig = -4.0;
    let u_ex = |x: [f64; 3]| (PI * x[0]).sin() * (PI * x[1]).sin();
    let f = |x: [f64; 3]| (2.0 * PI * PI + sig) * u_ex(x);
    let a = solve_helmholtz(&g, &bcs, sig, f, LinearSolver::Direct).unwrap();
    let b = solve_helmholtz(&g, &bcs, sig, f, LinearSolver::BiCgStab(1e-12)).unwrap();
    let c = solve_helmholtz(&g, &bcs, sig, f, LinearSolver::Gmres(1e-12)).unwrap();
    for p in 0..a.len() {
        assert!((a[p] - b[p]).abs() < 1e-8 && (a[p] - c[p]).abs() < 1e-8);
    }
    assert!(g.max_error(&a, u_ex) < 5e-3);
}

#[test]
fn pure_neumann_is_singular() {
    let bcs = FaceBcs::uniform(&Bc::neumann(0.0));
    let g = Grid::unit(1, 9, &bcs).unwrap();
    assert!(solve_helmholtz(&g, &bcs, 0.0, |_| 1.0, LinearSolver::Direct).is_err());
}

/// 1D heat with Robin data on the right and time-dependent Dirichlet on the left.
#[test]
fn heat_robin_time_dependent_second_order_in_time_and_space() {
    let u = |t: f64, x: f64| (-t).exp() * x.cos() + t * x;
    // u_t - u_xx = x
    let src = |_t: f64, x: [f64; 3]| x[0];
    let l = 1.5;
    let robin_g = move |t: f64| u(t, l) + 2.0 * ((-t).exp() * -(l.sin()) + t);
    let bcs = FaceBcs::dirichlet_zero()
        .with(Face::XMin, Bc::Dirichlet(BcValue::time(move |t| u(t, 0.0))))
        .with(Face::XMax, Bc::Robin { alpha: 1.0, beta: 2.0, g: BcValue::time(robin_g) });
    let run = |n: usize, steps: usize| {
        let g = Grid::new(&[n], &[0.0], &[l], &bcs).unwrap();
        let init = |x: [f64; 3]| u(0.0, x[0]);
        let prob = HeatProblem {
            grid: &g,
            bcs: &bcs,
            alpha: 1.0,
            init: &init,
            source: Some(&src),
            solver: LinearSolver::Auto,
        };
        let v = heat_fixed(&prob, 1.0, steps, 0.5).unwrap();
        g.max_error(&v, |x| u(1.0, x[0]))
    };
    let e1 = run(21, 20);
    let e2 = run(41, 40);
    let e3 = run(81, 80);
    assert!((e1 / e2).log2() > 1.8 && (e2 / e3).log2() > 1.8, "{e1} {e2} {e3}");
}

#[test]
fn heat_3d_periodic_decay_and_adaptive_error_control() {
    let bcs = FaceBcs::uniform(&Bc::Periodic);
    let g = Grid::new(&[8, 8, 8], &[0.0; 3], &[1.0; 3], &bcs).unwrap();
    let init = |x: [f64; 3]| (2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).cos() * (2.0 * PI * x[2]).sin();
    let prob = HeatProblem { grid: &g, bcs: &bcs, alpha: 0.05, init: &init, source: None, solver: LinearSolver::Auto };
    let h = g.h(0);
    let mu = 3.0 * (2.0 - 2.0 * (2.0 * PI * h).cos()) / (h * h) * 0.05;
    let (u, st) = heat_adaptive(&prob, 0.2, &AdaptiveOptions { tol: 1e-7, dt0: 1e-3, ..Default::default() }).unwrap();
    let ex = g.sample(|x| (-mu * 0.2).exp() * init(x));
    let err = u.iter().zip(&ex).fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    assert!(err < 2e-5, "err {err}");
    assert!(st.accepted > 2 && st.dt_max_used > st.dt_min_used);
}

#[test]
fn wave_1d_dirichlet_standing_wave_and_mur_absorbs() {
    let bcs = FaceBcs::dirichlet_zero();
    let g = Grid::new(&[201], &[0.0], &[1.0], &bcs).unwrap();
    let u0 = |x: [f64; 3]| (PI * x[0]).sin();
    let v0 = |_: [f64; 3]| 0.0;
    let prob = WaveProblem { grid: &g, bcs: &bcs, c: 1.0, u0: &u0, v0: &v0, source: None };
    let u = wave_leapfrog(&prob, 0.5, 100).unwrap();
    assert!(g.max_error(&u, |x| (PI * x[0]).sin() * (PI * 0.5).cos()) < 1e-4);
    assert!(wave_leapfrog(&prob, 0.5, 10).is_err());

    // right-going pulse leaves through absorbing ends
    let abs = FaceBcs::uniform(&Bc::Absorbing);
    let g = Grid::new(&[401], &[0.0], &[2.0], &abs).unwrap();
    let pulse = |x: f64| (-((x - 1.0) / 0.1).powi(2)).exp();
    let u0 = move |x: [f64; 3]| pulse(x[0]);
    let v0 = move |x: [f64; 3]| {
        let h = 1e-6;
        -(pulse(x[0] + h) - pulse(x[0] - h)) / (2.0 * h)
    };
    let prob = WaveProblem { grid: &g, bcs: &abs, c: 1.0, u0: &u0, v0: &v0, source: None };
    let u = wave_leapfrog(&prob, 2.5, 500).unwrap();
    let peak = u.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    assert!(peak < 1e-2, "residual {peak}");
}

#[test]
fn wave_2d_neumann_cosine_mode() {
    let bcs = FaceBcs::uniform(&Bc::neumann(0.0));
    let run = |n: usize| {
        let g = Grid::new(&[n, n], &[0.0; 2], &[1.0; 2], &bcs).unwrap();
        let u0 = |x: [f64; 3]| (PI * x[0]).cos() * (PI * x[1]).cos();
        let v0 = |_: [f64; 3]| 0.0;
        let prob = WaveProblem { grid: &g, bcs: &bcs, c: 1.0, u0: &u0, v0: &v0, source: None };
        let steps = (n - 1) * 2;
        let u = wave_leapfrog(&prob, 0.4, steps).unwrap();
        let w = PI * 2.0_f64.sqrt();
        g.max_error(&u, |x| (PI * x[0]).cos() * (PI * x[1]).cos() * (w * 0.4).cos())
    };
    let (e1, e2) = (run(21), run(41));
    assert!((e1 / e2).log2() > 1.8, "{e1} {e2}");
}

fn nl_exact(t: f64, x: [f64; 3], dim: usize) -> f64 {
    let mut s = (-t).exp();
    for a in 0..dim {
        s *= (PI * x[a]).sin();
    }
    1.0 + 0.5 * s
}

#[test]
fn nonlinear_diffusion_newton_krylov_manufactured_2d() {
    let d = |u: f64| 1.0 + u * u;
    let src = move |t: f64, x: [f64; 3], _u: f64| {
        let e = (-t).exp();
        let (sx, sy) = ((PI * x[0]).sin(), (PI * x[1]).sin());
        let (cx, cy) = ((PI * x[0]).cos(), (PI * x[1]).cos());
        let u = 1.0 + 0.5 * e * sx * sy;
        let ut = -0.5 * e * sx * sy;
        let ux = 0.5 * e * PI * cx * sy;
        let uy = 0.5 * e * PI * sx * cy;
        let lap = -0.5 * e * 2.0 * PI * PI * sx * sy;
        ut - ((1.0 + u * u) * lap + 2.0 * u * (ux * ux + uy * uy))
    };
    let bcs = FaceBcs::uniform(&Bc::dirichlet(1.0));
    let run = |n: usize, steps: usize| {
        let g = Grid::unit(2, n, &bcs).unwrap();
        let init = move |x: [f64; 3]| nl_exact(0.0, x, 2);
        let prob = NonlinearDiffusion { grid: &g, bcs: &bcs, d: &d, source: Some(&src), init: &init };
        let (u, st) = nonlinear_diffusion(&prob, 0.3, steps, 0.5).unwrap();
        assert!(st.newton_iters > 0);
        g.max_error(&u, |x| nl_exact(0.3, x, 2))
    };
    let (e1, e2) = (run(11, 10), run(21, 20));
    assert!((e1 / e2).log2() > 1.7, "{e1} {e2}");

    let g = Grid::unit(2, 15, &bcs).unwrap();
    let init = move |x: [f64; 3]| nl_exact(0.0, x, 2);
    let prob = NonlinearDiffusion { grid: &g, bcs: &bcs, d: &d, source: Some(&src), init: &init };
    let (u, stats, ns) =
        nonlinear_diffusion_adaptive(&prob, 0.3, &AdaptiveOptions { tol: 1e-6, dt0: 0.01, ..Default::default() })
            .unwrap();
    assert!(stats.accepted >= 1 && ns.krylov_iters > 0);
    assert!(g.max_error(&u, |x| nl_exact(0.3, x, 2)) < 2e-3);
}

#[test]
fn nonlinear_diffusion_neumann_conserves_mass_1d() {
    let bcs = FaceBcs::uniform(&Bc::neumann(0.0));
    let g = Grid::new(&[41], &[0.0], &[1.0], &bcs).unwrap();
    let d = |u: f64| u * u + 0.1;
    let init = |x: [f64; 3]| 1.0 + (-(x[0] - 0.5).powi(2) * 50.0).exp();
    let prob = NonlinearDiffusion { grid: &g, bcs: &bcs, d: &d, source: None, init: &init };
    let mass = |u: &[f64]| {
        let h = g.h(0);
        let n = u.len();
        h * (u.iter().sum::<f64>() - 0.5 * (u[0] + u[n - 1]))
    };
    let u0 = g.sample(init);
    let (u, _) = nonlinear_diffusion(&prob, 0.2, 20, 1.0).unwrap();
    assert!((mass(&u) - mass(&u0)).abs() < 5e-3 * mass(&u0));
    let spread = u.iter().fold(0.0_f64, |m, v| m.max(*v)) - u.iter().fold(f64::INFINITY, |m, v| m.min(*v));
    assert!(spread < 1.0);
}
