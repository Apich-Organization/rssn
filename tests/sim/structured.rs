//! Tests of the dimension-generic structured-grid layer.

use rssn::sim::boundary::BcValue;
use rssn::sim::boundary::BoundaryCondition as Bc;
use rssn::sim::boundary::Face;
use rssn::sim::boundary::FaceBcs;
use rssn::sim::structured::AdaptiveOptions;
use rssn::sim::structured::FaceAverage;
use rssn::sim::structured::Grid;
use rssn::sim::structured::HeatProblem;
use rssn::sim::structured::LinearSolver;
use rssn::sim::structured::NonlinearDiffusion;
use rssn::sim::structured::VariableHeatProblem;
use rssn::sim::structured::WaveProblem;
use rssn::sim::structured::heat_adaptive;
use rssn::sim::structured::heat_fixed;
use rssn::sim::structured::heat_variable_adaptive;
use rssn::sim::structured::heat_variable_fixed;
use rssn::sim::structured::nonlinear_diffusion;
use rssn::sim::structured::nonlinear_diffusion_adaptive;
use rssn::sim::structured::solve_helmholtz;
use rssn::sim::structured::solve_variable_helmholtz;
use rssn::sim::structured::wave_leapfrog;
use std::f64::consts::PI;

fn exact(
    dim: usize,
    x: [f64; 3],
) -> f64 {
    let e = (x[0] + if dim > 1 { x[1] } else { 0.0 }).exp();
    if dim == 3 {
        e * (2.0 * PI * x[2]).cos()
    } else {
        e
    }
}

/// `-Δ u` of [`exact`].
fn minus_lap(
    dim: usize,
    x: [f64; 3],
) -> f64 {
    let c = match dim {
        | 1 => 1.0,
        | 2 => 2.0,
        | _ => 2.0 - 4.0 * PI * PI,
    };
    -c * exact(dim, x)
}

fn mixed_bcs(dim: usize) -> FaceBcs {
    let u = move |x: [f64; 3]| exact(dim, x);
    if dim == 1 {
        return FaceBcs::dirichlet_zero()
            .with(Face::XMin, Bc::Neumann(BcValue::space(move |x| -u(x))))
            .with(
                Face::XMax,
                Bc::Robin {
                    alpha: 2.0,
                    beta: 1.0,
                    g: BcValue::space(move |x| 3.0 * u(x)),
                },
            );
    }
    let mut b = FaceBcs::dirichlet_zero()
        .with(Face::XMin, Bc::Dirichlet(BcValue::space(u)))
        .with(Face::XMax, Bc::Dirichlet(BcValue::space(u)))
        .with(Face::YMin, Bc::Neumann(BcValue::space(move |x| -u(x))))
        .with(
            Face::YMax,
            Bc::Robin {
                alpha: 2.0,
                beta: 1.0,
                g: BcValue::space(move |x| 3.0 * u(x)),
            },
        );
    if dim == 3 {
        b = b.with_axis(2, &Bc::Periodic);
    }
    b
}

fn poisson_error(
    dim: usize,
    n: usize,
) -> f64 {
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
        .with(
            Face::XMax,
            Bc::Robin {
                alpha: 1.0,
                beta: 2.0,
                g: BcValue::time(robin_g),
            },
        );
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
    assert!(
        (e1 / e2).log2() > 1.8 && (e2 / e3).log2() > 1.8,
        "{e1} {e2} {e3}"
    );
}

#[test]
fn heat_3d_periodic_decay_and_adaptive_error_control() {
    let bcs = FaceBcs::uniform(&Bc::Periodic);
    let g = Grid::new(&[8, 8, 8], &[0.0; 3], &[1.0; 3], &bcs).unwrap();
    let init =
        |x: [f64; 3]| (2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).cos() * (2.0 * PI * x[2]).sin();
    let prob = HeatProblem {
        grid: &g,
        bcs: &bcs,
        alpha: 0.05,
        init: &init,
        source: None,
        solver: LinearSolver::Auto,
    };
    let h = g.h(0);
    let mu = 3.0 * (2.0 - 2.0 * (2.0 * PI * h).cos()) / (h * h) * 0.05;
    let (u, st) = heat_adaptive(
        &prob,
        0.2,
        &AdaptiveOptions {
            tol: 1e-7,
            dt0: 1e-3,
            ..Default::default()
        },
    )
    .unwrap();
    let ex = g.sample(|x| (-mu * 0.2).exp() * init(x));
    let err = u
        .iter()
        .zip(&ex)
        .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    assert!(err < 2e-5, "err {err}");
    assert!(st.accepted > 2 && st.dt_max_used > st.dt_min_used);
}

#[test]
fn wave_1d_dirichlet_standing_wave_and_mur_absorbs() {
    let bcs = FaceBcs::dirichlet_zero();
    let g = Grid::new(&[201], &[0.0], &[1.0], &bcs).unwrap();
    let u0 = |x: [f64; 3]| (PI * x[0]).sin();
    let v0 = |_: [f64; 3]| 0.0;
    let prob = WaveProblem {
        grid: &g,
        bcs: &bcs,
        c: 1.0,
        u0: &u0,
        v0: &v0,
        source: None,
    };
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
    let prob = WaveProblem {
        grid: &g,
        bcs: &abs,
        c: 1.0,
        u0: &u0,
        v0: &v0,
        source: None,
    };
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
        let prob = WaveProblem {
            grid: &g,
            bcs: &bcs,
            c: 1.0,
            u0: &u0,
            v0: &v0,
            source: None,
        };
        let steps = (n - 1) * 2;
        let u = wave_leapfrog(&prob, 0.4, steps).unwrap();
        let w = PI * 2.0_f64.sqrt();
        g.max_error(&u, |x| {
            (PI * x[0]).cos() * (PI * x[1]).cos() * (w * 0.4).cos()
        })
    };
    let (e1, e2) = (run(21), run(41));
    assert!((e1 / e2).log2() > 1.8, "{e1} {e2}");
}

fn nl_exact(
    t: f64,
    x: [f64; 3],
    dim: usize,
) -> f64 {
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
        let prob = NonlinearDiffusion {
            grid: &g,
            bcs: &bcs,
            d: &d,
            source: Some(&src),
            init: &init,
        };
        let (u, st) = nonlinear_diffusion(&prob, 0.3, steps, 0.5).unwrap();
        assert!(st.newton_iters > 0);
        g.max_error(&u, |x| nl_exact(0.3, x, 2))
    };
    let (e1, e2) = (run(11, 10), run(21, 20));
    assert!((e1 / e2).log2() > 1.7, "{e1} {e2}");

    let g = Grid::unit(2, 15, &bcs).unwrap();
    let init = move |x: [f64; 3]| nl_exact(0.0, x, 2);
    let prob = NonlinearDiffusion {
        grid: &g,
        bcs: &bcs,
        d: &d,
        source: Some(&src),
        init: &init,
    };
    let (u, stats, ns) = nonlinear_diffusion_adaptive(
        &prob,
        0.3,
        &AdaptiveOptions {
            tol: 1e-6,
            dt0: 0.01,
            ..Default::default()
        },
    )
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
    let prob = NonlinearDiffusion {
        grid: &g,
        bcs: &bcs,
        d: &d,
        source: None,
        init: &init,
    };
    let mass = |u: &[f64]| {
        let h = g.h(0);
        let n = u.len();
        h * (u.iter().sum::<f64>() - 0.5 * (u[0] + u[n - 1]))
    };
    let u0 = g.sample(init);
    let (u, _) = nonlinear_diffusion(&prob, 0.2, 20, 1.0).unwrap();
    assert!((mass(&u) - mass(&u0)).abs() < 5e-3 * mass(&u0));
    let spread =
        u.iter().fold(0.0_f64, |m, v| m.max(*v)) - u.iter().fold(f64::INFINITY, |m, v| m.min(*v));
    assert!(spread < 1.0);
}

/// Constant `α` on the variable-coefficient path must reproduce the
/// constant-coefficient solver and operator to round-off.
#[test]
fn variable_diffusion_constant_alpha_matches_heat_fixed() {
    let bcs = mixed_bcs(1);
    let g = Grid::new(&[41], &[0.0], &[1.0], &bcs).unwrap();
    let src = |_t: f64, x: [f64; 3]| x[0] * x[0];
    let init = |x: [f64; 3]| 1.0 - x[0] * x[0];
    let fixed = HeatProblem {
        grid: &g,
        bcs: &bcs,
        alpha: 0.7,
        init: &init,
        source: Some(&src),
        solver: LinearSolver::Auto,
    };
    let alpha = vec![0.7; g.len()];
    let var = VariableHeatProblem {
        grid: &g,
        bcs: &bcs,
        alpha: &alpha,
        averaging: FaceAverage::Harmonic,
        init: &init,
        source: Some(&src),
        solver: LinearSolver::Auto,
    };
    let a = heat_fixed(&fixed, 0.2, 20, 0.5).unwrap();
    let b = heat_variable_fixed(&var, 0.2, 20, 0.5).unwrap();
    let d = a
        .iter()
        .zip(&b)
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
    assert!(d < 1e-11, "constant-alpha mismatch {d}");

    // operator level: -∇·(α∇u) + σu = f is α times the Laplacian solve
    let f = |x: [f64; 3]| x[0].exp();
    let ua = solve_helmholtz(&g, &bcs, 0.3, f, LinearSolver::Direct).unwrap();
    let ub = solve_variable_helmholtz(
        &g,
        &bcs,
        &alpha,
        FaceAverage::Harmonic,
        0.21,
        |x| 0.7 * f(x),
        LinearSolver::Direct,
    )
    .unwrap();
    let d = ua
        .iter()
        .zip(&ub)
        .fold(0.0_f64, |m, (x, y)| m.max((x - y).abs()));
    assert!(d < 1e-11, "constant-alpha operator mismatch {d}");
}

/// 1D layered medium: the harmonic face diffusivity represents a planar
/// layer interface exactly, so the piecewise-linear steady solution is
/// recovered (and reached by the time integration).
#[test]
fn variable_diffusion_layered_piecewise_linear_exact() {
    let (a1, a2) = (1.0, 100.0);
    let bcs = FaceBcs::dirichlet_zero()
        .with(Face::XMin, Bc::dirichlet(0.0))
        .with(Face::XMax, Bc::dirichlet(1.0));
    let flux = 2.0 * a1 * a2 / (a1 + a2);
    let (m1, m2) = (flux / a1, flux / a2);
    let exact = move |x: [f64; 3]| {
        if x[0] <= 0.5 {
            m1 * x[0]
        } else {
            m1 * 0.5 + m2 * (x[0] - 0.5)
        }
    };
    for n in [8, 16, 32] {
        let g = Grid::new(&[n], &[0.0], &[1.0], &bcs).unwrap();
        let alpha: Vec<f64> = (0..g.len())
            .map(|p| {
                if g.coord(p)[0] < 0.5 {
                    a1
                } else {
                    a2
                }
            })
            .collect();
        let u = solve_variable_helmholtz(
            &g,
            &bcs,
            &alpha,
            FaceAverage::Harmonic,
            0.0,
            |_| 0.0,
            LinearSolver::Direct,
        )
        .unwrap();
        let e = g.max_error(&u, exact);
        assert!(e < 1e-10, "n {n}: layered steady error {e}");

        // transient run to the same steady state (θ = 1)
        let init = |_x: [f64; 3]| 0.0;
        let prob = VariableHeatProblem {
            grid: &g,
            bcs: &bcs,
            alpha: &alpha,
            averaging: FaceAverage::Harmonic,
            init: &init,
            source: None,
            solver: LinearSolver::Auto,
        };
        let v = heat_variable_fixed(&prob, 5.0, 200, 1.0).unwrap();
        let e = g.max_error(&v, exact);
        assert!(e < 1e-5, "n {n}: layered transient error {e}");
    }
}

/// Manufactured solution with a smooth diffusivity: second order in space and
/// time, plus an adaptive-step smoke test.
#[test]
fn variable_diffusion_manufactured_second_order_2d() {
    let alpha = |x: [f64; 3]| 2.0 + (PI * x[0]).sin() * (PI * x[1]).cos();
    let exact = move |t: f64, x: [f64; 3]| (-t).exp() * (PI * x[0]).sin() * (PI * x[1]).sin();
    let src = move |t: f64, x: [f64; 3]| {
        let e = (-t).exp();
        let (sx, cx) = ((PI * x[0]).sin(), (PI * x[0]).cos());
        let (sy, cy) = ((PI * x[1]).sin(), (PI * x[1]).cos());
        let u = e * sx * sy;
        let lap = -2.0 * PI * PI * u;
        let gu = (e * PI * cx * sy, e * PI * sx * cy);
        let ga = (PI * cx * cy, -PI * sx * sy);
        -u - (alpha(x) * lap + ga.0 * gu.0 + ga.1 * gu.1)
    };
    let bcs = FaceBcs::dirichlet_zero()
        .with_axis(0, &Bc::Dirichlet(BcValue::spacetime(exact)))
        .with_axis(1, &Bc::Dirichlet(BcValue::spacetime(exact)));
    let run = |n: usize, steps: usize| {
        let g = Grid::unit(2, n, &bcs).unwrap();
        let nodes = g.sample(alpha);
        let init = move |x: [f64; 3]| exact(0.0, x);
        let prob = VariableHeatProblem {
            grid: &g,
            bcs: &bcs,
            alpha: &nodes,
            averaging: FaceAverage::Harmonic,
            init: &init,
            source: Some(&src),
            solver: LinearSolver::Auto,
        };
        let u = heat_variable_fixed(&prob, 0.3, steps, 0.5).unwrap();
        g.l2_error(&u, |x| exact(0.3, x))
    };
    let (e1, e2) = (run(11, 10), run(21, 20));
    assert!((e1 / e2).log2() > 1.8, "{e1} {e2}");

    let g = Grid::unit(2, 15, &bcs).unwrap();
    let nodes = g.sample(alpha);
    let init = move |x: [f64; 3]| exact(0.0, x);
    let prob = VariableHeatProblem {
        grid: &g,
        bcs: &bcs,
        alpha: &nodes,
        averaging: FaceAverage::Harmonic,
        init: &init,
        source: Some(&src),
        solver: LinearSolver::Auto,
    };
    let (u, st) = heat_variable_adaptive(
        &prob,
        0.3,
        &AdaptiveOptions {
            tol: 1e-6,
            dt0: 0.01,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(st.accepted >= 1);
    // the adaptive run must agree with the fixed-step Crank–Nicolson solution
    let fixed = heat_variable_fixed(&prob, 0.3, 20, 0.5).unwrap();
    let d = fixed
        .iter()
        .zip(&u)
        .fold(0.0_f64, |m, (a, b)| m.max((a - b).abs()));
    assert!(d < 5e-2, "adaptive vs fixed Crank–Nicolson mismatch {d}");
    assert!(g.l2_error(&u, |x| exact(0.3, x)) < 5e-2);
}

/// Periodic and homogeneous-Neumann domains: `∇·(α∇u)` is a discrete
/// divergence, so the mass is conserved (with the trapezoidal weights at a
/// Neumann wall) to round-off.
#[test]
fn variable_diffusion_mass_conservation_2d() {
    // periodic
    let bcs = FaceBcs::uniform(&Bc::Periodic);
    let g = Grid::unit(2, 16, &bcs).unwrap();
    let alpha = g.sample(|x| 1.0 + 0.5 * (2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).cos());
    let init = |x: [f64; 3]| 1.0 + (2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).sin();
    let prob = VariableHeatProblem {
        grid: &g,
        bcs: &bcs,
        alpha: &alpha,
        averaging: FaceAverage::Harmonic,
        init: &init,
        source: None,
        solver: LinearSolver::Auto,
    };
    let h = g.h(0) * g.h(1);
    let mass = |u: &[f64]| h * u.iter().sum::<f64>();
    let u0 = g.sample(init);
    for theta in [0.5, 1.0] {
        let u = heat_variable_fixed(&prob, 0.3, 12, theta).unwrap();
        let d = (mass(&u) - mass(&u0)).abs();
        assert!(d < 1e-10 * mass(&u0).abs(), "periodic mass drift {d}");
    }

    // homogeneous Neumann: trapezoidal (vertex) mass
    let bcs = FaceBcs::uniform(&Bc::neumann(0.0));
    let g = Grid::new(&[21, 21], &[0.0; 2], &[1.0; 2], &bcs).unwrap();
    let alpha = g.sample(|x| 1.0 + x[0] + 0.5 * (PI * x[1]).cos());
    let init = |x: [f64; 3]| 1.0 + (-40.0 * ((x[0] - 0.4).powi(2) + (x[1] - 0.6).powi(2))).exp();
    let prob = VariableHeatProblem {
        grid: &g,
        bcs: &bcs,
        alpha: &alpha,
        averaging: FaceAverage::Harmonic,
        init: &init,
        source: None,
        solver: LinearSolver::Auto,
    };
    let weight = |p: usize| {
        let ijk = g.ijk(p);
        let mut w = 1.0;
        for a in 0..2 {
            let n = g.n[a];
            if ijk[a] == 0 || ijk[a] == n - 1 {
                w *= 0.5;
            }
        }
        w
    };
    let h = g.h(0) * g.h(1);
    let mass = |u: &[f64]| {
        h * u
            .iter()
            .enumerate()
            .map(|(p, v)| weight(p) * v)
            .sum::<f64>()
    };
    let u0 = g.sample(init);
    for theta in [0.5, 1.0] {
        let u = heat_variable_fixed(&prob, 0.3, 12, theta).unwrap();
        let d = (mass(&u) - mass(&u0)).abs();
        assert!(d < 1e-10 * mass(&u0).abs(), "Neumann mass drift {d}");
    }
}

#[test]
fn variable_diffusion_rejects_bad_input() {
    let bcs = FaceBcs::uniform(&Bc::neumann(0.0));
    let g = Grid::unit(1, 9, &bcs).unwrap();
    let init = |_: [f64; 3]| 0.0;
    let alpha = vec![1.0; g.len() - 1];
    let prob = VariableHeatProblem {
        grid: &g,
        bcs: &bcs,
        alpha: &alpha,
        averaging: FaceAverage::Harmonic,
        init: &init,
        source: None,
        solver: LinearSolver::Auto,
    };
    assert!(heat_variable_fixed(&prob, 0.1, 2, 0.5).is_err());

    let alpha = vec![-1.0; g.len()];
    let prob = VariableHeatProblem {
        grid: &g,
        bcs: &bcs,
        alpha: &alpha,
        averaging: FaceAverage::Harmonic,
        init: &init,
        source: None,
        solver: LinearSolver::Auto,
    };
    assert!(heat_variable_fixed(&prob, 0.1, 2, 0.5).is_err());
}
