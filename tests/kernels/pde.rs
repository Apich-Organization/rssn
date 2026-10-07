//! Tests of the finite-difference PDE solvers against analytic solutions,
//! including observed convergence orders.

use rssn::kernels::ode_adaptive::OdeOptions;
use rssn::kernels::pde::{
    AdvectionScheme, Boundary, Grid1D, Grid2D, PdeError, PdeProblem, PdeSolution, PoissonSolver, StiffMethod,
    advection_1d, heat_1d, heat_2d, laplace_2d, method_of_lines_1d, pde_solver, poisson_2d, wave_1d, wave_2d,
};
use std::f64::consts::PI;

fn order(e_coarse: f64, e_fine: f64) -> f64 {
    (e_coarse / e_fine).log2()
}

fn dirichlet0() -> [Boundary; 2] {
    [Boundary::dirichlet_const(0.0), Boundary::dirichlet_const(0.0)]
}

#[test]
fn heat_1d_dirichlet_is_second_order() {
    let alpha = 0.7;
    let t = 0.2;
    let exact = |x: f64| (-PI * PI * alpha * t).exp() * (PI * x).sin();
    let mut errs = Vec::new();
    for &n in &[21usize, 41, 81] {
        let g = Grid1D { x0: 0.0, x1: 1.0, n };
        let f = heat_1d(alpha, &g, &dirichlet0(), &|x| (PI * x).sin(), t, 4 * (n - 1)).unwrap();
        errs.push(f.max_error(exact));
    }
    assert!(order(errs[0], errs[1]) > 1.9, "{errs:?}");
    assert!(order(errs[1], errs[2]) > 1.9, "{errs:?}");
    assert!(errs[2] < 1e-4);
}

#[test]
fn heat_1d_neumann_cosine_mode() {
    let t = 0.1;
    let exact = |x: f64| (-PI * PI * t).exp() * (PI * x).cos();
    let bc = [Boundary::neumann_const(0.0), Boundary::neumann_const(0.0)];
    let mut errs = Vec::new();
    for &n in &[21usize, 41] {
        let g = Grid1D { x0: 0.0, x1: 1.0, n };
        let f = heat_1d(1.0, &g, &bc, &|x| (PI * x).cos(), t, 4 * (n - 1)).unwrap();
        errs.push(f.max_error(exact));
    }
    assert!(order(errs[0], errs[1]) > 1.9, "{errs:?}");
}

#[test]
fn heat_1d_inhomogeneous_neumann_flux() {
    // u = x^2/2 + t is a solution of u_t = u_xx with u_x(1) = 1, u_x(0) = 0.
    let bc = [Boundary::neumann_const(0.0), Boundary::neumann_const(1.0)];
    let g = Grid1D { x0: 0.0, x1: 1.0, n: 41 };
    let f = heat_1d(1.0, &g, &bc, &|x| 0.5 * x * x, 0.3, 60).unwrap();
    assert!(f.max_error(|x| 0.5 * x * x + 0.3) < 1e-9);
}

#[test]
fn heat_1d_is_stable_for_huge_time_step() {
    let g = Grid1D { x0: 0.0, x1: 1.0, n: 51 };
    let f = heat_1d(1.0, &g, &dirichlet0(), &|x| (PI * x).sin(), 1.0, 2).unwrap();
    assert!(f.u.iter().all(|v| v.is_finite() && v.abs() < 1.0));
}

#[test]
fn heat_2d_mixed_bc_second_order() {
    // sin(pi x) sin(pi y) e^{-2 pi^2 t}, Dirichlet in x, Neumann through cos in y.
    let t = 0.05;
    let exact = |x: f64, y: f64| (-2.0 * PI * PI * t).exp() * (PI * x).sin() * (PI * y).cos();
    let bc = [
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::neumann_const(0.0),
        Boundary::neumann_const(0.0),
    ];
    let mut errs = Vec::new();
    for &n in &[11usize, 21] {
        let g = Grid2D::unit_square(n);
        let f = heat_2d(1.0, &g, &bc, &|x, y| (PI * x).sin() * (PI * y).cos(), t, 2 * (n - 1)).unwrap();
        errs.push(f.max_error(exact));
    }
    assert!(order(errs[0], errs[1]) > 1.9, "{errs:?}");
}

#[test]
fn heat_rejects_bad_input() {
    let g = Grid1D { x0: 0.0, x1: 1.0, n: 2 };
    assert!(matches!(heat_1d(1.0, &g, &dirichlet0(), &|_| 0.0, 1.0, 1), Err(PdeError::InvalidInput(_))));
    let g = Grid1D { x0: 0.0, x1: 1.0, n: 10 };
    assert!(heat_1d(-1.0, &g, &dirichlet0(), &|_| 0.0, 1.0, 1).is_err());
    assert!(heat_1d(1.0, &g, &dirichlet0(), &|_| 0.0, 1.0, 0).is_err());
}

#[test]
fn wave_1d_dirichlet_second_order() {
    let c = 1.5;
    let t = 0.4;
    let exact = |x: f64| (PI * x).sin() * (PI * c * t).cos();
    let mut errs = Vec::new();
    for &n in &[21usize, 41, 81] {
        let g = Grid1D { x0: 0.0, x1: 1.0, n };
        let dx = g.dx();
        // Courant number 0.5
        let steps = (t / (0.5 * dx / c)).round() as usize;
        let f = wave_1d(c, &g, &dirichlet0(), &|x| (PI * x).sin(), &|_| 0.0, t, steps).unwrap();
        errs.push(f.max_error(exact));
    }
    assert!(order(errs[0], errs[1]) > 1.9, "{errs:?}");
    assert!(order(errs[1], errs[2]) > 1.9, "{errs:?}");
}

#[test]
fn wave_1d_neumann_and_initial_velocity() {
    // u = cos(pi x) sin(pi t) / pi with u_t(x,0) = cos(pi x).
    let bc = [Boundary::neumann_const(0.0), Boundary::neumann_const(0.0)];
    let g = Grid1D { x0: 0.0, x1: 1.0, n: 101 };
    let t = 0.3;
    let f = wave_1d(1.0, &g, &bc, &|_| 0.0, &|x| (PI * x).cos(), t, 200).unwrap();
    assert!(f.max_error(|x| (PI * x).cos() * (PI * t).sin() / PI) < 2e-4);
}

#[test]
fn wave_cfl_violation_is_reported() {
    let g = Grid1D { x0: 0.0, x1: 1.0, n: 11 };
    let r = wave_1d(1.0, &g, &dirichlet0(), &|x| x, &|_| 0.0, 1.0, 5);
    assert!(matches!(r, Err(PdeError::Unstable { .. })));
    let g2 = Grid2D::unit_square(11);
    let bc = [
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
    ];
    // dt = 0.1/ ... Courant 1/sqrt(2) < c dt/h * sqrt(2) check: 1D-stable limit but 2D-unstable
    let r = wave_2d(1.0, &g2, &bc, &|_, _| 0.0, &|_, _| 0.0, 1.0, 10);
    assert!(matches!(r, Err(PdeError::Unstable { .. })));
}

#[test]
fn wave_2d_standing_mode_second_order() {
    let c = 1.0;
    let w = 2.0_f64.sqrt() * PI * c;
    let t = 0.3;
    let exact = |x: f64, y: f64| (PI * x).sin() * (PI * y).sin() * (w * t).cos();
    let bc = [
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
    ];
    let mut errs = Vec::new();
    for &n in &[21usize, 41] {
        let g = Grid2D::unit_square(n);
        // Courant (c dt sqrt(2)/h) = 0.5
        let steps = (t / (0.5 * g.dx() / (2.0_f64.sqrt() * c))).ceil() as usize;
        let f = wave_2d(c, &g, &bc, &|x, y| (PI * x).sin() * (PI * y).sin(), &|_, _| 0.0, t, steps).unwrap();
        errs.push(f.max_error(exact));
    }
    assert!(order(errs[0], errs[1]) > 1.9, "{errs:?}");
}

fn dirichlet4(f: fn(f64, f64) -> f64) -> [Boundary; 4] {
    [Boundary::dirichlet(f), Boundary::dirichlet(f), Boundary::dirichlet(f), Boundary::dirichlet(f)]
}

#[test]
fn laplace_2d_matches_harmonic_function_with_both_solvers() {
    let exact = |x: f64, y: f64| (PI * x).sin() * (PI * y).sinh() / PI.sinh();
    let bc = [
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet(|x, _| (PI * x).sin()),
    ];
    let mut errs = Vec::new();
    for &n in &[17usize, 33, 65] {
        let g = Grid2D::unit_square(n);
        let a = laplace_2d(&g, &bc, PoissonSolver::SparseLu).unwrap();
        let b = laplace_2d(&g, &bc, PoissonSolver::Cg { tol: 1e-12 }).unwrap();
        let diff = a.u.iter().zip(&b.u).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
        assert!(diff < 1e-9, "solvers disagree: {diff}");
        errs.push(a.max_error(exact));
    }
    assert!(order(errs[0], errs[1]) > 1.9, "{errs:?}");
    assert!(order(errs[1], errs[2]) > 1.9, "{errs:?}");
}

#[test]
fn poisson_2d_with_neumann_sides() {
    // u = cos(pi x) cos(pi y): u_x = 0 at x = 0, 1; Dirichlet in y.
    let exact = |x: f64, y: f64| (PI * x).cos() * (PI * y).cos();
    let bc = [
        Boundary::neumann_const(0.0),
        Boundary::neumann_const(0.0),
        Boundary::dirichlet(|x, _| (PI * x).cos()),
        Boundary::dirichlet(|x, _| -(PI * x).cos()),
    ];
    let mut errs = Vec::new();
    for &n in &[17usize, 33] {
        let g = Grid2D::unit_square(n);
        let u = poisson_2d(&g, &bc, &|x, y| 2.0 * PI * PI * exact(x, y), PoissonSolver::SparseLu).unwrap();
        errs.push(u.max_error(exact));
    }
    assert!(order(errs[0], errs[1]) > 1.9, "{errs:?}");
}

#[test]
fn poisson_rejects_singular_and_unsupported_combinations() {
    let g = Grid2D::unit_square(9);
    let neu = || {
        [
            Boundary::neumann_const(0.0),
            Boundary::neumann_const(0.0),
            Boundary::neumann_const(0.0),
            Boundary::neumann_const(0.0),
        ]
    };
    assert!(matches!(poisson_2d(&g, &neu(), &|_, _| 1.0, PoissonSolver::SparseLu), Err(PdeError::Unsupported(_))));
    let mixed = [
        Boundary::neumann_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
        Boundary::dirichlet_const(0.0),
    ];
    assert!(matches!(
        poisson_2d(&g, &mixed, &|_, _| 1.0, PoissonSolver::Cg { tol: 1e-10 }),
        Err(PdeError::Unsupported(_))
    ));
}

#[test]
fn poisson_2d_polynomial_is_exact_to_roundoff() {
    // quadratic solutions are reproduced exactly by the 5-point stencil
    fn q(x: f64, y: f64) -> f64 {
        x * x + 2.0 * y * y + x * y
    }
    let g = Grid2D::unit_square(12);
    let u = poisson_2d(&g, &dirichlet4(q), &|_, _| -6.0, PoissonSolver::SparseLu).unwrap();
    assert!(u.max_error(q) < 1e-11);
}

#[test]
fn advection_orders_and_cfl() {
    let a = 1.0;
    let t = 0.5;
    let exact = |x: f64| (2.0 * PI * (x - a * t)).sin();
    let u0 = |x: f64| (2.0 * PI * x).sin();
    let mut up = Vec::new();
    let mut lw = Vec::new();
    for &n in &[50usize, 100, 200] {
        let g = Grid1D { x0: 0.0, x1: 1.0, n };
        let steps = (t * n as f64 / 0.8).ceil() as usize;
        up.push(advection_1d(a, &g, AdvectionScheme::Upwind, &u0, t, steps).unwrap().max_error(exact));
        lw.push(advection_1d(a, &g, AdvectionScheme::LaxWendroff, &u0, t, steps).unwrap().max_error(exact));
    }
    assert!((order(up[1], up[2]) - 1.0).abs() < 0.2, "{up:?}");
    assert!(order(lw[1], lw[2]) > 1.8, "{lw:?}");
    assert!(lw[2] < up[2]);
    let g = Grid1D { x0: 0.0, x1: 1.0, n: 100 };
    assert!(matches!(
        advection_1d(1.0, &g, AdvectionScheme::Upwind, &u0, 1.0, 10),
        Err(PdeError::Unstable { .. })
    ));
    // negative speed
    let f = advection_1d(-1.0, &g, AdvectionScheme::LaxWendroff, &u0, 0.25, 50).unwrap();
    assert!(f.max_error(|x| (2.0 * PI * (x + 0.25)).sin()) < 5e-3);
}

#[test]
fn method_of_lines_heat_with_forcing_both_integrators() {
    // u = e^{-t} sin(pi x) solves u_t = u_xx + (pi^2 - 1) u.
    let t = 0.5_f64;
    let exact = |x: f64| (-t).exp() * (PI * x).sin();
    let f = |_t: f64, _x: f64, u: f64, _ux: f64, uxx: f64| uxx + (PI * PI - 1.0) * u;
    let opts = OdeOptions { rtol: 1e-9, atol: 1e-11, ..OdeOptions::default() };
    for method in [StiffMethod::Radau5, StiffMethod::Bdf(5)] {
        let mut errs = Vec::new();
        for &n in &[11usize, 21, 41] {
            let g = Grid1D { x0: 0.0, x1: 1.0, n };
            let s = method_of_lines_1d(&f, &g, &dirichlet0(), &|x| (PI * x).sin(), t, method, &opts).unwrap();
            assert!(s.nfev > 0 && s.steps > 0);
            errs.push(s.field.max_error(exact));
        }
        assert!(order(errs[0], errs[1]) > 1.8, "{method:?} {errs:?}");
        assert!(order(errs[1], errs[2]) > 1.8, "{method:?} {errs:?}");
    }
}

#[test]
fn method_of_lines_neumann_and_advection_diffusion() {
    // u_t = u_xx - a u_x with a cosine initial datum and insulated ends is
    // not a closed form, but total mass changes only by the advective flux;
    // check instead pure diffusion with Neumann ends conserves the mean.
    let g = Grid1D { x0: 0.0, x1: 1.0, n: 41 };
    let bc = [Boundary::neumann_const(0.0), Boundary::neumann_const(0.0)];
    let f = |_: f64, _: f64, _: f64, _: f64, uxx: f64| uxx;
    let u0 = |x: f64| (PI * x).cos() + 1.0;
    let s = method_of_lines_1d(&f, &g, &bc, &u0, 1.0, StiffMethod::Radau5, &OdeOptions::default()).unwrap();
    let h = g.dx();
    let mean = |u: &[f64]| {
        let n = u.len();
        h * (0.5 * u[0] + 0.5 * u[n - 1] + u[1..n - 1].iter().sum::<f64>())
    };
    let m0 = mean(&g.nodes().iter().map(|&x| u0(x)).collect::<Vec<_>>());
    assert!((mean(&s.field.u) - m0).abs() < 1e-7);
    // decayed towards the mean
    assert!(s.field.u.iter().all(|v| (v - 1.0).abs() < 0.01));
}

#[test]
fn dispatcher_routes_problems() {
    let u0 = |x: f64| (PI * x).sin();
    let p = PdeProblem::Heat1D {
        alpha: 1.0,
        grid: Grid1D { x0: 0.0, x1: 1.0, n: 41 },
        bc: dirichlet0(),
        u0: &u0,
        t_end: 0.1,
        steps: 40,
    };
    match pde_solver(&p).unwrap() {
        PdeSolution::D1(f) => assert!(f.max_error(|x| (-PI * PI * 0.1).exp() * (PI * x).sin()) < 1e-3),
        PdeSolution::D2(_) => panic!("expected a 1D field"),
    }
    let src = |_: f64, _: f64| -6.0;
    let p = PdeProblem::Poisson2D {
        grid: Grid2D::unit_square(9),
        bc: dirichlet4(|x, y| x * x + 2.0 * y * y + x * y),
        f: &src,
        solver: PoissonSolver::Cg { tol: 1e-12 },
    };
    match pde_solver(&p).unwrap() {
        PdeSolution::D2(f) => assert!(f.max_error(|x, y| x * x + 2.0 * y * y + x * y) < 1e-9),
        PdeSolution::D1(_) => panic!("expected a 2D field"),
    }
}
