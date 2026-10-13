//! P1/P2 simplex finite elements in 1D, 2D and 3D.

use rssn::sim::boundary::{BcValue, BoundaryCondition as Bc, MarkerBcs};
use rssn::sim::fem_simplex::{FemProblem, FemSpace, SimplexMesh, marker, simplex_rule, solve_heat, solve_static};
use rssn::sim::structured::LinearSolver;
use std::f64::consts::PI;

#[test]
fn grundmann_moeller_rule_integrates_monomials_exactly() {
    let fact = |m: usize| (1..=m).map(|x| x as f64).product::<f64>();
    for n in 0..=3usize {
        for s in 0..=3usize {
            let (pts, w) = simplex_rule(n, s);
            assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
            // all monomials of total degree <= 2s+1
            for a0 in 0..=(2 * s + 1) {
                for a1 in 0..=(2 * s + 1 - a0) {
                    for a2 in 0..=(2 * s + 1 - a0 - a1) {
                        let a = [a0, a1, a2];
                        if n < 2 && a2 > 0 || n < 1 && a1 > 0 {
                            continue;
                        }
                        let deg: usize = a.iter().take(n + 1).sum();
                        // monomial over barycentric coordinates 0..=n
                        let q: f64 = pts
                            .iter()
                            .zip(&w)
                            .map(|(p, wq)| wq * p.iter().zip(a.iter()).map(|(l, &e)| l.powi(e as i32)).product::<f64>())
                            .sum();
                        let exact = fact(n) * a.iter().take(n + 1).map(|&e| fact(e)).product::<f64>() / fact(n + deg);
                        assert!((q - exact).abs() < 1e-12, "n {n} s {s} a {a:?}: {q} vs {exact}");
                    }
                }
            }
        }
    }
}

fn exact(dim: usize, x: [f64; 3]) -> f64 {
    (0..dim).map(|a| x[a]).sum::<f64>().exp()
}

fn bcs_for(dim: usize) -> MarkerBcs<Bc> {
    let u = move |x: [f64; 3]| exact(dim, x);
    let dir = Bc::Dirichlet(BcValue::space(u));
    let mut b = MarkerBcs::new();
    if dim == 1 {
        return b
            .set(marker::XMIN, dir)
            .set(marker::XMAX, Bc::Robin { alpha: 2.0, beta: 1.0, g: BcValue::space(move |x| 3.0 * u(x)) });
    }
    b = b
        .set(marker::XMIN, dir.clone())
        .set(marker::XMAX, Bc::Neumann(BcValue::space(u)))
        .set(marker::YMIN, Bc::Robin { alpha: 2.0, beta: 1.0, g: BcValue::space(u) })
        .set(marker::YMAX, dir);
    if dim == 3 {
        b = b
            .set(marker::ZMIN, Bc::Neumann(BcValue::space(move |x| -u(x))))
            .set(marker::ZMAX, Bc::Neumann(BcValue::space(u)));
    }
    b
}

fn mesh(dim: usize, n: usize) -> SimplexMesh {
    match dim {
        1 => SimplexMesh::interval(n, 0.0, 1.0),
        2 => SimplexMesh::rectangle(n, n, [0.0; 2], [1.0; 2]),
        _ => SimplexMesh::boxed([n; 3], [0.0; 3], [1.0; 3]),
    }
}

fn poisson_l2(dim: usize, degree: usize, n: usize) -> f64 {
    let space = FemSpace::new(mesh(dim, n), degree).unwrap();
    let bcs = bcs_for(dim);
    // -Δu + c u = f with u = exp(sum x): f = (c - dim) u
    let c = 1.5;
    let f = move |_t: f64, x: [f64; 3]| (c - dim as f64) * exact(dim, x);
    let prob = FemProblem { k: 1.0, c, f: &f, bcs: &bcs };
    let u = solve_static(&space, &prob, 0.0, LinearSolver::Auto).unwrap();
    space.l2_error(&u, |x| exact(dim, x))
}

#[test]
fn mixed_marker_bcs_converge_with_optimal_order_p1_and_p2_in_1d_2d_3d() {
    let cases: [(usize, usize, [usize; 2]); 6] =
        [(1, 1, [8, 16]), (1, 2, [4, 8]), (2, 1, [8, 16]), (2, 2, [4, 8]), (3, 1, [4, 8]), (3, 2, [2, 4])];
    for (dim, deg, ns) in cases {
        let e1 = poisson_l2(dim, deg, ns[0]);
        let e2 = poisson_l2(dim, deg, ns[1]);
        let order = (e1 / e2).log2();
        let min = deg as f64 + 0.7;
        assert!(order > min, "dim {dim} P{deg}: order {order} ({e1}, {e2})");
    }
}

#[test]
fn bad_inputs_are_rejected() {
    assert!(FemSpace::new(SimplexMesh::interval(4, 0.0, 1.0), 3).is_err());
    let space = FemSpace::new(SimplexMesh::interval(4, 0.0, 1.0), 1).unwrap();
    let bcs = MarkerBcs::with_default(Bc::Periodic);
    let f = |_t: f64, _x: [f64; 3]| 0.0;
    let prob = FemProblem { k: 1.0, c: 1.0, f: &f, bcs: &bcs };
    assert!(solve_static(&space, &prob, 0.0, LinearSolver::Auto).is_err());
    let bcs = MarkerBcs::with_default(Bc::Robin { alpha: 1.0, beta: 0.0, g: BcValue::constant(0.0) });
    let prob = FemProblem { k: 1.0, c: 1.0, f: &f, bcs: &bcs };
    assert!(solve_static(&space, &prob, 0.0, LinearSolver::Auto).is_err());
}

#[test]
fn heat_equation_with_time_dependent_dirichlet_data_2d_p2_and_3d_p1() {
    // u = exp(-t) sin(pi x) + t y : u_t - Δu = -exp(-t) sin(pi x) + y + pi^2 exp(-t) sin(pi x) ... computed below
    let u = |t: f64, x: [f64; 3]| (-t).exp() * (PI * x[0]).sin() * (PI * x[1]).sin() + t * x[0];
    let f = move |t: f64, x: [f64; 3]| {
        // u_t = -e sin sin + x ; -Δu = 2 pi^2 e sin sin
        let e = (-t).exp() * (PI * x[0]).sin() * (PI * x[1]).sin();
        -e + x[0] + 2.0 * PI * PI * e
    };
    let bcs = MarkerBcs::with_default(Bc::Dirichlet(BcValue::spacetime(u)));
    let run = |n: usize, steps: usize| {
        let space = FemSpace::new(SimplexMesh::rectangle(n, n, [0.0; 2], [1.0; 2]), 2).unwrap();
        let prob = FemProblem { k: 1.0, c: 0.0, f: &f, bcs: &bcs };
        let v = solve_heat(&space, &prob, |x| u(0.0, x), 0.2, steps, 0.5, LinearSolver::Auto).unwrap();
        space.l2_error(&v, |x| u(0.2, x))
    };
    let (e1, e2) = (run(4, 8), run(8, 16));
    assert!((e1 / e2).log2() > 1.8 && e2 < 1e-3, "{e1} {e2}");

    // 3D: u = exp(-t) cos(pi x) with Neumann sides (zero flux in x, linear in y,z handled by symmetry)
    let u3 = |t: f64, x: [f64; 3]| (-PI * PI * t).exp() * (PI * x[0]).cos();
    let bcs3 = MarkerBcs::with_default(Bc::neumann(0.0));
    let f3 = |_t: f64, _x: [f64; 3]| 0.0;
    let space = FemSpace::new(SimplexMesh::boxed([16, 2, 2], [0.0; 3], [1.0; 3]), 1).unwrap();
    let prob = FemProblem { k: 1.0, c: 0.0, f: &f3, bcs: &bcs3 };
    let v = solve_heat(&space, &prob, |x| u3(0.0, x), 0.05, 20, 0.5, LinearSolver::Auto).unwrap();
    assert!(space.l2_error(&v, |x| u3(0.05, x)) < 6e-3);
    assert!(space.eval(&v, [0.5, 0.5, 0.5]).is_some());
}
