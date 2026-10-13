//! Nodal discontinuous Galerkin tests.

use rssn::sim::dg::{
    Burgers, DgBc, DgMesh, Euler1d, LinearAdvection, ModalFilter, lgl, lgl_diff_matrix, solve_dg,
};
use std::f64::consts::PI;
use std::sync::Arc;

fn periodic() -> [DgBc; 4] {
    [DgBc::Periodic, DgBc::Periodic, DgBc::Periodic, DgBc::Periodic]
}

#[test]
fn lgl_nodes_weights_and_differentiation() {
    let (x, w) = lgl(2);
    assert!((x[1]).abs() < 1e-14 && (w[0] - 1.0 / 3.0).abs() < 1e-14 && (w[1] - 4.0 / 3.0).abs() < 1e-14);
    for p in 1..=8 {
        let (x, w) = lgl(p);
        assert!((w.iter().sum::<f64>() - 2.0).abs() < 1e-13);
        // exact for degree 2p-1
        let q: f64 = x.iter().zip(&w).map(|(xi, wi)| wi * xi.powi((2 * p - 1) as i32 - 1)).sum();
        let exact = if (2 * p - 2) % 2 == 0 { 2.0 / (2 * p - 1) as f64 } else { 0.0 };
        assert!((q - exact).abs() < 1e-12, "p {p}");
        let d = lgl_diff_matrix(&x);
        for i in 0..=p {
            let dv: f64 = (0..=p).map(|j| d[i][j] * x[j].powi(p as i32)).sum();
            assert!((dv - p as f64 * x[i].powi(p as i32 - 1)).abs() < 1e-10);
        }
    }
}

fn advect_error(p: usize, ne: usize) -> f64 {
    let mesh = DgMesh::line(ne, 0.0, 1.0, p);
    let law = LinearAdvection { a: [1.0, 0.0] };
    let cfl = if p > 2 { 0.4 * (4.0 / ne as f64).powf((p as f64 - 2.0) / 3.0) } else { 0.4 };
    let s = solve_dg(&mesh, &law, &periodic(), &|x| vec![(2.0 * PI * x[0]).sin()], 1.0, cfl, None).unwrap();
    s.l2_error(0, |x| (2.0 * PI * x[0]).sin())
}

#[test]
fn linear_advection_1d_converges_with_order_p_plus_one() {
    for p in 1..=3 {
        let e1 = advect_error(p, 4);
        let e2 = advect_error(p, 8);
        let e3 = advect_error(p, 16);
        let o = (e2 / e3).log2();
        assert!(o > p as f64 + 0.5, "p {p}: order {o}, errors {e1} {e2} {e3}");
    }
}

#[test]
fn linear_advection_2d_converges() {
    let law = LinearAdvection { a: [1.0, 1.0] };
    let exact = |x: [f64; 2]| (2.0 * PI * (x[0] + x[1])).sin();
    let run = |ne: usize| {
        let mesh = DgMesh::rect([ne, ne], [0.0; 2], [1.0; 2], 2);
        let cfl = 0.4 * (3.0 / ne as f64).powf(0.0);
        let s = solve_dg(&mesh, &law, &periodic(), &|x| vec![exact(x)], 1.0, cfl, None).unwrap();
        s.l2_error(0, exact)
    };
    let (e1, e2) = (run(4), run(8));
    assert!((e1 / e2).log2() > 2.5, "{e1} {e2}");
}

fn burgers_exact(x: f64, t: f64) -> f64 {
    // u0 = 0.5 + 0.3 sin(2 pi x), characteristics x = x0 + u0 t (periodic)
    let u0 = |x0: f64| 0.5 + 0.3 * (2.0 * PI * x0).sin();
    let mut x0 = x - 0.5 * t;
    for _ in 0..100 {
        let f = x0 + u0(x0) * t - x;
        let df = 1.0 + 0.3 * 2.0 * PI * (2.0 * PI * x0).cos() * t;
        x0 -= f / df;
    }
    u0(x0)
}

#[test]
fn burgers_smooth_solution_converges_before_shock_formation() {
    let run = |ne: usize| {
        let mesh = DgMesh::line(ne, 0.0, 1.0, 2);
        let s = solve_dg(&mesh, &Burgers, &periodic(), &|x| vec![0.5 + 0.3 * (2.0 * PI * x[0]).sin()], 0.3, 0.3, None).unwrap();
        // the characteristic solution is periodic in x
        s.l2_error(0, |x| burgers_exact(x[0], 0.3))
    };
    let (e1, e2, e3) = (run(8), run(16), run(32));
    assert!((e2 / e3).log2() > 2.2 && (e1 / e2).log2() > 2.0, "{e1} {e2} {e3}");
}

#[test]
fn burgers_shock_with_filter_is_bounded_and_conservative() {
    let mesh = DgMesh::line(40, 0.0, 1.0, 3);
    let filt = Some(ModalFilter { alpha: 36.0, order: 4 });
    let s = solve_dg(&mesh, &Burgers, &periodic(), &|x| vec![(2.0 * PI * x[0]).sin()], 0.5, 0.3, filt).unwrap();
    let m0 = 0.0;
    assert!((s.integral(0) - m0).abs() < 1e-12);
    let mx = s.q.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    assert!(mx.is_finite() && mx < 1.2, "max {mx}");
}

#[test]
fn euler_density_wave_converges() {
    let law = Euler1d { gamma: 1.4 };
    let run = |ne: usize| {
        let mesh = DgMesh::line(ne, 0.0, 1.0, 2);
        let init = |x: [f64; 2]| {
            let rho = 1.0 + 0.2 * (2.0 * PI * x[0]).sin();
            vec![rho, rho, 1.0 / 0.4 + 0.5 * rho]
        };
        let s = solve_dg(&mesh, &law, &periodic(), &init, 1.0, 0.3, None).unwrap();
        s.l2_error(0, |x| 1.0 + 0.2 * (2.0 * PI * x[0]).sin())
    };
    let (e1, e2) = (run(8), run(16));
    assert!((e1 / e2).log2() > 2.4, "{e1} {e2}");
}

#[test]
fn fixed_inflow_and_outflow_boundaries() {
    let mesh = DgMesh::line(40, 0.0, 1.0, 3);
    let law = LinearAdvection { a: [1.0, 0.0] };
    let inflow = DgBc::Fixed(Arc::new(|t, _x| vec![(2.0 * t).sin()]));
    let bcs = [inflow, DgBc::Outflow, DgBc::Outflow, DgBc::Outflow];
    let t = 1.5;
    let s = solve_dg(&mesh, &law, &bcs, &|_| vec![0.0], t, 0.3, None).unwrap();
    // behind the front (x < t) u = sin(2 (t - x))
    let mut worst: f64 = 0.0;
    for e in 0..mesh.elements() {
        for nd in 0..mesh.nodes_per_element() {
            let x = s.coord(e, nd)[0];
            if x < t - 0.1 {
                worst = worst.max((s.value(e, nd, 0) - (2.0 * (t - x)).sin()).abs());
            }
        }
    }
    assert!(worst < 5e-4, "worst {worst}");
    let mut bad = bcs.clone();
    bad[0] = DgBc::Periodic;
    assert!(solve_dg(&mesh, &law, &bad, &|_| vec![0.0], 0.1, 0.3, None).is_err());
}
