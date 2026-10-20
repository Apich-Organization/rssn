//! DG: characteristic / hierarchical / WENO limiters, entropy-stable fluxes, 3D hexahedra.

use rssn::sim::dg::{DgBc, DgLaw, DgMesh, DgOptions, Euler1d, Euler2d, ShockLimiter, solve_dg_with};
use rssn::sim::dg3d::{Dg3Bc, Dg3Mesh, Dg3Options, Euler3d, LinearAdvection3d, euler_entropy, solve_dg3d};
use rssn::sim::euler_nd::{FvConfig, FvGrid, solve};
use std::f64::consts::PI;
use std::sync::Arc;

const GAMMA: f64 = 1.4;

#[test]
fn euler_two_point_flux_is_consistent_and_entropy_conservative() {
    let law = Euler1d { gamma: GAMMA };
    let to_cons = |rho: f64, u: f64, p: f64| [rho, rho * u, p / (GAMMA - 1.0) + 0.5 * rho * u * u];
    let (ql, qr) = (to_cons(1.3, 0.4, 2.0), to_cons(0.7, -0.9, 0.6));
    let mut f = [0.0; 3];
    let mut fe = [0.0; 3];
    law.two_point_flux(&ql, &ql, 0, &mut f);
    law.flux(&ql, 0, &mut fe);
    for c in 0..3 {
        assert!((f[c] - fe[c]).abs() < 1e-13, "consistency {c}");
    }
    let mut f_lr = [0.0; 3];
    let mut f_rl = [0.0; 3];
    law.two_point_flux(&ql, &qr, 0, &mut f_lr);
    law.two_point_flux(&qr, &ql, 0, &mut f_rl);
    for c in 0..3 {
        assert!((f_lr[c] - f_rl[c]).abs() < 1e-13, "symmetry {c}");
    }
    // Tadmor's condition: [[w]] . F# = [[psi]], psi = rho u, w the entropy variables of S = -rho s / (gamma - 1)
    let w = |rho: f64, u: f64, p: f64| {
        let s = (p / rho.powf(GAMMA)).ln();
        [(GAMMA - s) / (GAMMA - 1.0) - rho * u * u / (2.0 * p), rho * u / p, -rho / p]
    };
    let (wl, wr) = (w(1.3, 0.4, 2.0), w(0.7, -0.9, 0.6));
    let lhs: f64 = (0..3).map(|c| (wr[c] - wl[c]) * f_lr[c]).sum();
    let rhs = 0.7 * -0.9 - 1.3 * 0.4;
    assert!((lhs - rhs).abs() < 1e-12, "entropy conservation {lhs} vs {rhs}");
}

#[test]
fn characteristic_transform_diagonalises_the_euler_jacobian() {
    // L (dF/dq) R = diag(u - c, u, u + c) checked by finite differences of the flux
    let law = Euler2d { gamma: GAMMA };
    let (rho, u, v, p) = (1.2, 0.5, -0.3, 0.9);
    let q = [rho, rho * u, rho * v, p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)];
    let c = (GAMMA * p / rho).sqrt();
    for dir in 0..2 {
        let un = if dir == 0 { u } else { v };
        let expected = [un - c, un, un, un + c];
        for k in 0..4 {
            let mut w = [0.0; 4];
            w[k] = 1.0;
            let mut r = [0.0; 4];
            law.from_characteristic(&q, dir, &w, &mut r);
            // J r by differencing the flux
            let eps = 1e-6;
            let (mut qp, mut qm) = (q, q);
            for i in 0..4 {
                qp[i] += eps * r[i];
                qm[i] -= eps * r[i];
            }
            let (mut fp, mut fm) = ([0.0; 4], [0.0; 4]);
            law.flux(&qp, dir, &mut fp);
            law.flux(&qm, dir, &mut fm);
            let jr: Vec<f64> = (0..4).map(|i| (fp[i] - fm[i]) / (2.0 * eps)).collect();
            for i in 0..4 {
                assert!((jr[i] - expected[k] * r[i]).abs() < 1e-6, "dir {dir} eigenvector {k} comp {i}: {} vs {}", jr[i], expected[k] * r[i]);
            }
            let mut back = [0.0; 4];
            law.to_characteristic(&q, dir, &r, &mut back);
            for i in 0..4 {
                assert!((back[i] - w[i]).abs() < 1e-12, "L R = I");
            }
        }
    }
}

fn shu_osher_state(x: f64) -> Vec<f64> {
    let (rho, u, p) = if x < -4.0 { (3.857_143, 2.629_369, 10.333_333) } else { (1.0 + 0.2 * (5.0 * x).sin(), 0.0, 1.0) };
    vec![rho, rho * u, p / (GAMMA - 1.0) + 0.5 * rho * u * u]
}

fn shu_osher_reference(x: &[f64]) -> Vec<f64> {
    let grid = FvGrid::new(&[3200], &[-5.0], &[5.0]).unwrap();
    let cfg = FvConfig::euler(GAMMA);
    let init = |p: [f64; 3]| {
        let q = shu_osher_state(p[0]);
        let u = q[1] / q[0];
        [q[0], u, 0.0, 0.0, (GAMMA - 1.0) * (q[2] - 0.5 * q[0] * u * u)]
    };
    let sol = solve(&grid, &cfg, &init, 1.8).unwrap();
    x.iter()
        .map(|&xx| {
            let i = (((xx + 5.0) / 10.0 * 3200.0) as usize).min(3199);
            sol.primitive(i)[0]
        })
        .collect()
}

fn shu_osher_error(shock: ShockLimiter, p: usize, ne: usize) -> f64 {
    let mesh = DgMesh::line(ne, -5.0, 5.0, p);
    let left = shu_osher_state(-5.0);
    let bcs = [DgBc::Fixed(Arc::new(move |_, _| left.clone())), DgBc::Outflow, DgBc::Outflow, DgBc::Outflow];
    let opt = DgOptions { shock, positivity: true, ..DgOptions::default() };
    let s = solve_dg_with(&mesh, &Euler1d { gamma: GAMMA }, &bcs, &|x| shu_osher_state(x[0]), 1.8, 0.2, &opt).unwrap();
    // L1 distance of the density from the fine reference
    let mut xs = Vec::new();
    let mut rho = Vec::new();
    for e in 0..mesh.elements() {
        for nd in 0..mesh.nodes_per_element() {
            xs.push(s.coord(e, nd)[0]);
            rho.push(s.value(e, nd, 0));
        }
    }
    let reference = shu_osher_reference(&xs);
    rho.iter().zip(&reference).map(|(a, b)| (a - b).abs()).sum::<f64>() / rho.len() as f64
}

#[test]
fn shu_osher_with_characteristic_hierarchical_and_weno_limiters() {
    let tvb = shu_osher_error(ShockLimiter::Tvb { m: 50.0 }, 2, 160);
    let chr = shu_osher_error(ShockLimiter::TvbCharacteristic { m: 50.0 }, 2, 160);
    let hier = shu_osher_error(ShockLimiter::Hierarchical, 2, 160);
    let weno = shu_osher_error(ShockLimiter::Weno { m: 50.0 }, 2, 160);
    println!("shu-osher L1 vs reference: tvb {tvb:.4} char {chr:.4} hier {hier:.4} weno {weno:.4}");
    for (name, e) in [("tvb", tvb), ("characteristic", chr), ("hierarchical", hier), ("weno", weno)] {
        assert!(e < 0.2, "{name}: L1 error {e}");
    }
    // the oscillatory entropy wave is resolved much better with characteristic limiting than with component-wise TVB
    assert!(chr < 0.6 * tvb && chr < 0.08, "characteristic {chr} vs component-wise {tvb}");
    assert!(weno < tvb && hier < tvb, "weno {weno} / hierarchical {hier} vs tvb {tvb}");
    // refinement helps
    let fine = shu_osher_error(ShockLimiter::TvbCharacteristic { m: 50.0 }, 2, 320);
    assert!(fine < chr, "refinement {fine} vs {chr}");
}

#[test]
fn double_mach_reflection_sanity_on_a_coarse_mesh() {
    let law = Euler2d { gamma: GAMMA };
    let post = vec![8.0, 8.0 * 8.25 * (PI / 6.0).cos(), -8.0 * 8.25 * (PI / 6.0).sin(), 116.5 / (GAMMA - 1.0) + 0.5 * 8.0 * 8.25 * 8.25];
    let pre = vec![1.4, 0.0, 0.0, 1.0 / (GAMMA - 1.0)];
    let (post_i, pre_i) = (post.clone(), pre.clone());
    let init = move |x: [f64; 2]| if x[0] < 1.0 / 6.0 + x[1] / 3.0_f64.sqrt() { post_i.clone() } else { pre_i.clone() };
    let (post_b, post_t, pre_t) = (post.clone(), post.clone(), pre.clone());
    let bottom = DgBc::Custom(Arc::new(move |_t, x, qi| {
        if x[0] < 1.0 / 6.0 {
            post_b.clone()
        } else {
            let mut q = qi.to_vec();
            q[2] = -q[2];
            q
        }
    }));
    let top = DgBc::Custom(Arc::new(move |t, x, _| {
        if x[0] < 1.0 / 6.0 + (1.0 + 20.0 * t) / 3.0_f64.sqrt() { post_t.clone() } else { pre_t.clone() }
    }));
    let left = DgBc::Fixed(Arc::new(move |_, _| post.clone()));
    let bcs = [left, DgBc::Outflow, bottom, top];
    let mesh = DgMesh::rect([96, 24], [0.0, 0.0], [4.0, 1.0], 2);
    let opt = DgOptions { shock: ShockLimiter::TvbCharacteristic { m: 50.0 }, positivity: true, ..DgOptions::default() };
    let s = solve_dg_with(&mesh, &law, &bcs, &init, 0.2, 0.15, &opt).unwrap();
    let (mut rho_max, mut rho_min, mut p_min) = (0.0_f64, f64::INFINITY, f64::INFINITY);
    for e in 0..mesh.elements() {
        for nd in 0..mesh.nodes_per_element() {
            let q: Vec<f64> = (0..4).map(|v| s.value(e, nd, v)).collect();
            let p = 0.4 * (q[3] - 0.5 * (q[1] * q[1] + q[2] * q[2]) / q[0]);
            assert!(q.iter().all(|v| v.is_finite()) && p > 0.0 && q[0] > 0.0);
            rho_max = rho_max.max(q[0]);
            rho_min = rho_min.min(q[0]);
            p_min = p_min.min(p);
        }
    }
    // reference values: undisturbed gas 1.4, post-shock 8, peak density of the Mach stem ~ 20 on fine grids
    assert!(rho_min > 1.3 && rho_max > 8.0 && rho_max < 26.0, "density range {rho_min} {rho_max} (p_min {p_min})");
    let probe = |x: f64, y: f64| -> f64 {
        let (mut best, mut val) = (f64::INFINITY, 0.0);
        for e in 0..mesh.elements() {
            for nd in 0..mesh.nodes_per_element() {
                let c = s.coord(e, nd);
                let d = (c[0] - x).hypot(c[1] - y);
                if d < best {
                    best = d;
                    val = s.value(e, nd, 0);
                }
            }
        }
        val
    };
    assert!((probe(0.3, 0.5) - 8.0).abs() < 0.4, "post-shock density {}", probe(0.3, 0.5));
    assert!((probe(3.9, 0.5) - 1.4).abs() < 0.05, "pre-shock density {}", probe(3.9, 0.5));
}

#[test]
fn entropy_stable_dg_converges_and_never_creates_entropy() {
    let vortex = |x: [f64; 2], t: f64| {
        let beta = 5.0;
        let wrap = |d: f64| (d + 5.0).rem_euclid(10.0) - 5.0;
        let (dx, dy) = (wrap(x[0] - 5.0 - t), wrap(x[1] - 5.0 - t));
        let r2 = dx * dx + dy * dy;
        let e = ((1.0 - r2) / 2.0).exp();
        let u = 1.0 - beta / (2.0 * PI) * e * dy;
        let v = 1.0 + beta / (2.0 * PI) * e * dx;
        let temp = 1.0 - (GAMMA - 1.0) * beta * beta / (8.0 * GAMMA * PI * PI) * (1.0 - r2).exp();
        let rho = temp.powf(1.0 / (GAMMA - 1.0));
        let p = rho.powf(GAMMA);
        vec![rho, rho * u, rho * v, p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)]
    };
    let law = Euler2d { gamma: GAMMA };
    let per = [DgBc::Periodic, DgBc::Periodic, DgBc::Periodic, DgBc::Periodic];
    let opt = DgOptions { entropy_stable: true, ..DgOptions::default() };
    let run = |ne: usize| {
        let mesh = DgMesh::rect([ne, ne], [0.0; 2], [10.0; 2], 3);
        let cfl = 0.3 * (6.0 / ne as f64).powf(1.0 / 3.0);
        let s = solve_dg_with(&mesh, &law, &per, &|x| vortex(x, 0.0), 1.0, cfl, &opt).unwrap();
        s.l2_error(0, |x| vortex(x, 1.0)[0])
    };
    let (e1, e2) = (run(6), run(12));
    assert!((e1 / e2).log2() > 2.8, "entropy-stable vortex order: {e1} {e2}");
    // the discrete entropy of the exact vortex data does not grow
    let mesh = DgMesh::rect([8, 8], [0.0; 2], [10.0; 2], 3);
    let entropy = |s: &rssn::sim::dg::DgSolution| -> f64 {
        let n = s.mesh.nodes_per_element();
        let mut tot = 0.0;
        for e in 0..s.mesh.elements() {
            for nd in 0..n {
                let q: Vec<f64> = (0..4).map(|v| s.value(e, nd, v)).collect();
                let p = 0.4 * (q[3] - 0.5 * (q[1] * q[1] + q[2] * q[2]) / q[0]);
                let sp = (p / q[0].powf(GAMMA)).ln();
                let ij = [nd % (s.mesh.p + 1), nd / (s.mesh.p + 1)];
                let w = s.w[ij[0]] * s.w[ij[1]];
                tot += -w * q[0] * sp / (GAMMA - 1.0);
            }
        }
        tot
    };
    let s0 = solve_dg_with(&mesh, &law, &per, &|x| vortex(x, 0.0), 0.0, 0.2, &opt).unwrap();
    let s1 = solve_dg_with(&mesh, &law, &per, &|x| vortex(x, 0.0), 2.0, 0.2, &opt).unwrap();
    let (a, b) = (entropy(&s0), entropy(&s1));
    assert!(b <= a + 1e-9 * a.abs(), "entropy grew from {a} to {b}");
}

#[test]
fn dg_3d_advection_converges_with_order_p_plus_one() {
    let law = LinearAdvection3d { a: [1.0, 0.5, -0.5] };
    let per = [Dg3Bc::Periodic, Dg3Bc::Periodic, Dg3Bc::Periodic, Dg3Bc::Periodic, Dg3Bc::Periodic, Dg3Bc::Periodic];
    let f = |x: [f64; 3], t: f64| (2.0 * PI * (x[0] - t)).sin() * (2.0 * PI * (x[1] - 0.5 * t)).cos() * (2.0 * PI * (x[2] + 0.5 * t)).cos();
    let run = |ne: usize, p: usize| {
        let mesh = Dg3Mesh { ne: [ne; 3], lo: [0.0; 3], hi: [1.0; 3], p };
        let s = solve_dg3d(&mesh, &law, &per, &|x| vec![f(x, 0.0)], 0.5, 0.3, &Dg3Options::default()).unwrap();
        s.l2_error(0, |x| f(x, 0.5))
    };
    let (e1, e2) = (run(3, 2), run(6, 2));
    assert!((e1 / e2).log2() > 2.6, "p=2 order {}", (e1 / e2).log2());
    let (g1, g2) = (run(2, 3), run(4, 3));
    assert!((g1 / g2).log2() > 3.4, "p=3 order {}", (g1 / g2).log2());
}

#[test]
fn dg_3d_euler_vortex_converges_and_entropy_stable_variant_is_stable() {
    let vortex = |x: [f64; 3], t: f64| {
        let beta = 5.0;
        let wrap = |d: f64| (d + 5.0).rem_euclid(10.0) - 5.0;
        let (dx, dy) = (wrap(x[0] - 5.0 - t), wrap(x[1] - 5.0 - t));
        let r2 = dx * dx + dy * dy;
        let e = ((1.0 - r2) / 2.0).exp();
        let u = 1.0 - beta / (2.0 * PI) * e * dy;
        let v = 1.0 + beta / (2.0 * PI) * e * dx;
        let temp = 1.0 - (GAMMA - 1.0) * beta * beta / (8.0 * GAMMA * PI * PI) * (1.0 - r2).exp();
        let rho = temp.powf(1.0 / (GAMMA - 1.0));
        let p = rho.powf(GAMMA);
        vec![rho, rho * u, rho * v, 0.0, p / (GAMMA - 1.0) + 0.5 * rho * (u * u + v * v)]
    };
    let per = [Dg3Bc::Periodic, Dg3Bc::Periodic, Dg3Bc::Periodic, Dg3Bc::Periodic, Dg3Bc::Periodic, Dg3Bc::Periodic];
    let law = Euler3d { gamma: GAMMA };
    for es in [false, true] {
        let opt = Dg3Options { entropy_stable: es, positivity: false };
        let run = |ne: usize| {
            let mesh = Dg3Mesh { ne: [ne, ne, 2], lo: [0.0; 3], hi: [10.0, 10.0, 2.0], p: 2 };
            let s = solve_dg3d(&mesh, &law, &per, &|x| vortex(x, 0.0), 0.5, 0.25, &opt).unwrap();
            let err = s.l2_error(0, |x| vortex(x, 0.5)[0]);
            (err, s)
        };
        let ((e1, _), (e2, s2)) = (run(6), run(12));
        let order = (e1 / e2).log2();
        assert!(order > 2.0, "3D Euler vortex order {order} (entropy stable: {es}) errors {e1} {e2}");
        // z-invariance of the extruded flow
        let m = s2.mesh;
        let np = m.nodes_per_element();
        let mut dz: f64 = 0.0;
        for e in 0..m.elements() {
            let n = (s2.mesh.p + 1) * (s2.mesh.p + 1);
            for nd in 0..n {
                dz = dz.max((s2.value(e, nd, 0) - s2.value(e, nd + n, 0)).abs());
            }
            let _ = np;
        }
        assert!(dz < 1e-10, "extruded flow stays z-invariant: {dz}");
        // entropy bound: total entropy does not increase for the ES scheme
        if es {
            let ent = |s: &rssn::sim::dg3d::Dg3Solution| s.integral(|q| -q[0] * euler_entropy(GAMMA, 3, q) / (GAMMA - 1.0));
            let mesh = Dg3Mesh { ne: [6, 6, 2], lo: [0.0; 3], hi: [10.0, 10.0, 2.0], p: 2 };
            let s0 = solve_dg3d(&mesh, &law, &per, &|x| vortex(x, 0.0), 0.0, 0.25, &opt).unwrap();
            let s1 = solve_dg3d(&mesh, &law, &per, &|x| vortex(x, 0.0), 1.0, 0.25, &opt).unwrap();
            assert!(ent(&s1) <= ent(&s0) + 1e-9 * ent(&s0).abs(), "entropy grew");
        }
    }
}
