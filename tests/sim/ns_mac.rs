//! Incompressible Navier-Stokes projection solver tests.

use rssn::sim::boundary::{BcValue, Face, FluidBoundary, FluidFaceBcs};
use rssn::sim::ns_mac::{NsConfig, NsSolver};
use std::f64::consts::PI;

fn tg_exact(x: [f64; 3], t: f64, nu: f64) -> [f64; 3] {
    let d = (-2.0 * nu * t).exp();
    [-x[0].cos() * x[1].sin() * d, x[0].sin() * x[1].cos() * d, 0.0]
}

fn taylor_green_error(n: usize, nz: usize) -> (f64, f64) {
    let nu = 0.1;
    let (dim, len): (usize, Vec<f64>) = if nz > 1 { (3, vec![2.0 * PI, 2.0 * PI, 1.0]) } else { (2, vec![2.0 * PI, 2.0 * PI]) };
    let ns: Vec<usize> = if nz > 1 { vec![n, n, nz] } else { vec![n, n] };
    let mut cfg = NsConfig::new(&ns, &len, nu).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::Periodic);
    let h = 2.0 * PI / n as f64;
    let dt = 0.4 * h;
    let steps = (1.0 / dt).round() as usize;
    let dt = 1.0 / steps as f64;
    let mut s = NsSolver::new(cfg, dt, &|x| tg_exact(x, 0.0, nu)).unwrap();
    s.run(steps).unwrap();
    let mut worst: f64 = 0.0;
    let nn = s.n();
    for c in 0..2 {
        for j in 0..nn[1] + usize::from(c == 1) {
            for i in 0..nn[0] + usize::from(c == 0) {
                let x = [
                    (i as f64 + if c == 0 { 0.0 } else { 0.5 }) * h,
                    (j as f64 + if c == 1 { 0.0 } else { 0.5 }) * h,
                    0.0,
                ];
                let v = s.sample_velocity(c, x);
                worst = worst.max((v - tg_exact(x, 1.0, nu)[c]).abs());
            }
        }
    }
    let div = s.divergence().iter().fold(0.0_f64, |m, d| m.max(d.abs()));
    let _ = dim;
    (worst, div)
}

#[test]
fn taylor_green_vortex_second_order_in_2d_and_extruded_3d() {
    let (e1, d1) = taylor_green_error(16, 1);
    let (e2, d2) = taylor_green_error(32, 1);
    let (e3, _) = taylor_green_error(64, 1);
    assert!((e1 / e2).log2() > 1.8 && (e2 / e3).log2() > 1.8, "{e1} {e2} {e3}");
    assert!(d1 < 1e-8 && d2 < 1e-8, "divergence {d1} {d2}");
    // 3D extrusion reproduces the 2D result
    let (e3d, d3) = taylor_green_error(16, 4);
    assert!((e3d - e1).abs() < 1e-8 + 1e-6 * e1, "3D {e3d} vs 2D {e1}");
    assert!(d3 < 1e-8);
}

fn poiseuille_error(n: usize, dim3: bool) -> f64 {
    let (nu, g, hgt) = (0.1, 1.0, 1.0);
    let (ns, len): (Vec<usize>, Vec<f64>) =
        if dim3 { (vec![4, n, 4], vec![1.0, hgt, 1.0]) } else { (vec![4, n], vec![1.0, hgt]) };
    let mut cfg = NsConfig::new(&ns, &len, nu).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall()).with_axis(0, &FluidBoundary::Periodic);
    if dim3 {
        cfg.bcs = cfg.bcs.clone().with_axis(2, &FluidBoundary::Periodic);
    }
    cfg.force = Some(Box::new(move |_t, _x| [g, 0.0, 0.0]));
    let mut s = NsSolver::new(cfg, 0.1, &|_| [0.0; 3]).unwrap();
    s.run(150).unwrap();
    let mut worst: f64 = 0.0;
    for k in 0..n {
        let y = (k as f64 + 0.5) / n as f64;
        let ex = g / (2.0 * nu) * y * (hgt - y);
        worst = worst.max((s.sample_velocity(0, [0.5, y, 0.5]) - ex).abs());
    }
    worst
}

#[test]
fn poiseuille_flow_with_no_slip_walls_2d_and_3d() {
    let e1 = poiseuille_error(8, false);
    let e2 = poiseuille_error(16, false);
    assert!((e1 / e2).log2() > 1.8 && e2 < 8e-3, "{e1} {e2}");
    let e3 = poiseuille_error(16, true);
    assert!((e3 - e2).abs() < 1e-9, "3D duct with periodic z equals 2D: {e3} {e2}");
}

fn cavity(n: usize, re: f64, t_end: f64, dt: f64) -> NsSolver {
    let mut cfg = NsConfig::new(&[n, n], &[1.0, 1.0], 1.0 / re).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall()).with(Face::YMax, FluidBoundary::moving_wall([1.0, 0.0, 0.0]));
    let mut s = NsSolver::new(cfg, dt, &|_| [0.0; 3]).unwrap();
    s.run((t_end / dt).round() as usize).unwrap();
    s
}

#[test]
fn lid_driven_cavity_re100_matches_ghia_reference_points() {
    let s = cavity(64, 100.0, 30.0, 0.01);
    // Ghia, Ghia & Shin (1982), Re = 100: u along x = 0.5
    let u_ref = [
        (0.9766, 0.84123),
        (0.9688, 0.78871),
        (0.9609, 0.73722),
        (0.9531, 0.68717),
        (0.8516, 0.23151),
        (0.7344, 0.00332),
        (0.6172, -0.13641),
        (0.5000, -0.20581),
        (0.4531, -0.21090),
        (0.2813, -0.15662),
        (0.1719, -0.10150),
        (0.1016, -0.06434),
        (0.0703, -0.04775),
        (0.0625, -0.04192),
        (0.0547, -0.03717),
    ];
    for (y, u) in u_ref {
        let v = s.sample_velocity(0, [0.5, y, 0.0]);
        assert!((v - u).abs() < 0.015, "u({y}) = {v}, Ghia {u}");
    }
    // v along y = 0.5
    let v_ref = [
        (0.9688, -0.05906),
        (0.9531, -0.08864),
        (0.9063, -0.16914),
        (0.8594, -0.22445),
        (0.8047, -0.24533),
        (0.5, 0.05454),
        (0.2344, 0.17527),
        (0.1563, 0.16077),
        (0.0938, 0.12317),
        (0.0703, 0.10091),
        (0.0625, 0.09233),
    ];
    for (x, v) in v_ref {
        let val = s.sample_velocity(1, [x, 0.5, 0.0]);
        assert!((val - v).abs() < 0.015, "v({x}) = {val}, Ghia {v}");
    }
    assert!(s.divergence().iter().fold(0.0_f64, |m, d| m.max(d.abs())) < 1e-7);
}

#[test]
fn inflow_outflow_channel_keeps_parabolic_profile_and_pressure_drop() {
    let (n0, n1) = (48, 16);
    let nu = 0.1;
    let mut cfg = NsConfig::new(&[n0, n1], &[4.0, 1.0], nu).unwrap();
    let prof = BcValue::space(|x| 4.0 * x[1] * (1.0 - x[1]));
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::XMin, FluidBoundary::Inflow([prof, BcValue::constant(0.0), BcValue::constant(0.0)]))
        .with(Face::XMax, FluidBoundary::outflow());
    let mut s = NsSolver::new(cfg, 0.02, &|_| [0.0; 3]).unwrap();
    s.run(900).unwrap();
    for k in 0..n1 {
        let y = (k as f64 + 0.5) / n1 as f64;
        let ex = 4.0 * y * (1.0 - y);
        let v = s.sample_velocity(0, [3.9, y, 0.0]);
        assert!((v - ex).abs() < 5e-3, "outlet u({y}) = {v} vs {ex}");
        assert!(s.sample_velocity(1, [3.9, y, 0.0]).abs() < 5e-3);
    }
    // dp/dx = nu u'' = -0.8
    let i1 = s.cell_index([6, 8, 0]);
    let i2 = s.cell_index([36, 8, 0]);
    let dpdx = (s.p[i2] - s.p[i1]) / (30.0 * s.cell_size(0));
    assert!((dpdx + 0.8).abs() < 0.04, "dp/dx {dpdx}");
    // mass flux conserved: inflow equals outflow
    let (fin, fout) = (s.flux_through(0, 0), s.flux_through(0, n0));
    assert!((fin - fout).abs() < 1e-6 * fin.abs(), "{fin} {fout}");
}

#[test]
fn flow_past_square_obstacle_has_no_slip_wake_and_conserves_mass() {
    let (n0, n1) = (96, 32);
    let h = 1.0 / n1 as f64;
    let mut cfg = NsConfig::new(&[n0, n1], &[3.0, 1.0], 0.03).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::XMin, FluidBoundary::inflow([1.0, 0.0, 0.0]))
        .with(Face::XMax, FluidBoundary::outflow());
    let mut solid = vec![false; n0 * n1];
    for j in 0..n1 {
        for i in 0..n0 {
            let (x, y) = ((i as f64 + 0.5) * h, (j as f64 + 0.5) * h);
            if (0.75..1.05).contains(&x) && (0.3125..0.6875).contains(&y) {
                solid[i + n0 * j] = true;
            }
        }
    }
    cfg.solid = Some(solid.clone());
    let mut s = NsSolver::new(cfg, 0.01, &|_| [1.0, 0.0, 0.0]).unwrap();
    s.run(1500).unwrap();
    // mass conservation
    let (fin, fout) = (s.flux_through(0, 0), s.flux_through(0, n0));
    assert!((fin - fout).abs() < 1e-3 * fin.abs(), "{fin} {fout}");
    // zero velocity on faces of solid cells
    for j in 0..n1 {
        for i in 0..n0 {
            if solid[i + n0 * j] {
                let x = [(i as f64 + 0.5) * h, (j as f64 + 0.5) * h, 0.0];
                assert!(s.sample_velocity(0, x).abs() < 1e-12 && s.sample_velocity(1, x).abs() < 1e-12);
            }
        }
    }
    // recirculation behind the body
    let mut min_u = f64::INFINITY;
    for k in 0..40 {
        let x = 1.1 + 0.02 * k as f64;
        min_u = min_u.min(s.sample_velocity(0, [x, 0.5, 0.0]));
    }
    assert!(min_u < -1e-3, "no reverse flow found (min u {min_u})");
    assert!(s.max_speed() < 3.0);
}

#[test]
fn lid_driven_cavity_3d_is_divergence_free_and_symmetric() {
    let n = 12;
    let mut cfg = NsConfig::new(&[n, n, n], &[1.0; 3], 0.01).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall()).with(Face::YMax, FluidBoundary::moving_wall([1.0, 0.0, 0.0]));
    let mut s = NsSolver::new(cfg, 0.02, &|_| [0.0; 3]).unwrap();
    s.run(100).unwrap();
    assert!(s.divergence().iter().fold(0.0_f64, |m, d| m.max(d.abs())) < 1e-7);
    // reflection symmetry in z: u(x,y,z) = u(x,y,1-z), w antisymmetric
    for k in 0..n {
        let z = (k as f64 + 0.5) / n as f64;
        let (a, b) = (s.sample_velocity(0, [0.4, 0.8, z]), s.sample_velocity(0, [0.4, 0.8, 1.0 - z]));
        assert!((a - b).abs() < 1e-8, "u symmetry {a} {b}");
        let (c, d) = (s.sample_velocity(2, [0.4, 0.8, z]), s.sample_velocity(2, [0.4, 0.8, 1.0 - z]));
        assert!((c + d).abs() < 1e-8, "w antisymmetry {c} {d}");
    }
    assert!(s.sample_velocity(0, [0.5, 0.95, 0.5]) > 0.2);
    assert!(s.kinetic_energy() > 1e-3);
}
