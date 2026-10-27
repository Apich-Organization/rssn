//! Cut cells, convection schemes, LES and buoyancy in the MAC solver.

use rssn::sim::boundary::{Face, FluidBoundary, FluidFaceBcs};
use rssn::sim::cutcell::{open_fraction, open_length_fraction, simplex_positive_fraction};
use rssn::sim::ns_mac::{
    ConvLimiter, Convection, CutCellOptions, ImmersedBody, LesModel, NsConfig, NsSolver, ScalarBc, Thermal,
};
use std::f64::consts::PI;

#[test]
fn simplex_clipping_matches_closed_forms() {
    // triangle with one positive vertex: f^2 / ((f - g)(f - h))
    let f = simplex_positive_fraction(&[1.0, -1.0, -2.0]);
    assert!((f - 1.0 / (2.0 * 3.0)).abs() < 1e-9, "{f}");
    // tetrahedron cut by x1 + x2 + x3 + x4 type planes: symmetric two/two split is one half
    let g = simplex_positive_fraction(&[1.0, 1.0, -1.0, -1.0]);
    assert!((g - 0.5).abs() < 1e-9, "{g}");
    assert!((simplex_positive_fraction(&[1.0, 2.0, 3.0]) - 1.0).abs() < 1e-12);
    assert!(simplex_positive_fraction(&[-1.0, -2.0, -3.0]).abs() < 1e-12);
}

#[test]
fn open_fractions_are_exact_for_planes_and_converge_for_curved_bodies() {
    // plane x + 0.3 y + 0.1 z = 0.45 through the unit cube: compare with fine sampling-free formula
    let plane = |x: [f64; 3]| 0.45 - x[0] - 0.3 * x[1] - 0.1 * x[2];
    // area of the 2D square [0,1]^2 with x + 0.3 y < 0.45 (z fixed 0): integrate in y
    let mut exact = 0.0;
    let nq = 20_000;
    for k in 0..nq {
        let y = (k as f64 + 0.5) / nq as f64;
        exact += (0.45 - 0.3 * y).clamp(0.0, 1.0) / nq as f64;
    }
    let a = open_fraction(&plane, [0.0; 3], [1.0, 1.0, 0.0], &[0, 1], 0);
    assert!((a - exact).abs() < 1e-9, "face {a} vs {exact}");
    // sphere of radius 0.3 about the cell centre (fully inside the box): open fraction 1 - 4/3 pi r^3
    let sphere = |x: [f64; 3]| ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2) + (x[2] - 0.5).powi(2)).sqrt() - 0.3;
    let ve = 1.0 - 4.0 / 3.0 * PI * 0.027;
    let err_at = |depth: usize| (open_fraction(&sphere, [0.0; 3], [1.0; 3], &[0, 1, 2], depth) - ve).abs();
    let (e2, e3, e4) = (err_at(2), err_at(3), err_at(4));
    assert!(e4 < 3e-3 && e3 > 2.5 * e4 && e2 > 2.5 * e3, "sphere errors {e2} {e3} {e4} (second order in the leaf size)");
    // a disc: 2D area fraction and the length fraction on a chord
    let disc = |x: [f64; 3]| ((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2)).sqrt() - 0.4;
    let ar = open_fraction(&disc, [0.0; 3], [1.0, 1.0, 0.0], &[0, 1], 5);
    assert!((ar - (1.0 - PI * 0.16)).abs() < 1e-3, "disc {ar}");
    let l = open_length_fraction(&disc, [0.0, 0.5, 0.0], 0, 1.0);
    assert!((l - 0.2).abs() < 1e-12, "chord {l}");
}

fn cylinder_solver(re: f64, h: f64, opts: CutCellOptions) -> (NsSolver, [f64; 2]) {
    // diameter 1, centre at (14, 16) in a 40 x 32 box (7 diameters is not enough upstream at these Reynolds numbers)
    let (lx, ly) = (40.0, 32.0);
    let n = [(lx / h).round() as usize, (ly / h).round() as usize];
    let mut cfg = NsConfig::new(&n, &[lx, ly], 1.0 / re).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::XMin, FluidBoundary::inflow([1.0, 0.0, 0.0]))
        .with(Face::XMax, FluidBoundary::outflow())
        .with(Face::YMin, FluidBoundary::Slip)
        .with(Face::YMax, FluidBoundary::Slip);
    cfg.cut_cell = opts;
    cfg.solver = rssn::sim::structured::LinearSolver::Direct;
    cfg.body = Some(ImmersedBody {
        level_set: Box::new(|x| ((x[0] - 14.0).powi(2) + (x[1] - 16.0).powi(2)).sqrt() - 0.5),
        velocity: Box::new(|_| [0.0; 3]),
    });
    let dt = 0.25 * h;
    let s = NsSolver::new(cfg, dt, &|_| [1.0, 0.0, 0.0]).unwrap();
    (s, [dt, h])
}

#[test]
fn cut_cell_volume_fractions_and_apertures_integrate_to_the_body() {
    let (s, [_, h]) = cylinder_solver(20.0, 0.2, CutCellOptions::default());
    let n = s.n();
    let mut area = 0.0;
    for j in 0..n[1] {
        for i in 0..n[0] {
            area += (1.0 - s.volume_fraction([i, j, 0])) * h * h;
        }
    }
    assert!((area - PI * 0.25).abs() < 2e-3, "body area {area}");
    // vertical apertures through the centre line are exactly the chord complements
    let full = s.face_aperture(0, [0, 80, 0]);
    assert!((full - 1.0).abs() < 1e-12);
}

fn cylinder_drag(re: f64, h: f64, t_end: f64) -> f64 {
    let (mut s, [dt, _]) = cylinder_solver(re, h, CutCellOptions::default());
    s.run((t_end / dt).round() as usize).unwrap();
    let f = s.force_on_region([11.5, 13.5, 0.0], [16.5, 18.5, 0.0]);
    2.0 * f[0]
}

#[test]
fn cylinder_drag_matches_literature_re20_and_re40() {
    // Dennis & Chang (1970): Cd = 2.05 (Re = 20), 1.52 (Re = 40); other solvers 2.0-2.1 and 1.5-1.6.
    // The finite box (14 D upstream, 26 D downstream, free-slip sides at +-16 D) and h = D / 5 add a few per cent.
    let c20 = cylinder_drag(20.0, 0.2, 45.0);
    assert!((c20 - 2.05).abs() < 0.17, "Cd(20) = {c20}");
    let c40 = cylinder_drag(40.0, 0.2, 80.0);
    assert!((c40 - 1.55).abs() < 0.2, "Cd(40) = {c40}");
    assert!(c20 > c40);
}

const R_IN: f64 = 0.5;
const R_OUT: f64 = 1.0;

fn couette_error(n: usize, nu: f64, opts: CutCellOptions) -> f64 {
    let half = 1.2;
    let omega = 1.0;
    let mut cfg = NsConfig::new(&[n, n], &[2.0 * half, 2.0 * half], nu).unwrap();
    cfg.cut_cell = opts;
    let c = half;
    cfg.body = Some(ImmersedBody {
        level_set: Box::new(move |x| {
            let r = ((x[0] - c).powi(2) + (x[1] - c).powi(2)).sqrt();
            (r - R_IN).min(R_OUT - r)
        }),
        velocity: Box::new(move |x| {
            let (dx, dy) = (x[0] - c, x[1] - c);
            let r = (dx * dx + dy * dy).sqrt();
            if r < 0.5 * (R_IN + R_OUT) { [-omega * dy, omega * dx, 0.0] } else { [0.0; 3] }
        }),
    });
    let dt = 0.4 * (2.0 * half / n as f64);
    let mut s = NsSolver::new(cfg, dt, &|_| [0.0; 3]).unwrap();
    s.run(((3.0 / nu).max(12.0) / dt) as usize).unwrap();
    let b = omega * R_IN / (1.0 / R_IN - R_IN / (R_OUT * R_OUT));
    let a = -b / (R_OUT * R_OUT);
    let mut worst: f64 = 0.0;
    for ix in 0..60 {
        for iy in 0..60 {
            let x = [c - 0.95 + 1.9 * ix as f64 / 59.0, c - 0.95 + 1.9 * iy as f64 / 59.0, 0.0];
            let (dx, dy) = (x[0] - c, x[1] - c);
            let r = (dx * dx + dy * dy).sqrt();
            if r < R_IN + 0.08 || r > R_OUT - 0.08 {
                continue;
            }
            let ut = a * r + b / r;
            let (ue, ve) = (-ut * dy / r, ut * dx / r);
            worst = worst.max((s.sample_velocity(0, x) - ue).abs()).max((s.sample_velocity(1, x) - ve).abs());
        }
    }
    worst
}

#[test]
fn taylor_couette_is_second_order_with_conservative_cut_cell_convection() {
    // Re = Omega R^2 / nu = 10: the centripetal convective term is large
    let nu = 0.025;
    let cons = CutCellOptions::default();
    let plain = CutCellOptions { conservative: false, beta_min: 0.0 };
    let (e1, e2) = (couette_error(40, nu, cons), couette_error(80, nu, cons));
    let order = (e1 / e2).log2();
    assert!(order > 1.6 && e2 < 0.02, "conservative: order {order} ({e1}, {e2})");
    let (p1, p2) = (couette_error(40, nu, plain), couette_error(80, nu, plain));
    assert!(e2 <= 1.5 * p2 + 1e-4 && e1 <= 1.5 * p1 + 1e-4, "plain {p1} {p2} vs conservative {e1} {e2}");
}

fn tg_error(n: usize, conv: Convection) -> f64 {
    let nu = 0.01;
    let mut cfg = NsConfig::new(&[n, n], &[2.0 * PI, 2.0 * PI], nu).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::Periodic);
    cfg.convection = conv;
    let h = 2.0 * PI / n as f64;
    let steps = (1.0 / (0.3 * h)).round() as usize;
    let dt = 1.0 / steps as f64;
    let exact = |x: [f64; 3], t: f64| {
        let d = (-2.0 * nu * t).exp();
        [-x[0].cos() * x[1].sin() * d, x[0].sin() * x[1].cos() * d]
    };
    let mut s = NsSolver::new(cfg, dt, &|x| [exact(x, 0.0)[0], exact(x, 0.0)[1], 0.0]).unwrap();
    s.run(steps).unwrap();
    let mut worst: f64 = 0.0;
    for c in 0..2 {
        for j in 0..n + usize::from(c == 1) {
            for i in 0..n + usize::from(c == 0) {
                let x = [
                    (i as f64 + if c == 0 { 0.0 } else { 0.5 }) * h,
                    (j as f64 + if c == 1 { 0.0 } else { 0.5 }) * h,
                    0.0,
                ];
                worst = worst.max((s.sample_velocity(c, x) - exact(x, 1.0)[c]).abs());
            }
        }
    }
    worst
}

#[test]
fn quick_and_tvd_convection_converge_on_the_taylor_green_vortex() {
    for conv in [Convection::Quick, Convection::Tvd(ConvLimiter::VanLeer)] {
        let (e1, e2) = (tg_error(16, conv), tg_error(32, conv));
        let order = (e1 / e2).log2();
        assert!(order > 1.4 && e2 < 8e-3, "{conv:?}: order {order} ({e1}, {e2})");
    }
}

#[test]
fn tvd_convection_keeps_an_advected_front_monotone_where_central_overshoots() {
    // a uniform stream carries a step of temperature around a periodic box
    let run = |conv: Convection| {
        let n = 64;
        let mut cfg = NsConfig::new(&[n, 4], &[1.0, 1.0 / 16.0], 1e-3).unwrap();
        cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::Periodic);
        cfg.convection = conv;
        cfg.thermal = Some(Thermal {
            kappa: 0.0,
            gbeta: [0.0; 3],
            t_ref: 0.0,
            bcs: [ScalarBc::Periodic; 6],
            init: Box::new(|x| if (x[0] - 0.5).abs() < 0.25 { 1.0 } else { 0.0 }),
        });
        let mut s = NsSolver::new(cfg, 0.004, &|_| [1.0, 0.0, 0.0]).unwrap();
        s.run(100).unwrap();
        let t = s.temperature().unwrap();
        let mass: f64 = t.iter().sum();
        let mx = t.iter().fold(f64::NEG_INFINITY, |m, v| m.max(*v));
        let mn = t.iter().fold(f64::INFINITY, |m, v| m.min(*v));
        (mn, mx, mass)
    };
    let (cmn, cmx, _) = run(Convection::Central);
    let (tmn, tmx, mass) = run(Convection::Tvd(ConvLimiter::VanLeer));
    assert!(tmn > -1e-9 && tmx < 1.0 + 1e-9, "tvd range {tmn} {tmx}");
    assert!(cmn < -0.02 || cmx > 1.02, "central range {cmn} {cmx} should oscillate");
    assert!((mass - 32.0 * 4.0).abs() < 1e-9, "temperature is conserved: {mass}");
}

#[test]
fn smagorinsky_and_wale_eddy_viscosity_have_the_analytic_values() {
    // simple shear u = gamma y in a periodic-x box between moving walls: Smagorinsky nu_t = (Cs Delta)^2 gamma, WALE = 0
    let gamma = 3.0;
    let make = |les: LesModel| {
        let n = [16, 16];
        let mut cfg = NsConfig::new(&n, &[1.0, 1.0], 1e-2).unwrap();
        cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
            .with(Face::XMin, FluidBoundary::Periodic)
            .with(Face::XMax, FluidBoundary::Periodic)
            .with(Face::YMax, FluidBoundary::moving_wall([gamma, 0.0, 0.0]));
        cfg.les = les;
        NsSolver::new(cfg, 1e-3, &|x| [gamma * x[1], 0.0, 0.0]).unwrap()
    };
    let cs = 0.17;
    let s = make(LesModel::Smagorinsky { cs });
    let nut = s.eddy_viscosity();
    let delta = 1.0 / 16.0;
    // interior cells (away from the walls where the ghost layer enters)
    let v = nut[8 + 16 * 8];
    assert!((v - (cs * delta).powi(2) * gamma).abs() < 1e-9, "smagorinsky {v}");
    let w = make(LesModel::Wale { cw: 0.325 });
    assert!(w.eddy_viscosity()[8 + 16 * 8].abs() < 1e-12, "WALE must vanish in pure shear");
}

#[test]
fn smagorinsky_dissipates_an_under_resolved_taylor_green_vortex() {
    let run = |les: LesModel| {
        let n = 16;
        let mut cfg = NsConfig::new(&[n, n, n], &[2.0 * PI; 3], 1e-4).unwrap();
        cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::Periodic);
        cfg.les = les;
        cfg.convection = Convection::Quick;
        let init = |x: [f64; 3]| [x[0].sin() * x[1].cos() * x[2].cos(), -x[0].cos() * x[1].sin() * x[2].cos(), 0.0];
        let mut s = NsSolver::new(cfg, 0.02, &init).unwrap();
        let e0 = s.kinetic_energy();
        s.run(100).unwrap();
        s.kinetic_energy() / e0
    };
    let e_none = run(LesModel::None);
    let e_smag = run(LesModel::Smagorinsky { cs: 0.17 });
    let e_wale = run(LesModel::Wale { cw: 0.5 });
    assert!(e_smag < e_none - 1e-3, "smagorinsky {e_smag} vs none {e_none}");
    assert!(e_wale <= e_none + 1e-9 && e_wale > e_smag, "wale {e_wale}");
    assert!(e_smag > 0.0 && e_smag.is_finite());
}

fn rb_growth(ra: f64) -> f64 {
    let (pr, k) = (1.0, 3.117);
    let lx = 2.0 * PI / k;
    let (nx, ny) = (40, 20);
    let mut cfg = NsConfig::new(&[nx, ny], &[lx, 1.0], pr).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::XMin, FluidBoundary::Periodic)
        .with(Face::XMax, FluidBoundary::Periodic);
    cfg.thermal = Some(Thermal {
        kappa: 1.0,
        gbeta: [0.0, ra * pr, 0.0],
        t_ref: 0.5,
        bcs: [ScalarBc::Periodic, ScalarBc::Periodic, ScalarBc::Dirichlet(1.0), ScalarBc::Dirichlet(0.0), ScalarBc::Adiabatic, ScalarBc::Adiabatic],
        init: Box::new(move |x| 1.0 - x[1] + 1e-4 * (k * x[0]).cos() * (PI * x[1]).sin()),
    });
    let dt = 0.004;
    let mut s = NsSolver::new(cfg, dt, &|_| [0.0; 3]).unwrap();
    s.run((0.4 / dt) as usize).unwrap();
    let a0 = s.max_speed();
    let t0 = s.t;
    s.run((0.8 / dt) as usize).unwrap();
    let a1 = s.max_speed();
    (a1 / a0).ln() / (s.t - t0)
}

#[test]
fn rayleigh_benard_onset_is_at_the_classical_critical_rayleigh_number() {
    // rigid-rigid plates, Pr = 1: Ra_c = 1707.76; the growth rate crosses zero there
    let (r1, r2, r3) = (1500.0, 1800.0, 2100.0);
    let (g1, g2, g3) = (rb_growth(r1), rb_growth(r2), rb_growth(r3));
    assert!(g1 < 0.0 && g3 > 0.0, "growth rates {g1} {g2} {g3}");
    // quadratic interpolation of the zero crossing
    let rc = if g2 < 0.0 { r2 + (r3 - r2) * (-g2) / (g3 - g2) } else { r1 + (r2 - r1) * (-g1) / (g2 - g1) };
    assert!((rc - 1707.76).abs() < 0.05 * 1707.76, "estimated Ra_c = {rc} (growth {g1} {g2} {g3})");
}
