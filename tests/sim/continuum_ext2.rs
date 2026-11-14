//! Immersed boundaries in the MAC solver and the Re = 1000 cavity.

use rssn::sim::boundary::BcValue;
use rssn::sim::boundary::Face;
use rssn::sim::boundary::FluidBoundary;
use rssn::sim::boundary::FluidFaceBcs;
use rssn::sim::ns_mac::ImmersedBody;
use rssn::sim::ns_mac::NsConfig;
use rssn::sim::ns_mac::NsSolver;

const R_IN: f64 = 0.5;
const R_OUT: f64 = 1.0;
const OMEGA: f64 = 1.0;

fn couette(n: usize) -> f64 {
    let half = 1.2;
    let mut cfg = NsConfig::new(&[n, n], &[2.0 * half, 2.0 * half], 0.1).unwrap();
    // shift so that the cylinders are centred in the box
    let c = half;
    cfg.body = Some(ImmersedBody {
        level_set: Box::new(move |x| {
            let r = ((x[0] - c).powi(2) + (x[1] - c).powi(2)).sqrt();
            (r - R_IN).min(R_OUT - r)
        }),
        velocity: Box::new(move |x| {
            let (dx, dy) = (x[0] - c, x[1] - c);
            let r = (dx * dx + dy * dy).sqrt();
            if r < 0.5 * (R_IN + R_OUT) {
                [-OMEGA * dy, OMEGA * dx, 0.0]
            } else {
                [0.0; 3]
            }
        }),
    });
    let dt = 0.4 * (2.0 * half / n as f64);
    let mut s = NsSolver::new(cfg, dt, &|_| [0.0; 3]).unwrap();
    s.run((12.0 / dt) as usize).unwrap();
    // exact u_theta = A r + B / r with u_theta(R_IN) = OMEGA R_IN, u_theta(R_OUT) = 0
    let b = OMEGA * R_IN / (1.0 / R_IN - R_IN / (R_OUT * R_OUT));
    let a = -b / (R_OUT * R_OUT);
    let mut worst: f64 = 0.0;
    for ix in 0..60 {
        for iy in 0..60 {
            let x = [
                c - 0.95 + 1.9 * ix as f64 / 59.0,
                c - 0.95 + 1.9 * iy as f64 / 59.0,
                0.0,
            ];
            let (dx, dy) = (x[0] - c, x[1] - c);
            let r = (dx * dx + dy * dy).sqrt();
            if r < R_IN + 0.08 || r > R_OUT - 0.08 {
                continue;
            }
            let ut = a * r + b / r;
            let (ue, ve) = (-ut * dy / r, ut * dx / r);
            worst = worst
                .max((s.sample_velocity(0, x) - ue).abs())
                .max((s.sample_velocity(1, x) - ve).abs());
        }
    }
    worst
}

#[test]
fn taylor_couette_flow_between_immersed_cylinders_converges() {
    let e1 = couette(40);
    let e2 = couette(80);
    let order = (e1 / e2).log2();
    assert!(e2 < 0.02, "fine-grid error {e2}");
    assert!(order > 1.7, "observed order {order} ({e1}, {e2})");
}

#[test]
fn flow_around_a_circular_cylinder_conserves_flux() {
    let (n0, n1) = (96, 32);
    let mut cfg = NsConfig::new(&[n0, n1], &[3.0, 1.0], 0.05).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::XMin, FluidBoundary::inflow([1.0, 0.0, 0.0]))
        .with(Face::XMax, FluidBoundary::outflow())
        .with(Face::YMin, FluidBoundary::Slip)
        .with(Face::YMax, FluidBoundary::Slip);
    cfg.body = Some(ImmersedBody {
        level_set: Box::new(|x| ((x[0] - 1.0).powi(2) + (x[1] - 0.5).powi(2)).sqrt() - 0.2),
        velocity: Box::new(|_| [0.0; 3]),
    });
    let mut s = NsSolver::new(cfg, 0.01, &|_| [1.0, 0.0, 0.0]).unwrap();
    s.run(800).unwrap();
    let (fin, fout) = (s.flux_through(0, 0), s.flux_through(0, n0));
    assert!((fin - fout).abs() < 2e-3 * fin.abs(), "{fin} {fout}");
    // no flow inside the body, stagnation ahead of it, acceleration beside it
    assert!(s.sample_velocity(0, [1.0, 0.5, 0.0]).abs() < 1e-9);
    assert!(s.sample_velocity(0, [0.76, 0.5, 0.0]) < 0.5);
    assert!(s.sample_velocity(0, [1.0, 0.78, 0.0]) > 1.1);
    assert!(s.max_speed() < 3.0);
}

#[test]
fn lid_driven_cavity_re1000_matches_ghia() {
    let n = 96;
    let mut cfg = NsConfig::new(&[n, n], &[1.0, 1.0], 1.0e-3).unwrap();
    cfg.bcs = FluidFaceBcs::uniform(&FluidBoundary::wall())
        .with(Face::YMax, FluidBoundary::moving_wall([1.0, 0.0, 0.0]));
    let dt = 0.006;
    let mut s = NsSolver::new(cfg, dt, &|_| [0.0; 3]).unwrap();
    s.run((70.0 / dt) as usize).unwrap();
    // Ghia, Ghia & Shin (1982), Re = 1000
    let u_ref = [
        (0.9766, 0.65928),
        (0.9688, 0.57492),
        (0.9609, 0.51117),
        (0.9531, 0.46604),
        (0.8516, 0.33304),
        (0.7344, 0.18719),
        (0.6172, 0.05702),
        (0.5000, -0.06080),
        (0.4531, -0.10648),
        (0.2813, -0.27805),
        (0.1719, -0.38289),
        (0.1016, -0.29730),
        (0.0703, -0.22220),
        (0.0625, -0.20196),
        (0.0547, -0.18109),
    ];
    for (y, u) in u_ref {
        let v = s.sample_velocity(0, [0.5, y, 0.0]);
        assert!((v - u).abs() < 0.04, "u({y}) = {v}, Ghia {u}");
    }
    let v_ref = [
        (0.9688, -0.21388),
        (0.9531, -0.33714),
        (0.9063, -0.51550),
        (0.8594, -0.42665),
        (0.8047, -0.31966),
        (0.5000, 0.02526),
        (0.2344, 0.32235),
        (0.1563, 0.37095),
        (0.0938, 0.32627),
        (0.0703, 0.29012),
        (0.0625, 0.27485),
    ];
    for (x, v) in v_ref {
        let val = s.sample_velocity(1, [x, 0.5, 0.0]);
        assert!((val - v).abs() < 0.04, "v({x}) = {val}, Ghia {v}");
    }
    let _ = BcValue::constant(0.0);
}
