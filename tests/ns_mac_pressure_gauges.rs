//! Pressure references for disconnected MAC fluid components.

use std::f64::consts::TAU;

use rssn::sim::boundary::BcValue;
use rssn::sim::boundary::Face;
use rssn::sim::boundary::FluidBoundary;
use rssn::sim::ns_mac::Density;
use rssn::sim::ns_mac::ImmersedBody;
use rssn::sim::ns_mac::NsConfig;
use rssn::sim::ns_mac::NsSolver;

const NX: usize = 7;
const NY: usize = 4;

fn divided(
    n: &[usize],
    outlet: bool,
) -> NsConfig {
    let len: Vec<_> = n.iter().map(|&count| count as f64).collect();
    let mut cfg = NsConfig::new(n, &len, 0.1).unwrap();
    cfg.rotational = false;
    let solid: Vec<bool> = (0..n.iter().product()).map(|p| p % n[0] == 3).collect();
    cfg.solid = Some(solid);
    if outlet {
        cfg.bcs = cfg
            .bcs
            .with(Face::XMin, FluidBoundary::Outflow(BcValue::constant(2.0)));
    }
    cfg
}

fn initial([x, y, z]: [f64; 3]) -> [f64; 3] {
    [
        0.15 * (x + y).sin(),
        0.1 * (y - x).cos(),
        0.09 * (x + z).sin(),
    ]
}

fn assert_projected(solver: &NsSolver) {
    let summary = solver.divergence_summary().unwrap();
    assert!(summary.max_abs < 8.0e-13, "{summary:?}");
}

#[test]
fn each_closed_component_gets_an_independent_pressure_reference() {
    for n in [&[NX, NY][..], &[NX, NY, 3][..]] {
        let mut solver = NsSolver::new(divided(n, false), 0.001, &initial).unwrap();
        assert!(solver.divergence_summary().unwrap().rms > 0.05);
        solver.step().unwrap();
        assert_projected(&solver);
        assert_eq!(solver.p[0], 0.0);
        assert_eq!(solver.p[4], 0.0);
        assert!(solver.p[1].abs() > 0.1, "left pressure {:?}", solver.p);
        assert!(solver.p[5].abs() > 0.1, "right pressure {:?}", solver.p);
        solver.run(3).unwrap();
        assert_projected(&solver);
        assert_eq!(solver.p[0], 0.0);
        assert_eq!(solver.p[4], 0.0);
    }
}

#[test]
fn a_pressure_outlet_grounds_only_its_own_component() {
    let mut solver = NsSolver::new(divided(&[NX, NY], true), 0.001, &|x| {
        if x[0] < 3.0 {
            [0.0; 3]
        } else {
            initial(x)
        }
    })
    .unwrap();
    solver.step().unwrap();
    assert_projected(&solver);
    assert_eq!(solver.p[4], 0.0);
    assert!(solver.p[5].abs() > 0.1);
    solver.run(2).unwrap();
    assert_projected(&solver);
    assert!((solver.p[0] - 2.0).abs() < 5.0e-8, "{}", solver.p[0]);
    assert_eq!(solver.p[4], 0.0);
}

#[test]
fn an_outlet_behind_solid_cells_does_not_ground_a_pressure_component() {
    let mut cfg = divided(&[NX, NY], true);
    cfg.solid = Some((0..NX * NY).map(|p| p % NX == 0).collect());
    let mut solver = NsSolver::new(cfg, 0.001, &initial).unwrap();
    solver.step().unwrap();
    assert_projected(&solver);
    assert_eq!(solver.p[1], 0.0);
    assert!(solver.p[2].abs() > 0.1);
}

#[test]
fn a_periodic_seam_connects_the_two_sides_of_the_solid_column() {
    let mut cfg = divided(&[NX, NY], false);
    cfg.bcs = cfg.bcs.with_axis(0, &FluidBoundary::Periodic);
    let mut solver = NsSolver::new(cfg, 0.001, &|x| {
        [0.25 * (TAU * x[0] / NX as f64).sin(), 0.0, 0.0]
    })
    .unwrap();
    solver.step().unwrap();
    assert_projected(&solver);
    assert_eq!(solver.p[0], 0.0);
    assert!(solver.p[4].abs() > 0.1, "{}", solver.p[4]);
}

#[test]
fn a_one_way_cut_seam_inherits_the_downstream_pressure_reference() {
    let mut cfg = NsConfig::new(&[NX, NY], &[NX as f64, NY as f64], 0.1).unwrap();
    cfg.rotational = false;
    cfg.bcs = cfg.bcs.with_axis(0, &FluidBoundary::Periodic);
    cfg.body = Some(ImmersedBody {
        level_set: Box::new(|x| ((x[0] - 3.5).abs() - 0.55).min(6.8 - x[0])),
        velocity: Box::new(|_| [0.0; 3]),
    });
    let mut solver = NsSolver::new(cfg, 0.001, &|x| [0.12 + 0.02 * x[1], 0.0, 0.0]).unwrap();
    assert_eq!(solver.face_aperture(0, [0, 0, 0]), 1.0);
    assert_eq!(solver.face_aperture(0, [NX, 0, 0]), 0.0);
    solver.step().unwrap();
    assert_projected(&solver);
    assert_eq!(solver.p[4], 0.0);
    assert!(solver.p[0].abs() > 0.1);
}

fn density_config() -> NsConfig {
    let mut cfg = divided(&[NX, NY], true);
    cfg.density = Some(Density {
        init: Box::new(|x| 1.5 + 0.1 * x[0] + 0.03 * x[1]),
        transport: true,
        gravity: [0.0; 3],
    });
    cfg
}

fn assert_bits(
    a: &[f64],
    b: &[f64],
) {
    assert_eq!(a.len(), b.len());
    for (&x, &y) in a.iter().zip(b) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
}

#[test]
fn transport_rebuilds_pressure_with_the_current_density() {
    let mut automatic = NsSolver::new(density_config(), 0.001, &initial).unwrap();
    let mut explicit = NsSolver::new(density_config(), 0.001, &initial).unwrap();
    let original = automatic.density().unwrap().to_vec();
    automatic.step().unwrap();
    explicit.step().unwrap();
    let transported = automatic.density().unwrap().to_vec();
    assert_ne!(original, transported);
    explicit.set_density(&transported).unwrap();
    automatic.step().unwrap();
    explicit.step().unwrap();
    assert_projected(&automatic);
    assert_eq!(automatic.p[4], 0.0);
    assert_bits(&automatic.p, &explicit.p);
    for (a, b) in automatic.u.iter().zip(&explicit.u) {
        assert_bits(a, b);
    }
}

#[test]
fn pressure_rebuild_refreshes_changed_links_and_failed_builds_keep_density() {
    let mut cfg = density_config();
    cfg.bcs = cfg.bcs.with(Face::XMin, FluidBoundary::outflow());
    let mut solver = NsSolver::new(cfg, 0.001, &|_| [0.0; 3]).unwrap();
    let original = solver.density().unwrap().to_vec();
    let tiny = f64::MIN_POSITIVE / 64.0;
    assert!(solver.set_density(&[tiny; NX * NY]).is_err());
    assert_bits(solver.density().unwrap(), &original);
    solver.step().unwrap();
    assert_projected(&solver);

    // These coefficients fall below the pressure cutoff and isolate every cell.
    solver.set_density(&[1.0e13; NX * NY]).unwrap();
    solver.step().unwrap();
    assert_projected(&solver);

    // Restoring the coefficients must discard the temporary cellwise gauges.
    solver.set_density(&original).unwrap();
    solver.u[0].fill(0.04);
    solver.step().unwrap();
    assert_projected(&solver);
    assert_eq!(solver.p[4], 0.0);
    assert!(solver.p[5].abs() > 0.1);
}
