//! The normal-incidence Yee source with staggered auxiliary fields.

use rssn::sim::boundary::BoundaryCondition;
use rssn::sim::boundary::Face;
use rssn::sim::boundary::FaceBcs;
use rssn::sim::yee_fdtd::Comp;
use rssn::sim::yee_fdtd::PlaneWave;
use rssn::sim::yee_fdtd::Yee;
use rssn::sim::yee_fdtd::YeeConfig;

fn drive(time: f64) -> f64 {
    0.25 + 0.125 * time
}

fn source(
    ia: usize,
    ib: usize,
    psi: f64,
) -> PlaneWave {
    PlaneWave {
        ia,
        ib,
        psi,
        waveform: Box::new(drive),
    }
}

#[test]
fn auxiliary_probe_matches_full_staggered_updates_and_the_terminal_cell() {
    let hx = 0.25;
    let mut config = YeeConfig::new(&[12], &[hx]).unwrap();
    config.cfl = 0.5;
    let mut solver = Yee::new(config).unwrap();
    assert!(solver.plane_wave_sample(0).is_none());
    solver.set_plane_wave(source(3, 8, 0.7), 0).unwrap();
    let count = 32;
    let mut e = vec![0.0; count];
    let mut h = vec![0.0; count];
    let mut time = 0.0;
    let initial = solver.plane_wave_sample(0).unwrap();
    assert_eq!(initial.electric, 0.0);
    assert_eq!(initial.electric_time, 0.0);
    assert_eq!(initial.magnetic_time, -0.5 * solver.dt);
    assert!(solver.plane_wave_sample(count).is_none());

    for step in 0..48 {
        if step == 5 {
            solver.dt *= 0.5;
            let recorded = solver.plane_wave_sample(0).unwrap();
            assert_eq!(recorded.magnetic_time, time - solver.dt);
        }
        let ratio = solver.dt / hx;
        let t_half = time + 0.5 * solver.dt;
        for i in 0..count - 1 {
            h[i] -= ratio * (e[i + 1] - e[i]);
        }
        time += solver.dt;
        for i in 1..count {
            e[i] -= ratio * (h[i] - h[i - 1]);
        }
        e[0] = drive(time);
        solver.step();
        for i in 0..count {
            let actual = solver.plane_wave_sample(i).unwrap();
            assert_eq!(
                actual.electric.to_bits(),
                e[i].to_bits(),
                "E[{i}], step {step}"
            );
            assert_eq!(
                actual.magnetic.to_bits(),
                h[i].to_bits(),
                "H[{i}], step {step}"
            );
            assert_eq!(actual.electric_time, time);
            assert_eq!(actual.magnetic_time, t_half);
        }
    }
    assert!(e[count - 1].abs() > 0.0);

    solver.set_plane_wave(source(3, 8, 0.0), 0).unwrap();
    let reset = solver.plane_wave_sample(0).unwrap();
    assert_eq!(reset.electric, 0.0);
    assert_eq!(reset.magnetic, 0.0);
    assert_eq!(reset.electric_time, time);
    assert_eq!(reset.magnetic_time, time - 0.5 * solver.dt);
}

#[test]
fn manufactured_plane_wave_respects_electric_and_magnetic_slab_edges() {
    let [ia, ib] = [4, 12];
    let angle = 0.7_f64;
    let (s, c) = angle.sin_cos();
    let n = [20, 6, 4];
    let spacing = [0.23, 0.31, 0.27];
    for dimensions in 1..=3 {
        let mut config = YeeConfig::new(&n[..dimensions], &spacing[..dimensions]).unwrap();
        config.cfl = 0.7;
        config.bcs = FaceBcs::uniform(&BoundaryCondition::Periodic);
        let mut solver = Yee::new(config).unwrap();
        solver.set_plane_wave(source(ia, ib, angle), 45).unwrap();
        solver.run(30);
        let [nx, ny, nz] = solver.n();
        for i in 0..nx {
            let aux = solver.plane_wave_sample(i).unwrap();
            let e = if (ia..=ib).contains(&i) {
                aux.electric
            } else {
                0.0
            };
            let h = if (ia..ib).contains(&i) {
                aux.magnetic
            } else {
                0.0
            };
            for j in 0..ny {
                for k in 0..nz {
                    let p = solver.idx(i, j, k);
                    for (component, expected) in [
                        (Comp::Ex, 0.0),
                        (Comp::Ey, c * e),
                        (Comp::Ez, s * e),
                        (Comp::Hx, 0.0),
                        (Comp::Hy, -s * h),
                        (Comp::Hz, c * h),
                    ] {
                        let actual = solver.field(component)[p];
                        assert!(
                            (actual - expected).abs() < 4e-14,
                            "{dimensions}D {component:?} at ({i}, {j}, {k}): {actual} versus {expected}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn setup_checks_both_staggered_faces_and_keeps_the_previous_live_source() {
    let mut config = YeeConfig::new(&[24], &[0.25]).unwrap();
    config.bcs = FaceBcs::uniform(&BoundaryCondition::Absorbing);
    config.pml_cells = 5;
    let mut solver = Yee::new(config).unwrap();
    solver.set_plane_wave(source(6, 18, 0.0), 20).unwrap();
    solver.step();
    let before = solver.plane_wave_sample(0).unwrap();
    assert!(solver.set_plane_wave(source(5, 18, 0.0), 20).is_err());
    assert!(solver.set_plane_wave(source(6, 19, 0.0), 20).is_err());
    assert!(
        solver
            .set_plane_wave(source(6, usize::MAX, 0.0), 20)
            .is_err()
    );
    assert!(
        solver
            .set_plane_wave(source(6, 18, 0.0), usize::MAX)
            .is_err()
    );
    assert_eq!(solver.plane_wave_sample(0), Some(before));

    let mut config = YeeConfig::new(&[24, 20], &[0.25, 0.5]).unwrap();
    config.bcs = FaceBcs::dirichlet_zero().with(Face::YMin, BoundaryCondition::Absorbing);
    config.pml_cells = 5;
    let mut solver = Yee::new(config).unwrap();
    assert!(solver.set_plane_wave(source(3, 20, 0.0), 20).is_ok());
}
