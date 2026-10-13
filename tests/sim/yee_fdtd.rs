//! Yee FDTD tests: cavities (convergence), PML reflection, TF/SF, Fresnel.

use rssn::sim::boundary::{BoundaryCondition as Bc, FaceBcs};
use rssn::sim::yee_fdtd::{Comp, CurrentSource, PlaneWave, Yee, YeeConfig};
use std::f64::consts::PI;

fn cavity_1d_error(n: usize, eps: f64) -> f64 {
    let h = 1.0 / n as f64;
    let mut cfg = YeeConfig::new(&[n], &[h]).unwrap();
    cfg.cfl = 0.9;
    cfg.eps_r = Some(Box::new(move |_| eps));
    let mut y = Yee::new(cfg).unwrap();
    for i in 0..n {
        y.field_mut(Comp::Ey)[i] = (PI * i as f64 * h).sin();
    }
    y.half_step_h();
    let steps = (1.0 / y.dt).round() as usize;
    y.run(steps);
    let w = PI / eps.sqrt();
    (0..n).map(|i| (y.field(Comp::Ey)[i] - (PI * i as f64 * h).sin() * (w * y.t).cos()).abs()).fold(0.0, f64::max)
}

#[test]
fn cavity_1d_second_order_and_dielectric_slows_the_wave() {
    let (e1, e2) = (cavity_1d_error(50, 1.0), cavity_1d_error(100, 1.0));
    assert!((e1 / e2).log2() > 1.8 && e2 < 2e-3, "{e1} {e2}");
    let d = cavity_1d_error(100, 4.0);
    assert!(d < 2e-3, "dielectric cavity {d}");
}

fn cavity_2d_error(n: usize) -> f64 {
    let h = 1.0 / n as f64;
    let (nx, ny) = (n, 3 * n / 2);
    let cfg = YeeConfig::new(&[nx, ny], &[h, h]).unwrap();
    let mut y = Yee::new(cfg).unwrap();
    let (a, b) = (1.0, 1.5);
    for j in 0..ny {
        for i in 0..nx {
            let p = y.idx(i, j, 0);
            y.field_mut(Comp::Ez)[p] = (PI * i as f64 * h / a).sin() * (PI * j as f64 * h / b).sin();
        }
    }
    y.half_step_h();
    let steps = (1.0 / y.dt).round() as usize;
    y.run(steps);
    let w = PI * (1.0 / (a * a) + 1.0 / (b * b)).sqrt();
    let mut worst: f64 = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let ex = (PI * i as f64 * h / a).sin() * (PI * j as f64 * h / b).sin() * (w * y.t).cos();
            worst = worst.max((y.field(Comp::Ez)[y.idx(i, j, 0)] - ex).abs());
        }
    }
    worst
}

#[test]
fn cavity_2d_tm_mode_second_order() {
    let (e1, e2) = (cavity_2d_error(16), cavity_2d_error(32));
    assert!((e1 / e2).log2() > 1.8 && e2 < 1e-2, "{e1} {e2}");
}

fn cavity_3d_error(n: usize) -> f64 {
    let h = 1.0 / n as f64;
    let cfg = YeeConfig::new(&[n, n, 4], &[h, h, h]).unwrap();
    let mut y = Yee::new(cfg).unwrap();
    for k in 0..4 {
        for j in 0..n {
            for i in 0..n {
                let p = y.idx(i, j, k);
                y.field_mut(Comp::Ez)[p] = (PI * i as f64 * h).sin() * (PI * j as f64 * h).sin();
            }
        }
    }
    y.half_step_h();
    // omega T = 3 pi / 2: the phase error is measured where cos(omega T) is stationary, so rounding
    // the final time to a whole number of steps does not perturb the observed order
    let w = PI * 2.0_f64.sqrt();
    let steps = (1.5 * PI / w / y.dt).round() as usize;
    y.run(steps);
    let mut worst: f64 = 0.0;
    for k in 0..4 {
        for j in 0..n {
            for i in 0..n {
                let ex = (PI * i as f64 * h).sin() * (PI * j as f64 * h).sin() * (w * y.t).cos();
                worst = worst.max((y.field(Comp::Ez)[y.idx(i, j, k)] - ex).abs());
            }
        }
    }
    worst
}

#[test]
fn cavity_3d_mode_second_order() {
    let (e1, e2) = (cavity_3d_error(8), cavity_3d_error(16));
    assert!((e1 / e2).log2() > 1.8 && e2 < 2e-2, "{e1} {e2}");
}

fn pulse(t: f64) -> f64 {
    (-((t - 3.0) / 1.0).powi(2)).exp()
}

/// Returns (max |Ey| in the scattered-field region before the reflection can
/// arrive, max |Ey| in the scattered-field region overall).
fn tfsf_1d(eps: Option<f64>, total_time: f64, ib: usize) -> (f64, f64, f64) {
    let (n, h) = (600, 0.02);
    let mut cfg = YeeConfig::new(&[n], &[h]).unwrap();
    cfg.cfl = 0.9;
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = 20;
    if let Some(e) = eps {
        cfg.eps_r = Some(Box::new(move |x| if x[0] >= 6.0 { e } else { 1.0 }));
    }
    let mut y = Yee::new(cfg).unwrap();
    let steps = (total_time / y.dt) as usize;
    y.set_plane_wave(PlaneWave { ia: 60, ib, psi: 0.0, waveform: Box::new(pulse) }, steps).unwrap();
    let (mut early, mut all, mut tf_peak): (f64, f64, f64) = (0.0, 0.0, 0.0);
    // round trip to the interface at x = 6 and back to ia = x 1.2
    let arrive = 8.0;
    for _ in 0..steps {
        y.step();
        let sf = y.field(Comp::Ey)[22..57].iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        tf_peak = tf_peak.max(y.field(Comp::Ey)[70..190].iter().fold(0.0_f64, |m, v| m.max(v.abs())));
        all = all.max(sf);
        if y.t < arrive {
            early = early.max(sf);
        }
    }
    (early, all, tf_peak)
}

#[test]
fn tfsf_plane_wave_is_confined_and_pml_reflection_is_small() {
    // vacuum with PML at both ends: the scattered field region stays quiet
    let (early, all, tf) = tfsf_1d(None, 24.0, 200);
    assert!(tf > 0.9, "incident pulse reaches the TF region ({tf})");
    assert!(early < 1e-10, "TF/SF leakage {early}");
    assert!(all < 1e-3, "PML reflection {all}");
}

#[test]
fn fresnel_reflection_from_dielectric_half_space() {
    let (early, all, _) = tfsf_1d(Some(4.0), 24.0, 450);
    assert!(early < 1e-10);
    // amplitude reflection coefficient (1 - n)/(1 + n) with n = 2
    assert!((all - 1.0 / 3.0).abs() < 0.015, "reflected amplitude {all}");
}

/// Gaussian derivative: zero time integral, so the point source leaves no static charge field.
fn gauss_wave(t: f64) -> f64 {
    let z = (t - 0.45) / 0.12;
    -2.0 * z * (-z * z).exp()
}

fn decay_ratio(dim: usize, n: usize, pml: usize, steps: usize) -> f64 {
    let h = 1.0 / n as f64;
    let mut cfg = YeeConfig::new(&vec![n; dim], &vec![h; dim]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = pml;
    let mut y = Yee::new(cfg).unwrap();
    let c = n / 2;
    y.add_source(CurrentSource { cell: [c, if dim > 1 { c } else { 0 }, if dim > 2 { c } else { 0 }], comp: Comp::Ez, waveform: Box::new(gauss_wave) });
    let mut peak: f64 = 0.0;
    for _ in 0..steps {
        y.step();
        peak = peak.max(y.energy());
    }
    y.energy() / peak
}

#[test]
fn pml_absorbs_point_source_radiation_in_1d_2d_3d() {
    // 1D uses an Ez source; the radiated pulse leaves through both PMLs
    assert!(decay_ratio(1, 200, 20, 600) < 1e-5);
    let r2 = decay_ratio(2, 80, 10, 260);
    assert!(r2 < 3e-3, "2D energy ratio {r2}");
    let r3 = decay_ratio(3, 24, 6, 90);
    assert!(r3 < 2e-2, "3D energy ratio {r3}");
}

#[test]
fn pec_cavity_conserves_energy_and_rejects_unsupported_faces() {
    let n = 40;
    let h = 1.0 / n as f64;
    let mut y = Yee::new(YeeConfig::new(&[n, n], &[h, h]).unwrap()).unwrap();
    for j in 0..n {
        for i in 0..n {
            let p = y.idx(i, j, 0);
            y.field_mut(Comp::Ez)[p] = (PI * i as f64 * h).sin() * (2.0 * PI * j as f64 * h).sin();
        }
    }
    y.half_step_h();
    let e0 = y.energy();
    let mut lo = e0;
    let mut hi = e0;
    for _ in 0..400 {
        y.step();
        lo = lo.min(y.energy());
        hi = hi.max(y.energy());
    }
    assert!(hi / lo < 1.2 && (hi - e0).abs() < 0.1 * e0 && (lo - e0).abs() < 0.1 * e0, "energy range {lo} {hi} vs {e0}");

    let mut cfg = YeeConfig::new(&[10], &[0.1]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::neumann(0.0));
    assert!(Yee::new(cfg).is_err());
    let mut cfg = YeeConfig::new(&[10], &[0.1]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Periodic).with(rssn::sim::boundary::Face::XMax, Bc::dirichlet(0.0));
    assert!(Yee::new(cfg).is_err());
}

#[test]
fn periodic_plane_wave_in_2d_travels_one_period_with_tfsf_free_initial_data() {
    // Standing wave in a periodic-x, PEC-y strip: Ez = cos(2 pi x) sin(pi y) has omega = pi sqrt(5)
    let (nx, ny) = (40, 20);
    let h = 1.0 / 20.0;
    let mut cfg = YeeConfig::new(&[nx, ny], &[h, h]).unwrap();
    cfg.bcs = FaceBcs::dirichlet_zero().with_axis(0, &Bc::Periodic);
    let mut y = Yee::new(cfg).unwrap();
    let ly = ny as f64 * h;
    for j in 0..ny {
        for i in 0..nx {
            let p = y.idx(i, j, 0);
            y.field_mut(Comp::Ez)[p] = (PI * i as f64 * h).cos().powi(0) * (2.0 * PI * i as f64 * h / (nx as f64 * h)).cos() * (PI * j as f64 * h / ly).sin();
        }
    }
    y.half_step_h();
    let steps = (0.5 / y.dt).round() as usize;
    y.run(steps);
    let lx = nx as f64 * h;
    let w = PI * (4.0 / (lx * lx) + 1.0 / (ly * ly)).sqrt();
    let mut worst: f64 = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let ex = (2.0 * PI * i as f64 * h / lx).cos() * (PI * j as f64 * h / ly).sin() * (w * y.t).cos();
            worst = worst.max((y.field(Comp::Ez)[y.idx(i, j, 0)] - ex).abs());
        }
    }
    assert!(worst < 2e-2, "periodic standing wave error {worst}");
}
