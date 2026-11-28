//! Yee FDTD tests: cavities (convergence), PML reflection, TF/SF, Fresnel.

use rssn::sim::boundary::BoundaryCondition as Bc;
use rssn::sim::boundary::FaceBcs;
use rssn::sim::yee_fdtd::Comp;
use rssn::sim::yee_fdtd::CurrentSource;
use rssn::sim::yee_fdtd::PlaneWave;
use rssn::sim::yee_fdtd::Yee;
use rssn::sim::yee_fdtd::YeeConfig;
use std::f64::consts::PI;

fn cavity_1d_error(
    n: usize,
    eps: f64,
) -> f64 {
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
    (0..n)
        .map(|i| (y.field(Comp::Ey)[i] - (PI * i as f64 * h).sin() * (w * y.t).cos()).abs())
        .fold(0.0, f64::max)
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
            y.field_mut(Comp::Ez)[p] =
                (PI * i as f64 * h / a).sin() * (PI * j as f64 * h / b).sin();
        }
    }
    y.half_step_h();
    let steps = (1.0 / y.dt).round() as usize;
    y.run(steps);
    let w = PI * (1.0 / (a * a) + 1.0 / (b * b)).sqrt();
    let mut worst: f64 = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let ex =
                (PI * i as f64 * h / a).sin() * (PI * j as f64 * h / b).sin() * (w * y.t).cos();
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
fn tfsf_1d(
    eps: Option<f64>,
    total_time: f64,
    ib: usize,
) -> (f64, f64, f64) {
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
    y.set_plane_wave(
        PlaneWave {
            ia: 60,
            ib,
            psi: 0.0,
            waveform: Box::new(pulse),
        },
        steps,
    )
    .unwrap();
    let (mut early, mut all, mut tf_peak): (f64, f64, f64) = (0.0, 0.0, 0.0);
    // round trip to the interface at x = 6 and back to ia = x 1.2
    let arrive = 8.0;
    for _ in 0..steps {
        y.step();
        let sf = y.field(Comp::Ey)[22..57]
            .iter()
            .fold(0.0_f64, |m, v| m.max(v.abs()));
        tf_peak = tf_peak.max(
            y.field(Comp::Ey)[70..190]
                .iter()
                .fold(0.0_f64, |m, v| m.max(v.abs())),
        );
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

fn decay_ratio(
    dim: usize,
    n: usize,
    pml: usize,
    steps: usize,
) -> f64 {
    let h = 1.0 / n as f64;
    let mut cfg = YeeConfig::new(&vec![n; dim], &vec![h; dim]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = pml;
    let mut y = Yee::new(cfg).unwrap();
    let c = n / 2;
    y.add_source(CurrentSource {
        cell: [c, if dim > 1 { c } else { 0 }, if dim > 2 { c } else { 0 }],
        comp: Comp::Ez,
        waveform: Box::new(gauss_wave),
    });
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
    assert!(
        hi / lo < 1.2 && (hi - e0).abs() < 0.1 * e0 && (lo - e0).abs() < 0.1 * e0,
        "energy range {lo} {hi} vs {e0}"
    );

    // PMC walls are accepted; Robin faces are not supported
    let mut cfg = YeeConfig::new(&[10], &[0.1]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::neumann(0.0));
    assert!(Yee::new(cfg).is_ok());
    let mut cfg = YeeConfig::new(&[10], &[0.1]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::robin(1.0, 1.0, 0.0));
    assert!(Yee::new(cfg).is_err());
    let mut cfg = YeeConfig::new(&[10], &[0.1]).unwrap();
    cfg.bcs =
        FaceBcs::uniform(&Bc::Periodic).with(rssn::sim::boundary::Face::XMax, Bc::dirichlet(0.0));
    assert!(Yee::new(cfg).is_err());
}

/// A PMC wall sits on the tangential-H Yee plane, half a cell inside the
/// nominal face. A PEC-PMC cavity of `n` cells is therefore the exact left
/// half of a PEC-PEC cavity of `2n - 1` cells (the doubled PEC domain).
#[test]
fn pmc_wall_is_the_mirror_half_of_a_doubled_pec_cavity_1d() {
    let n = 60;
    let h = 1.0 / n as f64;
    let n2 = 2 * n - 1;
    let mut cfg = YeeConfig::new(&[n], &[h]).unwrap();
    cfg.bcs = FaceBcs::dirichlet_zero().with(rssn::sim::boundary::Face::XMax, Bc::neumann(0.0));
    let mut a = Yee::new(cfg).unwrap();
    let mut b = Yee::new(YeeConfig::new(&[n2], &[h]).unwrap()).unwrap();
    // fundamental mode of the doubled cavity: even about its midpoint,
    // where H vanishes (the PMC plane of the folded cavity)
    let mode = |i: usize| (PI * i as f64 / n2 as f64).sin();
    for i in 0..n {
        a.field_mut(Comp::Ey)[i] = mode(i);
    }
    for i in 0..n2 {
        b.field_mut(Comp::Ey)[i] = mode(i);
    }
    a.half_step_h();
    b.half_step_h();
    a.run(400);
    b.run(400);
    let worst = (0..n)
        .map(|i| (a.field(Comp::Ey)[i] - b.field(Comp::Ey)[i]).abs())
        .fold(0.0, f64::max);
    assert!(worst < 1e-10, "PMC vs doubled PEC cavity mismatch {worst}");
}

/// 2D PMC smoke/stability test: a cavity with magnetic walls on both `y`
/// faces and PEC `x` faces is lossless, finite and energy conserving. The
/// initial condition is a PMC-compatible TM mode (the tangential `H` vanishes
/// on both magnetic-wall planes).
#[test]
fn pmc_walls_are_stable_in_2d() {
    let (nx, ny) = (32, 24);
    let h = 1.0 / 24.0;
    let mut cfg = YeeConfig::new(&[nx, ny], &[h, h]).unwrap();
    cfg.bcs = FaceBcs::dirichlet_zero().with_axis(1, &Bc::neumann(0.0));
    let mut y = Yee::new(cfg).unwrap();
    for j in 0..ny {
        for i in 0..nx {
            let p = y.idx(i, j, 0);
            let yj = PI * (j as f64 - 0.5) / (ny as f64 - 1.0);
            y.field_mut(Comp::Ez)[p] = (PI * i as f64 / nx as f64).sin() * yj.cos();
        }
    }
    y.half_step_h();
    let e0 = y.energy();
    let (mut lo, mut hi) = (e0, e0);
    for _ in 0..400 {
        y.step();
        lo = lo.min(y.energy());
        hi = hi.max(y.energy());
    }
    assert!(y.f.iter().flatten().all(|v| v.is_finite()));
    assert!(hi / lo < 1.3, "2D PMC energy range {lo} {hi} vs {e0}");
}

/// 3D PMC stability: magnetic walls on two faces, a smooth initial field,
/// energy neither grows nor diverges.
#[test]
fn pmc_faces_stay_stable_in_3d() {
    let (nx, ny, nz) = (16, 12, 8);
    let h = 1.0 / 16.0;
    let mut cfg = YeeConfig::new(&[nx, ny, nz], &[h, h, h]).unwrap();
    cfg.bcs = FaceBcs::dirichlet_zero()
        .with(rssn::sim::boundary::Face::XMax, Bc::neumann(0.0))
        .with(rssn::sim::boundary::Face::ZMin, Bc::Symmetry);
    let mut y = Yee::new(cfg).unwrap();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let p = y.idx(i, j, k);
                let (x, yv, z) = (
                    (i as f64 + 0.5) * h,
                    (j as f64 + 0.5) * h,
                    (k as f64 + 0.5) * h,
                );
                // normal E vanishes on the magnetic-wall planes
                y.field_mut(Comp::Ex)[p] = (PI * (i as f64 + 0.5) / (nx as f64 - 0.5)).sin()
                    * (PI * yv).sin()
                    * (PI * z).cos();
                y.field_mut(Comp::Ey)[p] =
                    (PI * x).sin() * (PI * j as f64 / ny as f64).sin() * (PI * z).cos();
                y.field_mut(Comp::Ez)[p] =
                    (PI * x).sin() * (PI * yv).sin() * (PI * k as f64 / nz as f64).sin();
            }
        }
    }
    y.half_step_h();
    let e0 = y.energy();
    let mut hi = e0;
    for _ in 0..300 {
        y.step();
        hi = hi.max(y.energy());
    }
    assert!(hi < 1.5 * e0, "3D PMC energy grew from {e0} to {hi}");
    let finite = y.f.iter().flatten().all(|v| v.is_finite());
    assert!(finite, "3D PMC run produced a non-finite field");
    assert_eq!(y.ohmic_loss(), 0.0);
}

/// `σ = 0` (explicitly configured, per cell) must reproduce the lossless
/// code path bit for bit and dissipate nothing.
#[test]
fn zero_sigma_reproduces_the_lossless_path_bit_for_bit() {
    let n = 100;
    let h = 1.0 / n as f64;
    let mut cfg = YeeConfig::new(&[n], &[h]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = 10;
    let mut a = Yee::new(cfg).unwrap();
    let mut cfg = YeeConfig::new(&[n], &[h]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = 10;
    cfg.sigma = Some(Box::new(|_| 0.0));
    let mut b = Yee::new(cfg).unwrap();
    for y in [&mut a, &mut b] {
        for i in 0..n {
            y.field_mut(Comp::Ey)[i] = (-((i as f64 * h - 0.5) / 0.05).powi(2)).exp();
        }
        y.half_step_h();
    }
    a.run(300);
    b.run(300);
    for c in 0..6 {
        for i in 0..n {
            assert_eq!(
                a.f[c][i].to_bits(),
                b.f[c][i].to_bits(),
                "component {c} cell {i} differs"
            );
        }
    }
    assert_eq!(b.ohmic_loss(), 0.0);
    assert!(a.sigma().is_empty());
    assert_eq!(b.sigma().len(), n);
}

/// Analytic attenuation of a monochromatic wave in a lossy medium:
/// `α = ω sqrt(ε μ / 2) sqrt(sqrt(1 + (σ/(ω ε))²) - 1)`.
#[test]
fn ohmic_conductivity_attenuates_a_plane_wave_with_the_analytic_alpha() {
    let n = 600;
    let h = 1.0 / 150.0; // domain length 4
    let omega = 2.0 * PI * 2.0; // two cycles per unit length
    let sigma = 2.0;
    // ε = μ = c = 1
    let alpha = omega / 2.0_f64.sqrt() * ((1.0 + (sigma / omega).powi(2)).sqrt() - 1.0).sqrt();
    let mut cfg = YeeConfig::new(&[n], &[h]).unwrap();
    cfg.cfl = 0.9;
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = 40;
    cfg.sigma = Some(Box::new(move |_| sigma));
    let mut y = Yee::new(cfg).unwrap();
    y.add_source(CurrentSource {
        cell: [150, 0, 0],
        comp: Comp::Ey,
        waveform: Box::new(move |t| (omega * t).sin()),
    });
    // reach the steady state (transient decay time 2 ε / σ), then measure the
    // oscillation amplitude over one period at x = 1.5 and x = 2.5
    let steps = (12.0 / y.dt).round() as usize;
    y.run(steps);
    let period_steps = (2.0 * PI / omega / y.dt).round() as usize;
    let (i1, i2) = (225, 375);
    let (mut amp1, mut amp2) = (0.0_f64, 0.0_f64);
    for _ in 0..period_steps {
        y.step();
        amp1 = amp1.max(y.field(Comp::Ey)[i1].abs());
        amp2 = amp2.max(y.field(Comp::Ey)[i2].abs());
    }
    let measured = (amp1 / amp2).ln() / ((i2 - i1) as f64 * h);
    assert!(
        (measured - alpha).abs() < 0.04 * alpha,
        "measured α = {measured}, analytic α = {alpha} (amplitudes {amp1} {amp2})"
    );
}

/// The ohmic term dissipates energy at the rate `σ E²`: a lossy PEC cavity
/// decays as `e^{-σ t}`, and `energy() + ohmic_loss()` stays at its start
/// value.
#[test]
fn ohmic_loss_accounts_for_the_damped_energy() {
    let n = 100;
    let h = 1.0 / n as f64;
    let sigma = 0.3;
    let mut cfg = YeeConfig::new(&[n], &[h]).unwrap();
    cfg.cfl = 0.5;
    cfg.sigma = Some(Box::new(move |_| sigma));
    let mut y = Yee::new(cfg).unwrap();
    for i in 0..n {
        y.field_mut(Comp::Ey)[i] = (PI * i as f64 * h).sin();
    }
    y.half_step_h();
    let e0 = y.energy();
    let steps = (2.0 / y.dt).round() as usize;
    y.run(steps);
    let e1 = y.energy();
    let loss = y.ohmic_loss();
    let damped = e0 * (-sigma * y.t).exp();
    assert!(
        (e1 - damped).abs() < 0.04 * e0,
        "energy {e1} vs damped envelope {damped}"
    );
    assert!(
        (e0 - e1 - loss).abs() < 0.04 * e0,
        "energy budget: e0 {e0}, e1 {e1}, ohmic {loss}"
    );
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
            y.field_mut(Comp::Ez)[p] = (PI * i as f64 * h).cos().powi(0)
                * (2.0 * PI * i as f64 * h / (nx as f64 * h)).cos()
                * (PI * j as f64 * h / ly).sin();
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
            let ex = (2.0 * PI * i as f64 * h / lx).cos()
                * (PI * j as f64 * h / ly).sin()
                * (w * y.t).cos();
            worst = worst.max((y.field(Comp::Ez)[y.idx(i, j, 0)] - ex).abs());
        }
    }
    assert!(worst < 2e-2, "periodic standing wave error {worst}");
}
