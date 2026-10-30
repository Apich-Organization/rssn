//! FDTD: dispersive media, CFS-PML, Dey-Mittra subcell PEC, rotated TF/SF boxes.

use num_complex::Complex64;
use rssn::sim::boundary::BoundaryCondition as Bc;
use rssn::sim::boundary::FaceBcs;
use rssn::sim::yee_fdtd::Comp;
use rssn::sim::yee_fdtd::CurrentSource;
use rssn::sim::yee_fdtd::ObliquePlaneWave;
use rssn::sim::yee_fdtd::Pole;
use rssn::sim::yee_fdtd::Yee;
use rssn::sim::yee_fdtd::YeeConfig;
use std::f64::consts::PI;

// ------------------------------------------------------------ dispersive media

fn probe_series(
    poles: Option<Vec<Pole>>,
    dt_probe: &mut f64,
) -> Vec<f64> {
    let (h, n) = (0.01, 1300);
    let (xs, xp, xm) = (2.0, 2.5, 5.0);
    let mut cfg = YeeConfig::new(&[n], &[h]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = 25;
    cfg.cfl = 0.9;
    if let Some(p) = poles {
        cfg.dispersive = Some(Box::new(move |x| {
            if x[0] >= xm {
                p.clone()
            } else {
                Vec::new()
            }
        }));
    }
    let mut y = Yee::new(cfg).unwrap();
    let (fc, tau, t0) = (0.8, 0.9, 3.0);
    y.add_source(CurrentSource {
        cell: [(xs / h) as usize, 0, 0],
        comp: Comp::Ey,
        waveform: Box::new(move |t| {
            (2.0 * PI * fc * (t - t0)).sin() * (-((t - t0) / tau).powi(2)).exp()
        }),
    });
    let ip = (xp / h) as usize;
    let steps = (16.0 / y.dt) as usize;
    *dt_probe = y.dt;
    let mut out = Vec::with_capacity(steps);
    for _ in 0..steps {
        y.step();
        out.push(y.field(Comp::Ey)[ip]);
    }
    out
}

fn dft(
    series: &[f64],
    dt: f64,
    f: f64,
) -> Complex64 {
    series
        .iter()
        .enumerate()
        .map(|(i, &v)| Complex64::from_polar(v * dt, -2.0 * PI * f * (i + 1) as f64 * dt))
        .sum()
}

fn reflection_check(
    poles: Vec<Pole>,
    freqs: &[f64],
    eps_inf: f64,
) {
    let (mut dt, mut dt2) = (0.0, 0.0);
    let vac = probe_series(None, &mut dt);
    // the ADE medium uses eps_inf = 1 here
    assert!((eps_inf - 1.0).abs() < 1e-12);
    let mat = probe_series(Some(poles.clone()), &mut dt2);
    assert!((dt - dt2).abs() < 1e-15);
    let refl: Vec<f64> = mat.iter().zip(&vac).map(|(a, b)| a - b).collect();
    for &f in freqs {
        let omega = 2.0 * PI * f;
        let eps = Complex64::new(eps_inf, 0.0)
            + poles
                .iter()
                .map(|p| Complex64::new(p.susceptibility(omega).0, p.susceptibility(omega).1))
                .sum::<Complex64>();
        let n = eps.sqrt();
        let r_exact = ((Complex64::new(1.0, 0.0) - n) / (Complex64::new(1.0, 0.0) + n)).norm();
        let r_num = (dft(&refl, dt, f) / dft(&vac, dt, f)).norm();
        assert!(
            (r_num - r_exact).abs() < 0.02 + 0.03 * r_exact,
            "f = {f}: |R| = {r_num}, Fresnel {r_exact} (eps = {eps})"
        );
    }
}

#[test]
fn lorentz_medium_reflection_matches_the_fresnel_coefficient() {
    // resonance at f0 = 1 with damping 0.6 rad/s: strong dispersion across the pulse band
    reflection_check(
        vec![Pole::lorentz(1.5, 2.0 * PI * 1.0, 0.6)],
        &[0.35, 0.6, 0.85, 1.15, 1.3],
        1.0,
    );
}

#[test]
fn drude_medium_reflection_matches_the_fresnel_coefficient() {
    // plasma frequency 0.7, collision rate 0.4 rad/s: metal-like below and transparent above
    reflection_check(
        vec![Pole::drude(2.0 * PI * 0.7, 0.4)],
        &[0.3, 0.5, 0.75, 1.0, 1.3],
        1.0,
    );
}

#[test]
fn two_pole_medium_matches_its_susceptibility() {
    let poles = vec![
        Pole::lorentz(0.8, 2.0 * PI * 0.6, 0.5),
        Pole::lorentz(1.0, 2.0 * PI * 1.4, 0.8),
    ];
    reflection_check(poles, &[0.35, 0.5, 0.8, 1.0], 1.0);
}

// ------------------------------------------------------------------ CFS-PML

fn guide_energy(
    kappa: f64,
    alpha: f64,
) -> f64 {
    // parallel-plate waveguide (PEC walls y = 0, 1), TM_1 mode with cut-off f = 0.5, PML at both ends
    let (nx, ny) = (120, 20);
    let h = 1.0 / ny as f64;
    let mut cfg = YeeConfig::new(&[nx, ny], &[h, h]).unwrap();
    cfg.bcs = FaceBcs::dirichlet_zero()
        .with(rssn::sim::boundary::Face::XMin, Bc::Absorbing)
        .with(rssn::sim::boundary::Face::XMax, Bc::Absorbing);
    cfg.pml_cells = 12;
    cfg.cfl = 0.9;
    cfg.pml_kappa_max = kappa;
    cfg.pml_alpha_max = alpha;
    let mut y = Yee::new(cfg).unwrap();
    let (fc, tau, t0) = (0.5, 3.0, 8.0);
    for j in 1..ny {
        let amp = (PI * j as f64 * h).sin();
        y.add_source(CurrentSource {
            cell: [nx / 2, j, 0],
            comp: Comp::Ez,
            waveform: Box::new(move |t| {
                amp * (2.0 * PI * fc * (t - t0)).sin() * (-((t - t0) / tau).powi(2)).exp()
            }),
        });
    }
    let steps = (120.0 / y.dt) as usize;
    for _ in 0..steps {
        y.step();
    }
    y.energy()
}

#[test]
fn cfs_pml_absorbs_near_cutoff_waves_better_than_the_plain_pml() {
    let plain = guide_energy(1.0, 0.0);
    let cfs = guide_energy(5.0, 1.0);
    println!("residual guide energy: plain {plain:e}, cfs {cfs:e}");
    assert!(cfs.is_finite() && plain.is_finite());
    assert!(cfs < 0.6 * plain, "CFS {cfs:e} vs plain {plain:e}");
}

// ------------------------------------------------------------ Dey-Mittra PEC

/// Frequency of the lowest TM mode of a circular PEC cavity of radius 1 from the ring-down of `Ez`.
fn cavity_frequency(
    conformal: bool,
    cells: usize,
) -> f64 {
    let r_cav = 1.0;
    let half = 1.25;
    let h = 2.0 * half / cells as f64;
    let mut cfg = YeeConfig::new(&[cells, cells], &[h, h]).unwrap();
    cfg.cfl = 0.9;
    let c = half;
    if conformal {
        cfg.pec_body = Some(Box::new(move |x| {
            r_cav - ((x[0] - c).powi(2) + (x[1] - c).powi(2)).sqrt()
        }));
    }
    let mut y = Yee::new(cfg).unwrap();
    // initial Ez pulse (off-centre but symmetric modes dominate near the centre)
    let n = y.n();
    for j in 0..n[1] {
        for i in 0..n[0] {
            let (x, yy) = (i as f64 * h - c, j as f64 * h - c);
            let r = x.hypot(yy);
            if r < r_cav - 0.5 * h || conformal {
                let idx = y.idx(i, j, 0);
                let inside = r < r_cav;
                if inside {
                    y.field_mut(Comp::Ez)[idx] = (-(r / 0.3).powi(2)).exp();
                }
            }
        }
    }
    let nsteps = (400.0 / y.dt) as usize;
    let mut series = Vec::with_capacity(nsteps);
    let ip = y.idx((c / h) as usize + (0.4 / h) as usize, (c / h) as usize, 0);
    let mask: Vec<usize> = if conformal {
        Vec::new()
    } else {
        (0..y.size())
            .filter(|&p| {
                let (i, j) = (p % n[0], p / n[0]);
                (i as f64 * h - c).hypot(j as f64 * h - c) >= r_cav
            })
            .collect()
    };
    for _ in 0..nsteps {
        y.step();
        if !conformal {
            let f = y.field_mut(Comp::Ez);
            for &p in &mask {
                f[p] = 0.0;
            }
        }
        series.push(y.field(Comp::Ez)[ip]);
    }
    peak_frequency(&series, y.dt, 0.30, 0.50)
}

/// Frequency of the spectral maximum in `[f_lo, f_hi]` (Hann window, zero padding, parabolic interpolation).
fn peak_frequency(
    series: &[f64],
    dt: f64,
    f_lo: f64,
    f_hi: f64,
) -> f64 {
    let n = series.len();
    let mean = series.iter().sum::<f64>() / n as f64;
    let win: Vec<f64> = (0..n)
        .map(|i| (series[i] - mean) * 0.5 * (1.0 - (2.0 * PI * i as f64 / n as f64).cos()))
        .collect();
    let df = 1.0 / (n as f64 * dt) / 8.0;
    let spectrum = |f: f64| dft(&win, dt, f).norm();
    let mut best = (0.0, f_lo);
    let mut f = f_lo;
    while f <= f_hi {
        let a = spectrum(f);
        if a > best.0 {
            best = (a, f);
        }
        f += df;
    }
    let (fm, f0, fp) = (spectrum(best.1 - df), best.0, spectrum(best.1 + df));
    best.1 + 0.5 * df * (fm - fp) / (fm - 2.0 * f0 + fp)
}

#[test]
fn dey_mittra_cavity_resonance_converges_faster_than_the_staircase() {
    let exact = 2.404_825_557_695_773 / (2.0 * PI);
    let mut dm_err = Vec::new();
    let mut st_err = Vec::new();
    for cells in [40, 80] {
        dm_err.push(((cavity_frequency(true, cells) - exact) / exact).abs());
        st_err.push(((cavity_frequency(false, cells) - exact) / exact).abs());
    }
    println!("DM errors {dm_err:?}, staircase errors {st_err:?}");
    assert!(dm_err[1] < 1e-3, "Dey-Mittra error {}", dm_err[1]);
    assert!(
        dm_err[1] < 0.1 * st_err[1] && dm_err[0] < 0.1 * st_err[0],
        "DM {dm_err:?} vs staircase {st_err:?}"
    );
    // the staircase converges at best linearly; the conformal scheme improves with refinement too
    assert!(dm_err[1] < dm_err[0], "DM {dm_err:?} does not converge");
}

// ----------------------------------------------------------- rotated TF/SF

fn rotated_run(angle: f64) -> (f64, f64) {
    let (n, pml, width) = (120usize, 12usize, 0.15);
    let h = 2.0 / n as f64;
    let mut cfg = YeeConfig::new(&[n, n], &[h, h]).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Absorbing);
    cfg.pml_cells = pml;
    cfg.cfl = 0.9;
    let mut y = Yee::new(cfg).unwrap();
    let (lo, hi) = (0.6, 1.4);
    let t0 = 4.0 * width;
    let origin = [0.0; 3];
    let wave = ObliquePlaneWave {
        lo: [lo, lo, 0.0],
        hi: [hi, hi, 0.0],
        direction: [1.0, 0.0, 0.0],
        polarization: [0.0, 0.0, 1.0],
        origin,
        waveform: Box::new(move |tau| (-((tau - t0) / width).powi(2)).exp()),
    };
    let (s, c) = angle.sin_cos();
    y.set_oblique_plane_wave_rotated(wave, [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        .unwrap();
    let steps = (3.2 / y.dt) as usize;
    let (ic, jc) = (n / 2, n / 2);
    let probe = y.idx(ic, jc, 0);
    let (mut inside_err, mut leak): (f64, f64) = (0.0, 0.0);
    let nn = y.n();
    let centre = 1.0;
    for _ in 0..steps {
        y.step();
        let t = y.t;
        let pos = [ic as f64 * h, jc as f64 * h];
        let ex = (-((t - pos[0] - t0) / width).powi(2)).exp();
        inside_err = inside_err.max((y.field(Comp::Ez)[probe] - ex).abs());
        for p in 0..y.size() {
            let (i, j) = (p % nn[0], p / nn[0]);
            if i < pml + 2 || i + pml + 2 >= nn[0] || j < pml + 2 || j + pml + 2 >= nn[1] {
                continue;
            }
            let (dx, dy) = (i as f64 * h - centre, j as f64 * h - centre);
            // local coordinates in the frame of the rotated box
            let lx = c * dx + s * dy;
            let ly = -s * dx + c * dy;
            let half = 0.5 * (hi - lo);
            if lx.abs() > half + 0.1 || ly.abs() > half + 0.1 {
                leak = leak.max(y.field(Comp::Ez)[p].abs());
            }
        }
    }
    (inside_err, leak)
}

#[test]
fn tfsf_on_a_rotated_box_confines_the_incident_wave() {
    let (in0, leak0) = rotated_run(0.0);
    let (in1, leak1) = rotated_run(25.0_f64.to_radians());
    let (in2, leak2) = rotated_run(45.0_f64.to_radians());
    println!("inside errors {in0} {in1} {in2}, leaks {leak0} {leak1} {leak2}");
    assert!(in0 < 0.06 && leak0 < 0.06);
    assert!(
        in1 < 0.08 && leak1 < 0.08,
        "25 deg: inside {in1}, leak {leak1}"
    );
    assert!(
        in2 < 0.08 && leak2 < 0.08,
        "45 deg: inside {in2}, leak {leak2}"
    );
    // invalid rotations
    let mut y = Yee::new(YeeConfig::new(&[16, 16], &[0.1, 0.1]).unwrap()).unwrap();
    let mk = || ObliquePlaneWave {
        lo: [0.3, 0.3, 0.0],
        hi: [1.2, 1.2, 0.0],
        direction: [1.0, 0.0, 0.0],
        polarization: [0.0, 0.0, 1.0],
        origin: [0.0; 3],
        waveform: Box::new(|t| t),
    };
    assert!(
        y.set_oblique_plane_wave_rotated(mk(), [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .is_err()
    );
    assert!(
        y.set_oblique_plane_wave_rotated(mk(), [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
            .is_err()
    );
}
