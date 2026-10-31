//! `sim::acoustics`: the staggered-grid scalar wave solver. Direct checks of
//! the time step, the d'Alembert split of a pulse, the discrete conservation
//! laws (mass and energy), the domain of dependence, soft sources, symmetry,
//! wave-speed isotropy in 2D and the input checks.

use rssn::sim::acoustics::Acoustics;
use rssn::sim::acoustics::AcousticsConfig;
use rssn::sim::boundary::BoundaryCondition as Bc;
use rssn::sim::boundary::FaceBcs;
use std::f64::consts::PI;

fn solver(
    n: &[usize],
    c: f64,
    bc: Option<Bc>,
) -> Acoustics {
    let h: Vec<f64> = n.iter().map(|&m| 1.0 / m as f64).collect();
    let mut cfg = AcousticsConfig::new(n, &h, c).unwrap();
    if let Some(bc) = bc {
        cfg.bcs = FaceBcs::uniform(&bc);
    }
    Acoustics::new(cfg).unwrap()
}

#[test]
fn configuration_and_time_step() {
    assert!(AcousticsConfig::new(&[], &[], 1.0).is_err(), "no axes");
    assert!(
        AcousticsConfig::new(&[8, 8, 8, 8], &[0.1; 4], 1.0).is_err(),
        "four axes"
    );
    assert!(
        AcousticsConfig::new(&[3], &[0.1], 1.0).is_err(),
        "fewer than four cells"
    );
    assert!(
        AcousticsConfig::new(&[8, 8], &[0.1], 1.0).is_err(),
        "mismatched spacings"
    );
    assert!(
        AcousticsConfig::new(&[8], &[0.1], 0.0).is_err(),
        "zero wave speed"
    );
    assert!(
        AcousticsConfig::new(&[8], &[0.1], -1.0).is_err(),
        "negative wave speed"
    );
    for cfl in [0.0, -0.1, 1.2] {
        let mut cfg = AcousticsConfig::new(&[16], &[0.1], 1.0).unwrap();
        cfg.cfl = cfl;
        assert!(Acoustics::new(cfg).is_err(), "cfl {cfl}");
    }
    // dt = cfl / (c sqrt(sum 1/h²)) in every dimension.
    for (n, h, c) in [
        (vec![32], vec![0.05], 2.0),
        (vec![16, 16], vec![0.1, 0.05], 1.5),
        (vec![8, 8, 8], vec![0.1, 0.2, 0.25], 0.7),
    ] {
        let cfg = AcousticsConfig::new(&n, &h, c).unwrap();
        let cfl = cfg.cfl;
        let inv: f64 = h.iter().map(|x| 1.0 / (x * x)).sum();
        let a = Acoustics::new(cfg).unwrap();
        assert!(
            (a.dt - cfl / (c * inv.sqrt())).abs() < 1e-15,
            "{} vs {}",
            a.dt,
            cfl / (c * inv.sqrt())
        );
        assert_eq!(a.t, 0.0);
    }
    // A PML thicker than half the grid and Robin faces are rejected.
    let mut thick = AcousticsConfig::new(&[16], &[0.1], 1.0).unwrap();
    thick.bcs = FaceBcs::uniform(&Bc::Absorbing);
    thick.pml_cells = 8;
    assert!(Acoustics::new(thick).is_err());
    // Indexing is x-fastest.
    let a = solver(&[5, 6, 7], 1.0, None);
    assert_eq!(a.n(), [5, 6, 7]);
    assert_eq!(a.idx(0, 0, 0), 0);
    assert_eq!(a.idx(1, 0, 0), 1);
    assert_eq!(a.idx(0, 1, 0), 5);
    assert_eq!(a.idx(0, 0, 1), 30);
    assert_eq!(a.p.len(), 210);
}

#[test]
fn a_pulse_splits_into_two_halves_travelling_at_the_wave_speed() {
    // d'Alembert: p(x, t) = [f(x - ct) + f(x + ct)] / 2 before the walls are
    // reached; here c = 2, so the halves have moved 0.3 in t = 0.15.
    let (n, c, t_end) = (800, 2.0, 0.15);
    let mut a = solver(&[n], c, None);
    let sigma = 0.04;
    let f = |x: f64| (-(x - 0.5).powi(2) / (2.0 * sigma * sigma)).exp();
    a.set_pressure(|x| f(x[0]));
    a.start_from_rest();
    let steps = (t_end / a.dt).round() as usize;
    a.run(steps);
    let t = a.t;
    let h = 1.0 / n as f64;
    let mut worst: f64 = 0.0;
    for i in 0..n {
        let x = (i as f64 + 0.5) * h;
        let exact = 0.5 * (f(x - c * t) + f(x + c * t));
        worst = worst.max((a.p[a.idx(i, 0, 0)] - exact).abs());
    }
    assert!(worst < 5e-3, "d'Alembert error {worst}");
    // The two peaks have amplitude 1/2 at 0.5 ± ct.
    let argmax =
        |range: std::ops::Range<usize>| range.max_by(|&i, &j| a.p[i].total_cmp(&a.p[j])).unwrap();
    let (left, right) = (argmax(0..n / 2), argmax(n / 2..n));
    let position = |i: usize| (i as f64 + 0.5) * h;
    assert!(
        (position(left) - (0.5 - c * t)).abs() < 2.0 * h,
        "left peak at {}",
        position(left)
    );
    assert!(
        (position(right) - (0.5 + c * t)).abs() < 2.0 * h,
        "right peak at {}",
        position(right)
    );
    assert!(
        (a.p[left] - 0.5).abs() < 0.01 && (a.p[right] - 0.5).abs() < 0.01,
        "{} {}",
        a.p[left],
        a.p[right]
    );
}

#[test]
fn rigid_and_periodic_boxes_conserve_mass_and_energy() {
    for (name, bc) in [("rigid", None), ("periodic", Some(Bc::Periodic))] {
        for dim in [1_usize, 2] {
            let n = if dim == 1 { 128 } else { 40 };
            let mut a = solver(&vec![n; dim], 1.0, bc.clone());
            a.set_pressure(move |x| {
                let r2: f64 = (0..dim).map(|k| (x[k] - 0.4).powi(2)).sum();
                (-r2 / (2.0 * 0.07 * 0.07)).exp() + 0.3 * (2.0 * PI * x[0]).cos()
            });
            a.start_from_rest();
            let mass = |a: &Acoustics| a.p.iter().sum::<f64>();
            let (m0, e0) = (mass(&a), a.energy());
            let mut worst_energy: f64 = 0.0;
            for _ in 0..8 {
                a.run(25);
                // Σ p is exactly conserved by the flux form.
                assert!(
                    (mass(&a) - m0).abs() < 1e-9 * m0.abs().max(1.0),
                    "{name} {dim}D mass {} vs {m0}",
                    mass(&a)
                );
                worst_energy = worst_energy.max((a.energy() / e0 - 1.0).abs());
            }
            assert!(
                worst_energy < 0.05,
                "{name} {dim}D energy drift {worst_energy}"
            );
        }
    }
}

#[test]
fn information_travels_at_most_one_cell_per_step() {
    // A point disturbance in the middle: after m steps the pressure is
    // exactly zero beyond m cells (the stencil's domain of dependence).
    let n = 101;
    let mut a = solver(&[n], 1.0, None);
    let mid = n / 2;
    a.p[mid] = 1.0;
    for m in 1..=20 {
        a.step();
        for i in 0..n {
            if i.abs_diff(mid) > m {
                assert_eq!(a.p[i], 0.0, "step {m}, cell {i}");
            }
        }
    }
    assert!(
        a.p[mid + 10] != 0.0 && a.p[mid - 10] != 0.0,
        "the disturbance is felt ten cells away"
    );
    // The same in 2D along both axes and the diagonal (Manhattan reach).
    let n = 41;
    let mut b = solver(&[n, n], 1.0, None);
    let mid = n / 2;
    let centre = b.idx(mid, mid, 0);
    b.p[centre] = 1.0;
    for m in 1..=10 {
        b.step();
        for j in 0..n {
            for i in 0..n {
                if i.abs_diff(mid) + j.abs_diff(mid) > m {
                    assert_eq!(b.p[b.idx(i, j, 0)], 0.0, "step {m}, cell ({i}, {j})");
                }
            }
        }
    }
}

#[test]
fn a_soft_source_injects_dt_times_its_value_and_radiates_symmetrically() {
    let n = 100;
    let mut a = solver(&[n], 1.0, None);
    let cell = 50;
    a.add_source(cell, |t| (2.0 * PI * 5.0 * t).sin());
    a.step();
    // One step from rest: p = dt s(dt/2) in the source cell, nothing elsewhere.
    let want = a.dt * (2.0 * PI * 5.0 * 0.5 * a.dt).sin();
    assert!((a.p[cell] - want).abs() < 1e-15, "{} vs {want}", a.p[cell]);
    assert!(a.p.iter().enumerate().all(|(i, &p)| i == cell || p == 0.0));
    // The cell is not at the centre of the grid (cells 0..99, centre between
    // 49 and 50); mirror the source about the centre for the symmetry check.
    let mut b = solver(&[n], 1.0, None);
    b.add_source(49, |t| (2.0 * PI * 5.0 * t).sin());
    b.add_source(50, |t| (2.0 * PI * 5.0 * t).sin());
    b.run(60);
    for i in 0..n {
        assert!(
            (b.p[i] - b.p[n - 1 - i]).abs() < 1e-13,
            "asymmetry at {i}: {} vs {}",
            b.p[i],
            b.p[n - 1 - i]
        );
    }
    assert!(b.p.iter().any(|&p| p.abs() > 1e-4), "the source radiates");
    // Superposition: a source of twice the amplitude gives twice the field.
    let (mut c1, mut c2) = (solver(&[n], 1.0, None), solver(&[n], 1.0, None));
    c1.add_source(30, |t| t.sin());
    c2.add_source(30, |t| 2.0 * t.sin());
    c1.run(80);
    c2.run(80);
    for i in 0..n {
        assert!((c2.p[i] - 2.0 * c1.p[i]).abs() < 1e-13);
    }
}

#[test]
fn the_wave_speed_is_isotropic_in_two_dimensions() {
    // A Gaussian bump: the leading front travels at c along an axis and along
    // the diagonal alike (second-order scheme: a small anisotropy only).
    let n = 200;
    let mut a = solver(&[n, n], 1.0, None);
    let h = 1.0 / n as f64;
    let sigma = 0.025;
    a.set_pressure(|x| {
        (-((x[0] - 0.5).powi(2) + (x[1] - 0.5).powi(2)) / (2.0 * sigma * sigma)).exp()
    });
    a.start_from_rest();
    let steps = (0.3 / a.dt).round() as usize;
    a.run(steps);
    let t = a.t;
    // Radial position of the maximum along the +x axis and the diagonal.
    let along = |dx: f64, dy: f64| -> f64 {
        let mut best = (f64::NEG_INFINITY, 0.0);
        let mut s = 0.0;
        while s < 0.45 {
            let (x, y) = (0.5 + s * dx, 0.5 + s * dy);
            let (i, j) = ((x / h) as usize, (y / h) as usize);
            let p = a.p[a.idx(i.min(n - 1), j.min(n - 1), 0)];
            if p > best.0 {
                best = (p, s);
            }
            s += h / 2.0;
        }
        best.1
    };
    let axis = along(1.0, 0.0);
    let diagonal = along(
        std::f64::consts::FRAC_1_SQRT_2,
        std::f64::consts::FRAC_1_SQRT_2,
    );
    assert!(
        (axis - diagonal).abs() < 3.0 * h,
        "axis {axis}, diagonal {diagonal}"
    );
    // The front sits at ct, shifted slightly inwards by the 2D tail.
    assert!((axis - t).abs() < 0.04, "front at {axis}, ct = {t}");
}

#[test]
fn pressure_release_walls_conserve_energy_and_invert_the_reflection() {
    let n = 160;
    let mut a = solver(&[n], 1.0, Some(Bc::dirichlet(0.0)));
    a.set_pressure(|x| (-(x[0] - 0.5).powi(2) / (2.0 * 0.04_f64.powi(2))).exp());
    a.start_from_rest();
    let e0 = a.energy();
    // Right-moving half reaches the wall at 1 after t = 0.5 and returns
    // inverted: at t = 1.0 the original centre carries -1/2 (twice -1/4).
    a.run((1.0 / a.dt).round() as usize);
    let centre = a.p[a.idx(n / 2, 0, 0)];
    assert!(
        centre < -0.8,
        "after the double reflection the bump is inverted: {centre}"
    );
    assert!(
        (a.energy() / e0 - 1.0).abs() < 0.05,
        "energy {}",
        a.energy() / e0
    );
}
