//! `sim::acoustics`: the staggered-grid scalar wave solver. Direct checks of
//! the time step, the d'Alembert split of a pulse, the discrete conservation
//! laws (mass and energy), the domain of dependence, soft sources, symmetry,
//! wave-speed isotropy in 2D and the input checks; and of the variable medium,
//! the impedance reflection, Snell refraction, total internal reflection, the
//! energy of a heterogeneous box and the second-order convergence in `c(x)`.

use rssn::sim::acoustics::Acoustics;
use rssn::sim::acoustics::AcousticsConfig;
use rssn::sim::boundary::BcValue;
use rssn::sim::boundary::BoundaryCondition as Bc;
use rssn::sim::boundary::Face;
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

#[test]
fn the_dirichlet_datum_sets_the_boundary_gradient() {
    // One velocity update from a known pressure field. The face sits half a
    // cell from the first cell centre, so `v = -dt (p - g)/(h/2)` there; with
    // `g = 0` that is the old pressure-release rule exactly, which is what
    // keeps the pre-existing reflection test passing.
    let n = 16;
    let h = 1.0 / n as f64;
    for g in [0.0, 0.25, -1.5] {
        let mut cfg = AcousticsConfig::new(&[n], &[h], 1.0).unwrap();
        cfg.bcs = FaceBcs::uniform(&Bc::neumann(0.0)).with(Face::XMin, Bc::dirichlet(g));
        let mut a = Acoustics::new(cfg).unwrap();
        a.p.fill(1.0);
        let dt = a.dt;
        a.step();
        let want = -dt * (1.0 - g) / (0.5 * h);
        assert!(
            (a.v[0][0] - want).abs() < 1e-12,
            "g = {g}: {} vs {want}",
            a.v[0][0]
        );
    }
}

#[test]
fn the_dirichlet_datum_enters_linearly() {
    // The scheme is linear and starts from rest, so doubling the datum doubles
    // the field exactly. This is what would break if the datum were folded
    // into the state instead of into the boundary gradient.
    let run = |g: f64| {
        let n = 24;
        let h = 1.0 / n as f64;
        let mut cfg = AcousticsConfig::new(&[n], &[h], 1.0).unwrap();
        cfg.bcs = FaceBcs::uniform(&Bc::neumann(0.0)).with(Face::XMin, Bc::dirichlet(g));
        let mut a = Acoustics::new(cfg).unwrap();
        a.run(300);
        a.p
    };
    let one = run(1.0);
    let two = run(2.0);
    for (single, double) in one.iter().zip(&two) {
        assert!(
            (2.0 * single - double).abs() < 1e-12 * (1.0 + double.abs()),
            "{single} vs {double}"
        );
    }
}

#[test]
fn the_datum_is_evaluated_at_the_face_not_at_the_origin() {
    // In 2D the low `x` face runs along `y`, so a datum that varies with `y`
    // must produce a face velocity that varies with `y`. A datum evaluated at
    // the origin instead would give the same velocity on every row.
    let n = [8, 8];
    let h = 1.0 / 8.0;
    let mut cfg = AcousticsConfig::new(&n, &[h, h], 1.0).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::neumann(0.0))
        .with(Face::XMin, Bc::Dirichlet(BcValue::space(|x| x[1])));
    let mut a = Acoustics::new(cfg).unwrap();
    a.p.fill(1.0);
    let dt = a.dt;
    a.step();
    // `v[0]` is `(n[0] + 1) x n[1]`, x-fastest.
    let stride = n[0] + 1;
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for j in 0..n[1] {
        let y = (j as f64 + 0.5) * h;
        let want = -dt * (1.0 - y) / (0.5 * h);
        let got = a.v[0][j * stride];
        assert!((got - want).abs() < 1e-12, "row {j}: {got} vs {want}");
        lo = lo.min(got);
        hi = hi.max(got);
    }
    assert!(
        hi - lo > 0.5,
        "the face velocity must vary with y: {lo}..{hi}"
    );
}

#[test]
fn a_time_dependent_datum_drives_the_field() {
    // A half-sine pulse at the face; the box must be quiet before it starts,
    // move while it lasts, and not return to rest instantly afterwards.
    let n = 64;
    let h = 1.0 / n as f64;
    let mut cfg = AcousticsConfig::new(&[n], &[h], 1.0).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::neumann(0.0)).with(
        Face::XMin,
        Bc::Dirichlet(BcValue::time(|t| (PI * t).sin().max(0.0))),
    );
    let mut a = Acoustics::new(cfg).unwrap();
    assert_eq!(a.p.iter().fold(0.0_f64, |m, &p| m.max(p.abs())), 0.0);
    a.run(200);
    let during = a.p.iter().fold(0.0_f64, |m, &p| m.max(p.abs()));
    assert!(during > 0.2, "the pulse must reach the interior: {during}");
}

// ------------------------------------------------------- variable medium

/// A wave packet `cos(k d) e^{-(d² + s²)/2σ²}` travelling along `theta`, with
/// `d` the distance along the beam from `x0` and `s` the distance across it.
fn wave_packet(
    a: &mut Acoustics,
    x0: [f64; 2],
    theta: f64,
    k: f64,
    sigma: f64,
) {
    let (sn, cs) = theta.sin_cos();
    a.set_pressure(move |x| {
        let (dx, dy) = (x[0] - x0[0], x[1] - x0[1]);
        let d = dx * cs + dy * sn;
        let s = dy * cs - dx * sn;
        (k * d).cos() * (-(d * d + s * s) / (2.0 * sigma * sigma)).exp()
    });
    a.start_from_rest();
}

/// A value recorded before the variable medium existed. The recorded bits come
/// from one platform's `sin`/`exp`; the last bits of a field that has run for
/// hundreds of steps depend on the math library, so the comparison is to
/// round-off of the field amplitude, not to the bit. The scheme itself is
/// pinned bit for bit by `the_uniform_medium_is_the_varying_one_with_a_constant_field`.
fn frozen(
    got: f64,
    want_bits: u64,
    what: &str,
) {
    let want = f64::from_bits(want_bits);
    assert!(
        (got - want).abs() <= 1e-11 * (1.0 + want.abs()),
        "{what}: {got:e} vs {want:e}"
    );
}

#[test]
fn the_constant_medium_keeps_the_frozen_field_bits() {
    // The conservative form has to reduce to the old `p_t = -c²∇·v`,
    // `v_t = -∇p` scheme term for term. These are the exact bit patterns the
    // solver produced before the variable medium existed, so any change in the
    // arithmetic order moves at least one of them.
    let n = 64;
    let h = 1.0 / n as f64;
    let mut cfg = AcousticsConfig::new(&[n], &[h], 1.5).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::neumann(0.0)).with(Face::XMax, Bc::dirichlet(0.0));
    let mut a = Acoustics::new(cfg).unwrap();
    a.add_source(20, |t| (2.0 * PI * 3.0 * t).sin());
    a.run(300);
    assert_eq!(a.dt.to_bits(), 4_576_557_941_333_898_035, "1D dt");
    assert_eq!(a.t.to_bits(), 4_613_515_605_776_007_133, "1D t");
    for (i, want) in [
        (0usize, 0xbf16_b34a_9aa4_a0ea_u64),
        (1, 0x3ef1_4891_e9c2_6ed3),
        (20, 0x3f5f_d966_e7f6_17c0),
        (31, 0x3f64_61d0_cb03_6c08),
        (63, 0xbf03_01fc_f217_a5f1),
    ] {
        frozen(a.p[i], want, &format!("1D cell {i}"));
    }
    frozen(a.energy(), 0x3ec9_2f55_35e7_48eb, "1D energy");

    let n = [24, 20];
    let h = [1.0 / 24.0, 1.0 / 20.0];
    let mut cfg = AcousticsConfig::new(&n, &h, 2.0).unwrap();
    cfg.bcs = FaceBcs::uniform(&Bc::Periodic);
    let mut b = Acoustics::new(cfg).unwrap();
    b.set_pressure(|x| {
        (-((x[0] - 0.4).powi(2) + (x[1] - 0.6).powi(2)) / (2.0 * 0.05 * 0.05)).exp()
    });
    b.start_from_rest();
    b.run(150);
    assert_eq!(b.dt.to_bits(), 4_579_457_048_345_532_227, "2D dt");
    for (i, j, want) in [
        (0usize, 0usize, 0x3f9d_e532_c62d_aa26_u64),
        (5, 7, 0x3fa3_cbec_c52e_0c41),
        (23, 19, 0x3fbc_7b66_5e4d_33be),
    ] {
        frozen(b.p[b.idx(i, j, 0)], want, &format!("2D cell ({i}, {j})"));
    }
    frozen(b.energy(), 0x3f4e_9de5_19e2_a014, "2D energy");
}

#[test]
fn the_uniform_medium_is_the_varying_one_with_a_constant_field() {
    // A constant `c` and `ρ` must give the same bits, not merely the same
    // numbers: the conservative form has to collapse term for term onto the
    // old `p_t = -c²∇·v`, `v_t = -∇p` scheme.
    for (c, rho) in [(1.5, 1.0), (0.75, 2.0)] {
        let n = 128;
        let h = 1.0 / n as f64;
        let mut cfg = AcousticsConfig::new(&[n], &[h], c).unwrap();
        cfg.rho = rho;
        let mut plain = Acoustics::new(cfg).unwrap();
        let mut varied = Acoustics::new(AcousticsConfig::new(&[n], &[h], c).unwrap()).unwrap();
        varied.set_medium(move |_| (c, rho)).unwrap();
        assert_eq!(plain.dt.to_bits(), varied.dt.to_bits(), "c = {c}");
        for s in [&mut plain, &mut varied] {
            s.add_source(40, |t| (2.0 * PI * 4.0 * t).sin());
            s.run(400);
        }
        for i in 0..n {
            assert_eq!(
                plain.p[i].to_bits(),
                varied.p[i].to_bits(),
                "c = {c}, cell {i}"
            );
        }
        assert_eq!(plain.energy().to_bits(), varied.energy().to_bits());
    }
}

/// Peak of the bump in `0.8 < x < 1` of a run whose interface sits on the face
/// at `x = 1`, signed by the extremum it belongs to.
fn reflected_peak(
    a: &Acoustics,
    n: usize,
    h: f64,
) -> f64 {
    let mut best: f64 = 0.0;
    for i in 0..n {
        let x = (i as f64 + 0.5) * h;
        if (0.8..1.0).contains(&x) && a.p[i].abs() > best.abs() {
            best = a.p[i];
        }
    }
    best
}

/// Ratio of the reflected to the incident bump for a normal-incidence crossing.
fn normal_reflection(
    c2: f64,
    rho2: f64,
) -> f64 {
    let n = 1600;
    let h = 2.0 / n as f64;
    let mut a = Acoustics::new(AcousticsConfig::new(&[n], &[h], 1.0).unwrap()).unwrap();
    a.set_medium(move |x| {
        if x[0] < 1.0 {
            (1.0, 1.0)
        } else {
            (c2, rho2)
        }
    })
    .unwrap();
    a.set_pressure(|x| (-(x[0] - 0.8).powi(2) / (2.0 * 0.02 * 0.02)).exp());
    a.start_from_rest();
    // The bump is centred at 0.9 at t = 0.1 and, mirrored about the interface,
    // at 0.9 again at t = 0.3: both peaks are measured 0.1 from the interface,
    // so the numerical dispersion cancels in the ratio.
    a.run((0.1 / a.dt).round() as usize);
    let incident = reflected_peak(&a, n, h);
    a.run((0.2 / a.dt).round() as usize);
    reflected_peak(&a, n, h) / incident
}

#[test]
fn an_interface_reflects_with_the_impedance_ratio() {
    // R = (Z₂ - Z₁)/(Z₁ + Z₂) with Z = ρc, for a fast, a heavy and an
    // impedance-matched medium. The last one has R = 0 although both c and ρ
    // jump, which only a conservative interface treatment gets right.
    for (c2, rho2) in [(0.5, 1.0), (1.0, 4.0), (2.0, 0.5)] {
        let z2 = c2 * rho2;
        let want = (z2 - 1.0) / (z2 + 1.0);
        let got = normal_reflection(c2, rho2);
        assert!(
            (got - want).abs() < 0.02,
            "c = {c2}, rho = {rho2}: R = {got}, want {want}"
        );
    }
}

/// Direction of the transmitted beam: the energy centroid of the second
/// medium at two times gives its group velocity.
fn refraction_angle(
    c2: f64,
    theta1: f64,
) -> f64 {
    // h is small enough that the discrete group velocity points within a
    // quarter of a degree of the wavevector (the error is O((kh)²)).
    let n = 512;
    let h = 1.0 / n as f64;
    let mut a = Acoustics::new(AcousticsConfig::new(&[n, n], &[h, h], 1.0).unwrap()).unwrap();
    a.set_medium(move |x| {
        if x[0] < 0.5 {
            (1.0, 1.0)
        } else {
            (c2, 1.0)
        }
    })
    .unwrap();
    // The beam enters the interface at (0.5, 0.45), 0.2 back along its axis.
    wave_packet(
        &mut a,
        [0.3, 0.45 - 0.2 * theta1.tan()],
        theta1,
        100.0,
        0.05,
    );
    let centroid = |a: &Acoustics| -> [f64; 2] {
        let (mut w, mut sx, mut sy) = (0.0, 0.0, 0.0);
        for j in 0..n {
            for i in (n / 2)..n {
                let p = a.p[a.idx(i, j, 0)];
                let q = p * p;
                w += q;
                sx += q * (i as f64 + 0.5) * h;
                sy += q * (j as f64 + 0.5) * h;
            }
        }
        [sx / w, sy / w]
    };
    // Both samples are clear of the interface and of the walls.
    a.run((0.55 / a.dt).round() as usize);
    let before = centroid(&a);
    a.run((0.15 / a.dt).round() as usize);
    let after = centroid(&a);
    assert!(after[0] > before[0], "the beam must travel into medium 2");
    (after[1] - before[1]).atan2(after[0] - before[0])
}

#[test]
fn an_oblique_beam_refracts_by_snell_law() {
    for (c2, theta1) in [(0.5, 30.0_f64), (0.75, 45.0)] {
        let theta1 = theta1.to_radians();
        let got = refraction_angle(c2, theta1);
        let want = (c2 * theta1.sin()).asin();
        assert!(
            (got - want).abs() < 0.03,
            "c2 = {c2}: theta2 = {} deg, want {} deg",
            got.to_degrees(),
            want.to_degrees()
        );
        assert!(
            (theta1.sin() - got.sin() / c2).abs() < 0.03,
            "sin(theta1)/c1 = {} vs sin(theta2)/c2 = {}",
            theta1.sin(),
            got.sin() / c2
        );
    }
}

/// Peak pressure 0.1 deep in the fast medium, for a beam at `theta1` from a
/// slow medium into a fast one.
fn depth_amplitude(
    x0: [f64; 2],
    theta1: f64,
) -> f64 {
    let n = 256;
    let h = 1.0 / n as f64;
    let mut a = Acoustics::new(AcousticsConfig::new(&[n, n], &[h, h], 1.0).unwrap()).unwrap();
    a.set_medium(|x| {
        if x[0] < 0.5 {
            (1.0, 1.0)
        } else {
            (2.0, 1.0)
        }
    })
    .unwrap();
    wave_packet(&mut a, x0, theta1, 160.0, 0.04);
    a.run((0.45 / a.dt).round() as usize);
    let mut best: f64 = 0.0;
    for j in 0..n {
        for i in 0..n {
            if (i as f64 + 0.5) * h > 0.6 {
                best = best.max(a.p[a.idx(i, j, 0)].abs());
            }
        }
    }
    best
}

#[test]
fn past_the_critical_angle_nothing_propagates_into_the_fast_medium() {
    // c1 = 1 into c2 = 2: the critical angle is 30 degrees.
    let below = depth_amplitude([0.3, 0.25], 20.0_f64.to_radians());
    let beyond = depth_amplitude([0.3, 0.254], 60.0_f64.to_radians());
    assert!(
        below > 0.15,
        "below the critical angle the beam must cross: {below}"
    );
    assert!(
        beyond < 0.01 * below,
        "past the critical angle only the evanescent tail may remain: {beyond}"
    );
}

#[test]
fn a_closed_box_with_a_varying_medium_conserves_the_discrete_energy() {
    // ½ Σ (p²/K + ρv²) ΔV is the invariant; ½ Σ (p²/c² + v²) is not, because
    // the strain energy is p²/K with K = ρc² and the kinetic energy is ρv².
    type Medium = fn([f64; 3]) -> (f64, f64);
    let media: [(&str, Medium); 2] = [
        ("smooth", |x| {
            let c = 1.0 + 0.25 * (2.0 * PI * x[0]).sin() * (2.0 * PI * x[1]).cos();
            let r = 1.0 + 0.2 * (2.0 * PI * (x[0] + x[1])).cos();
            (c, r)
        }),
        ("step", |x| {
            if x[0] < 0.5 {
                (0.7, 1.0)
            } else {
                (1.4, 3.0)
            }
        }),
    ];
    for (name, medium) in media {
        let n = 96;
        let h = 1.0 / n as f64;
        let mut a = Acoustics::new(AcousticsConfig::new(&[n, n], &[h, h], 1.0).unwrap()).unwrap();
        a.set_medium(medium).unwrap();
        a.set_pressure(|x| {
            (-((x[0] - 0.4).powi(2) + (x[1] - 0.6).powi(2)) / (2.0 * 0.06 * 0.06)).exp()
        });
        a.start_from_rest();
        let e0 = a.energy();
        let mut worst: f64 = 0.0;
        for _ in 0..8 {
            a.run(40);
            worst = worst.max((a.energy() / e0 - 1.0).abs());
        }
        assert!(worst < 0.05, "{name}: energy drift {worst}");
    }
}

#[test]
fn a_smooth_medium_converges_at_second_order() {
    // The same physical problem on n, 2n and 4n cells with a common time step,
    // so the three runs differ in space only. The coarse centres fall halfway
    // between two fine cells, and the fine average is itself second order.
    let dt = 0.4 / (1.3 * 256.0);
    let run = |n: usize| -> Vec<f64> {
        let h = 1.0 / n as f64;
        let mut cfg = AcousticsConfig::new(&[n], &[h], 1.3).unwrap();
        cfg.cfl = 0.5;
        let mut a = Acoustics::new(cfg).unwrap();
        a.set_medium(|x| (1.0 + 0.3 * (2.0 * PI * x[0]).sin(), 1.0))
            .unwrap();
        a.dt = dt;
        a.set_pressure(|x| (-(x[0] - 0.5).powi(2) / (2.0 * 0.06 * 0.06)).exp());
        a.start_from_rest();
        a.run((0.25 / dt).round() as usize);
        a.p
    };
    let difference = |coarse: &[f64], fine: &[f64]| -> f64 {
        let mut s = 0.0;
        for (i, &c) in coarse.iter().enumerate() {
            let mid = 0.5 * (fine[2 * i] + fine[2 * i + 1]);
            s += (c - mid).powi(2);
        }
        (s / coarse.len() as f64).sqrt()
    };
    let (p64, p128, p256) = (run(64), run(128), run(256));
    let (e1, e2) = (difference(&p64, &p128), difference(&p128, &p256));
    let order = (e1 / e2).log2();
    assert!(e1 > 0.0, "the medium must change the field");
    assert!(order > 1.85, "measured order {order} ({e1} -> {e2})");
}
