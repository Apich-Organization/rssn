//! Electrostatic PIC: shape-function order, self-force, momentum, charge and cost.

use std::f64::consts::PI;
use std::time::Instant;

use rssn::kernels::fftconv::fft;
use rssn::sim::pic::*;

const MODES: i32 = 64;
const SHAPES: [Shape; 2] = [Shape::Linear, Shape::Cubic];

/// Periodic Gaussian density with unit integral.
fn gaussian_density(
    x: f64,
    sigma: f64,
    x0: f64,
    len: f64,
) -> f64 {
    let mut r = 1.0 / len;
    for m in 1..=MODES {
        let k = 2.0 * PI * m as f64 / len;
        r += 2.0 / len * (-0.5 * k * k * sigma * sigma).exp() * (k * (x - x0)).cos();
    }
    r
}

/// Exact periodic field of [`gaussian_density`].
fn gaussian_field(
    x: f64,
    sigma: f64,
    x0: f64,
    len: f64,
) -> f64 {
    let mut e = 0.0;
    for m in 1..=MODES {
        let k = 2.0 * PI * m as f64 / len;
        let amp = (-0.5 * k * k * sigma * sigma).exp() / len;
        e += 2.0 * amp / k * (k * (x - x0)).sin();
    }
    e
}

/// Particles that sample the periodic Gaussian by midpoint quadrature.
fn gaussian_particles(
    grid: &Grid,
    sigma: f64,
    x0: f64,
) -> (Vec<f64>, Vec<f64>) {
    let np = 1usize << 17;
    let dx = grid.length() / np as f64;
    let mut pos = Vec::with_capacity(np);
    let mut q = Vec::with_capacity(np);
    for j in 0..np {
        let x = (j as f64 + 0.5) * dx;
        pos.push(x);
        q.push(gaussian_density(x, sigma, x0, grid.length()) * dx);
    }
    (pos, q)
}

/// Relative `L2` error of the force seen by off-grid probe particles.
///
/// The mesh field is solved from the deposited Gaussian and then gathered at
/// points that do not coincide with cell centres; the residual is the
/// interpolation error of the shape, which is what a particle feels.
fn field_error(
    grid: &Grid,
    shape: Shape,
    sigma: f64,
    x0: f64,
) -> f64 {
    let (pos, q) = gaussian_particles(grid, sigma, x0);
    let rho = grid.deposit(&pos, &q, shape);
    let field = grid.solve(&rho, shape).expect("solve");
    let probes = 4096;
    let len = grid.length();
    let mut x = Vec::with_capacity(probes);
    for j in 0..probes {
        x.push((j as f64 + 0.37) * len / probes as f64);
    }
    let e = grid.gather(&field.e, &x, shape).expect("gather");
    let mut num = 0.0;
    let mut den = 0.0;
    for (j, &ej) in e.iter().enumerate() {
        let exact = gaussian_field(x[j], sigma, x0, len);
        let d = ej - exact;
        num += d * d;
        den += exact * exact;
    }
    (num / den).sqrt()
}

/// Least-squares slope of `ln(err)` against `ln(n)`.
fn slope(
    ns: &[usize],
    errs: &[f64],
) -> f64 {
    let x: Vec<f64> = ns.iter().map(|&n| (n as f64).ln()).collect();
    let y: Vec<f64> = errs.iter().map(|&e| e.ln()).collect();
    let mx = x.iter().sum::<f64>() / x.len() as f64;
    let my = y.iter().sum::<f64>() / y.len() as f64;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    for i in 0..x.len() {
        sxy += (x[i] - mx) * (y[i] - my);
        sxx += (x[i] - mx) * (x[i] - mx);
    }
    sxy / sxx
}

/// Spectral `d/dx` of a mesh quantity.
fn divergence(
    j: &[f64],
    len: f64,
) -> Vec<f64> {
    let n = j.len();
    let samples: Vec<(f64, f64)> = j.iter().map(|&v| (v, 0.0)).collect();
    let jk = fft(&samples, false);
    let mut dk = vec![(0.0, 0.0); n];
    for (m, d) in dk.iter_mut().enumerate() {
        let mf = if m <= n / 2 {
            m as f64
        } else {
            m as f64 - n as f64
        };
        let k = 2.0 * PI * mf / len;
        // i k (a + i b) = -k b + i k a
        *d = (-k * jk[m].1, k * jk[m].0);
    }
    fft(&dk, true).iter().map(|c| c.0).collect()
}

#[test]
fn pic_linear_deposit_is_the_cic_reference_bit_for_bit() {
    let grid = Grid::new(16, 1.0).unwrap();
    let dx = grid.spacing();
    let pos = [0.0, 0.0625, 0.265_625, 0.5, 0.999];
    let q = [1.0, -0.5, 0.25, 2.0, -1.5];
    let got = grid.deposit(&pos, &q, Shape::Linear);

    // Independent reference written from the CIC formula.
    let mut want = vec![0.0; 16];
    for (&x, &qi) in pos.iter().zip(q.iter()) {
        let u = x / dx;
        let fl = u.floor();
        let f = u - fl;
        let i0 = (fl as i64).rem_euclid(16) as usize;
        let i1 = (i0 + 1) % 16;
        want[i0] += qi * (1.0 - f);
        want[i1] += qi * f;
    }
    for w in &mut want {
        *w /= dx;
    }
    for (a, b) in got.iter().zip(want.iter()) {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "linear deposit is not bit-identical"
        );
    }

    // Literal guard: one particle at u = 4.25 has weights 0.75 / 0.25.
    let one = grid.deposit(&[0.265_625], &[1.0], Shape::Linear);
    assert_eq!(one[4], 12.0);
    assert_eq!(one[5], 4.0);

    let gathered = grid.gather(&got, &pos, Shape::Linear).unwrap();
    let mut want_gather = vec![0.0; pos.len()];
    for (g, &x) in want_gather.iter_mut().zip(pos.iter()) {
        let u = x / dx;
        let fl = u.floor();
        let f = u - fl;
        let i0 = (fl as i64).rem_euclid(16) as usize;
        let i1 = (i0 + 1) % 16;
        *g = (1.0 - f) * got[i0] + f * got[i1];
    }
    for (a, b) in gathered.iter().zip(want_gather.iter()) {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "linear gather is not bit-identical"
        );
    }
}

#[test]
fn pic_shapes_reproduce_the_zeroth_and_first_moments() {
    let grid = Grid::new(64, 1.0).unwrap();
    for shape in SHAPES {
        for step in 0..20 {
            let x = 0.31 + step as f64 * 0.031;
            let rho = grid.deposit(&[x], &[1.0], shape);
            let total: f64 = rho.iter().sum::<f64>() * grid.spacing();
            assert!(
                (total - 1.0).abs() < 1.0e-14,
                "{shape:?} charge {total} at {x}"
            );
            let mut moment = 0.0;
            for (i, &r) in rho.iter().enumerate() {
                moment += i as f64 * grid.spacing() * r * grid.spacing();
            }
            assert!(
                (moment - x).abs() < 1.0e-14,
                "{shape:?} first moment {moment} vs {x}"
            );
        }
    }
}

#[test]
fn pic_gaussian_field_converges_at_the_shape_order() {
    let ns = [32usize, 64, 128, 256];
    let sigma = 0.1;
    let x0 = 0.5;
    for shape in SHAPES {
        let errs: Vec<f64> = ns
            .iter()
            .map(|&n| field_error(&Grid::new(n, 1.0).unwrap(), shape, sigma, x0))
            .collect();
        let s = -slope(&ns, &errs);
        println!("{shape:?}: errors {errs:?} slope {s:.3}");
        assert!(
            (s - shape.order() as f64).abs() < 0.35,
            "{shape:?} slope {s} is not order {}",
            shape.order()
        );
        assert!(
            errs.windows(2).all(|w| w[1] < w[0]),
            "{shape:?} not decreasing"
        );
    }
}

#[test]
fn pic_single_particle_has_no_self_force() {
    let positions = [0.1, 0.234_567, 0.5, 0.876_543];
    let mut worst = 0.0_f64;
    for shape in SHAPES {
        for n in [16usize, 32, 64, 128] {
            let grid = Grid::new(n, 1.0).unwrap();
            for &x in &positions {
                let rho = grid.deposit(&[x], &[1.0], shape);
                let field = grid.solve(&rho, shape).unwrap();
                let e = grid.gather(&field.e, &[x], shape).unwrap();
                worst = worst.max(e[0].abs());
                assert!(e[0].abs() < 1.0e-11, "{shape:?} n={n} x={x} E={}", e[0]);
            }
        }
    }
    println!("self-force: max |E| = {worst:.3e}");

    // Sensitivity control: deposit and gather with different shapes. The
    // mismatch plays the role of an inconsistent shape derivative and gives a
    // force that is many orders larger than the same-shape roundoff. It does
    // not shrink with the mesh: a kernel mismatch at the particle is a fixed
    // fraction of the self-field, which is exactly why the same-shape test is
    // the sharp one.
    let mut mismatch = Vec::new();
    for n in [16usize, 32, 64, 128] {
        let grid = Grid::new(n, 1.0).unwrap();
        let x = 0.25 + 0.3 / n as f64;
        let rho = grid.deposit(&[x], &[1.0], Shape::Linear);
        let field = grid.solve(&rho, Shape::Linear).unwrap();
        let e = grid.gather(&field.e, &[x], Shape::Cubic).unwrap();
        mismatch.push(e[0].abs());
    }
    println!("mismatched-shape self-force: {mismatch:?}");
    assert!(
        mismatch[0] > 1.0e-8,
        "control is not sensitive: {mismatch:?}"
    );
    assert!(
        mismatch.iter().all(|&m| m > 1.0e-3),
        "control collapsed to roundoff: {mismatch:?}"
    );
}

#[test]
fn pic_drifting_plasma_conserves_momentum() {
    let np = 512;
    let n = 64;
    let v0 = 0.3;
    for shape in SHAPES {
        let cfg = PicConfig::new(n, 1.0, shape).unwrap();
        let mut pos = Vec::with_capacity(np);
        for j in 0..np {
            let base = (j as f64 + 0.5) / np as f64;
            pos.push(base + 0.4 / np as f64 * (7.0 * j as f64).sin());
        }
        let vel = vec![v0; np];
        let charge = vec![1.0; np];
        let mass = vec![1.0; np];
        let mut pic = Pic::new(&cfg, pos, vel, charge, mass).unwrap();
        let p0 = pic.total_momentum();
        let v_start = pic.velocities().to_vec();
        for _ in 0..100 {
            pic.step(1.0e-3).unwrap();
        }
        let p1 = pic.total_momentum();
        let relative = (p1 - p0).abs() / p0.abs();
        let kick = pic
            .velocities()
            .iter()
            .zip(v_start.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        println!("{shape:?}: momentum drift {relative:.3e}, max |dv| {kick:.3e}");
        assert!(relative < 1.0e-12, "{shape:?} momentum drift {relative}");
        assert!(
            kick > 1.0e-6,
            "{shape:?} test is trivial, max |dv| = {kick}"
        );
    }
}

#[test]
fn pic_charge_is_conserved() {
    let grid = Grid::new(64, 1.0).unwrap();
    let mut state = 12345u64;
    let mut pos = Vec::new();
    let mut q = Vec::new();
    for j in 0..2000 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        pos.push((state >> 11) as f64 / (1u64 << 53) as f64);
        q.push(if j % 2 == 0 { 1.0 } else { -0.5 });
    }
    let total: f64 = q.iter().sum();
    for shape in SHAPES {
        let rho = grid.deposit(&pos, &q, shape);
        let mesh: f64 = rho.iter().sum::<f64>() * grid.spacing();
        let relative = (mesh - total).abs() / total.abs();
        println!("{shape:?}: charge residual {relative:.3e}");
        assert!(relative < 1.0e-13, "{shape:?} charge residual {relative}");
    }
}

#[test]
fn pic_continuity_residual_is_small() {
    let n = 64;
    let grid = Grid::new(n, 1.0).unwrap();
    let (pos, q) = gaussian_particles(&grid, 0.15, 0.5);
    let v = 0.3;
    let dt = 1.0e-4;
    let len = grid.length();
    let back: Vec<f64> = pos.iter().map(|&x| (x - v * dt).rem_euclid(len)).collect();
    let fwd: Vec<f64> = pos.iter().map(|&x| (x + v * dt).rem_euclid(len)).collect();
    let jq: Vec<f64> = q.iter().map(|&qi| qi * v).collect();
    for shape in SHAPES {
        let rho0 = grid.deposit(&back, &q, shape);
        let rho1 = grid.deposit(&fwd, &q, shape);
        let j = grid.deposit(&pos, &jq, shape);
        let dj = divergence(&j, len);
        let mut num = 0.0;
        let mut den = 0.0;
        let mut peak = 0.0_f64;
        for i in 0..n {
            let dtrho = (rho1[i] - rho0[i]) / (2.0 * dt);
            let r = dtrho + dj[i];
            num += r * r;
            den += dtrho * dtrho;
            peak = peak.max(r.abs());
        }
        let relative = (num / den).sqrt();
        println!("{shape:?}: continuity L2 residual {relative:.3e}, peak {peak:.3e}");
        assert!(
            relative < 1.0e-2,
            "{shape:?} continuity residual {relative}"
        );
    }
}

#[test]
fn pic_cost_grows_with_the_stencil_width() {
    let grid = Grid::new(256, 1.0).unwrap();
    let np = 200_000;
    let mut state = 99u64;
    let mut pos = Vec::with_capacity(np);
    for _ in 0..np {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        pos.push((state >> 11) as f64 / (1u64 << 53) as f64);
    }
    let q = vec![1.0; np];
    let mut base = 0.0;
    for shape in SHAPES {
        let reps = 10;
        let start = Instant::now();
        for _ in 0..reps {
            let rho = grid.deposit(&pos, &q, shape);
            std::hint::black_box(&rho);
        }
        let ns = start.elapsed().as_nanos() as f64 / (reps as f64 * np as f64);
        if shape == Shape::Linear {
            base = ns;
        }
        println!(
            "{shape:?}: {ns:.2} ns/particle, ratio {:.2}, cells {}",
            ns / base,
            shape.cells()
        );
    }
    assert_eq!(Shape::Linear.cells(), 2);
    assert_eq!(Shape::Cubic.cells(), 4);
}
