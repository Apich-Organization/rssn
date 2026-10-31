//! Spectral method advection-diffusion solvers (ported from `physics_sm_test.rs`).

use assert_approx_eq::assert_approx_eq;
use num_complex::Complex;
use proptest::prelude::*;
use proptest::test_runner::RngSeed;
use rssn::sim::physics_sm::*;
use std::f64::consts::PI;

fn cfg() -> ProptestConfig {
    ProptestConfig {
        cases: 24,
        rng_seed: RngSeed::Fixed(0x5EED),
        failure_persistence: None,
        ..ProptestConfig::default()
    }
}

/// The solvers advance each Fourier mode with explicit Euler: u_hat *= 1 + dt (-i c k - d k^2).
fn euler_factor(
    c: f64,
    d: f64,
    k: f64,
    dt: f64,
    steps: i32,
) -> Complex<f64> {
    (Complex::new(1.0, 0.0) + Complex::new(-d * k * k, -c * k) * dt).powi(steps)
}

#[test]
fn one_dimensional_mode_follows_the_explicit_euler_amplification() {
    let n = 64;
    let dx = 2.0 * PI / n as f64;
    let (c, d, dt, steps) = (1.0, 0.05, 0.001, 100);
    let u0: Vec<f64> = (0..n).map(|i| (3.0 * i as f64 * dx).cos()).collect();
    let res = solve_advection_diffusion_1d(&u0, dx, c, d, dt, steps);
    // cos(3x) = Re(e^{3ix}); the mode with wavenumber +3 carries the amplification.
    let r = euler_factor(c, d, 3.0, dt, steps as i32);
    for (i, &v) in res.iter().enumerate() {
        let x = i as f64 * dx;
        let exact = (r * Complex::new(0.0, 3.0 * x).exp()).re;
        assert_approx_eq!(v, exact, 1e-9);
    }
    // Sanity against the continuous solution exp(-d k^2 t) cos(k (x - c t)).
    let t = dt * steps as f64;
    let x = 10.0 * dx;
    let cont = (-d * 9.0 * t).exp() * (3.0 * (x - c * t)).cos();
    assert!((res[10] - cont).abs() < 5e-3, "{} vs {cont}", res[10]);
}

#[test]
fn pure_diffusion_preserves_the_mean_and_damps_the_variance() {
    let n = 128;
    let dx = 2.0 * PI / n as f64;
    let u0: Vec<f64> = (0..n)
        .map(|i| 2.0 + (i as f64 * dx).sin() + 0.5 * (5.0 * i as f64 * dx).cos())
        .collect();
    let res = solve_advection_diffusion_1d(&u0, dx, 0.0, 0.1, 0.001, 200);
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    assert_approx_eq!(mean(&res), 2.0, 1e-10);
    let var = |v: &[f64]| v.iter().map(|x| (x - 2.0).powi(2)).sum::<f64>() / v.len() as f64;
    assert!(var(&res) < var(&u0));
    // The sine mode decays by (1 - d dt)^steps.
    let f = (1.0f64 - 0.1 * 0.001).powi(200);
    assert_approx_eq!(
        res[32],
        2.0 + f * 1.0 + 0.5 * (1.0f64 - 0.1 * 25.0 * 0.001).powi(200) * (5.0 * 32.0 * dx).cos(),
        1e-9
    );
}

#[test]
fn zero_steps_returns_the_initial_data() {
    let n = 32;
    let u0: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
    let res = solve_advection_diffusion_1d(&u0, 0.1, 1.0, 1.0, 0.01, 0);
    for (a, b) in res.iter().zip(&u0) {
        assert_approx_eq!(a, b, 1e-12);
    }
}

#[test]
fn scenario_1d_has_expected_size_mass_and_peak_position() {
    let res = simulate_1d_advection_diffusion_scenario();
    assert_eq!(res.len(), 128);
    let dx = 2.0 * PI / 128.0;
    // Gaussian of width^2 = 0.25, mass ~ sqrt(pi * 0.5) / dx.
    let mass: f64 = res.iter().sum();
    assert_approx_eq!(mass, (PI * 0.5).sqrt() / dx, 0.05);
    // After t = 2 the pulse has been advected by c t = 2 to the right of pi.
    let imax = res
        .iter()
        .enumerate()
        .fold(0, |b, (i, &v)| if v > res[b] { i } else { b });
    let x = imax as f64 * dx;
    assert!((x - (PI + 2.0)).abs() < 0.2, "peak at {x}");
}

#[test]
fn fft2d_round_trip_and_dc_component() {
    let (w, h) = (8, 4);
    let orig: Vec<Complex<f64>> = (0..w * h)
        .map(|i| Complex::new((i as f64).sin(), (0.3 * i as f64).cos()))
        .collect();
    let mut data = orig.clone();
    fft2d(&mut data, w, h);
    let sum: Complex<f64> = orig.iter().sum();
    assert_approx_eq!(data[0].re, sum.re, 1e-9);
    assert_approx_eq!(data[0].im, sum.im, 1e-9);
    ifft2d(&mut data, w, h);
    for (a, b) in data.iter().zip(&orig) {
        assert!((a - b).norm() < 1e-9);
    }
}

#[test]
fn fft3d_dc_component_is_the_sum() {
    let (w, h, d) = (4, 4, 8);
    let orig: Vec<Complex<f64>> = (0..w * h * d)
        .map(|i| Complex::new((0.7 * i as f64).sin(), 0.0))
        .collect();
    let mut data = orig.clone();
    fft3d(&mut data, w, h, d);
    assert_eq!(data.len(), orig.len());
    let sum: f64 = orig.iter().map(|c| c.re).sum();
    assert_approx_eq!(data[0].re, sum, 1e-9);
}

#[test]
fn fft3d_round_trip() {
    let (w, h, d) = (4, 4, 8);
    let orig: Vec<Complex<f64>> = (0..w * h * d)
        .map(|i| Complex::new((0.7 * i as f64).sin(), 0.0))
        .collect();
    let mut data = orig.clone();
    fft3d(&mut data, w, h, d);
    ifft3d(&mut data, w, h, d);
    for (a, b) in data.iter().zip(&orig) {
        assert!((a - b).norm() < 1e-9);
    }
}

#[test]
fn fft3d_puts_a_plane_wave_in_the_right_bin_for_non_cubic_grids() {
    // Layout is x fastest, then y, then z. exp(2 pi i (a x/w + b y/h + c z/d)) -> single spike at (a, b, c).
    let (w, h, d) = (4usize, 2usize, 8usize);
    let (a, b, c) = (1usize, 1usize, 3usize);
    let mut data: Vec<Complex<f64>> = (0..w * h * d)
        .map(|idx| {
            let (x, y, z) = (idx % w, (idx / w) % h, idx / (w * h));
            let ph = 2.0
                * PI
                * (a as f64 * x as f64 / w as f64
                    + b as f64 * y as f64 / h as f64
                    + c as f64 * z as f64 / d as f64);
            Complex::new(0.0, ph).exp()
        })
        .collect();
    fft3d(&mut data, w, h, d);
    for (idx, v) in data.iter().enumerate() {
        let (x, y, z) = (idx % w, (idx / w) % h, idx / (w * h));
        let want = if (x, y, z) == (a, b, c) {
            (w * h * d) as f64
        } else {
            0.0
        };
        assert!((v.norm() - want).abs() < 1e-9, "bin ({x},{y},{z}) = {v}");
    }
}

#[test]
fn two_dimensional_mode_follows_the_explicit_euler_amplification() {
    let (w, h) = (16usize, 8usize);
    let (dx, dy) = (2.0 * PI / w as f64, 2.0 * PI / h as f64);
    let cfg = AdvectionDiffusionConfig {
        width: w,
        height: h,
        dx,
        dy,
        c: (0.5, -0.25),
        d: 0.02,
        dt: 0.002,
        steps: 50,
    };
    let mut u0 = vec![0.0; w * h];
    for j in 0..h {
        for i in 0..w {
            u0[j * w + i] = (2.0 * i as f64 * dx + j as f64 * dy).cos();
        }
    }
    let res = solve_advection_diffusion_2d(&u0, &cfg);
    // Mode (kx, ky) = (2, 1): -i (c.x kx + c.y ky) - d (kx^2 + ky^2).
    let lam = Complex::new(-0.02 * 5.0, -(0.5 * 2.0 - 0.25 * 1.0));
    let r = (Complex::new(1.0, 0.0) + lam * 0.002).powi(50);
    for j in 0..h {
        for i in 0..w {
            let ph = Complex::new(0.0, 2.0 * i as f64 * dx + j as f64 * dy).exp();
            assert_approx_eq!(res[j * w + i], (r * ph).re, 1e-9);
        }
    }
}

#[test]
fn scenario_2d_has_expected_size_and_conserves_mass() {
    let res = simulate_2d_advection_diffusion_scenario();
    assert_eq!(res.len(), 64 * 64);
    let dx = 2.0 * PI / 64.0;
    let mass: f64 = res.iter().sum::<f64>() * dx * dx;
    assert_approx_eq!(mass, PI * 0.5, 0.02); // integral of exp(-r^2 / 0.5)
    assert!(res.iter().all(|v| v.is_finite() && *v < 1.05));
}

#[test]
fn three_dimensional_mode_follows_the_explicit_euler_amplification() {
    let n = 8usize;
    let dx = 2.0 * PI / n as f64;
    let cfg = AdvectionDiffusionConfig3d {
        width: n,
        height: n,
        depth: n,
        dx,
        dy: dx,
        dz: dx,
        c: (1.0, 0.0, 0.0),
        d: 0.05,
        dt: 0.001,
        steps: 40,
    };
    let mut u0 = vec![0.0; n * n * n];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                u0[(k * n + j) * n + i] = 1.0 + (i as f64 * dx).cos();
            }
        }
    }
    let res = solve_advection_diffusion_3d(&u0, &cfg);
    let r = euler_factor(1.0, 0.05, 1.0, 0.001, 40);
    for i in 0..n {
        let exact = 1.0 + (r * Complex::new(0.0, i as f64 * dx).exp()).re;
        assert_approx_eq!(res[i], exact, 1e-9);
        assert_approx_eq!(res[(3 * n + 2) * n + i], exact, 1e-9);
    }
}

#[test]
fn scenario_3d_has_expected_size_and_is_finite() {
    let s = simulate_3d_advection_diffusion_scenario();
    assert_eq!(s.len(), 16 * 16 * 16);
    assert!(s.iter().all(|v| v.is_finite()));
}

#[test]
fn scenario_3d_conserves_mass() {
    let s = simulate_3d_advection_diffusion_scenario();
    let d = 2.0 * PI / 16.0;
    let mass: f64 = s.iter().sum::<f64>() * d * d * d;
    assert_approx_eq!(mass, (PI * 0.5).powf(1.5), 0.05);
}

proptest! {
    #![proptest_config(cfg())]

    #[test]
    fn prop_spectral_solver_is_linear(scale in 0.1f64..10.0, pos in 0usize..64) {
        const N: usize = 64;
        let mut ic1 = vec![0.0; N];
        ic1[pos] = 1.0;
        let ic2: Vec<f64> = ic1.iter().map(|x| x * scale).collect();
        let r1 = solve_advection_diffusion_1d(&ic1, 1.0, 1.0, 0.01, 0.01, 10);
        let r2 = solve_advection_diffusion_1d(&ic2, 1.0, 1.0, 0.01, 0.01, 10);
        for i in 0..N {
            prop_assert!((r1[i] * scale - r2[i]).abs() < 1e-9 * scale);
        }
    }

    #[test]
    fn prop_total_mass_is_conserved(seed in 0u64..1000, c in -2.0f64..2.0, d in 0.0f64..0.1) {
        let n = 32;
        let mut s = seed;
        let u0: Vec<f64> = (0..n).map(|_| { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); ((s >> 40) % 100) as f64 / 10.0 }).collect();
        let res = solve_advection_diffusion_1d(&u0, 0.2, c, d, 0.001, 20);
        let (m0, m1): (f64, f64) = (u0.iter().sum(), res.iter().sum());
        prop_assert!((m0 - m1).abs() < 1e-9 * (1.0 + m0));
    }
}
