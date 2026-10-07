//! Split-step Fourier tests for the Schrodinger equation and the GPE.

use num_complex::Complex;
use rssn::sim::split_step::{SplitOrder, SplitStep};
use std::f64::consts::PI;

fn harmonic(dim: usize) -> impl Fn(f64, [f64; 3]) -> f64 {
    move |_t, x| 0.5 * (0..dim).map(|a| x[a] * x[a]).sum::<f64>()
}

fn gaussian_guess(s: &SplitStep, w: f64) -> Vec<Complex<f64>> {
    let dim = s.grid.dim;
    s.sample(|x| Complex::new((-(0..dim).map(|a| x[a] * x[a]).sum::<f64>() / (2.0 * w * w)).exp(), 0.0))
}

#[test]
fn imaginary_time_finds_harmonic_oscillator_ground_state_in_1d_2d_3d() {
    for (dim, n, len) in [(1, 64, 16.0), (2, 32, 12.0), (3, 16, 10.0)] {
        let s = SplitStep::centered(&vec![n; dim], &vec![len; dim], 1.0).unwrap();
        let v = harmonic(dim);
        let gs = s.ground_state(&gaussian_guess(&s, 2.0), &v, 0.0, 0.05, 2000, 1e-13);
        let exact = 0.5 * dim as f64;
        // splitting bias is O(dt^2)
        assert!((gs.energy - exact).abs() < 2e-3 * exact, "dim {dim}: E {} vs {exact}", gs.energy);
        assert!((gs.mu - gs.energy).abs() < 1e-10);
        // wavefunction close to the Gaussian
        let norm = PI.powf(-0.25 * dim as f64);
        let err = (0..s.grid.size())
            .map(|p| {
                let x = s.x(p);
                let r2: f64 = (0..dim).map(|a| x[a] * x[a]).sum();
                (gs.psi[p].norm() - norm * (-0.5 * r2).exp()).abs()
            })
            .fold(0.0, f64::max);
        assert!(err < 2e-3, "dim {dim}: psi error {err}");
    }
}

#[test]
fn imaginary_time_energy_bias_is_second_order_in_dt() {
    let s = SplitStep::centered(&[64], &[16.0], 1.0).unwrap();
    let v = harmonic(1);
    let bias = |dt: f64| {
        let gs = s.ground_state(&gaussian_guess(&s, 2.0), &v, 0.0, dt, 5000, 1e-14);
        (gs.energy - 0.5).abs()
    };
    let (b1, b2) = (bias(0.2), bias(0.1));
    let ratio = b1 / b2;
    assert!((3.0..5.5).contains(&ratio), "ratio {ratio} ({b1}, {b2})");
}

fn coherent(x: f64, t: f64, x0: f64) -> Complex<f64> {
    let xc = x0 * t.cos();
    let pc = -x0 * t.sin();
    let amp = PI.powf(-0.25);
    Complex::from_polar(amp * (-(x - xc).powi(2) / 2.0).exp(), pc * x + x0 * x0 * (2.0 * t).sin() / 4.0 - 0.5 * t)
}

#[test]
fn real_time_evolution_orders_strang_two_yoshida_four() {
    let s = SplitStep::centered(&[128], &[24.0], 1.0).unwrap();
    let v = harmonic(1);
    let psi0 = s.sample(|x| coherent(x[0], 0.0, 1.5));
    let t_end = 2.0;
    let err = |dt: f64, order| {
        let steps = (t_end / dt).round() as usize;
        let psi = s.evolve(&psi0, &v, 0.0, 0.0, dt, steps, order);
        (0..psi.len()).map(|p| (psi[p] - coherent(s.x(p)[0], t_end, 1.5)).norm()).fold(0.0, f64::max)
    };
    let (a1, a2) = (err(0.04, SplitOrder::Strang), err(0.02, SplitOrder::Strang));
    assert!((3.5..4.5).contains(&(a1 / a2)), "strang {a1} {a2}");
    let (b1, b2) = (err(0.1, SplitOrder::Yoshida4), err(0.05, SplitOrder::Yoshida4));
    assert!(b1 / b2 > 12.0, "yoshida {b1} {b2}");
    assert!(b2 < a2);
}

#[test]
fn free_gaussian_spreads_exactly_and_norm_is_conserved() {
    for (dim, n) in [(1, 64), (2, 32), (3, 16)] {
        let len = 12.0;
        let n = if dim == 3 { 32 } else { n };
        let s = SplitStep::centered(&vec![n; dim], &vec![len; dim], 1.0).unwrap();
        let v = |_t: f64, _x: [f64; 3]| 0.0;
        let psi0 = gaussian_guess(&s, 1.0);
        let psi = s.evolve(&psi0, &v, 0.0, 0.0, 0.1, 10, SplitOrder::Strang);
        let t = 1.0;
        let denom = Complex::new(1.0, t);
        let pre = denom.powf(-0.5 * dim as f64);
        // compare to the analytic solution of the (un-normalised) Gaussian initial datum
        let norm_free = s.norm2(&psi0);
        assert!((s.norm2(&psi) - norm_free).abs() < 1e-12 * norm_free);
        let mut worst: f64 = 0.0;
        for p in 0..s.grid.size() {
            let x = s.x(p);
            let r2: f64 = (0..dim).map(|a| x[a] * x[a]).sum();
            let ex = pre * (-Complex::new(r2, 0.0) / (2.0 * denom)).exp();
            worst = worst.max((psi[p] - ex).norm());
        }
        assert!(worst < 1e-8, "dim {dim}: {worst}");
    }
}

#[test]
fn gpe_conserves_norm_and_energy_and_bright_soliton_is_stationary() {
    let s = SplitStep::centered(&[64], &[20.0], 1.0).unwrap();
    let v = harmonic(1);
    let g = 5.0;
    let psi0 = s.sample(|x| Complex::new((-(x[0] - 0.5).powi(2) / 2.0).exp() * PI.powf(-0.25), 0.0));
    let e0 = s.energy(&psi0, &v, 0.0, g);
    let psi = s.evolve(&psi0, &v, g, 0.0, 0.005, 400, SplitOrder::Yoshida4);
    assert!((s.norm2(&psi) - s.norm2(&psi0)).abs() < 1e-12);
    assert!((s.energy(&psi, &v, 2.0, g) - e0).abs() < 1e-6 * e0.abs(), "energy drift");

    // attractive GPE (g = -1): psi = sech(x) exp(i t / 2)
    let s = SplitStep::centered(&[128], &[40.0], 1.0).unwrap();
    let zero = |_t: f64, _x: [f64; 3]| 0.0;
    let sech = |x: f64| 1.0 / x.cosh();
    let psi0 = s.sample(|x| Complex::new(sech(x[0]), 0.0));
    let t = 2.0;
    let psi = s.evolve(&psi0, &zero, -1.0, 0.0, 0.01, 200, SplitOrder::Yoshida4);
    let worst = (0..psi.len())
        .map(|p| (psi[p] - Complex::from_polar(sech(s.x(p)[0]), 0.5 * t)).norm())
        .fold(0.0, f64::max);
    assert!(worst < 1e-5, "soliton error {worst}");
}

#[test]
fn gpe_thomas_fermi_chemical_potential_in_1d() {
    let s = SplitStep::centered(&[256], &[40.0], 1.0).unwrap();
    let v = harmonic(1);
    let g = 200.0;
    let gs = s.ground_state(&gaussian_guess(&s, 3.0), &v, g, 0.005, 40_000, 1e-12);
    let mu_tf = (3.0 * g / (4.0 * 2.0_f64.sqrt())).powf(2.0 / 3.0);
    assert!(gs.converged);
    assert!((gs.mu - mu_tf).abs() < 0.04 * mu_tf, "mu {} vs TF {mu_tf}", gs.mu);
    assert!(gs.mu > gs.energy, "mu exceeds the energy per particle for repulsive g");
    assert!((s.norm2(&gs.psi) - 1.0).abs() < 1e-12);
}
