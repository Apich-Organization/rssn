//! Split-step Fourier tests for the Schrodinger equation and the GPE.

use num_complex::Complex;
use rssn::sim::split_step::SplitOrder;
use rssn::sim::split_step::SplitStep;
use std::f64::consts::PI;

fn harmonic(dim: usize) -> impl Fn(f64, [f64; 3]) -> f64 {
    move |_t, x| 0.5 * (0..dim).map(|a| x[a] * x[a]).sum::<f64>()
}

fn gaussian_guess(
    s: &SplitStep,
    w: f64,
) -> Vec<Complex<f64>> {
    let dim = s.grid.dim;
    s.sample(|x| {
        Complex::new(
            (-(0..dim).map(|a| x[a] * x[a]).sum::<f64>() / (2.0 * w * w)).exp(),
            0.0,
        )
    })
}

#[test]
fn imaginary_time_finds_harmonic_oscillator_ground_state_in_1d_2d_3d() {
    for (dim, n, len) in [(1, 64, 16.0), (2, 32, 12.0), (3, 16, 10.0)] {
        let s = SplitStep::centered(&vec![n; dim], &vec![len; dim], 1.0).unwrap();
        let v = harmonic(dim);
        let gs = s.ground_state(&gaussian_guess(&s, 2.0), &v, 0.0, 0.05, 2000, 1e-13);
        let exact = 0.5 * dim as f64;
        // splitting bias is O(dt^2)
        assert!(
            (gs.energy - exact).abs() < 2e-3 * exact,
            "dim {dim}: E {} vs {exact}",
            gs.energy
        );
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
fn imaginary_time_state_error_is_second_order_and_energy_error_fourth_order_in_dt() {
    let s = SplitStep::centered(&[64], &[16.0], 1.0).unwrap();
    let v = harmonic(1);
    let run = |dt: f64| {
        let gs = s.ground_state(&gaussian_guess(&s, 2.0), &v, 0.0, dt, 5000, 1e-15);
        let err = (0..s.grid.size())
            .map(|p| (gs.psi[p].norm() - PI.powf(-0.25) * (-0.5 * s.x(p)[0].powi(2)).exp()).abs())
            .fold(0.0, f64::max);
        (err, (gs.energy - 0.5).abs())
    };
    let ((a1, e1), (a2, e2)) = (run(0.4), run(0.2));
    // the energy is stationary at the ground state, so its error is the square of the state error
    assert!(
        (3.0..5.5).contains(&(a1 / a2)),
        "state ratio {} ({a1}, {a2})",
        a1 / a2
    );
    assert!(e1 / e2 > 10.0, "energy ratio {} ({e1}, {e2})", e1 / e2);
}

fn coherent(
    x: f64,
    t: f64,
    x0: f64,
) -> Complex<f64> {
    let xc = x0 * t.cos();
    let pc = -x0 * t.sin();
    let amp = PI.powf(-0.25);
    Complex::from_polar(
        amp * (-(x - xc).powi(2) / 2.0).exp(),
        pc * x + x0 * x0 * (2.0 * t).sin() / 4.0 - 0.5 * t,
    )
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
        (0..psi.len())
            .map(|p| (psi[p] - coherent(s.x(p)[0], t_end, 1.5)).norm())
            .fold(0.0, f64::max)
    };
    let (a1, a2) = (err(0.04, SplitOrder::Strang), err(0.02, SplitOrder::Strang));
    assert!((3.5..4.5).contains(&(a1 / a2)), "strang {a1} {a2}");
    let (b1, b2) = (
        err(0.1, SplitOrder::Yoshida4),
        err(0.05, SplitOrder::Yoshida4),
    );
    assert!(b1 / b2 > 12.0, "yoshida {b1} {b2}");
    assert!(b2 < a2);
}

#[test]
fn free_gaussian_spreads_exactly_and_norm_is_conserved() {
    for (dim, n) in [(1, 64), (2, 32), (3, 16)] {
        let len = 16.0;
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
        assert!(worst < 1e-5, "dim {dim}: {worst}");
    }
}

#[test]
fn gpe_conserves_norm_and_energy_and_bright_soliton_is_stationary() {
    let s = SplitStep::centered(&[64], &[20.0], 1.0).unwrap();
    let v = harmonic(1);
    let g = 5.0;
    let psi0 =
        s.sample(|x| Complex::new((-(x[0] - 0.5).powi(2) / 2.0).exp() * PI.powf(-0.25), 0.0));
    let e0 = s.energy(&psi0, &v, 0.0, g);
    let psi = s.evolve(&psi0, &v, g, 0.0, 0.005, 400, SplitOrder::Yoshida4);
    // 400 FFT round trips: the drift is round-off (~1e-12), not a scheme error.
    assert!(
        (s.norm2(&psi) - s.norm2(&psi0)).abs() < 1e-11,
        "norm drift {:e}",
        s.norm2(&psi) - s.norm2(&psi0)
    );
    assert!(
        (s.energy(&psi, &v, 2.0, g) - e0).abs() < 1e-6 * e0.abs(),
        "energy drift"
    );

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
    assert!(
        (gs.mu - mu_tf).abs() < 0.04 * mu_tf,
        "mu {} vs TF {mu_tf}",
        gs.mu
    );
    assert!(
        gs.mu > gs.energy,
        "mu exceeds the energy per particle for repulsive g"
    );
    assert!((s.norm2(&gs.psi) - 1.0).abs() < 1e-12);
}

mod rotating_frame {
    use super::*;

    fn packet(
        s: &SplitStep,
        centre: [f64; 2],
        widths: [f64; 2],
    ) -> Vec<Complex<f64>> {
        s.sample(|x| {
            Complex::new(
                (-((x[0] - centre[0]).powi(2) / (2.0 * widths[0].powi(2))
                    + (x[1] - centre[1]).powi(2) / (2.0 * widths[1].powi(2))))
                .exp(),
                0.0,
            )
        })
    }

    fn mean(
        s: &SplitStep,
        psi: &[Complex<f64>],
    ) -> [f64; 2] {
        let n = s.norm2(psi);
        let mut c = [0.0; 2];
        let cell = s.grid.h(0) * s.grid.h(1);
        for (p, v) in psi.iter().enumerate() {
            let x = s.x(p);
            c[0] += x[0] * v.norm_sqr() * cell / n;
            c[1] += x[1] * v.norm_sqr() * cell / n;
        }
        c
    }

    fn diff(
        a: &[Complex<f64>],
        b: &[Complex<f64>],
    ) -> f64 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).norm())
            .fold(0.0, f64::max)
    }

    #[test]
    fn a_rotation_by_shears_moves_the_pattern_and_conserves_the_norm() {
        let s = SplitStep::centered(&[64, 64], &[16.0, 16.0], 1.0).unwrap();
        let mut psi = packet(&s, [3.0, 0.0], [0.9, 0.9]);
        let n0 = s.norm2(&psi);
        let angle = 0.9_f64;
        s.rotate(&mut psi, angle);
        assert!((s.norm2(&psi) - n0).abs() < 1e-12 * n0);
        let c = mean(&s, &psi);
        assert!((c[0] - 3.0 * angle.cos()).abs() < 1e-6, "{c:?}");
        assert!((c[1] - 3.0 * angle.sin()).abs() < 1e-6, "{c:?}");
        // against the exact rotated Gaussian
        let exact = packet(&s, [3.0 * angle.cos(), 3.0 * angle.sin()], [0.9, 0.9]);
        assert!(diff(&psi, &exact) < 1e-5, "{}", diff(&psi, &exact));
        // rotating forth and back restores the data
        let mut back = psi.clone();
        s.rotate(&mut back, -angle);
        assert!(diff(&back, &packet(&s, [3.0, 0.0], [0.9, 0.9])) < 1e-5);
    }

    fn free_packet_difference(n: usize) -> f64 {
        // [L_z, Laplacian] = 0, so in a rotating frame the free evolution is the lab evolution
        // rotated clockwise by Omega t.
        let s = SplitStep::centered(&[n, n], &[24.0, 24.0], 1.0).unwrap();
        let psi0 = packet(&s, [1.5, -1.0], [1.0, 0.6]);
        let zero = |_t: f64, _x: [f64; 3]| 0.0;
        let (omega, dt, steps) = (0.7, 0.02, 50);
        let t = dt * steps as f64;
        let rot = s.evolve_rotating(&psi0, &zero, 0.0, omega, 0.0, dt, steps, SplitOrder::Strang);
        let mut lab = s.evolve(&psi0, &zero, 0.0, 0.0, dt, steps, SplitOrder::Strang);
        s.rotate(&mut lab, -omega * t);
        let n0 = s.norm2(&psi0);
        assert!(
            (s.norm2(&rot) - n0).abs() < 1e-10 * n0,
            "norm drift {}",
            (s.norm2(&rot) - n0).abs() / n0
        );
        diff(&rot, &lab)
    }

    #[test]
    fn a_free_anisotropic_packet_is_the_rotated_lab_solution() {
        let (d64, d128) = (free_packet_difference(64), free_packet_difference(128));
        println!("free packet, rotating vs rotated lab: n = 64 {d64:.2e}, n = 128 {d128:.2e}");
        assert!(d64 < 5e-3, "difference {d64}");
        assert!(
            d128 < 0.5 * d64,
            "no convergence with the grid: {d64} -> {d128}"
        );
    }

    #[test]
    fn a_packet_in_an_isotropic_trap_follows_the_rotated_oscillation() {
        // lab: <x> = x0 cos(t); rotating frame: that vector turned clockwise by Omega t
        let s = SplitStep::centered(&[64, 64], &[20.0, 20.0], 1.0).unwrap();
        let v = harmonic(2);
        let x0 = 2.0;
        let psi0 = packet(&s, [x0, 0.0], [1.0, 1.0]); // the ground state of the trap, displaced
        let (omega, dt, steps) = (0.35, 0.01, 300);
        let t = dt * steps as f64;
        let psi = s.evolve_rotating(&psi0, &v, 0.0, omega, 0.0, dt, steps, SplitOrder::Yoshida4);
        let c = mean(&s, &psi);
        let amp = x0 * t.cos();
        let want = [amp * (omega * t).cos(), -amp * (omega * t).sin()];
        assert!(
            (c[0] - want[0]).abs() < 2e-3 && (c[1] - want[1]).abs() < 2e-3,
            "{c:?} vs {want:?}"
        );
        let n0 = s.norm2(&psi0);
        assert!((s.norm2(&psi) - n0).abs() < 1e-11 * n0);
    }

    #[test]
    fn zero_rate_is_exactly_the_plain_evolution() {
        let s = SplitStep::centered(&[32, 32], &[12.0, 12.0], 1.0).unwrap();
        let v = harmonic(2);
        let psi0 = packet(&s, [1.0, 0.5], [1.0, 1.0]);
        let a = s.evolve_rotating(&psi0, &v, 0.0, 0.0, 0.0, 0.02, 20, SplitOrder::Strang);
        let b = s.evolve(&psi0, &v, 0.0, 0.0, 0.02, 20, SplitOrder::Strang);
        assert_eq!(a, b);
    }

    #[test]
    fn a_trap_with_gpe_nonlinearity_keeps_the_norm_in_a_rotating_frame() {
        let s = SplitStep::centered(&[32, 32], &[14.0, 14.0], 1.0).unwrap();
        let v = harmonic(2);
        let psi0 = packet(&s, [0.8, 0.0], [1.2, 0.9]);
        let n0 = s.norm2(&psi0);
        let psi = s.evolve_rotating(&psi0, &v, 2.0, 0.6, 0.0, 0.02, 200, SplitOrder::Strang);
        assert!((s.norm2(&psi) - n0).abs() < 1e-11 * n0);
    }
}
