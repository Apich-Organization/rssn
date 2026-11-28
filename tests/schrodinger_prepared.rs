//! Physical-unit and observer regressions for 2D Schrodinger propagation.

use std::f64::consts::TAU;
use std::num::NonZeroUsize;

use num_complex::Complex64;
use rssn::sim::models_schrodinger_quantum::SchrodingerParameters;
use rssn::sim::models_schrodinger_quantum::SchrodingerPlan2D;
use rssn::sim::models_schrodinger_quantum::run_schrodinger_simulation;

fn parameters(
    nx: usize,
    ny: usize,
) -> SchrodingerParameters {
    SchrodingerParameters {
        nx,
        ny,
        lx: 3.25,
        ly: 5.5,
        dt: 0.09,
        time_steps: 25,
        hbar: 2.4,
        mass: 1.7,
        potential: vec![0.3; nx * ny],
    }
}

fn mode(
    nx: usize,
    ny: usize,
    mx: i32,
    my: i32,
    amplitude: Complex64,
) -> Vec<Complex64> {
    (0..nx * ny)
        .map(|i| {
            let x = (i % nx) as f64 / nx as f64;
            let y = (i / nx) as f64 / ny as f64;
            amplitude * Complex64::from_polar(1.0, TAU * (f64::from(mx) * x + f64::from(my) * y))
        })
        .collect()
}

fn phase(
    params: &SchrodingerParameters,
    mx: i32,
    my: i32,
    step: usize,
) -> Complex64 {
    let kx = f64::from(mx) * TAU / params.lx;
    let ky = f64::from(my) * TAU / params.ly;
    let rate = params.hbar * (kx.powi(2) + ky.powi(2)) / (2.0 * params.mass)
        + params.potential[0] / params.hbar;
    Complex64::from_polar(1.0, -rate * params.dt * step as f64)
}

fn assert_rotation(
    input: &[Complex64],
    output: &[Complex64],
    multiplier: Complex64,
    tolerance: f64,
) {
    assert_eq!(input.len(), output.len());
    for (&sample, &computed) in input.iter().zip(output) {
        assert!(
            (computed - sample * multiplier).norm() < tolerance,
            "{sample:?} {computed:?} {multiplier:?}"
        );
    }
}

#[test]
fn rectangular_mixed_parity_mode_and_snapshot_schedule() {
    let params = parameters(7, 6);
    let amplitude = Complex64::new(0.2, 0.9);
    let initial = mode(7, 6, 3, -2, amplitude);
    let mut wave = initial.clone();
    let frames = run_schrodinger_simulation(&params, &mut wave).unwrap();
    assert_eq!(frames.len(), 3);
    for frame in frames {
        assert_eq!(frame.dim(), (6, 7));
        assert!(
            frame
                .iter()
                .all(|v| (v - amplitude.norm_sqr()).abs() < 2e-13)
        );
    }
    assert_rotation(&initial, &wave, phase(&params, 3, -2, 25), 6e-14);

    let plan = SchrodingerPlan2D::new(&params).unwrap();
    let mut recorded_steps = Vec::new();
    let observed = plan
        .evolve(&initial, 25, NonZeroUsize::new(7), |step, wave| {
            recorded_steps.push(step);
            assert_rotation(&initial, wave, phase(&params, 3, -2, step), 6e-14);
            Ok(())
        })
        .unwrap();
    assert_eq!(recorded_steps, [1, 8, 15, 22]);
    assert_rotation(&initial, &observed, phase(&params, 3, -2, 25), 6e-14);
}

#[test]
fn uniform_potential_advances_a_superposition_over_a_long_interval() {
    let mut params = parameters(11, 5);
    params.potential.fill(-0.35);
    let plan = SchrodingerPlan2D::new(&params).unwrap();
    let a = mode(11, 5, 1, 2, Complex64::new(0.5, -0.2));
    let b = mode(11, 5, -5, 0, Complex64::new(-0.3, 0.2));
    let initial: Vec<_> = a.iter().zip(&b).map(|(a, b)| a + b).collect();
    let steps = 200_000;
    let evolved = plan
        .evolve(&initial, steps, None, |_, _| {
            panic!("observer without a stride")
        })
        .unwrap();
    for ((computed, a), b) in evolved.iter().zip(a).zip(b) {
        let expected = a * phase(&params, 1, 2, steps) + b * phase(&params, -5, 0, steps);
        assert!(
            (computed - expected).norm() < 9e-10,
            "{computed:?} {expected:?}"
        );
    }
}

#[test]
fn a_spatial_potential_preserves_norm_and_reverses_signed_time() {
    let mut params = parameters(5, 8);
    params.potential = (0..40).map(|i| (i as f64 * 0.42).cos() * 0.3).collect();
    let initial: Vec<_> = (0..40)
        .map(|i| Complex64::new((i as f64 * 0.7).sin(), (i as f64 * 0.31).cos()))
        .collect();
    let norm = |wave: &[Complex64]| wave.iter().map(Complex64::norm_sqr).sum::<f64>();
    let forward = SchrodingerPlan2D::new(&params)
        .unwrap()
        .evolve(&initial, 35, None, |_, _| Ok(()))
        .unwrap();
    params.dt = -params.dt;
    let restored = SchrodingerPlan2D::new(&params)
        .unwrap()
        .evolve(&forward, 35, NonZeroUsize::new(4), |_, sample| {
            assert!((norm(sample) - norm(&initial)).abs() < 1e-11);
            Ok(())
        })
        .unwrap();
    assert_rotation(&initial, &restored, Complex64::ONE, 1e-13);
}

#[test]
fn finite_phases_survive_extreme_physical_unit_combinations() {
    let mut params = parameters(1, 1);
    params.hbar = 2.0_f64.powi(1023);
    params.potential[0] = params.hbar;
    params.dt = 1.0;
    params.time_steps = 1;
    let initial = [Complex64::ONE];
    let mut wave = initial;
    run_schrodinger_simulation(&params, &mut wave).unwrap();
    assert_rotation(&initial, &wave, Complex64::from_polar(1.0, -1.0), 4e-15);

    let mut params = parameters(3, 1);
    params.lx = 2.0_f64.powi(-650);
    params.hbar = params.lx;
    params.mass = 2.0_f64.powi(650);
    params.dt = 0.125;
    params.potential.fill(0.0);
    let initial = mode(3, 1, 1, 0, Complex64::ONE);
    let evolved = SchrodingerPlan2D::new(&params)
        .unwrap()
        .evolve(&initial, 1, None, |_, _| Ok(()))
        .unwrap();
    // hbar / (mass L_x^2) is one, though each factor separately underflows.
    let expected = Complex64::from_polar(1.0, -0.5 * TAU.powi(2) * params.dt);
    assert_rotation(&initial, &evolved, expected, 4e-15);
}

#[test]
fn kinetic_phases_preserve_the_smaller_axis_before_phase_reduction() {
    let mut params = parameters(2, 2);
    params.lx = 2.0_f64.powi(-40);
    params.ly = 2.0;
    params.hbar = 1.0;
    params.mass = 1.0;
    params.dt = 1.0;
    params.potential.fill(0.0);
    let initial = [
        Complex64::ONE,
        -Complex64::ONE,
        -Complex64::ONE,
        Complex64::ONE,
    ];
    let evolved = SchrodingerPlan2D::new(&params)
        .unwrap()
        .evolve(&initial, 1, None, |_, _| Ok(()))
        .unwrap();
    let expected = Complex64::from_polar(1.0, -0.5 * (TAU / params.lx).powi(2))
        * Complex64::from_polar(1.0, -0.5 * (TAU / params.ly).powi(2));
    assert_rotation(&initial, &evolved, expected, 4e-15);
}

#[test]
fn failed_shape_or_density_leaves_the_mutable_entry_unchanged() {
    let mut params = parameters(2, 1);
    params.potential.pop();
    let mut wave = [Complex64::ONE; 2];
    assert!(run_schrodinger_simulation(&params, &mut wave).is_err());
    assert_eq!(wave, [Complex64::ONE; 2]);
    params.potential.push(0.3);
    params.nx = 0;
    assert!(run_schrodinger_simulation(&params, &mut wave).is_err());
    params.nx = 2;

    let large = [Complex64::new(1e200, -1e200); 2];
    let mut wave = large;
    assert!(
        run_schrodinger_simulation(&params, &mut wave)
            .unwrap_err()
            .contains("density")
    );
    assert_eq!(wave, large);
    let plan = SchrodingerPlan2D::new(&params).unwrap();
    let interrupted = plan.evolve(&large, 5, NonZeroUsize::new(2), |step, state| {
        assert_eq!(step, 1);
        assert!((state[0].re / 1e200).abs() > 0.5);
        Err("observer stopped".into())
    });
    assert_eq!(interrupted.unwrap_err(), "observer stopped");
    assert_eq!(large, wave);
}

#[test]
fn zero_time_and_empty_step_ranges_keep_the_input_bits() {
    let mut params = parameters(2, 2);
    params.dt = 0.0;
    let initial = [
        Complex64::new(-0.0, 0.0),
        Complex64::new(0.2, -0.0),
        Complex64::ONE,
        Complex64::ZERO,
    ];
    let plan = SchrodingerPlan2D::new(&params).unwrap();
    let bits = |wave: &[Complex64]| {
        wave.iter()
            .map(|c| (c.re.to_bits(), c.im.to_bits()))
            .collect::<Vec<_>>()
    };
    let mut steps = Vec::new();
    let final_state = plan
        .evolve(&initial, 5, NonZeroUsize::new(2), |step, state| {
            steps.push(step);
            assert_eq!(bits(state), bits(&initial));
            Ok(())
        })
        .unwrap();
    assert_eq!(steps, [1, 3, 5]);
    assert_eq!(bits(&initial), bits(&final_state));
    let empty = plan
        .evolve(&initial, 0, NonZeroUsize::new(1), |_, _| {
            panic!("empty evolution sampled")
        })
        .unwrap();
    assert_eq!(bits(&initial), bits(&empty));
}

#[test]
fn the_prepared_step_is_deterministic_across_active_workers() {
    let mut params = parameters(257, 129);
    for (i, value) in params.potential.iter_mut().enumerate() {
        *value = (i as f64 * 0.021).sin();
    }
    let plan = SchrodingerPlan2D::new(&params).unwrap();
    let initial = mode(params.nx, params.ny, -25, 38, Complex64::new(0.35, -0.78));
    let run = |workers| {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .unwrap();
        pool.install(|| plan.evolve(&initial, 2, None, |_, _| Ok(())).unwrap())
    };
    let serial = run(1);
    let parallel = run(4);
    for (serial, parallel) in serial.iter().zip(parallel) {
        assert_eq!(serial.re.to_bits(), parallel.re.to_bits());
        assert_eq!(serial.im.to_bits(), parallel.im.to_bits());
    }
}

#[test]
fn odd_packed_axes_keep_the_physical_mode_in_both_rectangular_orders() {
    for (nx, ny) in [(127, 65), (65, 127)] {
        let params = parameters(nx, ny);
        let initial = mode(nx, ny, 2, -3, Complex64::new(0.12, 0.8));
        let plan = SchrodingerPlan2D::new(&params).unwrap();
        let mut steps = Vec::new();
        let result = plan
            .evolve(&initial, 3, NonZeroUsize::new(2), |step, wave| {
                steps.push(step);
                assert_rotation(&initial, wave, phase(&params, 2, -3, step), 2e-13);
                Ok(())
            })
            .unwrap();
        assert_eq!(steps, [1, 3]);
        assert_rotation(&initial, &result, phase(&params, 2, -3, 3), 2e-13);
    }
}

#[test]
fn spatial_history_matches_an_independent_dense_fourier_reference() {
    let (nx, ny) = (3, 4);
    let mut params = parameters(nx, ny);
    params.time_steps = 12;
    params.potential = (0..nx * ny).map(|i| (i as f64 * 0.6).cos()).collect();
    let initial: Vec<_> = (0..nx * ny)
        .map(|i| Complex64::new((i as f64 * 0.4).sin(), (i as f64 * 0.3).cos()))
        .collect();
    let root3 = 3.0_f64.sqrt() / 2.0;
    let roots_x = [
        Complex64::ONE,
        Complex64::new(-0.5, root3),
        Complex64::new(-0.5, -root3),
    ];
    let roots_y = [Complex64::ONE, Complex64::I, -Complex64::ONE, -Complex64::I];
    let transform = |input: &[Complex64], inverse| {
        let scale = if inverse {
            1.0 / (nx * ny) as f64
        } else {
            1.0
        };
        (0..nx * ny)
            .map(|p| {
                let (kx, ky) = (p % nx, p / nx);
                input
                    .iter()
                    .enumerate()
                    .map(|(i, value)| {
                        let phase = roots_x[kx * (i % nx) % nx] * roots_y[ky * (i / nx) % ny];
                        *value
                            * if inverse {
                                phase
                            } else {
                                phase.conj()
                            }
                    })
                    .sum::<Complex64>()
                    * scale
            })
            .collect::<Vec<_>>()
    };
    let mut actual = initial.clone();
    let frames = run_schrodinger_simulation(&params, &mut actual).unwrap();
    assert_eq!(frames.len(), 2);

    let potential: Vec<_> = params
        .potential
        .iter()
        .map(|v| Complex64::from_polar(1.0, -v * params.dt / (2.0 * params.hbar)))
        .collect();
    let k = |slot: usize, count: usize, span: f64| {
        let mode = if slot < count.div_ceil(2) {
            slot as f64
        } else {
            slot as f64 - count as f64
        };
        TAU * mode / span
    };
    let mut expected = initial;
    for step in 1..=params.time_steps {
        for (sample, v) in expected.iter_mut().zip(&potential) {
            *sample *= v;
        }
        let mut spec = transform(&expected, false);
        for (i, value) in spec.iter_mut().enumerate() {
            let ks = k(i % nx, nx, params.lx).powi(2) + k(i / nx, ny, params.ly).powi(2);
            *value *=
                Complex64::from_polar(1.0, -params.hbar * ks * params.dt / (2.0 * params.mass));
        }
        expected = transform(&spec, true);
        for (sample, v) in expected.iter_mut().zip(&potential) {
            *sample *= v;
        }
        if step == 1 || step == 11 {
            let frame = &frames[usize::from(step == 11)];
            for (value, sample) in frame.iter().zip(&expected) {
                let expected = sample.norm_sqr();
                assert!(
                    (value - expected).abs() < 3e-14 * (1.0 + expected),
                    "step {step}: density {value} vs {expected}"
                );
            }
        }
    }
    assert_rotation(&expected, &actual, Complex64::ONE, 4e-14);
}
