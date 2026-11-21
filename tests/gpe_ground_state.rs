//! Focused checks for the two-dimensional fixed-step imaginary-time solver.

use ndarray::Array2;
use num_complex::Complex64;
use rssn::sim::models_gpe_superfluidity::GpeParameters;
use rssn::sim::models_gpe_superfluidity::run_gpe_ground_state_finder;
use std::f64::consts::TAU;

fn parameters() -> GpeParameters {
    GpeParameters {
        nx: 3,
        ny: 5,
        lx: 6.0,
        ly: 7.5,
        d_tau: 0.07,
        time_steps: 3,
        g: 0.8,
        trap_strength: 1.4,
    }
}

fn close(
    got: f64,
    expected: f64,
    tolerance: f64,
) {
    assert!(
        (got - expected).abs() <= tolerance * expected.abs().max(1e-280),
        "{got} vs {expected}"
    );
}

fn potential(
    p: &GpeParameters,
    i: usize,
    j: usize,
) -> f64 {
    let x = (i as f64 / p.nx as f64 - 0.5) * p.lx;
    let y = (j as f64 / p.ny as f64 - 0.5) * p.ly;
    0.5 * ((p.trap_strength * x) * x + (p.trap_strength * y) * y)
}

fn area(p: &GpeParameters) -> f64 {
    p.lx / p.nx as f64 * (p.ly / p.ny as f64)
}

#[test]
fn no_steps_and_zero_increment_return_the_normalized_initial_profile() {
    let mut p = parameters();
    p.time_steps = 0;
    let no_steps = run_gpe_ground_state_finder(&p).unwrap();
    assert_eq!(no_steps.dim(), (5, 3));

    let guess_density =
        Array2::from_shape_fn((p.ny, p.nx), |(j, i)| (-0.2 * potential(&p, i, j)).exp());
    let total = guess_density.sum() * area(&p);
    for (&got, expected) in no_steps.iter().zip(guess_density.iter().map(|x| x / total)) {
        close(got, expected, 5e-15);
    }

    p.time_steps = 17;
    p.d_tau = 0.0;
    let zero_increment = run_gpe_ground_state_finder(&p).unwrap();
    assert!(
        zero_increment
            .iter()
            .zip(no_steps.iter())
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}

fn dft(
    data: &Array2<Complex64>,
    inverse: bool,
) -> Array2<Complex64> {
    let (ny, nx) = data.dim();
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut output = Array2::from_shape_fn(data.dim(), |(ky, kx)| {
        data.indexed_iter()
            .map(|((y, x), &val)| {
                let turns = (x * kx) as f64 / nx as f64 + (y * ky) as f64 / ny as f64;
                val * Complex64::from_polar(1.0, sign * TAU * turns)
            })
            .sum::<Complex64>()
    });
    if inverse {
        output.mapv_inplace(|val| val / (nx * ny) as f64);
    }
    output
}

fn direct_reference(p: &GpeParameters) -> Array2<f64> {
    let v = Array2::from_shape_fn((p.ny, p.nx), |(j, i)| potential(p, i, j));
    let potential_half = |psi: &mut Array2<Complex64>| {
        for (val, &trap) in psi.iter_mut().zip(&v) {
            *val *= (-0.5 * p.d_tau * (trap + p.g * val.norm_sqr())).exp();
        }
    };
    let mut psi = v.mapv(|val| Complex64::new((-0.1 * val).exp(), 0.0));
    for _ in 0..p.time_steps {
        potential_half(&mut psi);
        let mut spectrum = dft(&psi, false);
        for ((j, i), val) in spectrum.indexed_iter_mut() {
            let kx = TAU * i.min(p.nx - i) as f64 / p.lx;
            let ky = TAU * j.min(p.ny - j) as f64 / p.ly;
            *val *= (-0.5 * (kx * kx + ky * ky) * p.d_tau).exp();
        }
        psi = dft(&spectrum, true);
        potential_half(&mut psi);
        let norm = (psi.iter().map(Complex64::norm_sqr).sum::<f64>() * area(p)).sqrt();
        psi.mapv_inplace(|val| val / norm);
    }
    psi.mapv(|val| val.norm_sqr())
}

#[test]
fn rectangular_odd_grid_matches_a_direct_fixed_step_fourier_oracle() {
    for (nx, ny) in [(3, 5), (4, 3), (1, 5), (5, 1)] {
        for g in [-0.8, 0.0, 0.8] {
            let p = GpeParameters {
                nx,
                ny,
                g,
                ..parameters()
            };
            let result = run_gpe_ground_state_finder(&p).unwrap();
            let expected = direct_reference(&p);
            for (&got, &want) in result.iter().zip(&expected) {
                close(got, want, 4e-13);
            }
            close(result.sum() * area(&p), 1.0, 4e-15);
        }
    }
}

#[test]
fn large_area_small_area_and_anisotropic_free_states_keep_their_scale() {
    let mut p = GpeParameters {
        nx: 1,
        ny: 1,
        lx: 1e160,
        ly: 1e160,
        g: 0.0,
        trap_strength: 0.0,
        ..parameters()
    };
    for steps in [0, 1, 3] {
        p.time_steps = steps;
        let result = run_gpe_ground_state_finder(&p).unwrap();
        assert!(result[[0, 0]].to_bits().abs_diff(1e-320_f64.to_bits()) <= 1);
    }

    p.lx = 1e-154;
    p.ly = 1e-154;
    close(
        run_gpe_ground_state_finder(&p).unwrap()[[0, 0]],
        1e308,
        4e-15,
    );

    p.nx = 2;
    p.ny = 4;
    p.lx = 1e-309;
    p.ly = 1e305;
    let per_area = (1.0 / p.ly) / p.lx;
    for &got in &run_gpe_ground_state_finder(&p).unwrap() {
        close(got, per_area, 4e-14);
    }
}

#[test]
fn scales_the_initial_state_and_trap_before_squaring() {
    let p = GpeParameters {
        nx: 1,
        ny: 1,
        lx: 2.0,
        ly: 2.0,
        time_steps: 0,
        trap_strength: 5_500.0,
        ..parameters()
    };
    close(
        run_gpe_ground_state_finder(&p).unwrap()[[0, 0]],
        0.25,
        4e-15,
    );
    let tiny_amplitude = GpeParameters {
        g: 0.0,
        time_steps: 1,
        d_tau: 0.001,
        ..p
    };
    close(
        run_gpe_ground_state_finder(&tiny_amplitude).unwrap()[[0, 0]],
        0.25,
        4e-15,
    );
    let tiny_initial_weight = GpeParameters {
        trap_strength: 20_000.0,
        ..p
    };
    close(
        run_gpe_ground_state_finder(&tiny_initial_weight).unwrap()[[0, 0]],
        0.25,
        4e-15,
    );
    let p = GpeParameters {
        nx: 3,
        ny: 1,
        lx: 2e160,
        ly: 1e-118,
        trap_strength: 6e-321,
        ..p
    };
    let rho = run_gpe_ground_state_finder(&p).unwrap();
    let weight = |i| (-0.2 * potential(&p, i, 0)).exp();
    let sum = (0..p.nx).map(weight).sum::<f64>() * area(&p);
    for (i, &got) in rho.iter().enumerate() {
        close(got, weight(i) / sum, 4e-15);
    }
}

#[test]
fn invalid_grid_and_unrepresentable_density_are_errors() {
    let p = parameters();
    for bad in [
        GpeParameters { nx: 0, ..p },
        GpeParameters { ny: 0, ..p },
        GpeParameters { nx: usize::MAX, ..p },
        GpeParameters {
            nx: isize::MAX as usize,
            ny: 1,
            ..p
        },
        GpeParameters {
            lx: f64::from_bits(1),
            ..p
        },
        GpeParameters { lx: -1.0, ..p },
        GpeParameters {
            ly: f64::INFINITY,
            ..p
        },
        GpeParameters { d_tau: -0.01, ..p },
        GpeParameters {
            d_tau: f64::INFINITY,
            ..p
        },
        GpeParameters { g: f64::NAN, ..p },
        GpeParameters {
            trap_strength: f64::NEG_INFINITY,
            ..p
        },
        GpeParameters {
            lx: 1e-155,
            ly: 1e-155,
            g: 0.0,
            trap_strength: 0.0,
            ..p
        },
        GpeParameters {
            lx: 1e200,
            ly: 1e200,
            g: 0.0,
            trap_strength: 0.0,
            ..p
        },
    ] {
        assert!(run_gpe_ground_state_finder(&bad).is_err(), "{bad:?}");
    }
}

#[test]
fn fixed_step_results_match_in_one_and_two_worker_pools() {
    for (nx, ny) in [
        (65, 128),
        (128, 65),
        (63, 132),
        (128, 257),
        (128, 256),
        (512, 512),
    ] {
        let p = GpeParameters {
            nx,
            ny,
            lx: 14.0,
            ly: 9.0,
            time_steps: 2,
            ..parameters()
        };
        let run = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| run_gpe_ground_state_finder(&p))
                .unwrap()
        };
        let one = run(1);
        let two = run(2);
        assert!(
            one.iter()
                .zip(&two)
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
        close(one.sum() * area(&p), 1.0, 1e-13);
    }
}
