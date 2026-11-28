//! Fourier mode reuse for the 3D Forward Euler advection–diffusion solver.

use num_complex::Complex64;
use rssn::sim::physics_sm::AdvectionDiffusionConfig3d;
use rssn::sim::physics_sm::advance_advection_diffusion_3d_modes;
use rssn::sim::physics_sm::fft3d;
use rssn::sim::physics_sm::ifft3d;
use rssn::sim::physics_sm::solve_advection_diffusion_3d;
use std::f64::consts::TAU;

fn config() -> AdvectionDiffusionConfig3d {
    AdvectionDiffusionConfig3d {
        width: 7,
        height: 5,
        depth: 9,
        dx: TAU / 7.0,
        dy: 2.0 * TAU / 5.0,
        dz: 1.5 * TAU / 9.0,
        c: (0.4, -0.7, 0.125),
        d: 0.02,
        dt: 0.0005,
        steps: 113,
    }
}

fn phase(i: usize) -> f64 {
    TAU * ((i % 7) as f64 / 7.0 - 2.0 * ((i / 7) % 5) as f64 / 5.0 + 3.0 * (i / 35) as f64 / 9.0)
}

fn compare_complex(
    left: &[Complex64],
    right: &[Complex64],
) {
    assert_eq!(left.len(), right.len());
    for (a, b) in left.iter().zip(right) {
        assert_eq!(a.re.to_bits(), b.re.to_bits());
        assert_eq!(a.im.to_bits(), b.im.to_bits());
    }
}

#[test]
fn mixed_odd_grid_plane_wave_matches_euler_in_every_worker_pool() {
    let cfg = config();
    let input: Vec<_> = (0..315).map(|i| 1.2 + 0.75 * phase(i).cos()).collect();
    // The independent wave has physical k = (1, -1, 2), so |k|² = 6.
    let ck = cfg.c.0 - cfg.c.1 + 2.0 * cfg.c.2;
    let gain = Complex64::new(1.0 - cfg.dt * cfg.d * 6.0, -cfg.dt * ck).powu(113);
    let serial = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap()
        .install(|| solve_advection_diffusion_3d(&input, &cfg));
    for (i, &u) in serial.iter().enumerate() {
        let exact = 1.2 + (gain * Complex64::from_polar(0.75, phase(i))).re;
        assert!((u - exact).abs() < 1e-12);
    }
    for threads in [2, 4] {
        let output = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| solve_advection_diffusion_3d(&input, &cfg));
        for (a, b) in serial.iter().zip(&output) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}

#[test]
fn mode_reuse_composes_steps_and_uses_the_existing_fft_layout() {
    let cfg = config();
    let input: Vec<_> = (0..315).map(|i| 1.2 + 0.75 * phase(i).cos()).collect();
    let mut modes: Vec<_> = input.iter().map(|&u| Complex64::new(u, 0.0)).collect();
    fft3d(&mut modes, 7, 5, 9);
    let mut single = modes.clone();
    advance_advection_diffusion_3d_modes(&mut single, &cfg).unwrap();
    let first = AdvectionDiffusionConfig3d {
        steps: 19,
        ..cfg.clone()
    };
    let rest = AdvectionDiffusionConfig3d {
        steps: 94,
        ..cfg.clone()
    };
    advance_advection_diffusion_3d_modes(&mut modes, &first).unwrap();
    advance_advection_diffusion_3d_modes(&mut modes, &rest).unwrap();
    compare_complex(&single, &modes);

    ifft3d(&mut modes, 7, 5, 9);
    let physical = solve_advection_diffusion_3d(&input, &cfg);
    for (a, b) in modes.iter().zip(physical) {
        assert_eq!(a.re.to_bits(), b.to_bits());
    }
}

#[test]
fn zero_steps_and_zero_time_step_preserve_bits() {
    let fixed = AdvectionDiffusionConfig3d {
        width: 2,
        height: 1,
        depth: 1,
        ..config()
    };
    let original = [
        Complex64::new(-0.0, f64::MAX),
        Complex64::new(f64::from_bits(1), -0.0),
    ];
    for cfg in [
        AdvectionDiffusionConfig3d {
            steps: 0,
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            dt: 0.0,
            ..fixed.clone()
        },
    ] {
        let mut modes = original;
        advance_advection_diffusion_3d_modes(&mut modes, &cfg).unwrap();
        compare_complex(&modes, &original);
    }
}

#[test]
fn invalid_shape_and_nonfinite_updates_leave_modes_unchanged() {
    let fixed = AdvectionDiffusionConfig3d {
        width: 2,
        height: 1,
        depth: 1,
        ..config()
    };
    let initial = [Complex64::new(0.25, 0.75), Complex64::new(-2.0, -1.0)];
    for cfg in [
        AdvectionDiffusionConfig3d {
            width: usize::MAX,
            height: 2,
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            height: 0,
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            depth: 2,
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            dx: f64::MAX,
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            dz: 0.0,
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            d: -0.1,
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            dt: -0.0005,
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            c: (0.4, f64::INFINITY, 0.0),
            ..fixed.clone()
        },
        AdvectionDiffusionConfig3d {
            d: f64::MAX,
            ..fixed.clone()
        },
    ] {
        let mut modes = initial;
        assert!(advance_advection_diffusion_3d_modes(&mut modes, &cfg).is_err());
        compare_complex(&modes, &initial);
    }
    let invalid = [initial[0], Complex64::new(f64::NAN, 1.0)];
    let mut modes = invalid;
    assert!(advance_advection_diffusion_3d_modes(&mut modes, &fixed).is_err());
    compare_complex(&modes, &invalid);
}
