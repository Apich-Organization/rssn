//! Discrete time-start, energy and boundary tests for Dirichlet wave grids.

use std::f64::consts::PI;

use rssn::sim::physics_fdm::Dimensions;
use rssn::sim::physics_fdm::FdmGrid;
use rssn::sim::physics_fdm::FdmSolverConfig2D;
use rssn::sim::physics_fdm::FdmSolverConfig3D;
use rssn::sim::physics_fdm::solve_wave_dirichlet_sampled;
use rssn::sim::physics_fdm::solve_wave_equation_2d;
use rssn::sim::physics_fdm::solve_wave_equation_3d;

#[test]
fn zero_velocity_three_dimensional_start_has_half_acceleration() {
    let config = FdmSolverConfig3D {
        width: 5,
        height: 5,
        depth: 5,
        dx: 1.0,
        dy: 1.0,
        dz: 1.0,
        dt: 0.25,
        steps: 1,
    };
    let result = solve_wave_equation_3d(&config, 1.0, |x, y, z| {
        f64::from(u8::from((x, y, z) == (2, 2, 2)))
    });
    assert_eq!(result[(2, 2, 2)], 0.8125);
    for neighbor in [
        (1, 2, 2),
        (3, 2, 2),
        (2, 1, 2),
        (2, 3, 2),
        (2, 2, 1),
        (2, 2, 3),
    ] {
        assert_eq!(result[neighbor], 0.03125);
    }
    assert_eq!(result[(1, 1, 1)], 0.0);
}

#[test]
fn legacy_positive_steps_project_faces_before_reading_the_first_stencil() {
    let two = FdmSolverConfig2D {
        width: 4,
        height: 6,
        dx: 1.0,
        dy: 0.8,
        dt: 0.2,
        steps: 1,
    };
    let three = FdmSolverConfig3D {
        width: 4,
        height: 6,
        depth: 3,
        dx: 1.0,
        dy: 0.8,
        dz: 1.0,
        dt: 0.2,
        steps: 1,
    };
    let field = |x, y| {
        if x == 0 || x == 3 || y == 0 || y == 5 {
            11.0
        } else {
            0.0
        }
    };
    let field3 = |x, y, z| {
        if z == 0 || z == 2 {
            11.0
        } else {
            field(x, y)
        }
    };
    let initial = solve_wave_equation_2d(
        &FdmSolverConfig2D {
            steps: 0,
            ..two.clone()
        },
        1.0,
        field,
    );
    assert_eq!(initial[(0, 2)], 11.0);
    assert!(
        solve_wave_equation_2d(&two, 1.0, field)
            .as_slice()
            .iter()
            .all(|&v| v == 0.0)
    );
    assert!(
        solve_wave_equation_3d(&three, 1.0, field3)
            .as_slice()
            .iter()
            .all(|&v| v == 0.0)
    );
}

#[test]
fn signed_sampled_velocity_matches_the_discrete_standing_mode() {
    let [nx, ny, nz] = [9, 7, 5];
    let [dx, dy, dz] = [0.8, 1.1, 1.4];
    let c = 1.2;
    let dimensions = Dimensions::D3(nx, ny, nz);
    let mode: Vec<_> = (0..nx * ny * nz)
        .map(|p| {
            let (x, y, z) = (p % nx, p / nx % ny, p / (nx * ny));
            (PI * x as f64 / (nx - 1) as f64).sin()
                * (PI * y as f64 / (ny - 1) as f64).sin()
                * (PI * z as f64 / (nz - 1) as f64).sin()
        })
        .collect();
    let displacement = FdmGrid::from_data(mode.clone(), dimensions.clone());
    let velocity = FdmGrid::from_data(mode.iter().map(|&value| -0.7 * value).collect(), dimensions);
    for dt in [-0.2, 0.2] {
        let samples =
            solve_wave_dirichlet_sampled(&displacement, &velocity, &[dx, dy, dz], c, dt, 11, 4)
                .unwrap();
        assert_eq!(
            samples.iter().map(|sample| sample.step).collect::<Vec<_>>(),
            [0, 4, 8, 11]
        );
        let weighted = [nx, ny, nz]
            .into_iter()
            .zip([dx, dy, dz])
            .map(|(n, h)| {
                let s = c * dt / h;
                let k = (0.5 * PI / (n - 1) as f64).sin();
                s * s * k * k
            })
            .sum::<f64>();
        let theta = (1.0 - 2.0 * weighted).acos();
        for sample in samples {
            let phase = sample.step as f64 * theta;
            let amplitude = phase.cos() - 0.7 * dt * phase.sin() / theta.sin();
            for (obtained, phi) in sample.displacement.as_slice().iter().zip(&mode) {
                assert!((obtained - amplitude * phi).abs() < 3e-13);
            }
        }
    }
}

#[test]
fn leapfrog_keeps_the_discrete_dirichlet_energy() {
    let [nx, ny] = [13, 11];
    let [dx, dy] = [0.8, 1.1];
    let c = 1.2;
    let dt = 0.3;
    let dims = Dimensions::D2(nx, ny);
    let field = (0..nx * ny).map(|p| ((p as f64) * 0.23).sin()).collect();
    let rate = (0..nx * ny)
        .map(|p| 0.1 * ((p as f64) * 0.11).cos())
        .collect();
    let initial = FdmGrid::from_data(field, dims.clone());
    let velocity = FdmGrid::from_data(rate, dims);
    let state =
        solve_wave_dirichlet_sampled(&initial, &velocity, &[dx, dy], c, dt, 256, 1).unwrap();
    let energy = |older: &[f64], newer: &[f64]| {
        let mut energy = 0.0;
        for p in 0..older.len() {
            let rate = newer[p] - older[p];
            energy += rate * rate;
            for (stride, spacing, neighbor) in [(1, dx, p % nx + 1 < nx), (nx, dy, p / nx + 1 < ny)]
            {
                if neighbor {
                    let s = c * dt / spacing;
                    energy +=
                        s * s * (newer[p + stride] - newer[p]) * (older[p + stride] - older[p]);
                }
            }
        }
        0.5 * energy
    };
    let expected = energy(
        state[0].displacement.as_slice(),
        state[1].displacement.as_slice(),
    );
    for pair in state.windows(2) {
        let actual = energy(
            pair[0].displacement.as_slice(),
            pair[1].displacement.as_slice(),
        );
        assert!((actual - expected).abs() < 2e-12 * expected);
    }
}

#[test]
fn checked_sampling_handles_bounds_and_underflowed_travel_distance() {
    let dims = Dimensions::D2(3, 3);
    let initial = FdmGrid::from_data(
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        dims.clone(),
    );
    let velocity = FdmGrid::with_value(dims.clone(), 0.0);
    let tiny_spacing = f64::from_bits(1);
    let step = solve_wave_dirichlet_sampled(
        &initial,
        &velocity,
        &[tiny_spacing; 2],
        f64::MIN_POSITIVE,
        0.5 * f64::EPSILON,
        1,
        1,
    )
    .unwrap();
    assert_eq!(step.last().unwrap().displacement[(1, 1)], 0.5);

    let constant = solve_wave_dirichlet_sampled(
        &initial,
        &velocity,
        &[1.0; 2],
        1.0,
        -0.0,
        usize::MAX,
        usize::MAX,
    )
    .unwrap();
    assert_eq!(constant.len(), 2);
    assert_eq!(constant[1].step, usize::MAX);
    assert_eq!(
        constant[0].displacement.as_slice(),
        constant[1].displacement.as_slice()
    );

    let malformed = FdmGrid::from_data(vec![0.0; 8], dims);
    assert!(solve_wave_dirichlet_sampled(&initial, &malformed, &[1.0; 2], 1.0, 0.2, 1, 1).is_err());
    assert!(solve_wave_dirichlet_sampled(&initial, &velocity, &[1.0; 2], 1.0, 0.9, 1, 1).is_err());
    assert!(solve_wave_dirichlet_sampled(&initial, &velocity, &[1.0; 2], 1.0, 0.2, 1, 0).is_err());
}

#[test]
fn zero_wave_speed_advances_velocity_without_forming_a_laplacian() {
    let dims = Dimensions::D2(3, 3);
    let center = |value| {
        let mut grid = FdmGrid::with_value(dims.clone(), 0.0);
        grid[(1, 1)] = value;
        grid
    };
    let linear =
        solve_wave_dirichlet_sampled(&center(0.25), &center(0.125), &[1.0; 2], 0.0, -0.5, 3, 2)
            .unwrap();
    assert_eq!(linear.last().unwrap().displacement[(1, 1)], 0.0625);
    let rate = f64::from_bits(1);
    let huge_time =
        solve_wave_dirichlet_sampled(&center(1.0), &center(rate), &[1.0; 2], 0.0, f64::MAX, 2, 2)
            .unwrap();
    assert_eq!(
        huge_time.last().unwrap().displacement[(1, 1)],
        1.0 + 2.0 * (f64::MAX * rate)
    );
    let stationary = solve_wave_dirichlet_sampled(
        &center(f64::MAX),
        &center(0.0),
        &[1.0; 2],
        0.0,
        f64::MAX,
        2,
        2,
    )
    .unwrap();
    assert_eq!(stationary.last().unwrap().displacement[(1, 1)], f64::MAX);
}

#[test]
fn sampled_rows_are_bitwise_deterministic_across_worker_counts() {
    let [nx, ny] = [131, 257];
    let dims = Dimensions::D2(nx, ny);
    let initial = FdmGrid::from_data(
        (0..nx * ny).map(|p| (p as f64 * 0.017).sin()).collect(),
        dims.clone(),
    );
    let velocity = FdmGrid::with_value(dims, 0.25);
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                solve_wave_dirichlet_sampled(&initial, &velocity, &[0.6, 1.4], 1.0, 0.2, 7, 3)
                    .unwrap()
            })
    };
    let (serial, parallel) = (run(1), run(2));
    for (a, b) in serial.iter().zip(&parallel) {
        assert_eq!(a.step, b.step);
        assert!(
            a.displacement
                .as_slice()
                .iter()
                .zip(b.displacement.as_slice())
                .all(|(a, b)| a.to_bits() == b.to_bits())
        );
    }
}
