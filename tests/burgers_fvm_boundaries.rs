//! Discrete conservation, boundaries and scale limits for scalar Burgers FVM.

use rssn::sim::physics_fvm::BurgersBoundary;
use rssn::sim::physics_fvm::solve_burgers_nd;
use rssn::sim::physics_fvm::solve_burgers_nd_with_boundaries;

#[test]
fn periodic_one_step_matches_the_lax_friedrichs_stencil() {
    use BurgersBoundary::Periodic;
    let source = [-1.0, 0.25, 0.5, -0.75];
    let dt = 0.125;
    for (dims, spacing) in [
        (&[4][..], &[1.0][..]),
        (&[1, 4, 1][..], &[0.7, 1.0, 2.0][..]),
    ] {
        let actual = solve_burgers_nd_with_boundaries(
            &source,
            dims,
            spacing,
            dt,
            1,
            &[Periodic; 3][..dims.len()],
        )
        .unwrap();
        let viscosity = 0.5 / dims.len() as f64;
        for (i, value) in actual.iter().enumerate() {
            let left = source[(i + 3) % 4];
            let right = source[(i + 1) % 4];
            let expected = source[i] + viscosity * (left + right - 2.0 * source[i])
                - 0.25 * dt * (right * right - left * left);
            assert!((value - expected).abs() < 2e-15);
        }
    }
}

#[test]
fn mixed_boundaries_balance_the_physical_outflow_flux() {
    use BurgersBoundary::Outflow;
    use BurgersBoundary::Periodic;
    let [nx, ny, nz] = [5, 4, 3];
    let spacing = [0.7, 1.4, 2.1];
    let dt = 0.05;
    let source: Vec<_> = (0..nx * ny * nz)
        .map(|p: usize| {
            0.05 * (p % nx) as f64 + 0.15 * (p / nx % ny) as f64 - 0.17 * (p / (nx * ny)) as f64
        })
        .collect();
    let output = solve_burgers_nd_with_boundaries(
        &source,
        &[nx, ny, nz],
        &spacing,
        dt,
        1,
        &[Periodic, Outflow, Periodic],
    )
    .unwrap();
    let mut net_flux = 0.0;
    for z in 0..nz {
        for x in 0..nx {
            let low = source[x + nx * ny * z];
            let high = source[x + nx * (ny - 1) + nx * ny * z];
            net_flux += 0.5 * (high * high - low * low);
        }
    }
    let expected = source.iter().sum::<f64>() - dt / spacing[1] * net_flux;
    assert!((output.iter().sum::<f64>() - expected).abs() < 1e-13);
    let min = source.iter().copied().fold(f64::INFINITY, f64::min);
    let max = source.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert!(output.iter().all(|value| (min..=max).contains(value)));
}

#[test]
fn periodic_conservation_and_thread_determinism_hold_over_many_steps() {
    use BurgersBoundary::Periodic;
    let dims = [65, 137, 2];
    let source: Vec<_> = (0..dims.iter().product::<usize>())
        .map(|p| (p as f64 * 0.017).sin() * 0.27 + (p as f64 * 0.012).cos() * 0.31 + 0.12)
        .collect();
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                solve_burgers_nd_with_boundaries(
                    &source,
                    &dims,
                    &[0.8, 1.4, 2.0],
                    0.12,
                    16,
                    &[Periodic; 3],
                )
                .unwrap()
            })
    };
    let serial = run(1);
    let parallel = run(2);
    assert!(
        serial
            .iter()
            .zip(parallel)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    let initial_mass = source.iter().sum::<f64>();
    assert!((serial.iter().sum::<f64>() - initial_mass).abs() < 5e-10);
}

#[test]
fn extreme_scales_preserve_the_dimensionless_step() {
    use BurgersBoundary::Periodic;
    let mesh = [4];
    let spacing = [1.0];
    let boundaries = [Periodic];
    let zero_step_limit = solve_burgers_nd_with_boundaries(
        &[0.0, 1.0, 0.0, 0.0],
        &mesh,
        &spacing,
        f64::from_bits(1),
        1,
        &boundaries,
    )
    .unwrap();
    assert_eq!(zero_step_limit, [0.5, 0.0, 0.5, 0.0]);

    let source = [0.3, -0.7, 0.5, -0.2];
    for (amplitude, dt) in [(1e200, 1e-201), (1e-200, 1e199)] {
        let scaled = source.map(|value| value * amplitude);
        let actual =
            solve_burgers_nd_with_boundaries(&scaled, &mesh, &spacing, dt, 3, &boundaries).unwrap();
        let expected = solve_burgers_nd_with_boundaries(
            &source,
            &mesh,
            &spacing,
            amplitude * dt,
            3,
            &boundaries,
        )
        .unwrap();
        for (actual, expected) in actual.iter().zip(expected) {
            assert!((actual / amplitude - expected).abs() < 2e-15);
        }
    }
    let uniform = [f64::MAX; 4];
    assert_eq!(
        solve_burgers_nd(&uniform, &mesh, &spacing, f64::MAX, 3).unwrap(),
        uniform
    );
}

#[test]
fn invalid_shapes_are_rejected_and_default_outflow_agrees_with_explicit_faces() {
    use BurgersBoundary::Outflow;
    let source = [-0.1, 0.2, 0.25, -0.3, 0.15, 0.4];
    let default = solve_burgers_nd(&source, &[3, 2], &[0.9, 1.1], 0.15, 4).unwrap();
    let explicit =
        solve_burgers_nd_with_boundaries(&source, &[3, 2], &[0.9, 1.1], 0.15, 4, &[Outflow; 2])
            .unwrap();
    assert_eq!(default, explicit);
    assert_eq!(
        solve_burgers_nd(&source, &[6], &[1.0], 0.1, 0).unwrap(),
        source
    );

    assert!(solve_burgers_nd(&[], &[usize::MAX, 2], &[1.0; 2], 0.1, 1).is_err());
    assert!(solve_burgers_nd(&source, &[6], &[f64::INFINITY], 0.1, 1).is_err());
    assert!(solve_burgers_nd(&source, &[6], &[1.0], f64::INFINITY, 1).is_err());
    assert!(solve_burgers_nd(&[f64::NAN], &[1], &[1.0], 0.1, 0).is_err());
    assert!(
        solve_burgers_nd_with_boundaries(&source, &[3, 2], &[1.0; 2], 0.1, 1, &[Outflow]).is_err()
    );
}
