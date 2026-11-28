//! Separable Cartesian FEM solves and their load contracts.

use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

use rssn::sim::physics_fem::PreparedCartesianPoisson;
use rssn::sim::physics_fem::solve_poisson_2d;
use rssn::sim::physics_fem::solve_poisson_3d;

#[test]
fn the_single_free_node_in_q1_squares_and_cubes_has_the_exact_discrete_value() {
    for elements in [&[2, 2][..], &[2, 2, 2]] {
        let plan = PreparedCartesianPoisson::new(elements).unwrap();
        let calls = AtomicUsize::new(0);
        let result = plan
            .solve(|point| {
                calls.fetch_add(1, Ordering::Relaxed);
                assert!(point[..elements.len()].iter().all(|&x| x > 0.0 && x < 1.0));
                1.0
            })
            .unwrap();
        assert_eq!(plan.nodal_shape(), vec![3; elements.len()]);
        assert_eq!(
            calls.load(Ordering::Relaxed),
            4usize.pow(elements.len() as u32)
        );
        let centre = result.len() / 2;
        assert!((result[centre] - 3.0 / 32.0).abs() < 2e-16);
        assert!(
            result
                .iter()
                .enumerate()
                .all(|(i, x)| i == centre || x.to_bits() == 0)
        );
    }
}

#[test]
fn existing_entry_points_preserve_their_quadrature_callback_multiplicity() {
    let calls = AtomicUsize::new(0);
    let two = solve_poisson_2d(3, 5, |_, _| {
        calls.fetch_add(1, Ordering::Relaxed);
        1.0
    })
    .unwrap();
    assert_eq!(calls.swap(0, Ordering::Relaxed), 3 * 5 * 16);
    let three = solve_poisson_3d(3, |_, _, _| {
        calls.fetch_add(1, Ordering::Relaxed);
        1.0
    })
    .unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), 27 * 64);
    assert_eq!(two.len(), 24);
    assert_eq!(three.len(), 64);
}

#[test]
fn sources_can_change_and_can_span_extreme_float_scales() {
    let plan = PreparedCartesianPoisson::new(&[9, 6, 5]).unwrap();
    let source = |[x, y, z]: [f64; 3]| 1.0 + x - 2.0 * y + 3.0 * z;
    let normal = plan.solve(source).unwrap();
    for scale in [1e-300, -1e300] {
        let scaled = plan.solve(|p| scale * source(p)).unwrap();
        for (a, b) in scaled.iter().zip(&normal) {
            assert!((a / scale - b).abs() < 1e-14);
        }
    }
    assert!(normal.iter().any(|&x| x > 0.1));
}

#[test]
fn unequal_tensor_grids_satisfy_the_assembled_q1_equations() {
    let cells = [7, 6, 5];
    let plan = PreparedCartesianPoisson::new(&cells).unwrap();
    let force = |[x, y, z]: [f64; 3]| 1.0 + 2.0 * x - 3.0 * y + 4.0 * z;
    let values = plan.solve(force).unwrap();
    let width = cells[0] + 1;
    let height = cells[1] + 1;
    let spacing = cells.map(|n| 1.0 / n as f64);
    let load_scale = spacing.iter().product::<f64>();
    for z in 1..cells[2] {
        for y in 1..cells[1] {
            for x in 1..cells[0] {
                let mut applied = 0.0;
                for oz in 0..3 {
                    for oy in 0..3 {
                        for ox in 0..3 {
                            let offsets = [ox, oy, oz];
                            let mass = std::array::from_fn::<_, 3, _>(|a| {
                                spacing[a]
                                    * (if offsets[a] == 1 {
                                        2.0 / 3.0
                                    } else {
                                        1.0 / 6.0
                                    })
                            });
                            let stiffness = std::array::from_fn::<_, 3, _>(|a| {
                                (if offsets[a] == 1 { 2.0 } else { -1.0 }) / spacing[a]
                            });
                            let coefficient = stiffness[0] * mass[1] * mass[2]
                                + mass[0] * stiffness[1] * mass[2]
                                + mass[0] * mass[1] * stiffness[2];
                            let index = x + ox - 1 + width * (y + oy - 1 + height * (z + oz - 1));
                            applied += coefficient * values[index];
                        }
                    }
                }
                let point = [x, y, z];
                let load =
                    load_scale * force(std::array::from_fn(|a| point[a] as f64 * spacing[a]));
                assert!(
                    (applied - load).abs() < 3e-15,
                    "{point:?}: {applied} != {load}"
                );
            }
        }
    }
}

#[test]
fn invalid_shapes_and_sources_are_checked_and_boundary_only_grids_are_zero() {
    for elements in [&[3][..], &[3, 0], &[usize::MAX, 2], &[3, 3, 3, 3]] {
        assert!(PreparedCartesianPoisson::new(elements).is_err());
    }
    let plan = PreparedCartesianPoisson::new(&[4, 1, 3]).unwrap();
    assert!(plan.solve(|_| f64::NAN).is_err());
    assert_eq!(plan.solve(|_| 3.0).unwrap(), vec![0.0; 40]);
    assert!(solve_poisson_2d(0, 4, |_, _| panic!("invalid grid called the source")).is_err());
}

#[test]
fn the_strided_sine_passes_agree_across_worker_counts() {
    let plan = PreparedCartesianPoisson::new(&[27, 24, 17]).unwrap();
    let solve = |n| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build()
            .unwrap()
            .install(|| {
                plan.solve(|[x, y, z]| (x + y + z).sin() + x * y * z)
                    .unwrap()
            })
    };
    let serial = solve(1);
    let parallel = solve(2);
    for (a, b) in serial.into_iter().zip(parallel) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}
