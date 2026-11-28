//! Periodic scalar DG restart and prepared-state regression tests.

use rssn::sim::dg::DgAdvection1dPlan;
use rssn::sim::dg::DgMesh;
use rssn::sim::dg::advection_step;
use rssn::sim::dg::advection_step_nd;

#[test]
fn periodic_single_element_matches_the_exact_p1_ssprk_mode() {
    let mesh = DgMesh::line(1, 0.0, 1.0, 1);
    let initial = [0.25, 1.75];
    for velocity in [-2.0_f64, 2.0] {
        let plan = DgAdvection1dPlan::new(&mesh, velocity).unwrap();
        let mut state = initial;
        let mut output = [0.0; 2];
        let mut amplitude = -0.75;
        let mut time = -0.5;
        for dt in [0.0625, 0.03125, 0.125, 0.0625, 0.03125] {
            // The mass mode is constant; the P1 periodic difference eigenvalue is
            // -2 |a| / h. SSPRK3 applies its cubic stability polynomial.
            let z = -2.0 * velocity.abs() * dt;
            amplitude *= 1.0 + z + z * z / 2.0 + z * z * z / 6.0;
            plan.step_into(&state, dt, &mut output).unwrap();
            assert!((output[0] - (1.0 + amplitude)).abs() < 6e-15);
            assert!((output[1] - (1.0 - amplitude)).abs() < 6e-15);
            assert_eq!(plan.step(&state, dt).unwrap(), output);
            assert_eq!(
                advection_step(&mesh, velocity, &state, time, dt).unwrap(),
                output
            );
            assert_eq!(
                advection_step_nd(&mesh, [velocity, 5.0], &state, time, dt).unwrap(),
                output
            );
            state = output;
            time += dt;
        }
    }
}

#[test]
fn periodic_multielement_plan_conserves_the_quadrature_mean() {
    let mesh = DgMesh::line(8, -2.0, 2.0, 3);
    let initial: Vec<_> = (0..32)
        .map(|i| ((i * 11 % 29) as f64 - 9.0) / 16.0)
        .collect();
    let sum = |u: &[f64]| {
        u.as_chunks::<4>()
            .0
            .iter()
            .map(|cell| (cell[0] + 5.0 * cell[1] + 5.0 * cell[2] + cell[3]) / 24.0)
            .sum::<f64>()
    };
    let mass = sum(&initial);
    for velocity in [-0.5, 0.0, 0.5] {
        let plan = DgAdvection1dPlan::new(&mesh, velocity).unwrap();
        let evolve = || {
            let (mut state, mut output) = (initial.clone(), vec![0.0; initial.len()]);
            for _ in 0..64 {
                plan.step_into(&state, 0.0625, &mut output).unwrap();
                std::mem::swap(&mut state, &mut output);
            }
            state
        };
        let (left, right) = rayon::join(evolve, evolve);
        assert!((sum(&left) - mass).abs() < 2e-14);
        assert_eq!(
            left.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            right.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }
}

#[test]
fn wrong_inputs_and_unrepresentable_results_leave_the_output_untouched() {
    let mesh = DgMesh::line(1, 0.0, 1.0, 1);
    let plan = DgAdvection1dPlan::new(&mesh, 2.0).unwrap();
    let mut output = [-0.0, 123.45];
    let original = output.map(f64::to_bits);
    for (state, dt) in [
        (&[1.0][..], 0.0625),
        (&[1.0, f64::NAN][..], 0.0625),
        (&[1.0, 2.0][..], 0.0),
        (&[1.0, 2.0][..], f64::INFINITY),
        (&[1.0, 2.0][..], 0.25),
        (&[f64::MAX, f64::MAX][..], 0.0625),
    ] {
        assert!(plan.step_into(state, dt, &mut output).is_err());
        assert_eq!(output.map(f64::to_bits), original);
    }
    assert!(
        plan.step_into(&[1.0, 2.0], 0.0625, &mut output[..1])
            .is_err()
    );
    assert_eq!(output.map(f64::to_bits), original);
    assert!(DgAdvection1dPlan::new(&mesh, f64::NAN).is_err());
    assert!(
        DgAdvection1dPlan::new(
            &DgMesh {
                p: usize::MAX,
                ..mesh
            },
            2.0
        )
        .is_err()
    );
    assert!(DgAdvection1dPlan::new(&DgMesh { dim: 2, ..mesh }, 2.0).is_err());
    assert!(advection_step(&mesh, 2.0, &[1.0, 2.0], f64::INFINITY, 0.0625).is_err());
}
