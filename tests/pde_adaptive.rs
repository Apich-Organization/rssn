//! Analytic semidiscrete oracles distinguish temporal control from grid error.
use rssn::sim::fdm_nd::Boundary3D;
use rssn::sim::fdm_nd::BoundaryCondition;
use rssn::sim::fdm_nd::TimeStepping;
use rssn::sim::fdm_nd::adaptive::PdeAdaptiveOptions;
use rssn::sim::fdm_nd::adaptive::PdeTerminalEvent;
use rssn::sim::fdm_nd::adaptive::{
    self,
};

fn options(end: f64) -> PdeAdaptiveOptions {
    PdeAdaptiveOptions {
        t1: end,
        h0: end,
        h_max: end,
        rtol: 1e-6,
        atol: 1e-9,
        ..Default::default()
    }
}
fn mode(rank: usize) -> Vec<f64> {
    (0..4usize.pow(rank as u32))
        .map(|mut i| {
            let mut u = 1.0;
            for _ in 0..rank {
                u *= (std::f64::consts::FRAC_PI_2 * (i % 4) as f64).cos();
                i /= 4;
            }
            u
        })
        .collect()
}

#[test]
fn each_heat_method_controls_error_in_one_two_and_three_dimensions() {
    let bc = Boundary3D::uniform(BoundaryCondition::Periodic);
    for rank in 1..=3 {
        let extents = vec![4; rank];
        let initial = mode(rank);
        for method in [
            TimeStepping::Explicit,
            TimeStepping::BackwardEuler,
            TimeStepping::CrankNicolson,
            TimeStepping::Adi,
        ] {
            let result = adaptive::heat_adaptive(
                &extents,
                1.0,
                1.0,
                &initial,
                &bc,
                method,
                &options(0.2),
                None,
            )
            .unwrap();
            let exact = (-2.0 * rank as f64 * 0.2).exp();
            let error = (result.states.last().unwrap()[0] - exact).abs();
            assert!(error < 0.002, "{rank} {method:?}: {error}");
            assert_eq!(*result.times.last().unwrap(), 0.2);
            assert!(result.error_ratio <= 1.0);
            assert!(result.rejected > 0);
            assert_eq!(result.states.len(), result.accepted + 1);
        }
    }
}

#[test]
fn tighter_heat_tolerance_reduces_analytic_temporal_error() {
    let bc = Boundary3D::uniform(BoundaryCondition::Periodic);
    let initial = mode(1);
    let mut loose = options(0.4);
    loose.rtol = 1e-3;
    let mut tight = loose;
    tight.rtol = 1e-7;
    let a = adaptive::heat_adaptive(
        &[4],
        1.0,
        1.0,
        &initial,
        &bc,
        TimeStepping::CrankNicolson,
        &loose,
        None,
    )
    .unwrap();
    let b = adaptive::heat_adaptive(
        &[4],
        1.0,
        1.0,
        &initial,
        &bc,
        TimeStepping::CrankNicolson,
        &tight,
        None,
    )
    .unwrap();
    let exact = (-0.8_f64).exp();
    assert!(
        (b.states.last().unwrap()[0] - exact).abs() < (a.states.last().unwrap()[0] - exact).abs()
    );
    assert!(b.accepted > a.accepted);
}

#[test]
fn wave_keeps_velocity_state_and_advection_obeys_its_semidiscrete_solution() {
    let bc = Boundary3D::uniform(BoundaryCondition::Periodic);
    let initial = mode(1);
    let wave = adaptive::wave_adaptive(
        &[4],
        1.0,
        1.0,
        &initial,
        &[0.0; 4],
        1e-9,
        &bc,
        &options(0.2),
        None,
    )
    .unwrap();
    let last = wave.states.last().unwrap();
    let frequency = 2.0_f64.sqrt();
    assert!((last[0] - (frequency * 0.2).cos()).abs() < 1e-4);
    assert!((last[4] + frequency * (frequency * 0.2).sin()).abs() < 1e-4);
    assert!(wave.rejected > 0);
    for rank in 1..=3 {
        let mut initial = vec![0.0; 4usize.pow(rank as u32)];
        for (i, u) in initial.iter_mut().enumerate() {
            *u = (std::f64::consts::FRAC_PI_2 * (i % 4) as f64).cos();
        }
        let mut speeds = vec![0.0; rank];
        speeds[0] = 1.0;
        let result = adaptive::advection_adaptive(
            &vec![4; rank],
            1.0,
            &speeds,
            &initial,
            &bc,
            &options(0.2),
            None,
        )
        .unwrap();
        assert!(
            (result.states.last().unwrap()[0] - (-0.2_f64).exp() * 0.2_f64.cos()).abs() < 0.002
        );
    }
}

#[test]
fn terminal_events_reintegrate_and_budgets_never_force_acceptance() {
    let bc = Boundary3D::uniform(BoundaryCondition::Periodic);
    let initial = mode(1);
    let mut value = |_: f64, y: &[f64]| Ok(y[0] - 0.9);
    let mut event = PdeTerminalEvent {
        direction: -1,
        evaluate: &mut value,
    };
    let result = adaptive::heat_adaptive(
        &[4],
        1.0,
        1.0,
        &initial,
        &bc,
        TimeStepping::CrankNicolson,
        &options(0.3),
        Some(&mut event),
    )
    .unwrap();
    assert!((result.event_time.unwrap() + 0.9_f64.ln() / 2.0).abs() < 1e-4);
    assert!((result.states.last().unwrap()[0] - 0.9).abs() < 1e-10);
    let mut limited = options(1.0);
    limited.max_step_calls = 2;
    assert!(
        adaptive::heat_adaptive(
            &[4],
            1.0,
            1.0,
            &initial,
            &bc,
            TimeStepping::CrankNicolson,
            &limited,
            None
        )
        .is_err()
    );
    limited = options(1.0);
    limited.rtol = 0.0;
    limited.atol = 1e-30;
    limited.max_rejections = 1;
    assert!(
        adaptive::heat_adaptive(
            &[4],
            1.0,
            1.0,
            &initial,
            &bc,
            TimeStepping::Explicit,
            &limited,
            None
        )
        .is_err()
    );
    let mut conflicting = Boundary3D::uniform(BoundaryCondition::Dirichlet(0.0));
    conflicting.x_min = BoundaryCondition::Dirichlet(1.0);
    assert!(
        adaptive::heat_adaptive(
            &[4, 4],
            1.0,
            1.0,
            &mode(2),
            &conflicting,
            TimeStepping::Explicit,
            &options(1.0),
            None
        )
        .is_err()
    );
}

#[test]
fn advanced_evolvers_use_error_control_on_their_original_time_schemes() {
    use pde_api::adaptive as advanced;
    use rssn::sim::boundary;
    use rssn::sim::dg;
    use rssn::sim::fem_simplex;
    use rssn::sim::pde_api;
    let ctl = options(0.1);
    let spectral = advanced::spectral_etdrk4_adaptive(
        &[4],
        &[4.0],
        0.0,
        0.0,
        &[1.0; 4],
        |_, u| u.to_vec(),
        &ctl,
        None,
    )
    .unwrap();
    assert!((spectral.states.last().unwrap()[0] - 0.1_f64.exp()).abs() < 1e-5);
    let mesh = dg::DgMesh::line(2, 0.0, 1.0, 2);
    let dg = advanced::dg_advection_adaptive(&mesh, 1.0, &[1.0; 6], &ctl, None).unwrap();
    assert!(dg.states.iter().flatten().all(|u| (*u - 1.0).abs() < 1e-12));
    for scheme in [
        pde_api::FvmScheme::Upwind,
        pde_api::FvmScheme::LaxFriedrichs,
        pde_api::FvmScheme::MusclMinmod,
        pde_api::FvmScheme::MusclVanLeer,
    ] {
        let result =
            advanced::fvm_adaptive(1.0, &[1.0; 8], true, scheme, |u| u, |_| 1.0, &ctl, None)
                .unwrap();
        assert!(
            result
                .states
                .iter()
                .flatten()
                .all(|u| (*u - 1.0).abs() < 1e-12)
        );
    }
    let space =
        fem_simplex::FemSpace::new(fem_simplex::SimplexMesh::interval(4, 0.0, 1.0), 1).unwrap();
    let bcs = boundary::MarkerBcs::new()
        .set(
            fem_simplex::marker::XMIN,
            boundary::BoundaryCondition::neumann(0.0),
        )
        .set(
            fem_simplex::marker::XMAX,
            boundary::BoundaryCondition::neumann(0.0),
        );
    let source = |_: f64, _: [f64; 3]| 0.0;
    let problem = fem_simplex::FemProblem {
        k: 1.0,
        c: 1.0,
        f: &source,
        bcs: &bcs,
    };
    for theta in [0.5, 1.0] {
        let result =
            advanced::fem_heat_adaptive(&space, &problem, &[1.0; 5], theta, &ctl, None).unwrap();
        assert!((result.states.last().unwrap()[2] - (-0.1_f64).exp()).abs() < 0.001);
        assert!(result.rejected > 0);
    }
}

#[test]
fn native_two_dimensional_dg_uses_tolerance_control_and_fem_source_uses_physical_times() {
    use rssn::sim::boundary;
    use rssn::sim::dg;
    use rssn::sim::fem_simplex;
    use rssn::sim::pde_api;
    let mesh = dg::DgMesh::rect([2, 2], [0.0, 0.0], [1.0, 1.0], 2);
    let (xi, w) = dg::lgl(2);
    let layout = dg::DgSolution {
        mesh,
        xi,
        w,
        q: vec![0.0; 36],
        nvar: 1,
        t: 0.0,
        steps: 0,
    };
    let initial: Vec<_> = (0..4)
        .flat_map(|e| {
            let layout = &layout;
            (0..9).map(move |j| {
                let x = layout.coord(e, j);
                (2.0 * std::f64::consts::PI * x[0]).sin()
                    * (2.0 * std::f64::consts::PI * x[1]).cos()
            })
        })
        .collect();
    let mut reference = initial.clone();
    for i in 0..100 {
        reference = dg::advection_step_nd(&mesh, [1.0, 0.3], &reference, i as f64 * 0.0005, 0.0005)
            .unwrap();
    }
    let mut ctl = options(0.05);
    ctl.rtol = 1e-8;
    let result =
        pde_api::adaptive::dg_advection_nd_adaptive(&mesh, [1.0, 0.3], &initial, &ctl, None)
            .unwrap();
    assert!(result.rejected > 0);
    let error = result
        .states
        .last()
        .unwrap()
        .iter()
        .zip(reference)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0_f64, f64::max);
    assert!(error < 2e-5, "DG temporal error {error}");
    let space =
        fem_simplex::FemSpace::new(fem_simplex::SimplexMesh::interval(4, 0.0, 1.0), 1).unwrap();
    let bcs = boundary::MarkerBcs::new()
        .set(
            fem_simplex::marker::XMIN,
            boundary::BoundaryCondition::neumann(0.0),
        )
        .set(
            fem_simplex::marker::XMAX,
            boundary::BoundaryCondition::neumann(0.0),
        );
    let zero = |_: f64, _: [f64; 3]| 0.0;
    let problem = fem_simplex::FemProblem {
        k: 1.0,
        c: 0.0,
        f: &zero,
        bcs: &bcs,
    };
    let mut ctl = options(0.2);
    ctl.t0 = 2.0;
    ctl.t1 = 2.2;
    let result = pde_api::adaptive::fem_heat_source_adaptive(
        &space,
        &problem,
        &[0.0; 5],
        0.5,
        |t, coords| Ok(vec![1.0 + t * t; coords.len()]),
        &ctl,
        None,
    )
    .unwrap();
    let exact = 0.2 + (2.2_f64.powi(3) - 2.0_f64.powi(3)) / 3.0;
    assert!((result.states.last().unwrap()[2] - exact).abs() < 0.0001);
    assert!(result.rejected > 0);
}
