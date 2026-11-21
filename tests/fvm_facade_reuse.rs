//! Fixed and adaptive scalar finite-volume facades share the same step.

use std::cell::Cell;

use rssn::sim::fdm_nd::FdmError;
use rssn::sim::fdm_nd::adaptive::PdeAdaptiveOptions;
use rssn::sim::pde_api::FvmScheme;
use rssn::sim::pde_api::adaptive::fvm_adaptive;
use rssn::sim::pde_api::solve_fvm_conservation_1d;

#[test]
fn fixed_fvm_preserves_flux_visit_counts_and_periodic_mass() {
    let cells = 64;
    let steps = 11;
    let initial = |x: f64| 2.0 + (std::f64::consts::TAU * x).sin();
    for (scheme, flux_calls_per_cell) in [
        (FvmScheme::Upwind, 4),
        (FvmScheme::LaxFriedrichs, 4),
        (FvmScheme::MusclMinmod, 6),
        (FvmScheme::MusclVanLeer, 6),
    ] {
        let flux_calls = Cell::new(0);
        let speed_calls = Cell::new(0);
        let (x, result) = solve_fvm_conservation_1d(
            cells,
            1.0,
            0.003,
            steps,
            true,
            scheme,
            initial,
            |u| {
                flux_calls.set(flux_calls.get() + 1);
                u
            },
            |_| {
                speed_calls.set(speed_calls.get() + 1);
                1.0
            },
        )
        .unwrap();
        assert_eq!(flux_calls.get(), flux_calls_per_cell * steps * cells);
        assert_eq!(speed_calls.get(), 4 * steps * cells);
        assert!(
            (result.iter().sum::<f64>() - x.into_iter().map(initial).sum::<f64>()).abs() < 1e-10
        );
    }
}

#[test]
fn the_facades_reject_a_cell_width_rounded_to_zero() {
    let length = f64::from_bits(1);
    let fixed = solve_fvm_conservation_1d(
        3,
        length,
        1.0,
        0,
        false,
        FvmScheme::Upwind,
        |_| panic!("bad grid invoked initial condition"),
        |u| u,
        |_| 1.0,
    );
    assert!(matches!(fixed, Err(FdmError::Domain(_))));
    let adaptive = fvm_adaptive(
        length,
        &[1.0, 2.0, 3.0],
        false,
        FvmScheme::Upwind,
        |_| panic!("bad grid invoked flux"),
        |_| 1.0,
        &PdeAdaptiveOptions::default(),
        None,
    );
    assert!(matches!(adaptive, Err(FdmError::Domain(_))));
}
