//! `kernels::pde` currently exposes only the placeholder `pde_solver()`, which
//! prints a note and returns `()`. There is no numerical result to check, so the
//! tests only pin down the (deliberately trivial) contract: it is callable,
//! side-effect free apart from printing, repeatable and does not panic.

use rssn::kernels::pde::pde_solver;

#[test]
fn placeholder_returns_unit_and_does_not_panic() {
    let () = pde_solver();
}

#[test]
fn placeholder_is_repeatable() {
    for _ in 0..3 {
        pde_solver();
    }
}

#[test]
fn placeholder_is_callable_from_several_threads() {
    let handles: Vec<_> = (0..4).map(|_| std::thread::spawn(pde_solver)).collect();
    for h in handles {
        assert!(h.join().is_ok());
    }
}
