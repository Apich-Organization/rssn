//! Oriented generators, stable fiber tracking, and branch-lip extrapolation.
use num_complex::Complex64;
use rssn::kernels::complex_coords::algebraic_riemann_surface_superelliptic;
use rssn::kernels::complex_coords::branch_cut_jump_numeric;
use rssn::kernels::complex_coords::monodromy_around_loop;

#[test]
fn positive_local_infinity_generator_inverts_finite_cycles() {
    let surface = algebraic_riemann_surface_superelliptic(&[Complex64::ZERO, Complex64::ONE], 3);
    assert_eq!(surface.genus, 1);
    assert!(surface.has_branch_at_infinity);
    for sheet in 0..surface.degree {
        let end = surface.monodromy_generators.iter().fold(sheet, |i, p| p[i]);
        assert_eq!(end, sheet);
    }
}

#[test]
fn square_root_fibers_swap_after_winding_and_reject_ramification_path() {
    let zero = Complex64::ZERO;
    let one = Complex64::ONE;
    let coefficients = vec![vec![zero, -one], vec![zero, zero], vec![one, zero]];
    let path = [one, Complex64::I, -one, -Complex64::I];
    assert_eq!(
        monodromy_around_loop(&coefficients, &path),
        Some(vec![1, 0])
    );
    assert_eq!(monodromy_around_loop(&coefficients, &[one, -one]), None);
    assert_eq!(
        monodromy_around_loop(&coefficients, &[one, one * 2.0, one + Complex64::I]),
        Some(vec![0, 1])
    );
}

#[test]
fn logarithm_lip_jump_cancels_the_linear_normal_error() {
    let jump = branch_cut_jump_numeric(Complex64::ln, -Complex64::ONE);
    assert!(jump.re.abs() < 1e-14);
    assert!((jump.im - std::f64::consts::TAU).abs() < 1e-14);
}
