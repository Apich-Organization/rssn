//! Accuracy and portable-layout regressions for conformal kernels.
use num_complex::Complex64;
use rssn::kernels::conformal::SchwarzChristoffelMap;
use rssn::kernels::conformal::bergman_kernel_polynomial_basis;
use rssn::kernels::conformal::mobius_from_three_points;

#[test]
fn projective_fit_interpolates_nonaffine_triples_in_different_charts() {
    let one = Complex64::from(1.0);
    let two = Complex64::from(2.0);
    for source_scale in [1e-120, 1.0, 1e120] {
        for target_scale in [1e-120, 1.0, 1e120] {
            let f = mobius_from_three_points(
                Complex64::ZERO,
                one * source_scale,
                two * source_scale,
                Complex64::ZERO,
                0.5 * one * target_scale,
                (2.0 / 3.0) * one * target_scale,
            )
            .unwrap();
            let z = 3.0 * source_scale * one;
            let numerator = (f.a * z + f.b) / target_scale;
            let denominator = f.c * z + f.d;
            assert!((numerator.re / denominator.re - 0.75).abs() < 5e-14);
        }
    }
    assert!(mobius_from_three_points(one, one, two, Complex64::ZERO, one, two).is_err());
}

#[test]
fn bergman_basis_has_hermitian_symmetry_and_nonnegative_diagonal() {
    let basis = [
        vec![Complex64::from(1.0)],
        vec![Complex64::ZERO, Complex64::from(0.5)],
    ];
    let z = Complex64::new(0.3, 0.2);
    let w = Complex64::new(-0.25, 0.375);
    let value = bergman_kernel_polynomial_basis(&basis, z, w);
    assert!((value - (Complex64::from(1.0) + 0.25 * z * w.conj())).norm() < 1e-14);
    assert!((value - bergman_kernel_polynomial_basis(&basis, w, z).conj()).norm() < 1e-14);
    let diagonal = bergman_kernel_polynomial_basis(&basis, z, z);
    assert!((diagonal - Complex64::from(1.0 + 0.25 * z.norm_sqr())).norm() < 1e-14);
}

#[test]
fn schwarz_christoffel_identity_special_case_has_correct_anchor_and_boundary() {
    let map = SchwarzChristoffelMap {
        prevertices: vec![-1.0, 1.0],
        interior_angles_pi: vec![1.0, 1.0],
        scale: Complex64::new(2.0, -3.0),
        shift: Complex64::new(0.25, 0.125),
    };
    let z = Complex64::new(0.75, 0.25);
    assert!((map.eval(z).unwrap() - (map.shift + map.scale * (z - Complex64::I))).norm() < 1e-13);
    assert_eq!(map.eval(Complex64::I).unwrap(), map.shift);
    assert!(map.eval(Complex64::from(0.25)).is_err());
    let invalid = SchwarzChristoffelMap {
        prevertices: vec![1.0, -1.0],
        ..map
    };
    assert!(invalid.eval(z).is_err());
}

#[test]
fn fitted_square_retains_interior_evaluation() {
    use rssn::kernels::conformal::solve_sc_uhp;
    let square = [
        Complex64::ZERO,
        Complex64::ONE,
        Complex64::new(1.0, 1.0),
        Complex64::I,
    ];
    let map = solve_sc_uhp(&square, 40, 1e-7).unwrap();
    let value = map.eval(Complex64::new(0.5, 0.5)).unwrap();
    assert!(value.re > 0.0 && value.re < 1.0 && value.im > 0.0 && value.im < 1.0);
}

#[test]
fn shared_complex_tensor_codec_keeps_grunsky_square_preflight() {
    use rssn::api::Session;
    use rssn::kernel::CallCtx;
    use rssn::kernel::Value;
    use rssn::kernel::signature::TypeExpr;
    use rssn::kernel::value::TensorValue;

    let tensor = |shape| {
        Value::Tensor(TensorValue {
            shape,
            data: vec![],
            elem: TypeExpr::Complex,
        })
    };
    let session = Session::new();
    let context = CallCtx::new();
    let args = |shape| [tensor(shape), tensor(vec![0]), Value::Float(0.0)];
    assert!(
        session
            .call_value("num_univalent_grunsky_check", &args(vec![0, 3]), &context)
            .is_err()
    );
    assert_eq!(
        session.call_value("num_univalent_grunsky_check", &args(vec![0, 0]), &context),
        Ok(Value::Bool(true)),
    );
}
