//! Exact and high-dynamic-range checks for normalized univalent diagnostics.
use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::univalent::check_becker_univalence;
use rssn::kernels::univalent::check_bieberbach_de_branges;
use rssn::kernels::univalent::check_grunsky_inequality;
use rssn::kernels::univalent::grunsky_matrix;
use rssn::kernels::univalent::hankel_determinant;

#[test]
fn koebe_generates_the_identity_grunsky_matrix() {
    let order = 10;
    let series: Vec<_> = (0..2 * order + 2)
        .map(|n| Complex64::from(n as f64))
        .collect();
    let matrix = grunsky_matrix(&series, order).unwrap();
    for (i, row) in matrix.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            assert!((*value - Complex64::from(f64::from(i == j))).norm() < 1e-12);
        }
    }
    let test = vec![Complex64::new(0.5, -1.0); order];
    assert!(check_grunsky_inequality(&matrix, &test, 1e-12).unwrap());
    assert!(check_grunsky_inequality(&[vec![]], &[Complex64::from(1.0)], 0.0).is_err());
}

#[test]
fn coefficient_normalization_preserves_large_and_small_common_scales() {
    for scale in [1e-250, 1.0, 1e250] {
        let series = [0.0, 1.0, 2.0, 3.0].map(|y| Complex64::new(y * scale, y * scale));
        assert!(check_bieberbach_de_branges(&series, 1e-14));
        let matrix = grunsky_matrix(&series, 1).unwrap();
        assert!((matrix[0][0] - Complex64::from(1.0)).norm() < 1e-14);
        assert!(check_becker_univalence(
            Complex64::from(0.5),
            Complex64::new(scale, scale),
            Complex64::new(2.0 * scale, 2.0 * scale)
        ));
    }
    let x = Complex64::new(1e-20, 0.0);
    let determinant = hankel_determinant(&[x, x * 2.0, x * 3.0], 2, 0).unwrap();
    assert!(((determinant + x * x) / 1e-40).norm() < 1e-14);
    assert!(hankel_determinant(&[x], usize::MAX, usize::MAX).is_err());
}

#[test]
fn dispatch_checks_shapes_and_budgets_before_work() {
    let series = Value::Tensor(TensorValue {
        shape: vec![4],
        elem: TypeExpr::Complex,
        data: (0..4)
            .map(|v| {
                Value::Complex(ComplexValue::Float {
                    re: f64::from(v),
                    im: 0.0,
                })
            })
            .collect(),
    });
    let session = Session::new();
    assert_eq!(
        session.call_value(
            "num_univalent_grunsky_matrix",
            &[series.clone(), Value::int(93)],
            &CallCtx::new()
        ),
        Err(OpError::BudgetExceeded)
    );
    let result = session
        .call_value(
            "num_univalent_grunsky_matrix",
            &[series, Value::int(1)],
            &CallCtx::new(),
        )
        .unwrap();
    let Value::Tensor(tensor) = result else {
        panic!("expected a Grunsky tensor")
    };
    assert_eq!(tensor.shape, [1, 1]);
    assert_eq!(tensor.elem, TypeExpr::Complex);
    assert_eq!(
        tensor.data,
        [Value::Complex(ComplexValue::Float { re: 1.0, im: 0.0 })]
    );
}
