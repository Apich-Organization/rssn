//! Principal square roots through the measured numerical ABI.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;

#[test]
fn native_complex_roots_keep_small_components_and_halve_each_exponent() {
    let session = Session::new();
    let length = DimExpr::Base("Length".into());
    let area = length.clone().pow(Exp::Rat(Rational::int(2)));
    let root = |re, im| {
        let input = Value::Dyn(DynValue {
            dim: area.clone(),
            value: Box::new(Value::Complex(ComplexValue::Float { re, im })),
        });
        let result = session
            .call_value(
                "num_analysis_complex_scalar_method_sqrt",
                &[input],
                &CallCtx::new(),
            )
            .unwrap();
        let Value::Dyn(measured) = result else {
            panic!("expected a measured root");
        };
        assert_eq!(measured.dim, length);
        let Value::Complex(value) = *measured.value else {
            panic!("expected a complex root");
        };
        value.to_f64()
    };

    let (real, imaginary) = root(1e300, 1.0);
    assert!(real.to_bits().abs_diff(1e150_f64.to_bits()) <= 1);
    assert!(imaginary.to_bits().abs_diff(5e-151_f64.to_bits()) <= 1);

    let (real, imaginary) = root(-4.0, -0.0);
    assert_eq!(real.to_bits(), 0);
    assert_eq!(imaginary.to_bits(), (-2.0_f64).to_bits());

    let (real, imaginary) = root(f64::MAX, f64::MAX);
    assert!(real.is_finite() && real.is_sign_positive());
    assert!(imaginary.is_finite() && imaginary.is_sign_positive());
}
