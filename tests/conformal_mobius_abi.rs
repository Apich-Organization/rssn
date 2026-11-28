//! Projective scaling and typed composition of Möbius operations.
use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;

use num_complex::Complex64;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::value::ComplexValue;
use rssn::kernels::complex::Mobius;
use rssn::kernels::conformal::MobiusClass;
use rssn::kernels::conformal::classify_mobius;
use rssn::kernels::conformal::mobius_fixed_points;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

fn value(z: Complex64) -> Value {
    Value::Complex(ComplexValue::Float { re: z.re, im: z.im })
}

#[test]
fn projective_rescaling_does_not_change_class_or_fixed_points() {
    for scale in [1e-200, 1.0, 1e200] {
        let scaled = scale * (ONE + I);
        let rotation = Mobius::new(scaled * I, ZERO, ZERO, scaled);
        assert_eq!(classify_mobius(&rotation).unwrap(), MobiusClass::Elliptic);
        assert_eq!(mobius_fixed_points(&rotation), [ZERO]);
    }
    let map = Mobius::new(2.0 * ONE, ZERO, 1e-200 * ONE, ONE);
    let points = mobius_fixed_points(&map);
    assert_eq!(points.len(), 2);
    assert!((points[0].re / 1e200 - 1.0).abs() < 2e-15);
    assert_eq!(points[1], ZERO);
    let translation = Mobius::new(ONE, 1e-200 * ONE, ZERO, ONE);
    assert_eq!(
        classify_mobius(&translation).unwrap(),
        MobiusClass::Parabolic
    );
}

#[test]
fn typed_automorphisms_compose_with_class_and_fixed_points() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let mk = |angle| {
        session
            .call_value(
                "num_conformal_disk_automorphism",
                &[value(ZERO), Value::Float(angle)],
                &ctx,
            )
            .unwrap()
    };
    let identity = mk(0.0);
    let rotation = mk(FRAC_PI_2);
    let class = "num_conformal_classify_mobius";
    assert_eq!(
        session
            .call_value(class, std::slice::from_ref(&rotation), &ctx)
            .unwrap(),
        Value::String("elliptic".into())
    );
    assert_eq!(
        session
            .call_value(class, std::slice::from_ref(&identity), &ctx)
            .unwrap(),
        Value::String("identity".into())
    );
    let fixed = "num_conformal_mobius_fixed_points";
    for (map, all) in [(identity, true), (rotation, false)] {
        let Value::Record(result) = session.call_value(fixed, &[map], &ctx).unwrap() else {
            panic!("fixed-point record expected")
        };
        assert_eq!(result.get("all"), Some(&Value::Bool(all)));
        let Some(Value::Tensor(points)) = result.get("points") else {
            panic!("fixed-point vector expected")
        };
        if all {
            assert!(points.data.is_empty());
        } else {
            assert_eq!(points.data, [value(ZERO)]);
        }
    }
    let bergman = "num_conformal_bergman_kernel_disk";
    let got = session
        .call_value(bergman, &[value(ZERO), value(ZERO)], &ctx)
        .unwrap();
    let Value::Complex(got) = got else {
        panic!("complex Bergman kernel expected")
    };
    let (re, im) = got.to_f64();
    assert!((re - 1.0 / PI).abs() < 1e-14 && im == 0.0);
    assert!(
        session
            .call_value(bergman, &[value(ONE), value(ZERO)], &ctx)
            .is_err()
    );
}
