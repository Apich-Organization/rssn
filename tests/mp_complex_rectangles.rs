//! Checked complex rectangle operations through the numerical ABI.

use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::DimRuleInfo;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;

fn record(fields: impl IntoIterator<Item = (&'static str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .into_iter()
            .map(|(name, value)| (name.into(), value))
            .collect(),
    })
}

fn dyadic(
    mantissa: i64,
    exponent: i64,
    precision: u32,
) -> Value {
    record([
        ("mantissa", Value::int(mantissa)),
        ("exponent", Value::int(exponent)),
        ("precision", Value::int(precision)),
    ])
}

fn ball(
    midpoint: Value,
    radius: Value,
) -> Value {
    record([("midpoint", midpoint), ("radius", radius)])
}

fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected a record");
    };
    record.get(name).unwrap()
}

fn exact(value: &Value) -> BigRational {
    let Value::Int(mantissa) = field(value, "mantissa") else {
        panic!("expected an integer mantissa");
    };
    let Value::Int(exponent) = field(value, "exponent") else {
        panic!("expected an integer exponent");
    };
    let exponent = i64::try_from(exponent).unwrap();
    let shift = exponent.unsigned_abs();
    if exponent < 0 {
        BigRational::new(mantissa.clone(), BigInt::from(1_u32) << shift)
    } else {
        BigRational::from_integer(mantissa << shift)
    }
}

fn boundaries(value: &Value) -> (BigRational, BigRational) {
    let mid = exact(field(value, "midpoint"));
    let rad = exact(field(value, "radius"));
    (&mid - &rad, mid + rad)
}

fn encloses(
    actual: &Value,
    expected: &BigRational,
) {
    let (low, high) = boundaries(actual);
    assert!(
        &low <= expected && expected <= &high,
        "{expected} outside [{low}, {high}]"
    );
}

#[test]
fn squared_magnitude_and_reciprocal_keep_a_separated_rectangle() {
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context).unwrap();
    let real = ball(dyadic(1, 1, 1), dyadic(1, 0, 1));
    let rectangle = call("num_mp_complex_from_real", &[real]);
    let magnitude = call(
        "num_mp_complex_abs_sq",
        &[rectangle.clone(), Value::int(64)],
    );
    assert_eq!(
        boundaries(&magnitude),
        (
            BigRational::from_integer(1.into()),
            BigRational::from_integer(9.into())
        )
    );
    for precision in [1, 4, 64] {
        let inverse = call(
            "num_mp_complex_recip",
            &[rectangle.clone(), Value::int(precision)],
        );
        encloses(
            field(&inverse, "real"),
            &BigRational::new(1.into(), 3.into()),
        );
        encloses(
            field(&inverse, "real"),
            &BigRational::from_integer(1.into()),
        );
        assert_eq!(
            boundaries(field(&inverse, "imag")),
            (BigRational::default(), BigRational::default())
        );
    }
}

#[test]
fn reciprocal_accepts_a_separated_interval_with_a_large_dynamic_range() {
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context).unwrap();
    let radius = (1_i64 << 49) - 1;
    let real = ball(dyadic(1, 49, 1), dyadic(radius, 0, 49));
    let rectangle = call("num_mp_complex_from_real", &[real]);
    let inverse = call("num_mp_complex_recip", &[rectangle, Value::int(4)]);
    encloses(
        field(&inverse, "real"),
        &BigRational::new(1.into(), ((1_i64 << 50) - 1).into()),
    );
    encloses(
        field(&inverse, "real"),
        &BigRational::from_integer(1.into()),
    );
}

#[test]
fn products_and_squares_enclose_exact_corners() {
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context).unwrap();
    let positive = ball(dyadic(1, 1, 1), dyadic(1, 0, 1));
    let crosses_zero = ball(dyadic(1, 0, 1), dyadic(1, 1, 1));
    let left = call(
        "num_mp_complex_new",
        &[positive.clone(), crosses_zero.clone()],
    );
    let right = call("num_mp_complex_new", &[crosses_zero, positive]);
    let integer = |n: i64| BigRational::from_integer(n.into());
    for precision in [1, 3, 7] {
        let product = call(
            "num_mp_complex_mul",
            &[left.clone(), right.clone(), Value::int(precision)],
        );
        let norm = call(
            "num_mp_complex_abs_sq",
            &[left.clone(), Value::int(precision)],
        );
        for a in [1, 3] {
            for b in [-1, 0, 3] {
                encloses(&norm, &integer(a * a + b * b));
                for c in [-1, 3] {
                    for d in [1, 3] {
                        encloses(field(&product, "real"), &integer(a * c - b * d));
                        encloses(field(&product, "imag"), &integer(a * d + b * c));
                    }
                }
            }
        }
    }
}

#[test]
fn exact_points_use_composable_algebra_and_canonical_constants() {
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context).unwrap();
    let one = call("num_mp_complex_one", &[Value::int(8)]);
    let i = call("num_mp_complex_i", &[Value::int(8)]);
    let zero = call("num_mp_complex_zero", &[Value::int(8)]);
    let point = call("num_mp_complex_point", &[dyadic(3, 0, 2), dyadic(1, 2, 1)]);
    let squared = call(
        "num_mp_complex_powi",
        &[i.clone(), Value::int(2), Value::int(8)],
    );
    assert_eq!(
        squared,
        call("num_mp_complex_neg", std::slice::from_ref(&one))
    );
    assert_eq!(
        call(
            "num_mp_complex_div",
            &[point.clone(), point.clone(), Value::int(8)]
        ),
        one,
    );
    let conjugate = call("num_mp_complex_conj", std::slice::from_ref(&i));
    assert_eq!(
        call("num_mp_complex_add", &[i.clone(), conjugate, Value::int(8)]),
        zero
    );
    assert_eq!(
        call("num_mp_complex_sub", &[i.clone(), i, Value::int(8)]),
        zero
    );
}

#[test]
fn rejects_zero_overlap_malformed_nested_intervals_and_unbounded_precision() {
    let session = Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context);
    let zero = dyadic(0, 0, 1);
    let measured = Value::Dyn(DynValue {
        value: Box::new(dyadic(1, 0, 1)),
        dim: DimExpr::base("Length"),
    });
    let singular = call(
        "num_mp_complex_from_real",
        &[ball(zero.clone(), dyadic(1, 0, 1))],
    )
    .unwrap();
    for (name, args) in [
        (
            "num_mp_complex_recip",
            vec![singular.clone(), Value::int(16)],
        ),
        (
            "num_mp_complex_powi",
            vec![singular, Value::int(-1), Value::int(16)],
        ),
        (
            "num_mp_complex_from_real",
            vec![ball(zero.clone(), dyadic(-1, 0, 1))],
        ),
        ("num_mp_complex_zero", vec![Value::int(32_769)]),
    ] {
        assert!(
            matches!(call(name, &args), Err(OpError::Domain(_))),
            "{name}: {args:?}"
        );
    }
    let measured = ball(measured, zero);
    assert!(call("num_mp_complex_from_real", &[measured]).is_err());
}

#[test]
fn rectangle_manifest_uses_declared_nested_dimensionless_signatures() {
    let manifest = Session::new().manifest();
    let dyadic = TypeExpr::Record(
        ["mantissa", "exponent", "precision"]
            .map(|name| (name.into(), TypeExpr::Int))
            .into(),
    );
    let ball = TypeExpr::Record(vec![
        ("midpoint".into(), dyadic.clone()),
        ("radius".into(), dyadic.clone()),
    ]);
    let rectangle = TypeExpr::Record(vec![
        ("real".into(), ball.clone()),
        ("imag".into(), ball.clone()),
    ]);
    for (name, args, output) in [
        (
            "point",
            vec![dyadic.clone(), dyadic.clone()],
            rectangle.clone(),
        ),
        (
            "div",
            vec![rectangle.clone(), rectangle.clone(), TypeExpr::Int],
            rectangle.clone(),
        ),
        ("abs_sq", vec![rectangle.clone(), TypeExpr::Int], ball),
        (
            "powi",
            vec![rectangle.clone(), TypeExpr::Int, TypeExpr::Int],
            rectangle,
        ),
    ] {
        let op = manifest.op(&format!("num_mp_complex_{name}")).unwrap();
        let [signature] = op.overloads.as_slice() else {
            panic!("expected a declared signature for {name}");
        };
        assert_eq!(op.signature_source, SignatureSource::Declared);
        assert_eq!(op.dim_rule, Some(DimRuleInfo::Dimensionless));
        assert_eq!(
            op.rust_path.as_deref(),
            Some(format!("crate::kernels::mp::complex::MpComplex::{name}").as_str())
        );
        assert_eq!(signature.params, args);
        assert_eq!(signature.ret, output);
    }
}
