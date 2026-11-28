use std::sync::Arc;

use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;

fn complex(
    re: f64,
    dim: &DimExpr,
) -> Value {
    measured(Value::Complex(ComplexValue::Float { re, im: 0.0 }), dim)
}

fn measured(
    value: Value,
    dim: &DimExpr,
) -> Value {
    if *dim == DimExpr::Dimless {
        value
    } else {
        Value::Dyn(DynValue {
            value: Box::new(value),
            dim: dim.clone(),
        })
    }
}

fn complex_ty(dim: &DimExpr) -> TypeExpr {
    if *dim == DimExpr::Dimless {
        TypeExpr::Complex
    } else {
        TypeExpr::Dim(Box::new(TypeExpr::Complex), dim.clone())
    }
}

fn field<'a>(
    record: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = record else {
        panic!("expected record")
    };
    record.get(name).unwrap()
}

fn raw(value: &Value) -> &Value {
    if let Value::Dyn(value) = value {
        &value.value
    } else {
        value
    }
}

fn real(value: &Value) -> f64 {
    let Value::Complex(value) = raw(value) else {
        panic!("expected complex coordinate")
    };
    value.to_f64().0
}

#[test]
fn newton_preserves_measured_root_and_work_budget() {
    let session = rssn::api::Session::new();
    let length = DimExpr::Base("Length".into());
    let time = DimExpr::Base("Time".into());
    let rate = DimExpr::Mul(vec![
        time.clone(),
        DimExpr::Inverse(Box::new(length.clone())),
    ]);
    let mut ctx = CallCtx::new();
    let time_output = time.clone();
    let rate_output = rate.clone();
    let function = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |a: &[Value]| Ok(complex(real(&a[0]).powi(2) - 1.0, &time_output)),
            Signature::new(vec![complex_ty(&length)], complex_ty(&time)),
        )))
        .unwrap();
    let derivative = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |a: &[Value]| Ok(complex(2.0 * real(&a[0]), &rate_output)),
            Signature::new(vec![complex_ty(&length)], complex_ty(&rate)),
        )))
        .unwrap();
    let call = |seed, budget| {
        session.call_ho(
            "num_complex_dynamics_newton_basin_checked",
            &[
                complex(seed, &length),
                measured(Value::Float(2.0), &length),
                measured(Value::Float(4.0), &time),
                measured(Value::Float(1e-12), &time),
                Value::int(32),
                Value::int(budget),
            ],
            &[function, derivative],
            &ctx,
        )
    };
    let result = call(2.0, 64).unwrap();
    let Value::Variant(root) = field(&result, "root") else {
        panic!("expected optional root")
    };
    assert_eq!(root.tag.as_str(), "Some");
    let point = field(&root.payload[0], "point");
    assert!(
        rssn::kernel::overload::matches_type(&complex_ty(&length), &point.type_expr()).accepted()
    );
    assert!((real(point) - 1.0).abs() < 1e-10);
    let singular = call(0.0, 2).unwrap();
    let Value::Variant(root) = field(&singular, "root") else {
        panic!("expected optional root")
    };
    assert_eq!(root.tag.as_str(), "None");
    assert_eq!(field(&singular, "tested").as_f64(), Some(1.0));
    assert_eq!(field(&singular, "evaluations").as_f64(), Some(2.0));
    assert!(matches!(call(2.0, 1), Err(OpError::BudgetExceeded)));
}

#[test]
fn koenigs_linear_identity_and_fixed_point_validation() {
    let session = rssn::api::Session::new();
    let dim = DimExpr::Dimless;
    let mut ctx = CallCtx::new();
    let map = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |a: &[Value]| Ok(complex(0.5 * real(&a[0]), &DimExpr::Dimless)),
            Signature::new(vec![TypeExpr::Complex], TypeExpr::Complex),
        )))
        .unwrap();
    let call = |fixed, budget, callback| {
        session.call_ho(
            "num_complex_dynamics_koenigs_checked",
            &[
                complex(fixed, &dim),
                complex(0.5, &dim),
                complex(0.25, &dim),
                Value::Float(1.0),
                Value::Float(1e-12),
                Value::int(12),
                Value::int(budget),
            ],
            &[callback],
            &ctx,
        )
    };
    let result = call(0.0, 2, map).unwrap();
    assert!((real(field(&result, "coordinate")) - 0.25).abs() < 1e-14);
    assert_eq!(field(&result, "fixed_point_error").as_f64(), Some(0.0));
    assert_eq!(field(&result, "steps").as_f64(), Some(1.0));
    assert_eq!(field(&result, "evaluations").as_f64(), Some(2.0));
    assert!(
        call(0.25, 2, map)
            .unwrap_err()
            .to_string()
            .contains("fixed point")
    );
    assert!(matches!(call(0.0, 1, map), Err(OpError::BudgetExceeded)));
}
