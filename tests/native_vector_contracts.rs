//! ABI regressions for typed vector arithmetic and fixed-block reductions.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

fn measured(
    x: f64,
    dimension: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: dimension,
    })
}

fn vector(
    xs: &[f64],
    dimension: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![xs.len()],
        data: xs.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float.with_dim(dimension),
    })
}

fn assert_tensor(
    value: &Value,
    shape: &[usize],
    expected: &[f64],
    dimension: DimExpr,
) {
    let Value::Tensor(result) = value else {
        panic!("expected tensor, got {value:?}");
    };
    assert_eq!(result.shape, shape);
    assert!(matches_type(&TypeExpr::Float.with_dim(dimension), &result.elem).accepted());
    let numbers: Vec<_> = result
        .data
        .iter()
        .map(|cell| {
            match cell {
                | Value::Dyn(d) => d.value.as_f64().unwrap(),
                | other => other.as_f64().unwrap(),
            }
        })
        .collect();
    assert_eq!(numbers, expected);
}

#[test]
fn measured_vector_updates_keep_operands_and_apply_units_before_dispatch() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let mass = DimExpr::base("Mass");
    let velocity = DimExpr::Mul(vec![length.clone(), time.clone().pow(-1)]);
    let scale = measured(2.0, velocity);
    let times = vector(&[1.0, 3.0], time);
    let positions = vector(&[5.0, 7.0], length.clone());
    let updated = session
        .call_value(
            "num_vector_axpy",
            &[scale.clone(), times.clone(), positions.clone()],
            &ctx,
        )
        .unwrap();
    assert_tensor(&updated, &[2], &[7.0, 13.0], length.clone());
    assert_tensor(&positions, &[2], &[5.0, 7.0], length.clone());

    let masses = vector(&[5.0, 7.0], mass.clone());
    let meters_per_kg = DimExpr::Mul(vec![length.clone(), mass.pow(-1)]);
    let weighted = session
        .call_value(
            "num_vector_axpby",
            &[
                scale.clone(),
                times.clone(),
                measured(0.5, meters_per_kg),
                masses.clone(),
            ],
            &ctx,
        )
        .unwrap();
    assert_tensor(&weighted, &[2], &[4.5, 9.5], length);
    assert!(
        session
            .call_value("num_vector_axpy", &[scale, times, masses], &ctx)
            .is_err()
    );
}

#[test]
fn vector_products_quotients_and_fixed_scans_have_explicit_dimensions() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let product = DimExpr::Mul(vec![length.clone(), time.clone()]);
    let quotient_units = DimExpr::Mul(vec![length.clone(), time.clone().pow(-1)]);
    let positions = vector(&[2.0, 8.0], length.clone());
    let times = vector(&[1.0, 4.0], time.clone());
    let args = [positions.clone(), times];
    let multiply = session
        .call_value("num_vector_hadamard_mul", &args, &ctx)
        .unwrap();
    assert_tensor(&multiply, &[2], &[2.0, 32.0], product.clone());
    let quotient = session
        .call_value("num_vector_hadamard_div", &args, &ctx)
        .unwrap();
    assert_tensor(&quotient, &[2], &[2.0, 2.0], quotient_units);
    let outer = session
        .call_value("num_vector_outer_product", &args, &ctx)
        .unwrap();
    assert_tensor(&outer, &[2, 2], &[2.0, 8.0, 8.0, 32.0], product);
    for grain in [0, usize::MAX] {
        let prefix = session
            .call_value(
                "num_par_prefix_sum_f64_fixed",
                &[positions.clone(), Value::int(1), Value::int(grain)],
                &ctx,
            )
            .unwrap();
        assert_tensor(&prefix, &[2], &[2.0, 10.0], length.clone());
    }
    let singular = [positions, vector(&[1.0, 0.0], time)];
    let err = session
        .call_value("num_vector_hadamard_div", &singular, &ctx)
        .unwrap_err();
    assert!(matches!(err, OpError::Domain(_)));
}
