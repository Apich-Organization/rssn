//! The native grid borrows its dimensions from the validated session object.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

fn call(
    session: &Session,
    name: &str,
    values: &[Value],
) -> Result<Value, OpError> {
    session.call_value(name, values, &CallCtx::new())
}

fn measured(
    value: f64,
    dim: &DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim: dim.clone(),
    })
}

fn tensor(
    shape: &[usize],
    values: &[f64],
    dim: &DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        elem: TypeExpr::Float.with_dim(dim.clone()),
        data: values.iter().map(|&value| measured(value, dim)).collect(),
    })
}

#[test]
fn grid_checks_shape_dimensions_and_session_ownership() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let mass = DimExpr::base("Mass");
    let xs = tensor(&[3], &[0.0, 2.0, 3.0], &length);
    let ys = tensor(&[2], &[1.0, 5.0], &time);
    let cells = [0.0, 4.0, 6.0, 8.0, 12.0, 14.0];
    let prepare = "num_interp_bilinear_grid_prepare";
    let evaluate = "num_interp_bilinear_grid_evaluate";
    let grid = call(
        &session,
        prepare,
        &[xs.clone(), ys.clone(), tensor(&[2, 3], &cells, &mass)],
    )
    .unwrap();
    assert!(
        call(
            &session,
            prepare,
            &[xs.clone(), ys.clone(), tensor(&[3, 2], &cells, &mass)]
        )
        .is_err()
    );
    let duplicates = tensor(&[3], &[0.0, 0.0, 3.0], &length);
    assert!(
        call(
            &session,
            prepare,
            &[duplicates, ys, tensor(&[2, 3], &cells, &mass)]
        )
        .is_err()
    );

    let sample = |x, y, x_dim: &DimExpr, cell_dim: &DimExpr| {
        [
            grid.clone(),
            measured(x, x_dim),
            measured(y, &time),
            measured(1.0, cell_dim),
        ]
    };
    let query = sample(1.0, 3.0, &length, &mass);
    assert_eq!(call(&session, evaluate, &query), Ok(measured(6.0, &mass)));
    assert_eq!(
        call(&session, evaluate, &sample(-10.0, 7.0, &length, &mass)),
        Ok(measured(8.0, &mass))
    );
    assert!(call(&session, evaluate, &sample(1.0, 3.0, &time, &mass)).is_err());
    assert!(call(&session, evaluate, &sample(1.0, 3.0, &length, &length)).is_err());
    assert!(call(&session, evaluate, &sample(f64::NAN, 3.0, &length, &mass)).is_err());
    assert!(call(&Session::new(), evaluate, &query).is_err());
    assert_eq!(
        call(&session, "num_native_release", &[grid]),
        Ok(Value::Bool(true))
    );
    assert!(call(&session, evaluate, &query).is_err());
}
