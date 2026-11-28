use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::optim::LpStatus;
use rssn::kernels::optim::Pricing;
use rssn::kernels::optim::Relation;
use rssn::kernels::optim::simplex_with_max_iter;

fn measured(
    value: Value,
    dim: &DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(value),
        dim: dim.clone(),
    })
}

fn numbers(
    shape: Vec<usize>,
    values: &[f64],
) -> Value {
    Value::Tensor(TensorValue::of(
        shape,
        values.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float,
    ))
}

fn measured_numbers(
    shape: Vec<usize>,
    values: &[f64],
    dim: &DimExpr,
) -> Value {
    Value::Tensor(TensorValue::of(
        shape,
        values
            .iter()
            .map(|&value| measured(Value::Float(value), dim))
            .collect(),
        TypeExpr::Dim(Box::new(TypeExpr::Float), dim.clone()),
    ))
}

fn scalar(value: &Value) -> f64 {
    match value {
        | Value::Float(x) => *x,
        | Value::Dyn(x) => scalar(&x.value),
        | _ => panic!("expected a floating scalar, got {value:?}"),
    }
}

fn record_field<'a>(
    record: &'a RecordValue,
    name: &str,
) -> &'a Value {
    &record
        .fields
        .iter()
        .find(|(field, _)| field == name)
        .unwrap()
        .1
}

#[test]
fn pure_cleanup_uses_the_exact_pivot_budget() {
    let rows = [vec![-1.0]];
    let relation = [Relation::Eq];
    let rejected = simplex_with_max_iter(&[0.0], &rows, &relation, &[0.0], Pricing::Bland, 0);
    assert_eq!(rejected.unwrap_err(), "iteration limit exceeded");
    let accepted =
        simplex_with_max_iter(&[0.0], &rows, &relation, &[0.0], Pricing::Bland, 1).unwrap();
    assert_eq!(accepted.iterations, 1);
    assert!(matches!(accepted.status, LpStatus::Optimal { .. }));
}

#[test]
fn tableau_reports_unbounded_without_consuming_a_pivot_and_rejects_bad_data() {
    let r = simplex_with_max_iter(&[-1.0], &[], &[], &[], Pricing::SteepestEdge, 0).unwrap();
    assert!(matches!(r.status, LpStatus::Unbounded));
    assert_eq!(r.iterations, 0);
    assert!(simplex_with_max_iter(&[1.0], &[vec![0.0]], &[], &[], Pricing::Bland, 5).is_err());
    assert!(simplex_with_max_iter(&[f64::NAN], &[], &[], &[], Pricing::Bland, 5).is_err());
}

#[test]
fn measured_bounds_and_outputs_use_independent_reference_dimensions() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let length = DimExpr::Base("Length".into());
    let mass = DimExpr::Base("Mass".into());
    let constraint = DimExpr::Base("Time".into());
    let objective_coefficient = DimExpr::Mul(vec![
        mass.clone(),
        DimExpr::Inverse(Box::new(length.clone())),
    ]);
    let row_coefficient = DimExpr::Mul(vec![
        constraint.clone(),
        DimExpr::Inverse(Box::new(length.clone())),
    ]);
    let mut args = [
        measured_numbers(vec![2], &[-1.0, -2.0], &objective_coefficient),
        measured_numbers(vec![1, 2], &[1.0, 1.0], &row_coefficient),
        Value::Tensor(TensorValue::of(
            vec![1],
            vec![Value::String("Le".into())],
            TypeExpr::String,
        )),
        measured_numbers(vec![1], &[2.0], &constraint),
        measured(Value::Float(2.0), &length),
        measured(Value::Float(5.0), &mass),
        measured(Value::Float(4.0), &constraint),
        measured_numbers(vec![2], &[-2.0, 0.0], &length),
        measured_numbers(vec![2], &[1.0, 3.0], &length),
        Value::Tensor(TensorValue::of(
            vec![2],
            vec![Value::Bool(true); 2],
            TypeExpr::Bool,
        )),
        Value::int(20),
    ];
    let result = session
        .call_value("num_convex_solve_lp", &args, &ctx)
        .unwrap();
    let Value::Record(output) = result else {
        panic!("expected the LP result record")
    };
    assert_eq!(
        record_field(&output, "status"),
        &Value::String("Optimal".into())
    );
    let Value::Tensor(coordinates) = record_field(&output, "x") else {
        panic!("expected measured coordinates, got {output:?}")
    };
    for (coordinate, expected) in coordinates.data.iter().zip([-1.0, 3.0]) {
        assert!(matches!(coordinate, Value::Dyn(v) if v.dim == length));
        assert!((scalar(coordinate) - expected).abs() < 1e-12);
    }
    let Value::Variant(value) = record_field(&output, "value") else {
        panic!("expected an optional objective")
    };
    assert!(matches!(&value.payload[..], [Value::Dyn(v)] if v.dim == mass));
    assert!((scalar(&value.payload[0]) + 5.0).abs() < 1e-12);
    args[7] = measured_numbers(vec![2], &[-2.0, 0.0], &constraint);
    assert!(
        session
            .call_value("num_convex_solve_lp", &args, &ctx)
            .is_err()
    );
}

#[test]
fn budgeted_operator_calls_through_the_checked_budget() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let args = [
        numbers(vec![1], &[0.0]),
        numbers(vec![1, 1], &[-1.0]),
        Value::Tensor(TensorValue::of(
            vec![1],
            vec![Value::String("Eq".into())],
            TypeExpr::String,
        )),
        numbers(vec![1], &[0.0]),
        Value::String("Bland".into()),
        Value::Float(1.0),
        Value::Float(1.0),
        Value::Float(1.0),
        Value::int(0),
    ];
    assert!(matches!(
        session.call_value("num_analysis_simplex_budgeted", &args, &ctx),
        Err(OpError::BudgetExceeded)
    ));
}
