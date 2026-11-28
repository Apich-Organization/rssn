//! Checked native objective ownership, solver work and scoped root callbacks.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::plugin::NativeCallback;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use std::sync::Arc;
use std::sync::atomic::AtomicUsize;
use std::sync::atomic::Ordering;

fn named(
    session: &Session,
    name: &str,
    args: &[Value],
    ctx: &CallCtx<'_>,
) -> Result<Value, OpError> {
    session.call_value(name, args, ctx)
}

fn run(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    named(
        session,
        &format!("num_optimizer_{name}"),
        args,
        &CallCtx::new(),
    )
}

fn field<'a>(
    record: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = record else {
        panic!("expected a record")
    };
    &record
        .fields
        .iter()
        .find(|(field, _)| field == name)
        .unwrap()
        .1
}

fn tensor(
    shape: &[usize],
    elem: TypeExpr,
    values: Vec<Value>,
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        elem,
        data: values,
    })
}

fn vector(values: &[f64]) -> Value {
    tensor(
        &[values.len()],
        TypeExpr::Float,
        values.iter().copied().map(Value::Float).collect(),
    )
}

fn strings(values: &[&str]) -> Value {
    tensor(
        &[values.len()],
        TypeExpr::String,
        values.iter().map(|v| Value::String((*v).into())).collect(),
    )
}

fn measures(
    values: &[f64],
    shape: &[usize],
    dimension: DimExpr,
) -> Value {
    let data = values
        .iter()
        .map(|&v| measured(v, dimension.clone()))
        .collect();
    tensor(shape, TypeExpr::Float.with_dim(dimension), data)
}

fn measured(
    value: f64,
    dim: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim,
    })
}

fn record(fields: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: fields
            .into_iter()
            .map(|(key, value)| (key.into(), value))
            .collect(),
    })
}

fn options(evaluations: usize) -> Value {
    record(vec![
        ("iterations", Value::int(128)),
        ("evaluations", Value::int(evaluations)),
        ("tolerance", Value::Float(1e-14)),
        ("gradient_tolerance", Value::Float(1e-12)),
        ("seed", Value::int(7)),
    ])
}

fn root_options(
    solver: &str,
    lower: f64,
    upper: f64,
    evaluations: usize,
) -> Value {
    record(vec![
        ("solver", Value::String(solver.into())),
        ("lower", Value::Float(lower)),
        ("upper", Value::Float(upper)),
        ("tolerance", Value::Float(1e-12)),
        ("iterations", Value::int(128)),
        ("evaluations", Value::int(evaluations)),
    ])
}

#[test]
fn built_in_optima_and_analytic_gradients_agree() {
    let session = Session::new();
    for (name, point, optimum) in [
        ("rosenbrock", [1.0, 1.0], 0.0),
        ("sphere", [3.0, 4.0], 25.0),
        ("rastrigin", [0.0, 1.0], 1.0),
    ] {
        let coordinates = vector(&point);
        let plan = run(&session, name, &[Value::int(2)]).unwrap();
        let value = run(
            &session,
            &format!("{name}_cost"),
            std::slice::from_ref(&coordinates),
        )
        .unwrap();
        assert_eq!(value, Value::Float(optimum));
        let evaluation = run(&session, "value_and_gradient", &[plan.clone(), coordinates]).unwrap();
        assert_eq!(field(&evaluation, "value"), &value);
        assert_eq!(
            run(&session, "has_gradient", std::slice::from_ref(&plan)).unwrap(),
            Value::Bool(true)
        );
        assert!(run(&Session::new(), "has_gradient", &[plan]).is_err());
    }
    assert!(run(&session, "rosenbrock", &[Value::int(1)]).is_err());
}

#[test]
fn prepared_minimization_reuses_the_box_and_exhausts_budget() {
    let session = Session::new();
    let sphere = run(&session, "sphere", &[Value::int(2)]).unwrap();
    let bounded = run(
        &session,
        "with_bounds",
        &[sphere.clone(), vector(&[1.0, 1.0]), vector(&[4.0, 4.0])],
    )
    .unwrap();
    let method = Value::String("lbfgsb".into());
    assert_eq!(
        run(&session, "needs_bounds", std::slice::from_ref(&method)).unwrap(),
        Value::Bool(true)
    );
    assert_eq!(
        run(&session, "needs_gradient", std::slice::from_ref(&method)).unwrap(),
        Value::Bool(true)
    );
    let start = vector(&[2.0, 2.0]);
    let result = run(
        &session,
        "minimize",
        &[bounded, start.clone(), method.clone(), options(1024)],
    )
    .unwrap();
    assert_eq!(field(&result, "point"), &vector(&[1.0, 1.0]));
    assert_eq!(field(&result, "cost"), &Value::Float(2.0));
    assert!(
        run(
            &session,
            "minimize",
            &[sphere.clone(), start.clone(), method, options(1024)]
        )
        .is_err()
    );
    assert!(matches!(
        run(
            &session,
            "minimize",
            &[sphere, start, Value::String("bfgs".into()), options(1)]
        ),
        Err(OpError::BudgetExceeded)
    ));
}

#[test]
fn immutable_tape_can_be_released_after_preparing_a_function() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let tape = named(
        &session,
        "num_native_tape_lower",
        &[
            strings(&["(x-p)^2"]),
            strings(&["x"]),
            strings(&["p"]),
            vector(&[1.0]),
            vector(&[1.0]),
            vector(&[1.0]),
        ],
        &ctx,
    )
    .unwrap();
    let input = [
        tape.clone(),
        Value::int(0),
        vector(&[3.0]),
        Value::Bool(true),
    ];
    let a = run(&session, "problem_new", &input).unwrap();
    let b = run(&session, "problem_try_new", &input).unwrap();
    assert_eq!(
        named(&session, "num_native_release", &[tape], &ctx).unwrap(),
        Value::Bool(true)
    );
    for plan in [a, b] {
        let evaluation = run(
            &session,
            "value_and_gradient",
            &[plan.clone(), vector(&[1.0])],
        )
        .unwrap();
        assert_eq!(field(&evaluation, "value"), &Value::Float(4.0));
        assert_eq!(field(&evaluation, "gradient"), &vector(&[-4.0]));
        let best = run(
            &session,
            "minimize",
            &[
                plan,
                vector(&[1.0]),
                Value::String("bfgs".into()),
                options(256),
            ],
        )
        .unwrap();
        let Value::Tensor(point) = field(&best, "point") else {
            panic!("normalized best point")
        };
        assert!((point.data[0].as_f64().unwrap() - 3.0).abs() < 1e-7);
    }
}

#[test]
fn regression_residuals_preserve_response_units() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let lx = DimExpr::base("Length");
    let ty = DimExpr::base("Time");
    let slope = DimExpr::Mul(vec![ty.clone(), DimExpr::Inverse(Box::new(lx.clone()))]);
    let y_reference = measured(4.0, ty.clone());
    let model = named(
        &session,
        "num_analysis_linear_regression_new",
        &[
            measures(&[0.0, 1.0, 2.0], &[3, 1], lx.clone()),
            measures(&[0.0, 4.0, 8.0], &[3], ty.clone()),
            measured(2.0, lx),
            y_reference.clone(),
        ],
        &ctx,
    )
    .unwrap();
    assert_eq!(
        run(
            &session,
            "regression_coefficient_count",
            std::slice::from_ref(&model)
        )
        .unwrap(),
        Value::int(2)
    );
    let beta = tensor(
        &[2],
        TypeExpr::Float,
        vec![measured(1.0, ty.clone()), measured(2.0, slope)],
    );
    let output = run(
        &session,
        "regression_residuals",
        &[model.clone(), beta, y_reference.clone()],
    )
    .unwrap();
    assert_eq!(output, measures(&[1.0, -1.0, -3.0], &[3], ty));
    let plan = run(&session, "regression_problem", &[model]).unwrap();
    let v = run(&session, "value_and_gradient", &[plan, vector(&[0.0, 2.0])]).unwrap();
    assert_eq!(field(&v, "value"), &Value::Float(0.0));
}

#[test]
fn scalar_root_callbacks_validate_brackets_and_preserve_failures() {
    let session = Session::new();
    let mut ctx = CallCtx::new();
    let count = Arc::new(AtomicUsize::new(0));
    let calls = Arc::clone(&count);
    let function = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |args: &[Value]| {
                calls.fetch_add(1, Ordering::Relaxed);
                Ok(Value::Float(args[0].as_f64().unwrap().powi(2) - 2.0))
            },
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    let derivative = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| Ok(Value::Float(2.0 * args[0].as_f64().unwrap())),
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    let root = |config| {
        named(
            &session,
            "num_optimizer_find_root",
            &[Value::Function(function), config],
            &ctx,
        )
    };
    assert!(root(root_options("brent", 2.0, 0.0, 128)).is_err());
    assert_eq!(count.load(Ordering::Relaxed), 0);
    assert!(matches!(
        root(root_options("brent", 0.0, 2.0, 1)),
        Err(OpError::BudgetExceeded)
    ));
    for solver in [
        "bisection",
        "brent",
        "secant",
        "toms748",
        "newton",
        "halley",
    ] {
        let result = root(root_options(solver, 0.0, 2.0, 1024)).unwrap();
        let got = field(&result, "root").as_f64().unwrap();
        assert!((got - 2.0_f64.sqrt()).abs() < 1e-9, "{solver}: {got}");
    }
    let result = named(
        &session,
        "num_optimizer_find_root_with_derivative",
        &[
            Value::Function(function),
            Value::Function(derivative),
            root_options("newton", 0.0, 2.0, 128),
        ],
        &ctx,
    )
    .unwrap();
    assert_eq!(field(&result, "converged"), &Value::Bool(true));
    let failure = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |_: &[Value]| Err(OpError::Host("precise root callback failure".into())),
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    let error = named(
        &session,
        "num_optimizer_find_root",
        &[
            Value::Function(failure),
            root_options("brent", 0.0, 2.0, 128),
        ],
        &ctx,
    )
    .unwrap_err();
    assert!(error.to_string().contains("precise root callback failure"));
}

#[test]
fn manifest_pins_eighteen_optimizer_sources() {
    let manifest = Session::new().manifest();
    let bound = manifest
        .ops
        .iter()
        .filter(|entry| entry.name.starts_with("num_optimizer_"))
        .collect::<Vec<_>>();
    assert_eq!(bound.len(), 19);
    assert_eq!(
        bound
            .iter()
            .filter(|entry| entry.rust_path.is_some())
            .count(),
        18
    );
    assert!(bound.iter().all(|entry| {
        entry.signature_source == rssn::kernel::manifest::SignatureSource::Declared
    }));
}
