//! Measured fixed-phase integration and fused polynomial contracts.
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

fn tensor(
    shape: Vec<usize>,
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue { shape, data, elem })
}

fn measured(
    value: Value,
    dimension: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        dim: dimension,
        value: Box::new(value),
    })
}

fn scalar(
    value: f64,
    dimension: DimExpr,
) -> Value {
    measured(Value::Float(value), dimension)
}

fn components(value: &Value) -> (f64, f64) {
    match value {
        | Value::Dyn(value) => components(&value.value),
        | Value::Complex(ComplexValue::Float { re, im }) => (*re, *im),
        | _ => panic!("expected complex output: {value:?}"),
    }
}

fn accepts(
    actual: &TypeExpr,
    expected: &TypeExpr,
) {
    assert!(
        rssn::kernel::overload::matches_type(expected, actual).accepted(),
        "actual {actual:?}, expected {expected:?}"
    );
}

#[test]
fn measured_polynomial_antidifference_preserves_tiny_coordinates_and_widths() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &ctx).unwrap();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let coefficients = |values: &[f64]| {
        tensor(
            vec![values.len()],
            values.iter().copied().map(Value::Float).collect(),
            TypeExpr::Float.with_dim(time.clone()),
        )
    };
    let axis = |value| scalar(value, length.clone());
    let build = "num_analysis_polynomial_antidifference_new";
    let eval = "num_analysis_polynomial_antidifference_eval";
    let delta = "num_analysis_polynomial_antidifference_difference";
    let scalar_result = |value| {
        let Value::Dyn(dynamic) = value else {
            panic!("expected a measured antidifference")
        };
        assert_eq!(dynamic.dim, time);
        let Value::Float(value) = *dynamic.value else {
            panic!("expected a Float magnitude")
        };
        value
    };
    let plan = call(build, &[coefficients(&[0.0, 1.0]), axis(0.0), axis(2.0)]);
    assert_eq!(
        scalar_result(call(eval, &[plan.clone(), axis(2.0_f64.powi(-59))])),
        -2.0_f64.powi(-61)
    );
    let start = 3.0_f64;
    let end = start.next_up();
    assert_eq!(
        scalar_result(call(delta, &[plan.clone(), axis(start), axis(end)])),
        (end - start) * 0.5
    );
    assert_eq!(
        scalar_result(call(delta, &[plan.clone(), axis(end), axis(start)])),
        (start - end) * 0.5
    );
    assert!(
        session
            .call_value(delta, &[plan, axis(3.0), scalar(4.0, time.clone())], &ctx)
            .is_err()
    );

    let constant = call(build, &[coefficients(&[3.0]), axis(0.0), axis(f64::MAX)]);
    assert_eq!(
        scalar_result(call(delta, &[constant, axis(-f64::MAX), axis(f64::MAX)])),
        6.0
    );
}

#[test]
fn levin_plan_solves_independent_amplitude_fields_after_releasing_phase_programs() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name, args: &[Value]| {
        session
            .call_value(name, args, &ctx)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"))
    };
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let mass = DimExpr::base("Mass");
    let reference = |dimension| call("num_backend_reference", &[scalar(1.0, dimension)]);
    let coordinate_ref = reference(length.clone());
    let compile = |body: &str, output| {
        call(
            "num_backend_compile_current",
            &[
                Value::String(body.into()),
                tensor(vec![1], vec![Value::String("x".into())], TypeExpr::String),
                tensor(
                    vec![1],
                    vec![coordinate_ref.clone()],
                    coordinate_ref.type_expr(),
                ),
                reference(output),
                Value::int(1_000),
            ],
        )
    };
    let phase = compile("x", DimExpr::Dimless);
    let derivative = compile("1", length.clone().pow(-1));
    let build = "num_analysis_levin_plan_new";
    let node_count = 6;
    let plan = call(
        build,
        &[
            scalar(1.0, length.clone()),
            scalar(4.0, length.clone()),
            Value::Float(0.0),
            Value::int(node_count),
            Value::int(1),
            phase.clone(),
            derivative.clone(),
        ],
    );
    call("num_backend_release", &[phase]);
    call("num_backend_release", &[derivative]);
    let nodes = call("num_analysis_levin_plan_nodes", &[plan.clone()]);
    let Value::Tensor(nodes) = nodes else {
        panic!("expected measured Levin nodes")
    };
    assert_eq!(nodes.shape, [node_count]);
    accepts(&nodes.elem, &TypeExpr::Float.with_dim(length.clone()));
    let amplitude = compile("2", mass.clone());
    let result = call(
        "num_analysis_levin_plan_integrate",
        &[plan.clone(), amplitude, Value::int(node_count)],
    );
    accepts(
        &result.type_expr(),
        &TypeExpr::Complex.with_dim(DimExpr::Mul(vec![length.clone(), mass])),
    );
    let (re, im) = components(&result);
    assert!((re - 6.0).abs() < 1e-12 && im == 0.0);
    let times = tensor(
        vec![node_count],
        vec![Value::Float(5.0); node_count],
        TypeExpr::Float.with_dim(time.clone()),
    );
    let result = call("num_analysis_levin_plan_integrate_samples", &[plan, times]);
    accepts(
        &result.type_expr(),
        &TypeExpr::Complex.with_dim(DimExpr::Mul(vec![length, time])),
    );
    let (re, im) = components(&result);
    assert!((re - 15.0).abs() < 1e-12 && im == 0.0);
}

#[test]
fn fused_polynomial_evaluation_keeps_value_and_jacobian_shape_and_dimensions() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name, args: &[Value]| {
        session
            .call_value(name, args, &ctx)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"))
    };
    let length = DimExpr::base("Length");
    let mass = DimExpr::base("Mass");
    let slope = DimExpr::Mul(vec![mass.clone(), length.clone().pow(-1)]);
    let system = call(
        "num_native_polynomial_system_new",
        &[
            tensor(vec![1], vec![Value::Float(1.0)], TypeExpr::Float),
            tensor(vec![1, 1], vec![Value::int(2)], TypeExpr::Int),
            tensor(vec![2], vec![Value::int(0), Value::int(1)], TypeExpr::Int),
            scalar(2.0, length.clone()),
            scalar(3.0, mass.clone()),
        ],
    );
    let exact = |re, dimension| {
        measured(
            Value::Complex(ComplexValue::Float { re, im: 0.0 }),
            dimension,
        )
    };
    let point = tensor(
        vec![1],
        vec![Value::Complex(ComplexValue::Float { re: 4.0, im: 0.0 })],
        TypeExpr::Complex.with_dim(length.clone()),
    );
    let result = call(
        "num_native_ball_values_and_jacobian_at_terms",
        &[system.clone(), point],
    );
    let Value::Record(point_result) = result else {
        panic!("expected polynomial data record")
    };
    for (name, shape, expected, dimension) in [
        ("values", vec![1], 12.0, mass.clone()),
        ("jacobian", vec![1, 1], 6.0, slope.clone()),
    ] {
        let Some(Value::Tensor(out)) = point_result.get(name) else {
            panic!("missing {name}")
        };
        assert_eq!(out.shape, shape);
        accepts(&out.elem, &TypeExpr::Complex.with_dim(dimension));
        assert_eq!(components(&out.data[0]), (expected, 0.0));
    }
    let bound = call(
        "num_native_cball_point",
        &[scalar(4.0, length.clone()), scalar(0.0, length)],
    );
    let elem = bound.type_expr();
    let result = call(
        "num_native_ball_eval_value_and_jacobian_terms",
        &[system, tensor(vec![1], vec![bound], elem)],
    );
    let Value::Record(box_result) = result else {
        panic!("expected boxed polynomial data record")
    };
    for (name, shape, expected, dimension) in [
        ("values", vec![1], 12.0, mass),
        ("jacobian", vec![1, 1], 6.0, slope),
    ] {
        let Some(Value::Tensor(out)) = box_result.get(name) else {
            panic!("missing {name}")
        };
        assert_eq!(out.shape, shape);
        assert_eq!(
            call(
                "num_native_cball_contains",
                &[out.data[0].clone(), exact(expected, dimension)]
            ),
            Value::Bool(true)
        );
    }
}

#[test]
fn sampled_series_reports_measured_estimates_and_exact_consumption() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let length = DimExpr::base("Length");
    let terms = tensor(
        vec![3],
        vec![Value::Float(3.0), Value::Float(4.0), Value::Float(20.0)],
        TypeExpr::Float.with_dim(length.clone()),
    );
    let call = |budget| {
        session
            .call_value(
                "num_analysis_series_sum_terms",
                &[
                    terms.clone(),
                    scalar(1e-10, length.clone()),
                    Value::int(budget),
                ],
                &ctx,
            )
            .unwrap()
    };
    let Value::Record(short) = call(2) else {
        panic!("expected a series estimate record")
    };
    assert_eq!(short.get("terms"), Some(&Value::int(2)));
    let Some(Value::Variant(value)) = short.get("value") else {
        panic!("expected an estimate option")
    };
    assert_eq!(value.tag, "Some");
    let [Value::Dyn(estimate)] = value.payload.as_slice() else {
        panic!("expected one measured estimate")
    };
    assert_eq!(estimate.dim, length);
    assert_eq!(*estimate.value, Value::Float(7.0));
    let Some(Value::Variant(error)) = short.get("error") else {
        panic!("expected an error option")
    };
    assert_eq!(error.tag, "None");
    let Value::Record(empty) = call(0) else {
        panic!("expected an empty series estimate record")
    };
    assert_eq!(empty.get("terms"), Some(&Value::int(0)));
    assert!(matches!(empty.get("value"), Some(Value::Variant(v)) if v.tag == "None"));
    assert!(
        session
            .call_value(
                "num_analysis_series_sum_terms",
                &[terms, Value::Float(0.0), Value::int(2)],
                &ctx
            )
            .is_err()
    );
}
