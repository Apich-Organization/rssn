//! Measured state-space discretisation and simulation contracts.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::overload::matches_type;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

fn tensor(
    shape: &[usize],
    data: &[f64],
    dim: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        elem: TypeExpr::Float.with_dim(dim),
        data: data.iter().copied().map(Value::Float).collect(),
    })
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

fn product(parts: &[DimExpr]) -> DimExpr {
    DimExpr::Mul(parts.to_vec())
}

fn scalar(value: &Value) -> f64 {
    match value {
        | Value::Dyn(measured) => scalar(&measured.value),
        | Value::Float(value) => *value,
        | _ => panic!("expected a real state-space component"),
    }
}

fn checked_tensor<'a>(
    value: &'a Value,
    shape: &[usize],
    dim: DimExpr,
) -> &'a TensorValue {
    let Value::Tensor(tensor) = value else {
        panic!("expected a tensor")
    };
    assert_eq!(tensor.shape, shape);
    assert!(matches_type(&TypeExpr::Float.with_dim(dim), &tensor.elem).accepted());
    tensor
}

fn call(
    session: &Session,
    operator: &str,
    values: &[Value],
) -> Value {
    session
        .call_value(operator, values, &CallCtx::new())
        .unwrap()
}

#[test]
fn zoh_and_measured_trajectories_agree_with_an_analytic_rc_model() {
    let session = Session::new();
    let x = DimExpr::base("Length");
    let u = DimExpr::base("Current");
    let y = DimExpr::base("Mass");
    let time = DimExpr::base("Time");
    let dt = 0.125;
    let b_discrete = product(&[x.clone(), u.clone().pow(-1)]);
    let c_dim = product(&[y.clone(), x.clone().pow(-1)]);
    let d_dim = product(&[y.clone(), u.clone().pow(-1)]);
    let matrices = [
        tensor(&[1, 1], &[-2.0], time.clone().pow(-1)),
        tensor(
            &[1, 1],
            &[3.0],
            product(&[b_discrete.clone(), time.clone().pow(-1)]),
        ),
        tensor(&[2, 1], &[4.0, 2.0], c_dim.clone()),
        tensor(&[2, 1], &[5.0, -1.0], d_dim.clone()),
    ];
    let mut discretise = matrices.to_vec();
    discretise.push(measured(dt, time.clone()));
    let model = call(&session, "num_control_c2d_zoh", &discretise);
    let Value::Record(fields) = &model else {
        panic!("expected the discrete model record")
    };
    let a = checked_tensor(fields.get("a").unwrap(), &[1, 1], DimExpr::Dimless);
    let b = checked_tensor(fields.get("b").unwrap(), &[1, 1], b_discrete);
    checked_tensor(fields.get("c").unwrap(), &[2, 1], c_dim);
    checked_tensor(fields.get("d").unwrap(), &[2, 1], d_dim);
    let period = fields.get("dt").unwrap();
    assert!(matches_type(&TypeExpr::Float.with_dim(time.clone()), &period.type_expr()).accepted());
    assert_eq!(scalar(period).to_bits(), dt.to_bits());
    assert!((scalar(&a.data[0]) - (-2.0 * dt).exp()).abs() < 2e-14);
    assert!((scalar(&b.data[0]) - 1.5 * (1.0 - (-2.0 * dt).exp())).abs() < 2e-14);

    let initial = tensor(&[1], &[0.0], x);
    let inputs = tensor(&[3, 1], &[1.0; 3], u);
    let times = [0.0, dt, 2.0 * dt];
    let mut continuous = matrices.to_vec();
    continuous.extend([initial.clone(), inputs.clone(), tensor(&[3], &times, time)]);
    let smooth = call(&session, "num_control_simulate_continuous", &continuous);
    let discrete = call(
        &session,
        "num_control_simulate_discrete",
        &[model, initial, inputs],
    );
    for result in [&smooth, &discrete] {
        let output = checked_tensor(result, &[3, 2], y.clone());
        let (pairs, tail) = output.data.as_chunks::<2>();
        assert!(tail.is_empty());
        assert_eq!(pairs.len(), times.len());
        for (&time, [first, second]) in times.iter().zip(pairs) {
            let state = 1.5 * (1.0 - (-2.0 * time).exp());
            assert!((scalar(first) - (4.0 * state + 5.0)).abs() < 2e-13);
            assert!((scalar(second) - (2.0 * state - 1.0)).abs() < 2e-13);
        }
    }
}

#[test]
fn zero_state_feedthrough_and_empty_output_keep_matrix_extents() {
    let session = Session::new();
    let x = DimExpr::base("Length");
    let u = DimExpr::base("Current");
    let y = DimExpr::base("Mass");
    let time = DimExpr::base("Time");
    let matrices = [
        tensor(&[0, 0], &[], time.clone().pow(-1)),
        tensor(
            &[0, 1],
            &[],
            product(&[x.clone(), u.clone().pow(-1), time.clone().pow(-1)]),
        ),
        tensor(&[2, 0], &[], product(&[y.clone(), x.clone().pow(-1)])),
        tensor(
            &[2, 1],
            &[2.0, -3.0],
            product(&[y.clone(), u.clone().pow(-1)]),
        ),
    ];
    let mut discretise = matrices.to_vec();
    discretise.push(measured(1.0, time.clone()));
    let model = call(&session, "num_control_c2d_zoh", &discretise);
    let Value::Record(fields) = &model else {
        panic!("expected a discrete model")
    };
    checked_tensor(
        fields.get("b").unwrap(),
        &[0, 1],
        product(&[x.clone(), u.clone().pow(-1)]),
    );
    let initial = tensor(&[0], &[], x);
    let inputs = tensor(&[0, 1], &[], u.clone());
    let discrete = call(
        &session,
        "num_control_simulate_discrete",
        &[model.clone(), initial.clone(), inputs.clone()],
    );
    checked_tensor(&discrete, &[0, 2], y.clone());
    let mut continuous = matrices.to_vec();
    continuous.extend([initial.clone(), inputs, tensor(&[0], &[], time)]);
    let continuous = call(&session, "num_control_simulate_continuous", &continuous);
    checked_tensor(&continuous, &[0, 2], y.clone());

    let visible = call(
        &session,
        "num_control_simulate_discrete",
        &[model, initial, tensor(&[1, 1], &[4.0], u)],
    );
    let output = checked_tensor(&visible, &[1, 2], y);
    assert_eq!(
        output.data.iter().map(scalar).collect::<Vec<_>>(),
        [8.0, -12.0]
    );
}

#[test]
fn rejects_incompatible_units_extents_and_sample_times() {
    let session = Session::new();
    let time = DimExpr::base("Time");
    let a = tensor(&[1, 1], &[-2.0], time.clone().pow(-1));
    let b = tensor(&[1, 1], &[1.0], time.clone().pow(-1));
    let c = tensor(&[1, 1], &[1.0], DimExpr::Dimless);
    let d = tensor(&[1, 1], &[0.0], DimExpr::Dimless);
    let mut args = [a, b, c, d, measured(0.25, time.clone())];
    let good = call(&session, "num_control_c2d_zoh", &args);
    let rejects_zoh = |args: &[Value]| {
        assert!(
            session
                .call_value("num_control_c2d_zoh", args, &CallCtx::new())
                .is_err()
        );
    };
    args[4] = Value::Float(0.25);
    rejects_zoh(&args);
    args[4] = measured(-0.25, time);
    rejects_zoh(&args);
    args[4] = measured(0.25, DimExpr::base("Time"));
    args[3] = tensor(&[1, 1], &[0.0], DimExpr::base("Length"));
    rejects_zoh(&args);
    args[3] = tensor(&[1, 2], &[0.0, 0.0], DimExpr::Dimless);
    rejects_zoh(&args);

    assert!(
        session
            .call_value(
                "num_control_simulate_discrete",
                &[
                    good,
                    tensor(&[1], &[0.0], DimExpr::Dimless),
                    tensor(&[0, 0], &[], DimExpr::Dimless)
                ],
                &CallCtx::new(),
            )
            .is_err()
    );
}
