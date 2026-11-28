//! Native simplex residual plans own their geometry and preserve physical references.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

const PREPARE: &str = "num_grid_fem_residual_estimator_new";
const ESTIMATE: &str = "num_grid_fem_residual_estimator_indicators";
const ESTIMATE_VAR: &str = "num_grid_fem_residual_estimator_indicators_var";

fn call(
    s: &Session,
    op: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    s.call_value(op, args, &CallCtx::new())
}

fn must_call(
    s: &Session,
    op: &str,
    args: &[Value],
) -> Value {
    call(s, op, args).unwrap_or_else(|error| panic!("{op}: {error:?}"))
}

fn vector(
    data: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data,
        elem,
    })
}

fn floats(values: &[f64]) -> Value {
    vector(
        values.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float,
    )
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

fn power(
    base: &str,
    numer: i64,
    denom: i64,
) -> DimExpr {
    DimExpr::base(base).pow(Exp::Rat(Rational::new(numer, denom)))
}

fn frame(
    s: &Session,
    [length, time, value]: [f64; 3],
) -> Value {
    must_call(
        s,
        "num_grid_frame_new",
        &[
            measured(length, DimExpr::base("Length")),
            measured(time, DimExpr::base("Time")),
            measured(value, DimExpr::base("Temperature")),
        ],
    )
}

fn markers(
    s: &Session,
    frame: &Value,
) -> Value {
    let datum = must_call(s, "num_grid_dirichlet", &[frame.clone(), Value::Float(0.0)]);
    must_call(s, "num_grid_marker_with_default", &[datum])
}

fn compile(
    s: &Session,
    expression: &str,
    inputs: &[(&str, f64, &str)],
    output: Value,
) -> Value {
    let backend = must_call(s, "num_backend_interpreter", &[]);
    let references = inputs
        .iter()
        .map(|(_, scale, dim)| {
            must_call(
                s,
                "num_backend_reference",
                &[measured(*scale, DimExpr::base(dim))],
            )
        })
        .collect();
    let output = must_call(s, "num_backend_reference", &[output]);
    must_call(
        s,
        "num_backend_compile",
        &[
            backend,
            Value::String(expression.into()),
            vector(
                inputs
                    .iter()
                    .map(|(name, _, _)| Value::String((*name).into()))
                    .collect(),
                TypeExpr::String,
            ),
            vector(
                references,
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            output,
            Value::int(10000),
        ],
    )
}

fn source(
    s: &Session,
    expression: &str,
    [length, time, value]: [f64; 3],
) -> Value {
    compile(
        s,
        expression,
        &[
            ("t", time, "Time"),
            ("x", length, "Length"),
            ("y", length, "Length"),
            ("z", length, "Length"),
        ],
        measured(
            value / time,
            DimExpr::Mul(vec![DimExpr::base("Temperature"), power("Time", -1, 1)]),
        ),
    )
}

fn constant_problem(
    s: &Session,
    frame: &Value,
    source: Value,
) -> Value {
    must_call(
        s,
        "num_grid_fem_constant_problem",
        &[
            markers(s, frame),
            Value::Float(1.0),
            Value::Float(0.0),
            source,
            Value::int(10000),
        ],
    )
}

fn rectangle(
    s: &Session,
    frame: &Value,
) -> Value {
    must_call(
        s,
        "num_grid_simplex_rectangle",
        &[
            frame.clone(),
            Value::int(1),
            Value::int(1),
            floats(&[0.0; 2]),
            floats(&[1.0; 2]),
        ],
    )
}

fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected residual record");
    };
    &record.fields.iter().find(|(key, _)| key == name).unwrap().1
}

fn assert_indicator(
    result: &Value,
    expected: f64,
    cells: usize,
    length_power: Rational,
    scale: f64,
) {
    let Value::Tensor(indicators) = field(result, "indicators") else {
        panic!("expected indicator vector");
    };
    assert_eq!(indicators.shape, vec![cells]);
    assert_eq!(indicators.elem, TypeExpr::Float);
    for entry in &indicators.data {
        let value = entry.as_f64().unwrap();
        assert!((value - expected).abs() <= 3e-9, "{value} != {expected}");
    }
    let Value::Dyn(reference) = field(result, "reference") else {
        panic!("expected a dimensioned indicator reference");
    };
    assert_eq!(
        reference.dim,
        DimExpr::Mul(vec![
            DimExpr::base("Length").pow(Exp::Rat(length_power)),
            DimExpr::base("Temperature"),
            power("Time", -1, 1),
        ]),
    );
    assert!((reference.value.as_f64().unwrap() / scale - 1.0).abs() < 2e-15);
}

#[test]
fn constant_estimator_survives_source_space_release_and_checks_frames() {
    let s = Session::new();
    let scales = [4.0, 2.0, 3.0];
    let grid_frame = frame(&s, scales);
    let source = source(&s, "1+t", scales);
    let problem = constant_problem(&s, &grid_frame, source.clone());
    let mismatched = constant_problem(&s, &frame(&s, [8.0, 2.0, 3.0]), source);
    let mesh = rectangle(&s, &grid_frame);
    let space = must_call(&s, "num_grid_fem_space_new", &[mesh.clone(), Value::int(1)]);
    // Separate tiny plans can coexist without being charged each one's maximum bound.
    let plans: Vec<_> = (0..12)
        .map(|_| must_call(&s, PREPARE, std::slice::from_ref(&space)))
        .collect();
    let plan = &plans[0];
    let state = floats(&[0.0; 4]);
    let previous = must_call(
        &s,
        "num_grid_fem_error_indicators",
        &[
            space.clone(),
            problem.clone(),
            Value::Float(0.25),
            state.clone(),
        ],
    );
    must_call(&s, "num_backend_release", &[space]);
    must_call(&s, "num_backend_release", &[mesh]);
    for time in [0.25, 0.75] {
        let output = must_call(
            &s,
            ESTIMATE,
            &[
                plan.clone(),
                problem.clone(),
                Value::Float(time),
                state.clone(),
            ],
        );
        assert_indicator(&output, 1.0 + time, 2, Rational::int(2), 24.0);
        if time < 0.5 {
            assert_eq!(output, previous);
        }
    }
    assert!(
        call(
            &s,
            ESTIMATE,
            &[plan.clone(), mismatched, Value::Float(0.0), state.clone()]
        )
        .is_err()
    );
    assert!(
        call(
            &s,
            ESTIMATE,
            &[
                plan.clone(),
                problem.clone(),
                Value::Float(0.0),
                floats(&[0.0])
            ]
        )
        .is_err()
    );
    assert!(
        call(
            &s,
            ESTIMATE,
            &[
                plan.clone(),
                problem.clone(),
                measured(1.0, DimExpr::base("Time")),
                state.clone()
            ]
        )
        .is_err()
    );
    assert!(
        call(
            &Session::new(),
            ESTIMATE,
            &[
                plan.clone(),
                problem.clone(),
                Value::Float(0.0),
                state.clone()
            ]
        )
        .is_err()
    );
    must_call(&s, "num_backend_release", std::slice::from_ref(plan));
    assert!(
        call(
            &s,
            ESTIMATE,
            &[plan.clone(), problem, Value::Float(0.0), state]
        )
        .is_err()
    );
}

#[test]
fn variable_diffusivity_contributes_to_residual_and_must_stay_positive() {
    let s = Session::new();
    let scales = [4.0, 2.0, 3.0];
    let grid_frame = frame(&s, scales);
    let markers = markers(&s, &grid_frame);
    let spatial = [
        ("x", 4.0, "Length"),
        ("y", 4.0, "Length"),
        ("z", 4.0, "Length"),
    ];
    let diffusion_dim = DimExpr::Mul(vec![power("Length", 2, 1), power("Time", -1, 1)]);
    let reaction = compile(&s, "0", &spatial, measured(0.5, power("Time", -1, 1)));
    let source = source(&s, "0", scales);
    let make_problem = |expression| {
        let k = compile(
            &s,
            expression,
            &spatial,
            measured(8.0, diffusion_dim.clone()),
        );
        must_call(
            &s,
            "num_grid_fem_variable_problem",
            &[
                markers.clone(),
                k,
                reaction.clone(),
                source.clone(),
                Value::int(10000),
            ],
        )
    };
    let problem = make_problem("1+x");
    let invalid = make_problem("-1");
    let mesh = must_call(
        &s,
        "num_grid_simplex_interval",
        &[
            grid_frame,
            Value::int(1),
            Value::Float(0.0),
            Value::Float(1.0),
        ],
    );
    let space = must_call(&s, "num_grid_fem_space_new", &[mesh, Value::int(1)]);
    let plan = must_call(&s, PREPARE, std::slice::from_ref(&space));
    let state = floats(&[0.0, 1.0]);
    let output = must_call(
        &s,
        ESTIMATE_VAR,
        &[
            plan.clone(),
            problem.clone(),
            Value::Float(0.0),
            state.clone(),
        ],
    );
    // u=x, k=1+x and f=c=0 in normalized variables, with a unit interval.
    assert_indicator(&output, 1.0, 1, Rational::new(3, 2), 12.0);
    let previous = must_call(
        &s,
        "num_grid_fem_error_indicators_var",
        &[space, problem, Value::Float(0.0), state.clone()],
    );
    assert_eq!(output, previous);
    assert!(call(&s, ESTIMATE_VAR, &[plan, invalid, Value::Float(0.0), state]).is_err());
}

#[test]
fn residual_reference_combines_large_scales_without_intermediate_overflow() {
    let s = Session::new();
    let scales = [1e150, 1e150, 1e150];
    let grid_frame = frame(&s, scales);
    let mesh = rectangle(&s, &grid_frame);
    let space = must_call(&s, "num_grid_fem_space_new", &[mesh, Value::int(1)]);
    let plan = must_call(&s, PREPARE, &[space]);
    let problem = constant_problem(&s, &grid_frame, source(&s, "1", scales));
    let result = must_call(
        &s,
        ESTIMATE,
        &[plan, problem, Value::Float(0.0), floats(&[0.0; 4])],
    );
    assert_indicator(&result, 1.0, 2, Rational::int(2), 1e300);
}
