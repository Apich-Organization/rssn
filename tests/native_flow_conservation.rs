//! Actual native DG/FV methods retain physical component domains and enforce work budgets.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    n: &str,
    a: &[Value],
) -> Value {
    s.call_value(n, a, &CallCtx::new()).unwrap()
}
fn base(n: &str) -> DimExpr {
    DimExpr::Base(n.into())
}
fn measured(
    x: f64,
    d: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: d,
    })
}
fn num(x: f64) -> Value {
    Value::Float(x)
}
fn vector(
    v: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![v.len()],
        data: v,
        elem,
    })
}
fn numbers(v: &[f64]) -> Value {
    vector(v.iter().copied().map(num).collect(), TypeExpr::Float)
}
fn sizes(v: &[usize]) -> Value {
    vector(v.iter().map(|v| Value::int(*v)).collect(), TypeExpr::Int)
}
fn strings(v: &[&str]) -> Value {
    vector(
        v.iter().map(|s| Value::String((*s).into())).collect(),
        TypeExpr::String,
    )
}
fn record(v: Vec<(&str, Value)>) -> Value {
    Value::Record(RecordValue {
        fields: v.into_iter().map(|(n, v)| (n.into(), v)).collect(),
    })
}
fn field<'a>(
    v: &'a Value,
    n: &str,
) -> &'a Value {
    let Value::Record(v) = v else {
        panic!("record expected")
    };
    &v.fields.iter().find(|(k, _)| k == n).unwrap().1
}
fn number(v: &Value) -> f64 {
    match v {
        | Value::Dyn(v) => number(&v.value),
        | _ => v.as_f64().unwrap(),
    }
}
fn compile(
    s: &Session,
    expression: &str,
    inputs: &[(f64, DimExpr)],
    output: (f64, DimExpr),
) -> Value {
    let b = call(s, "num_backend_interpreter", &[]);
    let refs = inputs
        .iter()
        .map(|(v, d)| call(s, "num_backend_reference", &[measured(*v, d.clone())]))
        .collect();
    let out = call(s, "num_backend_reference", &[measured(output.0, output.1)]);
    let names = (0..inputs.len())
        .map(|i| Value::String(format!("x{i}")))
        .collect();
    call(
        s,
        "num_backend_compile",
        &[
            b,
            Value::String(expression.into()),
            vector(names, TypeExpr::String),
            vector(
                refs,
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            out,
            Value::int(1000000),
        ],
    )
}
#[test]
fn native_dg_constant_integral_and_flux_keep_units_and_small_budget_rejects() {
    let s = Session::new();
    let f = call(
        &s,
        "num_grid_frame_new",
        &[
            measured(2.0, base("Length")),
            measured(4.0, base("Time")),
            measured(3.0, base("Mass")),
        ],
    );
    let law = call(&s, "num_dg_advection", &[f.clone(), numbers(&[1.0])]);
    let flux = call(
        &s,
        "num_dg_flux",
        &[law.clone(), numbers(&[2.0]), Value::int(0)],
    );
    assert!((number(field(&flux, "component0")) - 3.0).abs() < 1e-12);
    assert!(matches!(field(&flux, "component0"), Value::Dyn(_)));
    let mesh = call(
        &s,
        "num_dg_mesh_line",
        &[f, Value::int(2), num(0.0), num(1.0), Value::int(2)],
    );
    let bcs = call(
        &s,
        "num_dg_bcs",
        &[
            law,
            strings(&["Periodic", "Periodic", "Outflow", "Outflow"]),
        ],
    );
    let initial = compile(
        &s,
        "2",
        &[(1.0, base("Length")), (1.0, base("Length"))],
        (3.0, base("Mass")),
    );
    let initial = vector(
        vec![initial],
        TypeExpr::Opaque("rssn/native/compiled/v1".into()),
    );
    let opts = record(vec![
        ("filter", Value::Bool(false)),
        ("filter_alpha", num(36.0)),
        ("filter_order", Value::int(8)),
    ]);
    let mut args = vec![
        mesh,
        bcs,
        initial,
        num(0.25),
        num(0.2),
        opts,
        Value::int(10000),
        Value::int(100000),
    ];
    let solution = call(&s, "num_solve_dg", &args);
    let total = call(&s, "num_dg_integral", &[solution, Value::int(0)]);
    assert!((number(&total) - 12.0).abs() < 1e-10);
    assert!(matches!(total, Value::Dyn(_)));
    args[7] = Value::int(1);
    let err = s
        .call_value("num_solve_dg", &args, &CallCtx::new())
        .unwrap_err();
    assert!(format!("{err:?}").contains("work budget exhausted"));
}
#[test]
fn native_fv_constant_state_and_amr_conserve_mass_per_omitted_area() {
    let s = Session::new();
    let r = call(
        &s,
        "num_grid_mechanical_references",
        &[
            measured(2.0, base("Length")),
            measured(4.0, base("Time")),
            measured(8.0, base("Mass")),
            measured(1.0, base("Temperature")),
        ],
    );
    let grid = call(
        &s,
        "num_fv_grid_new",
        &[r.clone(), sizes(&[8]), numbers(&[0.0]), numbers(&[1.0])],
    );
    let opts = record(vec![
        ("cfl", num(0.4)),
        ("reconstruction", Value::String("FirstOrder".into())),
        ("riemann", Value::String("Hllc".into())),
    ]);
    let cfg = call(&s, "num_fv_config_euler", &[r, num(1.4), opts]);
    let density = DimExpr::Mul(vec![
        base("Mass"),
        DimExpr::Inverse(Box::new(DimExpr::Pow(
            Box::new(base("Length")),
            rssn::kernel::signature::Exp::Rat(rssn::kernel::signature::Rational::new(3, 1)),
        ))),
    ]);
    let velocity = DimExpr::Mul(vec![
        base("Length"),
        DimExpr::Inverse(Box::new(base("Time"))),
    ]);
    let pressure = DimExpr::Mul(vec![
        base("Mass"),
        DimExpr::Inverse(Box::new(base("Length"))),
        DimExpr::Pow(
            Box::new(base("Time")),
            rssn::kernel::signature::Exp::Rat(rssn::kernel::signature::Rational::new(-2, 1)),
        ),
    ]);
    let inputs = vec![(1.0, base("Length")); 3];
    let refs = [
        (1.0, density),
        (0.5, velocity.clone()),
        (0.5, velocity.clone()),
        (0.5, velocity),
        (0.25, pressure),
    ];
    let values = ["1", "0.2", "0", "0", "1"]
        .into_iter()
        .zip(refs)
        .map(|(e, r)| compile(&s, e, &inputs, r))
        .collect();
    let initial = vector(values, TypeExpr::Opaque("rssn/native/compiled/v1".into()));
    let sol = call(
        &s,
        "num_fv_solve",
        &[
            grid.clone(),
            cfg.clone(),
            initial.clone(),
            num(0.1),
            Value::int(256),
            Value::int(10000),
        ],
    );
    let totals = call(&s, "num_fv_totals", &[sol]);
    assert!((number(field(&totals, "component0")) - 2.0).abs() < 1e-11);
    let a = call(
        &s,
        "num_fv_amr_config",
        &[
            cfg,
            Value::int(1),
            Value::int(4),
            num(0.05),
            Value::int(1),
            num(0.6),
            Value::Bool(true),
        ],
    );
    let sol = call(
        &s,
        "num_fv_solve_amr",
        &[
            grid,
            a,
            initial,
            num(0.1),
            Value::int(256),
            Value::int(10000),
        ],
    );
    let totals = call(&s, "num_fv_amr_totals", &[sol]);
    assert!((number(field(&totals, "component0")) - 2.0).abs() < 1e-11);
}

#[test]
fn continuum_derived_references_fail_before_callback_and_cancel_extreme_factors() {
    let s = Session::new();
    let frame = call(
        &s,
        "num_grid_frame_new",
        &[
            measured(1.0, base("Length")),
            measured(1e-308, base("Time")),
            measured(1e308, base("Mass")),
        ],
    );
    let boundary = call(&s, "num_grid_dirichlet", &[frame, num(0.0)]);
    let faces = call(&s, "num_grid_faces_uniform", &[boundary]);
    let grid = call(
        &s,
        "num_grid_structured_new",
        &[faces, sizes(&[3]), numbers(&[0.0]), numbers(&[1.0])],
    );
    let initial = compile(
        &s,
        "0",
        &vec![(1.0, base("Length")); 3],
        (1.0, base("Mass")),
    );
    let output = DimExpr::Mul(vec![base("Mass"), DimExpr::Inverse(Box::new(base("Time")))]);
    let forcing = compile(
        &s,
        "0",
        &[
            (1.0, base("Time")),
            (1.0, base("Length")),
            (1.0, base("Length")),
            (1.0, base("Length")),
        ],
        (1.0, output),
    );
    let error = s
        .call_value(
            "num_grid_heat_fixed",
            &[
                grid,
                initial,
                forcing,
                Value::int(100),
                num(1.0),
                num(0.1),
                Value::int(1),
                num(0.5),
            ],
            &CallCtx::new(),
        )
        .unwrap_err();
    assert!(format!("{error:?}").contains("overflows Float"));
    // Energy M*L^2/T^2 is representable although L^2 overflows and T^-2 underflows.
    let mechanical = call(
        &s,
        "num_grid_mechanical_references",
        &[
            measured(1e200, base("Length")),
            measured(1e200, base("Time")),
            measured(1.0, base("Mass")),
            measured(1.0, base("Temperature")),
        ],
    );
    let beam = call(
        &s,
        "num_grid_beam_new",
        &[mechanical, num(1.0), num(1.0), num(1.0), num(1.0), num(0.0)],
    );
    let stiffness = call(&s, "num_grid_beam_local_stiffness_matrix", &[beam]);
    assert!((number(field(&stiffness, "quadratic_form_reference")) - 1.0).abs() < 1e-14);
}

#[test]
fn native_mac_state_is_callable_and_keeps_velocity_pressure_and_length_domains() {
    let s = Session::new();
    let refs = call(
        &s,
        "num_grid_mechanical_references",
        &[
            measured(2.0, base("Length")),
            measured(4.0, base("Time")),
            measured(8.0, base("Mass")),
            measured(1.0, base("Temperature")),
        ],
    );
    let velocity = DimExpr::Mul(vec![
        base("Length"),
        DimExpr::Inverse(Box::new(base("Time"))),
    ]);
    let f = call(
        &s,
        "num_grid_frame_new",
        &[
            measured(2.0, base("Length")),
            measured(4.0, base("Time")),
            measured(0.5, velocity.clone()),
        ],
    );
    let wall = call(&s, "num_grid_fluid_wall", &[f]);
    let faces = call(&s, "num_grid_fluid_faces_uniform", &[wall]);
    let options = record(vec![
        ("rotational", Value::Bool(true)),
        ("solver", Value::String("Direct".into())),
        ("solver_tolerance", num(1e-12)),
        ("convection", Value::String("Central".into())),
        ("les_constant", num(0.1)),
        ("les", Value::String("None".into())),
        ("solid", vector(vec![], TypeExpr::Bool)),
        ("beta_min", num(0.5)),
        ("conservative", Value::Bool(true)),
    ]);
    let none = Value::Variant(rssn::kernel::value::VariantValue {
        tag: "None".into(),
        payload: vec![],
    });
    let callbacks = record(
        [
            "force_x",
            "force_y",
            "force_z",
            "level_set",
            "body_vx",
            "body_vy",
            "body_vz",
            "temperature",
            "density",
        ]
        .into_iter()
        .map(|n| (n, none.clone()))
        .collect(),
    );
    let config = call(
        &s,
        "num_ns_config_new",
        &[
            refs,
            faces,
            sizes(&[3, 3]),
            numbers(&[1.0, 1.0]),
            num(1.0),
            options,
            callbacks,
            Value::int(10000),
        ],
    );
    let initial = compile(&s, "0", &vec![(1.0, base("Length")); 3], (0.5, velocity));
    let initial = vector(
        vec![initial.clone(), initial.clone(), initial],
        TypeExpr::Opaque("rssn/native/compiled/v1".into()),
    );
    let state = call(
        &s,
        "num_ns_new",
        &[
            config.clone(),
            num(0.01),
            initial.clone(),
            Value::int(10000),
        ],
    );
    assert!(
        s.call_value(
            "num_ns_new",
            &[config, num(0.01), initial, Value::int(10000)],
            &CallCtx::new()
        )
        .is_err()
    );
    assert!(matches!(
        call(&s, "num_ns_step", std::slice::from_ref(&state)),
        Value::Bool(true)
    ));
    let cell_size = call(&s, "num_ns_cell_size", &[state.clone(), Value::int(0)]);
    assert!((number(&cell_size) - 2.0 / 3.0).abs() < 1e-12);
    assert!(matches!(cell_size, Value::Dyn(_)));
    let pressure = call(
        &s,
        "num_ns_sample_pressure",
        &[state.clone(), numbers(&[0.5, 0.5, 0.0])],
    );
    let speed = call(
        &s,
        "num_ns_sample_velocity",
        &[state.clone(), numbers(&[0.5, 0.5, 0.0]), Value::int(0)],
    );
    assert_eq!(number(&pressure), 0.0);
    assert_eq!(number(&speed), 0.0);
    let (Value::Dyn(p), Value::Dyn(u)) = (pressure, speed) else {
        panic!("physical witnesses expected")
    };
    assert_ne!(p.dim, u.dim);
    assert!(
        s.call_value(
            "num_ns_sample_velocity",
            &[state, numbers(&[0.5, 0.5, 0.0]), Value::int(2)],
            &CallCtx::new()
        )
        .is_err()
    );
}
