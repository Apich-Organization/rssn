//! Fast physical/native engine, Monte Carlo budget and wave contracts.
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use std::sync::Arc;
fn call(
    s: &rssn::Session,
    name: &str,
    args: &[Value],
) -> Value {
    s.call_value(name, args, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{name}: {e:?}"))
}
fn tensor(
    shape: Vec<usize>,
    values: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape,
        data: values,
        elem,
    })
}
fn floats(
    shape: Vec<usize>,
    values: &[f64],
    dim: DimExpr,
) -> Value {
    tensor(
        shape,
        values.iter().copied().map(Value::Float).collect(),
        TypeExpr::Float.with_dim(dim),
    )
}
fn ints(values: &[i64]) -> Value {
    tensor(
        vec![values.len()],
        values.iter().copied().map(Value::int).collect(),
        TypeExpr::Int,
    )
}
fn dim(name: &str) -> DimExpr {
    DimExpr::base(name)
}
fn quantity(
    x: f64,
    d: DimExpr,
) -> Value {
    if d == DimExpr::Dimless {
        Value::Float(x)
    } else {
        Value::Dyn(DynValue {
            dim: d,
            value: Box::new(Value::Float(x)),
        })
    }
}
fn number(x: &Value) -> f64 {
    match x {
        | Value::Dyn(x) => number(&x.value),
        | _ => x.as_f64().unwrap(),
    }
}
fn field<'a>(
    x: &'a Value,
    k: &str,
) -> &'a Value {
    let Value::Record(r) = x else {
        panic!("record")
    };
    r.get(k).unwrap()
}
#[test]
fn native_md_transfers_forcefield_and_preserves_free_particle_energy() {
    let s = rssn::Session::new();
    let f = call(&s, "num_native_md_normalized_frame", &[]);
    let b = call(
        &s,
        "num_native_engine_box_periodic_d3",
        &[
            f.clone(),
            floats(vec![3], &[10., 10., 10.], DimExpr::Dimless),
        ],
    );
    let system = call(
        &s,
        "num_native_engine_system_new_d3",
        &[
            b,
            floats(vec![2, 3], &[1., 1., 1., 3., 1., 1.], DimExpr::Dimless),
        ],
    );
    call(
        &s,
        "num_native_engine_system_set_state_d3",
        &[
            system.clone(),
            floats(vec![2, 3], &[1., 0., 0., -1., 0., 0.], DimExpr::Dimless),
            floats(vec![2], &[1., 1.], DimExpr::Dimless),
            floats(vec![2], &[0., 0.], DimExpr::Dimless),
        ],
    );
    let ff = call(
        &s,
        "num_native_engine_forcefield_new_d3",
        &[
            f,
            tensor(
                vec![0],
                vec![],
                TypeExpr::Opaque("rssn/native/sim/engine-pair/v1".into()),
            ),
            Value::Float(2.),
        ],
    );
    let sim = call(
        &s,
        "num_native_engine_simulation_new_d3",
        &[
            system,
            ff.clone(),
            Value::Float(0.2),
            Value::Float(0.01),
            Value::int(1),
        ],
    );
    assert!(
        s.call_value(
            "num_native_engine_forcefield_validate_d3",
            &[ff],
            &CallCtx::new()
        )
        .is_err()
    );
    let before = number(&call(
        &s,
        "num_native_engine_simulation_total_energy_d3",
        std::slice::from_ref(&sim),
    ));
    call(
        &s,
        "num_native_engine_simulation_run_d3",
        &[sim.clone(), Value::int(10)],
    );
    let after = number(&call(
        &s,
        "num_native_engine_simulation_total_energy_d3",
        std::slice::from_ref(&sim),
    ));
    assert!((before - 1.).abs() < 1e-14);
    assert!((after - before).abs() < 1e-14);
    let state = call(&s, "num_native_engine_simulation_state_d3", &[sim]);
    assert!((number(field(&state, "time")) - 0.1).abs() < 1e-14);
}
#[test]
fn physical_nbody_generic_source_field_restores_source_over_length_squared() {
    let s = rssn::Session::new();
    let charge = DimExpr::Mul(vec![dim("Current"), dim("Time")]);
    let f = call(
        &s,
        "num_native_md_frame",
        &[
            quantity(1., dim("Mass")),
            quantity(1., dim("Length")),
            quantity(1., dim("Time")),
            quantity(1., charge),
        ],
    );
    let out = call(
        &s,
        "num_native_nbody_direct_field_d3",
        &[
            f,
            floats(vec![2, 3], &[0., 0., 0., 2., 0., 0.], dim("Length")),
            floats(vec![2], &[1., 1.], dim("Mass")),
            quantity(1., dim("Mass")),
            quantity(0., dim("Length")),
        ],
    );
    let Value::Tensor(out) = out else {
        panic!("field")
    };
    assert!((number(&out.data[0]) - 0.25).abs() < 1e-14);
    assert!((number(&out.data[3]) + 0.25).abs() < 1e-14);
    let expected = TypeExpr::Float.with_dim(DimExpr::Mul(vec![dim("Mass"), dim("Length").pow(-2)]));
    assert!(rssn::kernel::overload::matches_type(&expected, &out.elem).accepted());
}
#[test]
fn native_wang_landau_budget_returns_error_without_an_unfinished_estimate() {
    let s = rssn::Session::new();
    let f = call(&s, "num_native_md_normalized_frame", &[]);
    let result = s.call_value(
        "num_native_spin_wang_landau_ising_d1",
        &[
            f,
            ints(&[3]),
            Value::Bool(true),
            Value::Float(0.5),
            Value::Float(0.8),
            Value::int(2),
            Value::int(1),
        ],
        &CallCtx::new(),
    );
    assert!(result.unwrap_err().to_string().contains("proposal budget"));
}
#[test]
fn native_sse_budget_invalidates_partial_mutable_state() {
    let s = rssn::Session::new();
    let f = call(&s, "num_native_md_normalized_frame", &[]);
    let bonds = call(
        &s,
        "num_native_sse_chain_bonds",
        &[Value::int(4), Value::Bool(true)],
    );
    let state = call(
        &s,
        "num_native_sse_new",
        &[
            f,
            Value::int(4),
            bonds,
            Value::Float(1.),
            Value::Float(1.),
            Value::int(1),
        ],
    );
    let err = s
        .call_value(
            "num_native_sse_sweep",
            &[state.clone(), Value::int(1)],
            &CallCtx::new(),
        )
        .unwrap_err();
    assert!(err.to_string().contains("budget"));
    assert!(
        s.call_value("num_native_sse_sample", &[state], &CallCtx::new())
            .is_err()
    );
}
#[test]
fn native_multigrid_fixed_domain_references_preserve_solution_units() {
    let s = rssn::Session::new();
    let temperature = dim("Temperature");
    let rhs_dim = DimExpr::Mul(vec![temperature.clone(), dim("Length").pow(-2)]);
    let out = call(
        &s,
        "num_native_mtm_solve_poisson_1d_multigrid",
        &[
            Value::int(3),
            floats(vec![3], &[2., 2., 2.], rhs_dim),
            Value::int(20),
            quantity(2., dim("Length")),
            quantity(3., temperature.clone()),
        ],
    );
    let Value::Tensor(solution) = field(&out, "solution") else {
        panic!("solution")
    };
    assert_eq!(solution.data.len(), 5);
    // -u''=2 K/m² on [0,2 m] gives u=x(2-x) K/m².
    // The centre is 1 K; the independent 3 K solution reference cancels.
    for (value, expected) in solution.data.iter().zip([0.0, 0.75, 1.0, 0.75, 0.0]) {
        assert!(
            (number(value) - expected).abs() < 1e-8,
            "physical Poisson node: {value:?}, expected {expected}"
        );
    }
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(temperature),
            &solution.elem
        )
        .accepted()
    );
}
#[test]
fn physical_split_step_constant_wave_has_dimensionless_norm_and_unitful_amplitude() {
    let s = rssn::Session::new();
    let solver = call(
        &s,
        "num_native_split_centered",
        &[
            ints(&[4]),
            floats(vec![1], &[2.], dim("Length")),
            quantity(1., dim("Mass")),
            quantity(1., dim("Length")),
            quantity(1., dim("Time")),
        ],
    );
    let psi_dim = dim("Length").pow(Exp::Rat(Rational::new(-1, 2)));
    let mut ctx = CallCtx::new();
    let p = psi_dim.clone();
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            move |_: &[Value]| {
                Ok(Value::Dyn(DynValue {
                    dim: p.clone(),
                    value: Box::new(Value::Complex(ComplexValue::Float { re: 1., im: 0. })),
                }))
            },
            Signature::new(
                vec![TypeExpr::Tensor {
                    rank: NatExpr::Lit(1),
                    shape: None,
                    elem: Box::new(TypeExpr::Float.with_dim(dim("Length"))),
                }],
                TypeExpr::Complex.with_dim(psi_dim.clone()),
            ),
        )))
        .unwrap();
    let psi = s
        .call_value(
            "num_native_split_sample",
            &[solver.clone(), Value::Function(callback), Value::int(4)],
            &ctx,
        )
        .unwrap();
    let norm = call(&s, "num_native_split_norm2", &[solver.clone(), psi.clone()]);
    assert!((number(&norm) - 2.).abs() < 1e-12);
    let kinetic = call(&s, "num_native_split_kinetic_energy", &[solver, psi]);
    assert!(number(&kinetic).abs() < 1e-30);
}

#[test]
fn native_gpe_density_normalizes_once_and_schrodinger_returns_mutated_final_wave() {
    use rssn::sim::models_gpe_superfluidity as gpe;
    use rssn::sim::models_schrodinger_quantum as quantum;
    let density = gpe::run_gpe_ground_state_finder(&gpe::GpeParameters {
        nx: 4,
        ny: 4,
        lx: 2.,
        ly: 2.,
        d_tau: 0.01,
        time_steps: 2,
        g: 0.,
        trap_strength: 0.,
    })
    .unwrap();
    let integral = density.iter().sum::<f64>() * 0.5 * 0.5;
    assert!((integral - 1.).abs() < 1e-12, "density integral={integral}");
    let mut psi = vec![rssn::kernels::ball::Cball::point(1., 0.).centre(); 16];
    let frames = quantum::run_schrodinger_simulation(
        &quantum::SchrodingerParameters {
            nx: 4,
            ny: 4,
            lx: 2.,
            ly: 2.,
            dt: 0.1,
            time_steps: 1,
            hbar: 1.,
            mass: 1.,
            potential: vec![1.; 16],
        },
        &mut psi,
    )
    .unwrap();
    assert_eq!(frames.len(), 1);
    for p in psi {
        assert!((p.re - 0.1f64.cos()).abs() < 1e-12);
        assert!((p.im + 0.1f64.sin()).abs() < 1e-12);
    }
}

#[test]
fn native_scenario_dispatch_executes_and_preserves_ordered_per_request_errors() {
    use rssn::kernel::value::RecordValue;
    let s = rssn::Session::new();
    let request = |name: &str, parameters: &str| {
        Value::Record(RecordValue {
            fields: vec![
                ("name".into(), Value::String(name.into())),
                ("parameters".into(), Value::String(parameters.into())),
            ],
        })
    };
    let valid = r#"{"cells":16,"t":0.01,"limiter":"van_leer"}"#;
    let value = call(
        &s,
        "num_native_scenario_run",
        &[
            Value::String("sod_shock_tube".into()),
            Value::String(valid.into()),
            Value::String("native-schema".into()),
            Value::int(100000),
        ],
    );
    let Value::String(value) = value else {
        panic!("JSON result")
    };
    let json: serde_json::Value = serde_json::from_str(&value).unwrap();
    assert_eq!(json["x"].as_array().unwrap().len(), 16);
    let requests = tensor(
        vec![3],
        vec![
            request("sod_shock_tube", valid),
            request(
                "sod_shock_tube",
                r#"{"cells":18446744073709551615,"t":0.1}"#,
            ),
            request("sod_shock_tube", valid),
        ],
        TypeExpr::Dyn(rssn::kernel::signature::RepExpr {
            repr: "term".into(),
            dim: None,
        }),
    );
    let result = call(
        &s,
        "num_native_scenario_run_many",
        &[
            requests,
            Value::String("native-schema".into()),
            Value::int(100000),
        ],
    );
    let Value::Tensor(t) = result else {
        panic!("results")
    };
    assert_eq!(t.data.len(), 3);
    for (i, expected) in ["Ok", "Err", "Ok"].iter().enumerate() {
        let Value::Variant(v) = &t.data[i] else {
            panic!("result envelope")
        };
        assert_eq!(v.tag, *expected);
    }
}
