//! Targeted execution tests for physical tape frames and stateful adapters.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::RepExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    n: &str,
    a: &[Value],
) -> Value {
    s.call_value(n, a, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{n}: {e:?}"))
}
fn base(n: &str) -> DimExpr {
    DimExpr::base(n)
}
fn ratio(
    a: DimExpr,
    b: DimExpr,
) -> DimExpr {
    DimExpr::Mul(vec![a, DimExpr::Inverse(Box::new(b))])
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
fn vector(v: Vec<Value>) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![v.len()],
        elem: TypeExpr::Dyn(RepExpr {
            repr: "term".into(),
            dim: None,
        }),
        data: v,
    })
}
fn strings(v: &[&str]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![v.len()],
        elem: TypeExpr::String,
        data: v.iter().map(|s| Value::String((*s).into())).collect(),
    })
}
fn data(v: &Value) -> &[Value] {
    let Value::Tensor(t) = v else {
        panic!("tensor expected")
    };
    &t.data
}
fn number(v: &Value) -> f64 {
    match v {
        | Value::Dyn(d) => number(&d.value),
        | _ => v.as_f64().unwrap(),
    }
}
fn unit(
    v: &Value,
    d: DimExpr,
) {
    let expected = TypeExpr::Float.with_dim(d);
    let actual = match v {
        | Value::Dyn(v) => TypeExpr::Float.with_dim(v.dim.clone()),
        | Value::Quantity(v) => TypeExpr::Float.with_dim(v.unit.dimension.clone()),
        | _ => v.type_expr(),
    };
    assert!(
        rssn::kernel::overload::matches_type(&expected, &actual).accepted(),
        "{expected:?} vs {actual:?}"
    );
}
fn tape(
    s: &Session,
    expression: &str,
    refs: Vec<Value>,
    out: Value,
) -> Value {
    let names = if refs.len() == 1 {
        vec!["x"]
    } else {
        vec!["x", "y"]
    };
    call(
        s,
        "num_native_tape_lower",
        &[
            strings(&[expression]),
            strings(&names),
            strings(&[]),
            vector(refs),
            vector(vec![]),
            vector(vec![out]),
        ],
    )
}
#[test]
fn tape_mixed_coordinates_restore_value_gradient_hessian_and_vjp_units() {
    let s = Session::new();
    let l = base("Length");
    let t = base("Time");
    let velocity = ratio(l.clone(), t.clone());
    let tape = tape(
        &s,
        "x/y",
        vec![measured(1., l.clone()), measured(1., t.clone())],
        measured(1., velocity.clone()),
    );
    let x = vector(vec![measured(4., l.clone()), measured(2., t.clone())]);
    let p = vector(vec![]);
    let out = call(
        &s,
        "num_native_tape_eval",
        &[tape.clone(), x.clone(), p.clone()],
    );
    assert_eq!(number(&data(&out)[0]), 2.);
    unit(&data(&out)[0], velocity.clone());
    for n in [
        "num_native_tape_gradient",
        "num_native_tape_gradient_reverse",
    ] {
        let g = call(&s, n, &[tape.clone(), Value::int(0), x.clone(), p.clone()]);
        assert_eq!(number(&data(&g)[0]), 0.5);
        assert_eq!(number(&data(&g)[1]), -1.);
        unit(&data(&g)[0], ratio(velocity.clone(), l.clone()));
        unit(&data(&g)[1], ratio(velocity.clone(), t.clone()));
    }
    let h = call(
        &s,
        "num_native_tape_hessian",
        &[tape.clone(), Value::int(0), x.clone(), p.clone()],
    );
    let hs = data(&h);
    assert_eq!(number(&hs[0]), 0.);
    assert_eq!(number(&hs[1]), -0.25);
    assert_eq!(number(&hs[2]), -0.25);
    assert_eq!(number(&hs[3]), 1.);
    unit(&hs[1], ratio(ratio(velocity.clone(), l.clone()), t.clone()));
    let energy = DimExpr::Mul(vec![
        base("Mass"),
        l.clone().pow(2),
        DimExpr::Inverse(Box::new(t.clone().pow(2))),
    ]);
    let weights = vector(vec![measured(3., ratio(energy.clone(), velocity))]);
    let g = call(
        &s,
        "num_native_tape_vjp",
        &[
            tape.clone(),
            x.clone(),
            p.clone(),
            weights,
            measured(1., energy.clone()),
        ],
    );
    assert_eq!(number(&data(&g)[0]), 1.5);
    assert_eq!(number(&data(&g)[1]), -3.);
    unit(&data(&g)[0], ratio(energy, l));
    let bad = vector(vec![measured(4., t.clone()), measured(2., t)]);
    assert!(
        s.call_value("num_native_tape_eval", &[tape, bad, p], &CallCtx::new())
            .is_err()
    );
}
#[test]
fn derivative_reference_products_do_not_overflow_before_cancellation() {
    let s = Session::new();
    let tape = tape(
        &s,
        "1e308*x",
        vec![measured(1e300, base("Length"))],
        measured(1e300, base("Time")),
    );
    let g = call(
        &s,
        "num_native_tape_gradient",
        &[
            tape,
            Value::int(0),
            vector(vec![measured(1., base("Length"))]),
            vector(vec![]),
        ],
    );
    assert_eq!(number(&data(&g)[0]), 1e308);
}
#[test]
fn particle_frame_keeps_energy_momentum_and_kelvin_distinct() {
    let s = Session::new();
    let mass = base("Mass");
    let l = base("Length");
    let t = base("Time");
    let q = DimExpr::Mul(vec![base("Current"), t.clone()]);
    let v = ratio(l.clone(), t.clone());
    let frame = call(
        &s,
        "num_native_md_frame",
        &[
            measured(1., mass.clone()),
            measured(1., l.clone()),
            measured(1., t.clone()),
            measured(1., q),
        ],
    );
    let positions = Value::Tensor(TensorValue {
        shape: vec![1],
        elem: TypeExpr::Float.with_dim(l.clone()),
        data: vec![Value::Float(1.)],
    });
    let velocities = Value::Tensor(TensorValue {
        shape: vec![1],
        elem: TypeExpr::Float.with_dim(v.clone()),
        data: vec![Value::Float(3.)],
    });
    let p = call(
        &s,
        "num_native_md_particle_new",
        &[
            frame.clone(),
            Value::int(0),
            measured(2., mass.clone()),
            positions,
            velocities,
        ],
    );
    let k = call(
        &s,
        "num_native_md_particle_kinetic_energy",
        std::slice::from_ref(&p),
    );
    assert_eq!(number(&k), 9.);
    unit(&k, DimExpr::Mul(vec![mass.clone(), v.clone().pow(2)]));
    let mom = call(
        &s,
        "num_native_md_particle_momentum",
        std::slice::from_ref(&p),
    );
    assert_eq!(number(&data(&mom)[0]), 6.);
    unit(&data(&mom)[0], DimExpr::Mul(vec![mass, v]));
    let cloud = call(
        &s,
        "num_native_md_cloud",
        &[frame, Value::int(1), vector(vec![p])],
    );
    let temperature = call(&s, "num_native_md_temperature", &[cloud]);
    unit(&temperature, base("Temperature"));
    assert!((number(&temperature) * 1.380649e-23 - 18.).abs() < 1e-12);
}
#[test]
fn malformed_prefix_extent_is_rejected_before_empty_nested_allocation() {
    let s = Session::new();
    let forged = Value::Tensor(TensorValue {
        shape: vec![usize::MAX, 0],
        elem: TypeExpr::Float,
        data: vec![],
    });
    let result = s.call_value(
        "num_native_tape_lower",
        &[
            forged,
            strings(&[]),
            strings(&[]),
            vector(vec![]),
            vector(vec![]),
            vector(vec![]),
        ],
        &CallCtx::new(),
    );
    assert!(result.is_err());
}
#[test]
fn physical_device_columns_and_loop_ir_execute_real_native_programs() {
    let s = Session::new();
    let t = tape(
        &s,
        "2*x",
        vec![measured(1., base("Length"))],
        measured(1., base("Time")),
    );
    let ir = call(&s, "num_native_device_kernel_from_tape", &[t]);
    let cpu = call(&s, "num_native_device_cpu", &[]);
    let xs = vector(vec![
        measured(4., base("Length")),
        measured(8., base("Length")),
    ]);
    let map = |column| {
        s.call_value(
            "num_native_device_map_physical",
            &[
                cpu.clone(),
                ir.clone(),
                vector(vec![column]),
                vector(vec![]),
            ],
            &CallCtx::new(),
        )
    };
    let out = map(xs).unwrap();
    let typed = |second| {
        Value::Tensor(TensorValue {
            shape: vec![2],
            elem: TypeExpr::Float.with_dim(base("Length")),
            data: vec![Value::Float(4.0), second],
        })
    };
    assert_eq!(map(typed(measured(8.0, base("Length")))).unwrap(), out);
    assert!(map(typed(measured(8.0, base("Time")))).is_err());
    let ys = data(&data(&out)[0]);
    assert_eq!(number(&ys[0]), 8.);
    assert_eq!(number(&ys[1]), 16.);
    unit(&ys[0], base("Time"));
    let b = call(
        &s,
        "num_native_lir_builder_new",
        &[
            Value::String("double".into()),
            strings(&["I64"]),
            strings(&["I64"]),
        ],
    );
    let p = call(&s, "num_native_lir_param", &[b.clone(), Value::int(0)]);
    let inst = call(
        &s,
        "num_native_lir_instruction",
        &[
            Value::String("IAdd".into()),
            vector(vec![p.clone(), p]),
            Value::Float(0.),
            Value::int(0),
            Value::String("".into()),
        ],
    );
    let out_reg = call(&s, "num_native_lir_ins", &[b.clone(), inst]);
    call(
        &s,
        "num_native_lir_ret",
        &[b.clone(), vector(vec![out_reg])],
    );
    let f = call(&s, "num_native_lir_finish", std::slice::from_ref(&b));
    assert!(
        s.call_value("num_native_lir_current", &[b], &CallCtx::new())
            .is_err()
    );
    let m = call(&s, "num_native_lir_module", &[vector(vec![f])]);
    let answer = call(
        &s,
        "num_native_lir_interpret",
        &[
            m,
            Value::int(0),
            vector(vec![Value::int(7)]),
            Value::int(100),
        ],
    );
    let Value::Record(a) = answer else { panic!() };
    let Value::Record(o) = a.get("outcome").unwrap() else {
        panic!()
    };
    assert_eq!(data(o.get("values").unwrap())[0], Value::int(14));
}
#[test]
fn relativistic_metric_restores_geometrized_lengths_and_physical_proper_time() {
    let s = Session::new();
    let f = call(
        &s,
        "num_native_relativity_frame",
        &[measured(2., base("Length"))],
    );
    let metric = call(
        &s,
        "num_native_relativity_schwarzschild",
        &[f, measured(1., base("Length"))],
    );
    let isco = call(
        &s,
        "num_native_relativity_isco",
        &[metric.clone(), Value::Bool(true)],
    );
    let Value::Variant(v) = isco else { panic!() };
    assert_eq!(v.tag, "Some");
    assert!((number(&v.payload[0]) - 6.).abs() < 1e-5);
    unit(&v.payload[0], base("Length"));
    let g = call(
        &s,
        "num_native_relativity_geodesic_circular",
        &[metric, measured(10., base("Length")), Value::Bool(true)],
    );
    let state = call(
        &s,
        "num_native_relativity_geodesic_step",
        &[g, measured(1e-10, base("Time")), Value::Bool(true)],
    );
    let Value::Record(state) = state else {
        panic!()
    };
    assert!((number(state.get("radius").unwrap()) - 10.).abs() < 1e-8);
    unit(state.get("proper_time").unwrap(), base("Time"));
}

#[test]
fn relativistic_horizon_restores_a_finite_length_across_normalized_overflow() {
    let session = Session::new();
    let length = base("Length");
    let frame = call(
        &session,
        "num_native_relativity_frame",
        &[measured(0.5, length.clone())],
    );
    let mass = f64::MAX * 0.375;
    let metric = call(
        &session,
        "num_native_relativity_schwarzschild",
        &[frame.clone(), measured(mass, length.clone())],
    );
    let horizon = call(&session, "num_native_relativity_horizon", &[metric]);
    let Value::Variant(horizon) = horizon else {
        panic!("expected an optional horizon")
    };
    assert_eq!(horizon.tag, "Some");
    unit(&horizon.payload[0], length.clone());
    assert_eq!(number(&horizon.payload[0]), 2.0 * mass);

    let parameters = |spin: f64, charge: f64| {
        call(
            &session,
            "num_native_relativity_metric",
            &[
                frame.clone(),
                measured(mass, length.clone()),
                measured(spin, length.clone()),
                measured(charge, length.clone()),
            ],
        )
    };
    let horizon = call(
        &session,
        "num_native_relativity_horizon",
        &[parameters(0.4 * mass, 0.6 * mass)],
    );
    let Value::Variant(horizon) = horizon else {
        panic!("expected an optional charged horizon")
    };
    assert_eq!(horizon.tag, "Some");
    let expected = mass * (1.0 + 0.48_f64.sqrt());
    assert!((number(&horizon.payload[0]) / expected - 1.0).abs() < 8.0 * f64::EPSILON);

    let naked = call(
        &session,
        "num_native_relativity_horizon",
        &[parameters(mass, f64::from_bits(1))],
    );
    let Value::Variant(naked) = naked else {
        panic!("expected an optional horizon")
    };
    assert_eq!(naked.tag, "None");
}

#[test]
fn relativistic_geodesic_encoding_failure_keeps_the_retained_state() {
    const C: f64 = 299_792_458.0;
    const REFERENCE: f64 = 1.0e307;
    let session = Session::new();
    let length = base("Length");
    let time = base("Time");
    let frame = call(
        &session,
        "num_native_relativity_frame",
        &[measured(REFERENCE, length.clone())],
    );
    let metric = call(
        &session,
        "num_native_relativity_schwarzschild",
        &[frame, measured(1.0e294, length.clone())],
    );
    let new = || {
        call(
            &session,
            "num_native_relativity_geodesic_new",
            &[
                metric.clone(),
                measured(REFERENCE, length.clone()),
                Value::Float(std::f64::consts::FRAC_PI_2),
                Value::Float(0.0),
                Value::Float(1.0),
                measured(0.0, length.clone()),
                measured(0.0, length.clone()),
                Value::Float(1.0),
                Value::Float(0.0),
            ],
        )
    };
    let step = "num_native_relativity_geodesic_step";
    let integrate = "num_native_relativity_geodesic_integrate";
    let state = "num_native_relativity_geodesic_state";
    let h = measured(32.0 * (REFERENCE / C), time.clone());
    let mut rollback = [false; 2];
    for (index, name) in [step, integrate].into_iter().enumerate() {
        let geodesic = new();
        let before = call(&session, state, std::slice::from_ref(&geodesic));
        let args = if name == step {
            vec![geodesic.clone(), h.clone(), Value::Bool(false)]
        } else {
            vec![
                geodesic.clone(),
                h.clone(),
                Value::int(1),
                Value::Bool(false),
                measured(0.2 * REFERENCE, length.clone()),
                measured(f64::MAX, length.clone()),
            ]
        };
        let error = session
            .call_value(name, &args, &CallCtx::new())
            .unwrap_err();
        assert!(
            error.to_string().contains("overflows Float"),
            "{name}: {error}"
        );
        rollback[index] = session
            .call_value(state, std::slice::from_ref(&geodesic), &CallCtx::new())
            .is_ok_and(|after| after == before);
        if !rollback[index] {
            continue;
        }
        let after = call(
            &session,
            step,
            &[
                geodesic,
                measured(0.03125 * (REFERENCE / C), time.clone()),
                Value::Bool(false),
            ],
        );
        let Value::Record(after) = after else {
            panic!("expected measured state")
        };
        let radius = after.get("radius").unwrap();
        assert!((number(radius) / REFERENCE - 1.03125).abs() < 1e-11);
        unit(radius, length.clone());
        unit(after.get("proper_time").unwrap(), time.clone());
    }
    assert_eq!(
        rollback,
        [true, true],
        "step and integrate both preserve the handle"
    );
}

#[test]
fn native_geodesic_effective_potential_keeps_coherent_equal_scale_ratios() {
    let session = Session::new();
    let length = base("Length");
    let time = base("Time");
    let angular_rate = time.clone().pow(-1);
    let radial_velocity = DimExpr::Mul(vec![length.clone(), angular_rate.clone()]);
    let frame = call(
        &session,
        "num_native_relativity_frame",
        &[measured(1.0, length.clone())],
    );
    let name = "num_native_geodesic_model_effective_potential";
    assert_eq!(
        session.manifest().op(name).unwrap().overloads[0].ret,
        TypeExpr::Float
    );
    let parameters = |mass: f64, radius: f64| {
        call(
            &session,
            "num_native_geodesic_model_parameters",
            &[
                frame.clone(),
                measured(mass, length.clone()),
                measured(radius, length.clone()),
                measured(0.0, radial_velocity.clone()),
                Value::Float(0.0),
                measured(0.0, angular_rate.clone()),
                measured(0.0, time.clone()),
                measured(1e-9, time.clone()),
            ],
        )
    };
    for scale in [1.0, 1e-200, 1e200] {
        let parameters = parameters(0.1 * scale, scale);
        for (angular, expected) in [(0.3, -0.064), (-0.3, -0.064), (0.0, -0.1)] {
            let value = call(
                &session,
                name,
                &[
                    parameters.clone(),
                    measured(scale, length.clone()),
                    measured(angular * scale, length.clone()),
                ],
            );
            assert!(
                (number(&value) - expected).abs() < 1e-15,
                "scale {scale}: {value:?}"
            );
        }
    }
    // r/2 - M would be half the minimum positive Float and round to zero.
    let subnormal = f64::from_bits(1);
    let radius = 5.0 * subnormal;
    let result = call(
        &session,
        name,
        &[
            parameters(2.0 * subnormal, radius),
            measured(radius, length.clone()),
            measured(radius, length.clone()),
        ],
    );
    assert!((number(&result) + 0.3).abs() < 1e-15);
    let horizon = call(
        &session,
        name,
        &[
            parameters(1.0, 2.0),
            measured(2.0, length.clone()),
            measured(f64::MAX, length),
        ],
    );
    assert_eq!(horizon, Value::Float(-0.5));
}

#[test]
fn numerical_shadows_keep_order_shape_and_geometrized_length_units() {
    let session = Session::new();
    let length = base("Length");
    let frame = call(
        &session,
        "num_native_relativity_frame",
        &[measured(2.0, length.clone())],
    );
    let arguments = [
        frame,
        measured(1.0, length.clone()),
        measured(0.0, length.clone()),
        Value::Float(1.2),
    ];
    let observer = measured(1e4, length.clone());
    let expected_radius = 27.0_f64.sqrt();
    let mut args = arguments.to_vec();
    args.push(Value::int(6));
    let analytic = call(&session, "num_native_relativity_shadow_analytic", &args);
    let Value::Tensor(analytic) = analytic else {
        panic!("expected image-plane coordinates")
    };
    assert_eq!(analytic.shape, [12, 2]);
    args.pop();
    args.extend([observer, Value::int(6)]);
    let trace = call(&session, "num_native_relativity_trace_shadow", &args);
    let Value::Tensor(trace) = trace else {
        panic!("expected image-plane coordinates")
    };
    assert_eq!(trace.shape, [6, 2]);
    for (data, tolerance) in [(analytic.data, 1e-12), (trace.data, 1e-6)] {
        for [x, y] in data.as_chunks::<2>().0 {
            unit(x, length.clone());
            unit(y, length.clone());
            assert!((number(x).hypot(number(y)) - expected_radius).abs() < tolerance);
        }
    }
    args[5] = Value::int(65);
    assert!(
        session
            .call_value("num_native_relativity_trace_shadow", &args, &CallCtx::new())
            .is_err()
    );
}

#[test]
fn native_equatorial_model_keeps_kerr_horizon_checks_and_sampled_units() {
    const LIGHT_SPEED: f64 = 299_792_458.0;
    const REFERENCE: f64 = 1_000.0;
    let s = Session::new();
    let ctx = CallCtx::new();
    let length = base("Length");
    let time = base("Time");
    let rate = time.clone().pow(-1);
    let velocity = DimExpr::Mul(vec![length.clone(), rate.clone()]);
    let frame = call(
        &s,
        "num_native_relativity_frame",
        &[measured(REFERENCE, length.clone())],
    );
    let parameter = |r: f64, rd: f64, phid: f64| {
        call(
            &s,
            "num_native_geodesic_model_parameters",
            &[
                frame.clone(),
                measured(REFERENCE, length.clone()),
                measured(r * REFERENCE, length.clone()),
                measured(rd * LIGHT_SPEED, velocity.clone()),
                Value::Float(0.0),
                measured(phid * LIGHT_SPEED / REFERENCE, rate.clone()),
                measured(0.25 * REFERENCE / LIGHT_SPEED, time.clone()),
                measured(0.0078125 * REFERENCE / LIGHT_SPEED, time.clone()),
            ],
        )
    };

    let (r, a) = (1.9_f64, 0.99);
    let root = r.sqrt();
    let r32 = r * root;
    let denominator = root * (root * (r32 - 3.0 * root + 2.0 * a)).sqrt();
    let params = parameter(r, 0.0, denominator.recip());
    let spin = measured(a * REFERENCE, length.clone());
    let charge = measured(0.0, length.clone());
    let state_of = |handle| call(&s, "num_native_relativity_geodesic_state", &[handle]);
    let handle = call(
        &s,
        "num_native_geodesic_model_equatorial_timelike",
        &[
            params.clone(),
            spin.clone(),
            charge.clone(),
            Value::Bool(false),
        ],
    );
    let Value::Record(state) = state_of(handle) else {
        panic!("expected geodesic state")
    };
    let angular = state.get("angular_momentum_per_rest_mass").unwrap();
    unit(angular, length.clone());
    let expected_angular = REFERENCE * (r * r - 2.0 * a * root + a * a) / denominator;
    assert!((number(angular) - expected_angular).abs() < 1e-9);
    assert!(
        (number(state.get("energy_per_rest_mass").unwrap()) - (r32 - 2.0 * root + a) / denominator)
            .abs()
            < 1e-13
    );

    let fixed = call(
        &s,
        "num_native_geodesic_model_run_kerr_newman_geodesic",
        &[
            params.clone(),
            spin.clone(),
            charge.clone(),
            Value::int(150),
        ],
    );
    let sampled_args = |stride, budget| {
        [
            params.clone(),
            spin.clone(),
            charge.clone(),
            Value::Bool(false),
            Value::int(stride),
            Value::int(budget),
        ]
    };
    let sampled = call(
        &s,
        "num_native_geodesic_model_run_kerr_newman_geodesic_sampled",
        &sampled_args(7, 150),
    );
    let (Value::Tensor(fixed), Value::Tensor(sampled)) = (fixed, sampled) else {
        panic!("expected measured Cartesian tensors")
    };
    let rows = fixed.shape[0];
    assert!(rows > 28);
    assert_eq!(sampled.shape, vec![(rows - 2) / 7 + 2, 2]);
    let expected = fixed
        .data
        .as_chunks::<2>()
        .0
        .iter()
        .enumerate()
        .filter(|(i, _)| i % 7 == 0 || i + 1 == rows)
        .flat_map(|(_, row)| row.iter().cloned())
        .collect::<Vec<_>>();
    assert_eq!(sampled.data, expected);
    unit(&sampled.data[0], length.clone());
    for args in [sampled_args(0, 150), sampled_args(7, 1)] {
        assert!(
            s.call_value(
                "num_native_geodesic_model_run_kerr_newman_geodesic_sampled",
                &args,
                &ctx
            )
            .is_err()
        );
    }
    for name in [
        "num_native_geodesic_model_run_geodesic_simulation",
        "num_native_geodesic_model_run_geodesic_simulation_rk45",
    ] {
        assert!(
            s.call_value(name, &[params.clone(), Value::int(150)], &ctx)
                .is_err()
        );
    }
    let falling = parameter(1.5, 0.1, 2.0);
    let g = call(
        &s,
        "num_native_geodesic_model_equatorial_timelike",
        &[
            falling,
            measured(REFERENCE, length),
            charge,
            Value::Bool(true),
        ],
    );
    let Value::Record(negative) = state_of(g) else {
        panic!("expected geodesic state")
    };
    assert!(number(negative.get("energy_per_rest_mass").unwrap()) < 0.0);
}
