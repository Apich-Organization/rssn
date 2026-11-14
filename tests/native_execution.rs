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
    let k = call(&s, "num_native_md_particle_kinetic_energy", &[p.clone()]);
    assert_eq!(number(&k), 9.);
    unit(&k, DimExpr::Mul(vec![mass.clone(), v.clone().pow(2)]));
    let mom = call(&s, "num_native_md_particle_momentum", &[p.clone()]);
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
    let out = call(
        &s,
        "num_native_device_map_physical",
        &[cpu, ir, vector(vec![xs]), vector(vec![])],
    );
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
    let f = call(&s, "num_native_lir_finish", &[b.clone()]);
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
