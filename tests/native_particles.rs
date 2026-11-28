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
fn legacy_sph_native_combined_update_matches_its_two_component_phases() {
    let session = rssn::Session::new();
    let make = || {
        call(
            &session,
            "num_native_mm_system_new",
            &[
                floats(vec![2, 2], &[0.75, 1.0, 1.25, 1.0], DimExpr::Dimless),
                Value::Float(1.0),
                Value::Float(1.0),
                floats(vec![2], &[0.0, 0.0], DimExpr::Dimless),
                Value::Float(0.1),
                Value::Float(1.0),
                Value::Float(1.0),
                floats(vec![2], &[2.0, 2.0], DimExpr::Dimless),
            ],
        )
    };
    let separate = make();
    let combined = make();
    for name in [
        "num_native_mm_compute_density_pressure",
        "num_native_mm_compute_forces",
    ] {
        assert_eq!(
            call(&session, name, std::slice::from_ref(&separate)),
            Value::Bool(true)
        );
    }
    assert_eq!(
        call(
            &session,
            "num_native_mm_compute_density_and_forces",
            std::slice::from_ref(&combined)
        ),
        Value::Bool(true)
    );
    let separate = call(&session, "num_native_mm_state", &[separate]);
    let combined = call(&session, "num_native_mm_state", &[combined]);
    assert_eq!(separate, combined);
    let Value::Tensor(position) = field(&combined, "positions") else {
        panic!("position tensor")
    };
    assert_eq!(position.shape, [2, 2]);
    assert_eq!(
        position.data,
        vec![
            Value::Float(0.75),
            Value::Float(1.0),
            Value::Float(1.25),
            Value::Float(1.0)
        ]
    );
    let Value::Tensor(forces) = field(&combined, "forces") else {
        panic!("force tensor")
    };
    assert_eq!(forces.shape, [2, 2]);
    assert!(number(&forces.data[0]).abs() > 0.0);
    assert_eq!(number(&forces.data[0]), -number(&forces.data[2]));
}

#[test]
fn native_rdf_merges_distinct_frames_and_preserves_state_on_incompatible_bins() {
    let session = rssn::Session::new();
    let charge = DimExpr::Mul(vec![dim("Current"), dim("Time")]);
    let make_frame = |length_ref| {
        call(
            &session,
            "num_native_md_frame",
            &[
                quantity(1.0, dim("Mass")),
                quantity(length_ref, dim("Length")),
                quantity(1.0, dim("Time")),
                quantity(1.0, charge.clone()),
            ],
        )
    };
    let frame = make_frame(4.0);
    let make_rdf = |frame: &Value, bins| {
        call(
            &session,
            "num_native_engine_rdf_new_d1",
            &[
                frame.clone(),
                Value::int(bins),
                quantity(1.0, dim("Length")),
            ],
        )
    };
    let make_system = |length, positions: &[f64]| {
        let cell = call(
            &session,
            "num_native_engine_box_periodic_d1",
            &[frame.clone(), floats(vec![1], &[length], dim("Length"))],
        );
        call(
            &session,
            "num_native_engine_system_new_d1",
            &[
                cell,
                floats(vec![positions.len(), 1], positions, dim("Length")),
            ],
        )
    };
    let first_system = make_system(4.0, &[0.0, 0.5]);
    let second_system = make_system(8.0, &[0.0, 0.5, 3.0]);
    let first = make_rdf(&frame, 1);
    let second = make_rdf(&frame, 1);
    for (rdf, system) in [(&first, &first_system), (&second, &second_system)] {
        assert_eq!(
            call(
                &session,
                "num_native_engine_rdf_add_frame_d1",
                &[rdf.clone(), system.clone()]
            ),
            Value::Bool(true)
        );
    }
    let result = |rdf: &Value| {
        call(
            &session,
            "num_native_engine_rdf_result_d1",
            std::slice::from_ref(rdf),
        )
    };
    let ratio = |result: &Value| {
        let Value::Tensor(values) = field(result, "density_ratio") else {
            panic!("density-ratio tensor");
        };
        assert_eq!(values.shape, [1]);
        assert_eq!(values.elem, TypeExpr::Float);
        number(values.data.first().unwrap())
    };
    let merge = |a: &Value, b: &Value| {
        session.call_value(
            "num_native_engine_rdf_merge_d1",
            &[a.clone(), b.clone()],
            &CallCtx::new(),
        )
    };
    let original_source = result(&second);
    assert_eq!(merge(&first, &second).unwrap(), Value::Bool(true));
    assert_eq!(result(&second), original_source);
    let combined = result(&first);
    let Value::Tensor(radius) = field(&combined, "radius") else {
        panic!("radius tensor");
    };
    assert_eq!(radius.elem, TypeExpr::Float.with_dim(dim("Length")));
    assert_eq!(number(radius.data.first().unwrap()), 0.5);
    assert!((ratio(&combined) - 16.0 / 17.0).abs() < 1e-14);
    assert!(merge(&first, &make_rdf(&frame, 2)).is_err());
    assert!(merge(&first, &make_rdf(&make_frame(2.0), 1)).is_err());
    assert_eq!(result(&first), combined);

    assert_eq!(merge(&first, &first).unwrap(), Value::Bool(true));
    assert_eq!(result(&first), combined);
    assert_eq!(merge(&first, &second).unwrap(), Value::Bool(true));
    assert!((ratio(&result(&first)) - 40.0 / 43.0).abs() < 1e-14);
}

#[test]
fn measured_morse_radial_resolves_near_equilibrium_energy_and_force() {
    let session = rssn::Session::new();
    let energy_dim = DimExpr::Mul(vec![dim("Length").pow(2), dim("Mass"), dim("Time").pow(-2)]);
    let force_dim = DimExpr::Mul(vec![dim("Length"), dim("Mass"), dim("Time").pow(-2)]);
    let frame = call(
        &session,
        "num_native_md_frame",
        &[
            quantity(2.0, dim("Mass")),
            quantity(4.0, dim("Length")),
            quantity(2.0, dim("Time")),
            quantity(1.0, DimExpr::Mul(vec![dim("Current"), dim("Time")])),
        ],
    );
    let distance = 1.0_f64.next_up();
    let mut args = [
        frame,
        quantity(distance, dim("Length")),
        quantity(2.0, energy_dim.clone()),
        quantity(1e-7, dim("Length").pow(-1)),
        quantity(1.0, dim("Length")),
    ];
    let result = call(&session, "num_native_md_morse_radial", &args);
    let delta = (distance - 1.0) * 1e-7;
    for (name, dim, expected) in [
        ("energy", energy_dim, 2.0 * delta * delta),
        ("radial_force", force_dim, -4e-7 * delta),
    ] {
        let Value::Dyn(value) = field(&result, name) else {
            panic!("{name} should carry a physical unit");
        };
        assert_eq!(value.dim, dim);
        assert!((number(&value.value) / expected - 1.0).abs() < 1e-13);
    }
    args[1] = quantity(distance, dim("Time"));
    assert!(
        session
            .call_value("num_native_md_morse_radial", &args, &CallCtx::new())
            .is_err()
    );
}

#[test]
fn sse_native_budget_counts_free_spins_and_invalidates_partial_updates() {
    let session = rssn::Session::new();
    let energy = DimExpr::Mul(vec![dim("Mass"), dim("Length").pow(2), dim("Time").pow(-2)]);
    let frame = call(
        &session,
        "num_native_md_frame",
        &[
            quantity(2.0, dim("Mass")),
            quantity(1.0, dim("Length")),
            quantity(1.0, dim("Time")),
            quantity(1.0, DimExpr::Mul(vec![dim("Current"), dim("Time")])),
        ],
    );
    let make = || {
        call(
            &session,
            "num_native_sse_new",
            &[
                frame.clone(),
                Value::int(2),
                tensor(
                    vec![1, 2],
                    vec![Value::int(0), Value::int(1)],
                    TypeExpr::Int,
                ),
                quantity(0.0, energy.clone()),
                quantity(1e-12, energy.clone().pow(-1)),
                Value::int(31),
            ],
        )
    };
    let sweep = |state: &Value, budget: usize| {
        session.call_value(
            "num_native_sse_sweep",
            &[state.clone(), Value::int(budget)],
            &CallCtx::new(),
        )
    };
    let sample = |state: &Value| {
        session.call_value(
            "num_native_sse_sample",
            std::slice::from_ref(state),
            &CallCtx::new(),
        )
    };
    let short = make();
    assert!(sweep(&short, 0).is_err());
    assert!(sample(&short).is_ok());
    assert!(sweep(&short, 49).is_err());
    assert!(sample(&short).is_err());
    assert!(sweep(&short, 100).is_err());

    let exact = make();
    let spacious = make();
    assert_eq!(sweep(&exact, 50).unwrap(), Value::Bool(true));
    assert_eq!(sweep(&spacious, 100).unwrap(), Value::Bool(true));
    let snapshot = sample(&exact).unwrap();
    assert_eq!(field(&snapshot, "operators"), &Value::int(0));
    assert_eq!(snapshot, sample(&spacious).unwrap());
    let shift = call(&session, "num_native_sse_constant", &[exact]);
    assert_eq!(number(&shift), 0.1);
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(energy), &shift.type_expr())
            .accepted()
    );
}

#[test]
fn prepared_fmm_keeps_source_dimensions_and_original_rows() {
    use rssn::kernel::value::RecordValue;
    let s = rssn::Session::new();
    let length = dim("Length");
    let charge = DimExpr::Mul(vec![dim("Current"), dim("Time")]);
    let frame = call(
        &s,
        "num_native_md_frame",
        &[
            quantity(2.0, dim("Mass")),
            quantity(4.0, length.clone()),
            quantity(2.0, dim("Time")),
            quantity(1.0, charge.clone()),
        ],
    );
    let options = Value::Record(RecordValue {
        fields: vec![
            ("order".into(), Value::int(6)),
            ("leaf_size".into(), Value::int(2)),
            ("theta".into(), Value::Float(0.5)),
            ("softening".into(), quantity(0.0, length.clone())),
        ],
    });
    let positions = floats(vec![2, 2], &[0.0, 0.0, 2.0, 0.0], length.clone());
    let plan = call(&s, "num_native_fmm_prepare2d", &[frame, positions, options]);
    assert_eq!(
        call(
            &s,
            "num_native_fmm_plan_count2d",
            std::slice::from_ref(&plan)
        ),
        Value::int(2)
    );
    let name = "num_native_fmm_evaluate2d";
    let reference = quantity(2.0, charge.clone());
    let expected_dimension =
        TypeExpr::Float.with_dim(DimExpr::Mul(vec![charge.clone(), length.pow(-1)]));
    for (strength, expected) in [([2.0, 4.0], [2.0, -1.0]), ([-6.0, 2.0], [1.0, 3.0])] {
        let out = call(
            &s,
            name,
            &[
                plan.clone(),
                floats(vec![2], &strength, charge.clone()),
                reference.clone(),
            ],
        );
        let Value::Tensor(out) = out else {
            panic!("FMM field tensor")
        };
        assert_eq!(out.shape, [2, 2]);
        assert!(rssn::kernel::overload::matches_type(&expected_dimension, &out.elem).accepted());
        assert_eq!(
            out.data.iter().map(number).collect::<Vec<_>>(),
            [expected[0], 0.0, expected[1], 0.0]
        );
    }
    assert!(
        s.call_value(
            name,
            &[plan, floats(vec![1], &[2.0], charge), reference],
            &CallCtx::new()
        )
        .is_err()
    );
}

#[test]
fn prepared_sph_density_supports_independent_fields_at_fixed_positions() {
    let s = rssn::Session::new();
    let mass = dim("Mass");
    let length = dim("Length");
    let time = dim("Time");
    let charge = DimExpr::Mul(vec![dim("Current"), time.clone()]);
    let internal_energy = DimExpr::Mul(vec![length.clone().pow(2), time.clone().pow(-2)]);
    let frame = call(
        &s,
        "num_native_md_frame",
        &[
            quantity(2.0, mass.clone()),
            quantity(4.0, length.clone()),
            quantity(3.0, time),
            quantity(1.0, charge),
        ],
    );
    let masses = floats(vec![2], &[2.0, 3.0], mass.clone());
    let gas = call(
        &s,
        "num_native_sph_new_d1",
        &[
            frame,
            floats(vec![2, 1], &[1.0, 5.0], length.clone()),
            masses.clone(),
            floats(vec![2], &[1.0, 1.0], internal_energy),
            floats(vec![1], &[0.0], length.clone()),
            floats(vec![1], &[10.0], length.clone()),
            tensor(
                vec![1],
                vec![Value::String("Open".into())],
                TypeExpr::String,
            ),
        ],
    );
    let plan = call(
        &s,
        "num_native_sph_density_plan_new_d1",
        &[gas, quantity(1.0, length.clone())],
    );
    let smoothing = floats(vec![2], &[0.5, 1.0], length.clone());
    let kernel = Value::String("CubicSpline".into());
    let evaluate = "num_native_sph_density_plan_evaluate_d1";
    let density =
        TypeExpr::Float.with_dim(DimExpr::Mul(vec![mass.clone(), length.clone().pow(-1)]));
    for (masses, expected) in [
        (masses, [8.0 / 3.0, 2.0]),
        (
            floats(vec![2], &[4.0, 6.0], mass.clone()),
            [16.0 / 3.0, 4.0],
        ),
    ] {
        let value = call(
            &s,
            evaluate,
            &[plan.clone(), masses, smoothing.clone(), kernel.clone()],
        );
        let Value::Tensor(rho) = value else {
            panic!("expected measured mass density");
        };
        assert_eq!(rho.shape, [2]);
        assert!(rssn::kernel::overload::matches_type(&density, &rho.elem).accepted());
        for (actual, expected) in rho.data.iter().zip(expected) {
            assert!((number(actual) - expected).abs() < 1e-13);
        }
    }
    let excessive = floats(vec![2], &[0.5, 1.1], length);
    let masses = floats(vec![2], &[2.0, 3.0], mass);
    assert!(
        s.call_value(
            evaluate,
            &[plan, masses, excessive, kernel],
            &CallCtx::new()
        )
        .is_err()
    );
}

#[test]
fn prepared_reciprocal_ewald_reuses_fixed_cells_and_restores_physical_fields() {
    let s = rssn::Session::new();
    let charge = DimExpr::Mul(vec![dim("Current"), dim("Time")]);
    let length = dim("Length");
    let time = dim("Time");
    let mass = dim("Mass");
    let energy = DimExpr::Mul(vec![
        mass.clone(),
        length.clone().pow(2),
        time.clone().pow(-2),
    ]);
    let coupling = DimExpr::Mul(vec![energy.clone(), length.clone(), charge.clone().pow(-2)]);
    let velocity = DimExpr::Mul(vec![length.clone(), time.clone().pow(-1)]);
    let f = call(
        &s,
        "num_native_md_frame",
        &[
            quantity(2.0, mass.clone()),
            quantity(3.0, length.clone()),
            quantity(4.0, time),
            quantity(5.0, charge.clone()),
        ],
    );
    let system = |side, offset| {
        let cell = call(
            &s,
            "num_native_engine_box_periodic_d3",
            &[f.clone(), floats(vec![3], &[side; 3], length.clone())],
        );
        let system = call(
            &s,
            "num_native_engine_system_new_d3",
            &[
                cell,
                floats(
                    vec![2, 3],
                    &[1.0, 1.0, 1.0, offset, 1.0, 1.0],
                    length.clone(),
                ),
            ],
        );
        call(
            &s,
            "num_native_engine_system_set_state_d3",
            &[
                system.clone(),
                floats(vec![2, 3], &[0.0; 6], velocity.clone()),
                floats(vec![2], &[1.0, 1.0], mass.clone()),
                floats(vec![2], &[1.0, -1.0], charge.clone()),
            ],
        );
        system
    };
    let first = system(12.0, 3.0);
    let second = system(12.0, 4.0);
    let different_cell = system(13.0, 3.0);
    let alpha = quantity(0.75, length.pow(-1));
    let constant = quantity(0.5, coupling);
    let cutoff = Value::int(2);
    let plan = call(
        &s,
        "num_native_ewald_reciprocal_plan_new",
        &[first.clone(), alpha.clone(), cutoff.clone()],
    );
    let evaluate = "num_native_ewald_reciprocal_plan_evaluate";
    for system in [first, second] {
        let expected = call(
            &s,
            "num_native_ewald_ewald_reciprocal_cell",
            &[
                system.clone(),
                alpha.clone(),
                constant.clone(),
                cutoff.clone(),
            ],
        );
        let result = call(&s, evaluate, &[plan.clone(), system, constant.clone()]);
        assert_eq!(result, expected);
        let result = field(&result, "energy");
        assert_ne!(number(result), 0.0);
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Float.with_dim(energy.clone()),
                &result.type_expr()
            )
            .accepted()
        );
    }
    assert!(
        s.call_value(evaluate, &[plan, different_cell, constant], &CallCtx::new())
            .is_err()
    );
}

#[test]
fn native_rigid_md_checks_partition_and_preserves_measured_free_motion() {
    use rssn::kernel::value::RecordValue;

    let session = rssn::Session::new();
    let manifest = session.manifest();
    let constructor = manifest
        .ops
        .iter()
        .find(|op| op.name == "num_native_rigid_md_new")
        .unwrap();
    assert_eq!(constructor.category, "numerical_vector");
    assert!(constructor.doc.contains("CSR"));
    assert!(constructor.overloads.iter().any(|signature| {
        signature.params.len() == 5
            && signature.params[0] == TypeExpr::Opaque("rssn/native/sim/engine-system/d3/v1".into())
            && signature.params[1]
                == TypeExpr::Opaque("rssn/native/sim/engine-force-field/d3/v1".into())
            && signature.ret == TypeExpr::Opaque("rssn/native/sim/rigid-md/v1".into())
    }));
    let length = dim("Length");
    let mass = dim("Mass");
    let time = dim("Time");
    let charge = DimExpr::Mul(vec![dim("Current"), time.clone()]);
    let inverse_time = time.clone().pow(-1);
    let velocity = DimExpr::Mul(vec![length.clone(), inverse_time.clone()]);
    let momentum = DimExpr::Mul(vec![mass.clone(), length.clone(), inverse_time]);
    let frame = call(
        &session,
        "num_native_md_frame",
        &[
            quantity(2.0, mass.clone()),
            quantity(4.0, length.clone()),
            quantity(2.0, time.clone()),
            quantity(1.0, charge.clone()),
        ],
    );
    let cell = call(
        &session,
        "num_native_engine_box_periodic_d3",
        &[frame.clone(), floats(vec![3], &[16.0; 3], length.clone())],
    );
    let system = call(
        &session,
        "num_native_engine_system_new_d3",
        &[
            cell,
            floats(vec![2, 3], &[4.0, 4.0, 4.0, 6.0, 4.0, 4.0], length.clone()),
        ],
    );
    assert_eq!(
        call(
            &session,
            "num_native_engine_system_set_state_d3",
            &[
                system.clone(),
                floats(vec![2, 3], &[0.0, 0.5, 0.0, 0.0, 0.5, 0.0], velocity),
                floats(vec![2], &[2.0; 2], mass),
                floats(vec![2], &[0.0; 2], charge),
            ],
        ),
        Value::Bool(true)
    );
    let force_field = call(
        &session,
        "num_native_engine_forcefield_new_d3",
        &[
            frame,
            tensor(
                vec![0],
                vec![],
                TypeExpr::Opaque("rssn/native/sim/engine-pair/v1".into()),
            ),
            quantity(2.0, length.clone()),
        ],
    );
    let arguments = |indices: &[i64]| {
        [
            system.clone(),
            force_field.clone(),
            quantity(0.5, length.clone()),
            quantity(0.25, time.clone()),
            Value::Record(RecordValue {
                fields: vec![
                    ("offsets".into(), ints(&[0, 2])),
                    ("indices".into(), ints(indices)),
                ],
            }),
        ]
    };
    let check_tensor = |value: &Value, shape: &[usize], numbers: &[f64], dimension| {
        let Value::Tensor(tensor) = value else {
            panic!("expected measured tensor");
        };
        assert_eq!(tensor.shape, shape);
        assert_eq!(tensor.data.iter().map(number).collect::<Vec<_>>(), numbers);
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Float.with_dim(dimension),
                &tensor.elem,
            )
            .accepted()
        );
    };
    for invalid in [[0, 0], [0, 2]] {
        assert!(
            session
                .call_value(
                    "num_native_rigid_md_new",
                    &arguments(&invalid),
                    &CallCtx::new()
                )
                .is_err()
        );
        assert_eq!(
            call(
                &session,
                "num_native_engine_forcefield_validate_d3",
                std::slice::from_ref(&force_field)
            ),
            Value::Bool(true)
        );
    }
    let handle = call(&session, "num_native_rigid_md_new", &arguments(&[1, 0]));
    let object = std::slice::from_ref(&handle);
    let before = call(&session, "num_native_rigid_md_state", object);
    assert!(
        session
            .call_value(
                "num_native_rigid_md_run",
                &[handle.clone(), Value::int(65_537)],
                &CallCtx::new()
            )
            .is_err()
    );
    assert_eq!(call(&session, "num_native_rigid_md_state", object), before);
    assert_eq!(
        number(&call(&session, "num_native_rigid_md_dof", object)),
        2.0
    );
    check_tensor(
        &call(&session, "num_native_rigid_md_momentum", object),
        &[3],
        &[0.0, 2.0, 0.0],
        momentum,
    );
    assert_eq!(
        number(&call(&session, "num_native_rigid_md_total_energy", object)),
        0.5
    );
    assert_eq!(
        call(
            &session,
            "num_native_rigid_md_run",
            &[handle.clone(), Value::int(4)]
        ),
        Value::Bool(true)
    );
    let after = call(&session, "num_native_rigid_md_state", object);
    assert_eq!(field(&after, "time"), &quantity(1.0, time));
    check_tensor(
        field(&after, "positions"),
        &[2, 3],
        &[4.0, 4.5, 4.0, 6.0, 4.5, 4.0],
        length,
    );
    assert!(
        number(&call(
            &session,
            "num_native_rigid_md_max_geometry_error",
            object
        )) < 1e-14
    );
    assert!(
        rssn::Session::new()
            .call_value("num_native_rigid_md_state", object, &CallCtx::new())
            .is_err()
    );
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
fn native_md_sparse_neighbours_preserve_a_periodic_interaction_in_a_large_box() {
    use rssn::kernel::value::RecordValue;

    let session = rssn::Session::new();
    let frame = call(&session, "num_native_md_normalized_frame", &[]);
    let parameters = Value::Record(RecordValue {
        fields: vec![
            ("depth".into(), Value::Float(1.0)),
            ("alpha".into(), Value::Float(1.0)),
            ("r0".into(), Value::Float(0.25)),
        ],
    });
    let term = call(
        &session,
        "num_native_engine_pair_new",
        &[frame.clone(), Value::String("Morse".into()), parameters],
    );
    let field_handle = call(
        &session,
        "num_native_engine_forcefield_new_d3",
        &[
            frame.clone(),
            tensor(
                vec![1],
                vec![term],
                TypeExpr::Opaque("rssn/native/sim/engine-pair/v1".into()),
            ),
            Value::Float(1.0),
        ],
    );
    let run = |side| {
        let periodic = call(
            &session,
            "num_native_engine_box_periodic_d3",
            &[frame.clone(), floats(vec![3], &[side; 3], DimExpr::Dimless)],
        );
        let sys = call(
            &session,
            "num_native_engine_system_new_d3",
            &[
                periodic,
                floats(
                    vec![2, 3],
                    &[side - 0.125, 0.0, 0.0, 0.125, 0.0, 0.0],
                    DimExpr::Dimless,
                ),
            ],
        );
        let list = call(
            &session,
            "num_native_engine_neighbours_new_d3",
            &[frame.clone(), Value::Float(1.0), Value::Float(0.0)],
        );
        assert_eq!(
            call(
                &session,
                "num_native_engine_neighbours_build_d3",
                &[list.clone(), sys.clone()]
            ),
            Value::Bool(true)
        );
        let forces = call(
            &session,
            "num_native_engine_compute_forces_d3",
            &[sys, field_handle.clone(), list],
        );
        call(&session, "num_native_engine_forces_data_d3", &[forces])
    };
    let large = run(1e7);
    assert_eq!(number(field(&large, "energy")), -1.0);
    assert_eq!(large, run(10.0));
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
fn native_gpe_returns_a_measured_normalized_density_without_evolution() {
    let session = rssn::Session::new();
    let time = dim("Time");
    let length = dim("Length");
    let mass = dim("Mass");
    let args = |n: [i64; 2], dt, steps| {
        [
            ints(&n),
            floats(vec![2], &[2.0, 4.0], length.clone()),
            quantity(1.0, mass.clone()),
            quantity(2.0, length.clone()),
            quantity(dt, time.clone()),
            Value::int(steps),
            quantity(
                0.0,
                DimExpr::Mul(vec![
                    mass.clone(),
                    length.clone().pow(4),
                    time.clone().pow(-2),
                ]),
            ),
            quantity(0.0, DimExpr::Mul(vec![mass.clone(), time.clone().pow(-2)])),
        ]
    };
    for (n, dt, steps) in [([1, 1], 1e32, 0), ([2, 1], 0.0, 9)] {
        let result = call(
            &session,
            "num_native_gpe_run_ground_state",
            &args(n, dt, steps),
        );
        let Value::Tensor(density) = result else {
            panic!("expected density tensor")
        };
        assert_eq!(density.shape, vec![n[1] as usize, n[0] as usize]);
        assert!(
            density
                .data
                .iter()
                .all(|v| (number(v) - 0.125).abs() < 1e-15)
        );
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Float.with_dim(length.clone().pow(-2)),
                &density.elem,
            )
            .accepted()
        );
    }
    for invalid in [args([1, 0], 0.0, 0), args([1, 1], -1.0, 0)] {
        assert!(
            session
                .call_value("num_native_gpe_run_ground_state", &invalid, &CallCtx::new())
                .is_err()
        );
    }
}

#[test]
fn native_potts_conditional_probabilities_are_normalized_and_read_only() {
    let s = rssn::Session::new();
    let f = call(&s, "num_native_md_normalized_frame", &[]);
    let context = CallCtx::new();
    for (dimensionality, shape) in [(1, vec![3]), (2, vec![3, 4]), (3, vec![3, 4, 3])] {
        let create = || {
            call(
                &s,
                &format!("num_native_spin_potts_new_d{dimensionality}"),
                &[
                    f.clone(),
                    ints(&shape),
                    Value::Bool(true),
                    Value::int(64),
                    Value::Float(0.0),
                    Value::int(884),
                ],
            )
        };
        let first = create();
        let second = create();
        let name = format!("num_native_spin_potts_heat_bath_probabilities_d{dimensionality}");
        let values = call(&s, &name, &[first.clone(), Value::int(0)]);
        let Value::Tensor(values) = values else {
            panic!("expected conditional probabilities")
        };
        assert_eq!(values.shape, vec![64]);
        assert_eq!(values.elem, TypeExpr::Float);
        assert_eq!(values.data, vec![Value::Float(1.0 / 64.0); 64]);
        assert!(
            s.call_value(&name, &[first.clone(), Value::int(4096)], &context)
                .is_err()
        );
        for model in [&first, &second] {
            let result = call(
                &s,
                &format!("num_native_spin_potts_heat_bath_sweep_d{dimensionality}"),
                std::slice::from_ref(model),
            );
            assert_eq!(result, Value::Bool(true));
        }
        let energy_name = format!("num_native_spin_potts_energy_d{dimensionality}");
        assert_eq!(
            call(&s, &energy_name, &[first]),
            call(&s, &energy_name, &[second])
        );
    }
}

#[test]
fn native_ising_online_measurements_keep_units_and_seeded_final_state() {
    let s = rssn::Session::new();
    let energy = DimExpr::Mul(vec![dim("Mass"), dim("Length").pow(2), dim("Time").pow(-2)]);
    let f = call(
        &s,
        "num_native_md_frame",
        &[
            quantity(2.0, dim("Mass")),
            quantity(1.0, dim("Length")),
            quantity(1.0, dim("Time")),
            quantity(1.0, DimExpr::Mul(vec![dim("Current"), dim("Time")])),
        ],
    );
    let temperature = quantity(6.0 / 1.380_649e-23, dim("Temperature"));
    let args = |width, height, warmup, sweeps, stride, seed| {
        [
            f.clone(),
            Value::int(width),
            Value::int(height),
            temperature.clone(),
            Value::int(warmup),
            Value::int(sweeps),
            Value::int(stride),
            Value::int(seed),
        ]
    };
    let name = "num_native_ising_measure_seeded";
    let single = call(&s, name, &args(1, 1, 7_usize, 8, 3, 31));
    assert_eq!(field(&single, "samples"), &Value::int(2));
    let observables = field(&single, "observables");
    assert_eq!(number(field(observables, "energy")), -4.0);
    assert_eq!(number(field(observables, "magnetization")), 1.0);
    assert_eq!(number(field(observables, "heat_capacity_over_kB")), 0.0);
    let Value::Dyn(measured) = field(observables, "susceptibility") else {
        panic!("expected inverse-energy susceptibility")
    };
    let expected = TypeExpr::Float.with_dim(energy.pow(-1));
    let actual = TypeExpr::Float.with_dim(measured.dim.clone());
    assert!(rssn::kernel::overload::matches_type(&expected, &actual).accepted());
    let Value::Variant(binder) = field(observables, "binder") else {
        panic!("expected optional Binder")
    };
    assert_eq!(binder.tag, "Some");
    assert!((number(&binder.payload[0]) - 2.0 / 3.0).abs() < 1e-14);

    let ensemble = call(&s, name, &args(7, 5, 9, 13, 4, 777));
    let plain = call(
        &s,
        "num_native_ising_run_ising_simulation_seeded",
        &[
            f.clone(),
            Value::int(7),
            Value::int(5),
            temperature.clone(),
            Value::int(22),
            Value::int(777),
        ],
    );
    assert_eq!(field(&ensemble, "spins"), field(&plain, "spins"));
    assert_eq!(field(&ensemble, "samples"), &Value::int(3));
    assert!((0..64).any(|seed| {
        let narrow = call(&s, name, &args(1, 2, 0, 1, 1, seed));
        let Value::Variant(binder) = field(field(&narrow, "observables"), "binder") else {
            panic!("expected optional Binder")
        };
        binder.tag == "None"
    }));
    for invalid in [
        args(1, 1, 0, 8, 0, 31),
        args(1, 1, 0, 8, 9, 31),
        args(1, 1, usize::MAX, 8, 3, 31),
    ] {
        assert!(s.call_value(name, &invalid, &CallCtx::new()).is_err());
    }
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
fn prepared_measured_schrodinger_reuses_a_single_cell_and_signed_time() {
    let session = rssn::Session::new();
    let length = dim("Length");
    let time = dim("Time");
    let mass = dim("Mass");
    let energy = DimExpr::Mul(vec![
        mass.clone(),
        length.clone().pow(2),
        time.clone().pow(-2),
    ]);
    let action = DimExpr::Mul(vec![energy.clone(), time.clone()]);
    let amplitude = TypeExpr::Complex.with_dim(length.clone().pow(-1));
    let new = "num_native_quantum_schrodinger_plan_new";
    let evolve = "num_native_quantum_schrodinger_plan_evolve";
    let arguments = [
        ints(&[1, 1]),
        floats(vec![2], &[2.0, 3.0], length),
        quantity(-0.25, time.clone()),
        quantity(1.0, action),
        quantity(1.0, mass),
        floats(vec![1], &[1.0], energy),
    ];
    let plan = call(&session, new, &arguments);
    let sample = |re| {
        tensor(
            vec![1],
            vec![Value::Complex(ComplexValue::Float { re, im: 0.0 })],
            amplitude.clone(),
        )
    };
    for re in [1.0, -2.0] {
        let wave = sample(re);
        let result = call(
            &session,
            evolve,
            &[plan.clone(), wave.clone(), Value::int(4)],
        );
        let result_time = field(&result, "time");
        assert_eq!(result_time.dimension(), Some(&time));
        assert_eq!(number(result_time), -1.0);
        let Value::Tensor(psi) = field(&result, "final_psi") else {
            panic!("final wave tensor")
        };
        assert_eq!(psi.shape, [1]);
        assert!(rssn::kernel::overload::matches_type(&amplitude, &psi.elem).accepted());
        let Value::Complex(ComplexValue::Float {
            re: actual_re,
            im: actual_im,
        }) = psi.data[0]
        else {
            panic!("complex final amplitude")
        };
        assert!((actual_re - re * 1.0_f64.cos()).abs() < 2e-15);
        assert!((actual_im - re * 1.0_f64.sin()).abs() < 2e-15);
        let unchanged = call(
            &session,
            evolve,
            &[plan.clone(), wave.clone(), Value::int(0)],
        );
        assert_eq!(number(field(&unchanged, "time")), 0.0);
        assert_eq!(field(&unchanged, "final_psi"), &wave);
    }

    let args = [plan.clone(), sample(1.0), Value::int(1)];
    assert!(
        rssn::Session::new()
            .call_value(evolve, &args, &CallCtx::new())
            .is_err()
    );
    let invalid = tensor(vec![1], vec![Value::Float(1.0)], TypeExpr::Complex);
    assert!(
        session
            .call_value(evolve, &[plan, invalid, Value::int(1)], &CallCtx::new())
            .is_err()
    );

    let mut stopped = arguments;
    stopped[2] = quantity(0.0, time);
    let stopped = call(&session, new, &stopped);
    let initial = sample(1.0);
    let result = call(&session, evolve, &[stopped, initial.clone(), Value::int(7)]);
    assert_eq!(number(field(&result, "time")), 0.0);
    assert_eq!(field(&result, "final_psi"), &initial);
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

#[test]
fn prepared_linear_split_step_supports_measured_odd_grids_and_inverse_time() {
    let s = rssn::Session::new();
    let hbar = 1.054_571_817e-34;
    let length = dim("Length");
    let time = dim("Time");
    let energy = DimExpr::Mul(vec![
        dim("Mass"),
        length.clone().pow(2),
        time.clone().pow(-2),
    ]);
    let amplitude = length.clone().pow(Exp::Rat(Rational::new(-1, 2)));
    let psi = tensor(
        vec![5],
        vec![Value::Complex(ComplexValue::Float { re: 0.5, im: 0.0 }); 5],
        TypeExpr::Complex.with_dim(amplitude.clone()),
    );
    let args = |dt| {
        [
            ints(&[5]),
            floats(vec![1], &[7.0], length.clone()),
            quantity(2.0, dim("Mass")),
            floats(vec![5], &[hbar; 5], energy.clone()),
            quantity(dt, time.clone()),
            quantity(4.0, length.clone()),
            quantity(2.0, time.clone()),
            Value::String("Yoshida4".into()),
        ]
    };
    let new = "num_native_split_linear_plan_new";
    let evolve = "num_native_split_linear_plan_evolve";
    let forward = call(&s, new, &args(0.125));
    let backward = call(&s, new, &args(-0.125));
    let stationary = call(&s, new, &args(0.0));
    let size = call(
        &s,
        "num_native_split_linear_plan_size",
        std::slice::from_ref(&forward),
    );
    assert_eq!(size, Value::int(5));
    let evolved = call(&s, evolve, &[forward.clone(), psi.clone(), Value::int(4)]);
    let round_trip = call(&s, evolve, &[backward, evolved.clone(), Value::int(4)]);
    let zero_time = call(&s, evolve, &[stationary, psi.clone(), Value::int(4)]);
    for (value, phase) in [(evolved, -0.5_f64), (round_trip, 0.0), (zero_time, 0.0)] {
        let Value::Tensor(t) = value else {
            panic!("expected a wavefunction tensor")
        };
        assert_eq!(t.shape, [5]);
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Complex.with_dim(amplitude.clone()),
                &t.elem
            )
            .accepted()
        );
        for value in t.data {
            let Value::Complex(ComplexValue::Float { re, im }) = value else {
                panic!("expected a complex sample")
            };
            assert!((re - 0.5 * phase.cos()).abs() < 5e-14);
            assert!((im - 0.5 * phase.sin()).abs() < 5e-14);
        }
    }
    let mut bad_units = args(0.125);
    bad_units[3] = floats(vec![5], &[hbar; 5], length);
    assert!(s.call_value(new, &bad_units, &CallCtx::new()).is_err());
    assert!(
        s.call_value(
            evolve,
            &[forward, psi, Value::int(100_001)],
            &CallCtx::new()
        )
        .is_err()
    );
}

#[test]
fn linear_split_plan_reuses_the_normalized_centered_solver_metadata() {
    let s = rssn::Session::new();
    let solver = call(
        &s,
        "num_native_split_centered_normalized",
        &[
            ints(&[4]),
            floats(vec![1], &[2.0], DimExpr::Dimless),
            Value::Float(1.0),
        ],
    );
    let plan = call(
        &s,
        "num_native_split_prepare_linear",
        &[
            solver,
            floats(vec![4], &[2.0; 4], DimExpr::Dimless),
            Value::Float(-0.25),
            Value::String("Strang".into()),
        ],
    );
    let psi = tensor(
        vec![4],
        vec![Value::Complex(ComplexValue::Float { re: 1.0, im: 0.0 }); 4],
        TypeExpr::Complex,
    );
    let out = call(
        &s,
        "num_native_split_linear_plan_evolve",
        &[plan.clone(), psi.clone(), Value::int(1)],
    );
    let Value::Tensor(out) = out else {
        panic!("expected a wavefunction tensor")
    };
    assert_eq!(out.elem, TypeExpr::Complex.with_dim(DimExpr::Dimless));
    for value in out.data {
        let Value::Complex(ComplexValue::Float { re, im }) = value else {
            panic!("expected a complex sample")
        };
        assert!((re - 0.5_f64.cos()).abs() < 5e-14);
        assert!((im - 0.5_f64.sin()).abs() < 5e-14);
    }
    let fresh = rssn::Session::new();
    assert!(
        fresh
            .call_value(
                "num_native_split_linear_plan_evolve",
                &[plan, psi, Value::int(1)],
                &CallCtx::new()
            )
            .is_err()
    );
}

#[test]
fn classical_nbody_retains_measured_snapshots_and_preflights_output_budget() {
    let session = rssn::Session::new();
    let length = dim("Length");
    let velocity = DimExpr::Mul(vec![length.clone(), dim("Time").pow(-1)]);
    let gravity = DimExpr::Mul(vec![
        length.clone().pow(3),
        dim("Mass").pow(-1),
        dim("Time").pow(-2),
    ]);
    let particle = call(
        &session,
        "num_native_classical_particle_new",
        &[
            quantity(3.0, dim("Mass")),
            floats(vec![3], &[-1.0, 2.0, 3.0], length.clone()),
            floats(vec![3], &[2.0, -4.0, 0.0], velocity.clone()),
        ],
    );
    let particle_kind = particle.type_expr();
    let mut args = [
        tensor(vec![1], vec![particle.clone()], particle_kind.clone()),
        quantity(0.5, dim("Time")),
        Value::int(2),
        quantity(0.0, gravity),
    ];
    let name = "num_native_classical_simulate_n_body";
    let manifest = session.manifest();
    let signature = &manifest.op(name).unwrap().overloads[0];
    assert_eq!(signature.params.len(), 4);
    assert_eq!(signature.params[2], TypeExpr::Int);
    assert!(
        signature
            .effects
            .effects
            .contains(&rssn::kernel::signature::Effect::Stateful)
    );
    assert!(
        matches!(&signature.params[0], TypeExpr::Tensor { rank: NatExpr::Lit(1), elem, .. } if **elem == particle_kind)
    );
    let history = call(&session, name, &args);
    let Value::Tensor(history) = history else {
        panic!("expected an outer snapshot tensor")
    };
    assert_eq!(history.shape, [3]);
    for (index, snapshot) in history.data.iter().enumerate() {
        let Value::Tensor(snapshot) = snapshot else {
            panic!("expected an inner particle tensor")
        };
        assert_eq!(snapshot.shape, [1]);
        let state = &snapshot.data[0];
        assert_eq!(field(state, "mass"), &quantity(3.0, dim("Mass")));
        for (key, units, components) in [
            (
                "position",
                length.clone(),
                [-1.0 + index as f64, 2.0 - 2.0 * index as f64, 3.0],
            ),
            ("velocity", velocity.clone(), [2.0, -4.0, 0.0]),
        ] {
            let Value::Tensor(vector) = field(state, key) else {
                panic!("expected a measured vector")
            };
            assert_eq!(vector.shape, [3]);
            assert!(
                rssn::kernel::overload::matches_type(
                    &TypeExpr::Float.with_dim(units),
                    &vector.elem
                )
                .accepted()
            );
            assert_eq!(
                vector.data.iter().map(number).collect::<Vec<_>>(),
                components
            );
        }
        if index == 0 {
            assert_eq!(
                call(
                    &session,
                    "num_native_classical_particle_data",
                    std::slice::from_ref(&particle)
                ),
                *state
            );
        }
    }
    args[2] = Value::int(0);
    let Value::Tensor(zero_steps) = call(&session, name, &args) else {
        panic!("expected a zero-step history")
    };
    assert_eq!(zero_steps.shape, [1]);
    assert_eq!(zero_steps.data, history.data[..1]);
    args[2] = Value::int(14563);
    let error = session
        .call_value(name, &args, &CallCtx::new())
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("classical N-body history exceeds the native value budget"),
        "{error}"
    );
    let Value::Tensor(initial) = &history.data[0] else {
        unreachable!()
    };
    assert_eq!(
        call(&session, "num_native_classical_particle_data", &[particle]),
        initial.data[0]
    );
}

#[test]
fn classical_nbody_rejects_nonfinite_evolution_without_mutating_input() {
    let session = rssn::Session::new();
    let length = dim("Length");
    let velocity = DimExpr::Mul(vec![length.clone(), dim("Time").pow(-1)]);
    let mut input = [
        quantity(1.0, dim("Mass")),
        floats(vec![3], &[f64::MAX, 0.0, 0.0], length.clone()),
        floats(vec![3], &[f64::MAX, 0.0, 0.0], velocity),
    ];
    let particle = call(&session, "num_native_classical_particle_new", &input);
    let before = call(
        &session,
        "num_native_classical_particle_data",
        std::slice::from_ref(&particle),
    );
    let gravity = DimExpr::Mul(vec![
        length.pow(3),
        dim("Mass").pow(-1),
        dim("Time").pow(-2),
    ]);
    let name = "num_native_classical_simulate_n_body";
    let error = session
        .call_value(
            name,
            &[
                tensor(vec![1], vec![particle.clone()], particle.type_expr()),
                quantity(1.0, dim("Time")),
                Value::int(1),
                quantity(0.0, gravity),
            ],
            &CallCtx::new(),
        )
        .unwrap_err();
    assert!(error.to_string().contains("non-finite result"), "{error}");
    assert_eq!(
        call(&session, "num_native_classical_particle_data", &[particle]),
        before
    );
    input[0] = quantity(0.0, dim("Mass"));
    assert!(
        session
            .call_value("num_native_classical_particle_new", &input, &CallCtx::new())
            .is_err()
    );
}
