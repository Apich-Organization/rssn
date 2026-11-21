//! Native legacy continuum contracts, including finite scale extremes.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    n: &str,
    a: &[Value],
) -> Value {
    s.call_value(n, a, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{n}:{e:?}"))
}
fn q(
    x: f64,
    n: &str,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: DimExpr::base(n),
    })
}
fn v(
    xs: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![xs.len()],
        elem,
        data: xs,
    })
}
fn measured(
    values: Vec<f64>,
    shape: Vec<usize>,
    d: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape,
        elem: TypeExpr::Float.with_dim(d),
        data: values.into_iter().map(Value::Float).collect(),
    })
}
fn number(v: &Value) -> f64 {
    if let Value::Dyn(v) = v {
        number(&v.value)
    } else {
        v.as_f64().unwrap()
    }
}
#[test]
fn measured_finite_differences_and_extreme_native_norms_are_correct() {
    let s = Session::new();
    let f = measured(
        vec![0., 1., 4., 1., 2., 5., 4., 5., 8.],
        vec![3, 3],
        DimExpr::base("Temperature"),
    );
    let lap = call(
        &s,
        "num_grid_legacy_compute_laplacian",
        &[f, q(1., "Length"), q(1., "Length")],
    );
    let Value::Tensor(lap) = lap else { panic!() };
    assert_eq!(number(&lap.data[4]), 4.);
    let expected = TypeExpr::Float.with_dim(DimExpr::Mul(vec![
        DimExpr::base("Temperature"),
        DimExpr::base("Length").pow(-2),
    ]));
    assert!(rssn::kernel::overload::matches_type(&expected, &lap.elem).accepted());
    let huge = rssn::kernels::matrix::Matrix::new(1, 2, vec![3e300, 4e300]);
    let rms = rssn::sim::physics_cfd::l2_norm(&huge);
    assert!(
        rssn::sim::physics_cfd::l2_norm(&rssn::kernels::matrix::Matrix::new(1, 1, vec![f64::NAN]))
            .is_nan()
    );
    assert!((rms / 1e300 - (12.5f64).sqrt()).abs() < 1e-14);
    let x = rssn::kernels::matrix::Matrix::new(1, 1, vec![3e300]);
    let y = rssn::kernels::matrix::Matrix::new(1, 1, vec![4e300]);
    assert_eq!(
        rssn::sim::physics_cfd::max_velocity_magnitude(&x, &y),
        5e300
    );
}
#[test]
fn explicit_fdm_grid_specialization_handles_empty_and_mutable_physical_data() {
    let s = Session::new();
    let frame = call(
        &s,
        "num_grid_frame_new",
        &[q(2., "Length"), q(3., "Time"), q(5., "Temperature")],
    );
    let empty = call(
        &s,
        "num_grid_legacy_fdm_new",
        &[frame.clone(), v(vec![Value::int(0)], TypeExpr::Int)],
    );
    assert_eq!(
        call(&s, "num_grid_legacy_fdm_is_empty", &[empty]),
        Value::Bool(true)
    );
    let grid = call(
        &s,
        "num_grid_legacy_fdm_with_value",
        &[
            frame,
            v(vec![Value::int(3)], TypeExpr::Int),
            q(10., "Temperature"),
        ],
    );
    let before = call(
        &s,
        "num_grid_legacy_fdm_as_slice",
        std::slice::from_ref(&grid),
    );
    let Value::Tensor(before) = before else {
        panic!()
    };
    assert_eq!(number(&before.data[1]), 10.);
    call(
        &s,
        "num_grid_legacy_fdm_as_mut_slice",
        &[
            grid.clone(),
            measured(vec![2., 4., 6.], vec![3], DimExpr::base("Temperature")),
        ],
    );
    let after = call(&s, "num_grid_legacy_fdm_as_slice", &[grid]);
    let Value::Tensor(after) = after else {
        panic!()
    };
    assert_eq!(number(&after.data[1]), 4.);
}
#[test]
fn fvm_advection1_preserves_measured_boundaries_and_rejects_without_updating() {
    let session = Session::new();
    let frame = call(
        &session,
        "num_grid_frame_new",
        &[q(4.0, "Length"), q(8.0, "Time"), q(2.0, "Temperature")],
    );
    let field = DimExpr::base("Temperature");
    let mesh = call(
        &session,
        "num_grid_legacy_fvm_mesh1_new",
        &[
            frame,
            v(vec![Value::int(3)], TypeExpr::Int),
            measured(vec![6.0], vec![1], DimExpr::base("Length")),
            measured(vec![2.0, 4.0, 8.0], vec![3], field.clone()),
        ],
    );
    let velocity = |value| {
        measured(
            vec![value],
            vec![1],
            DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time").pow(-1)]),
        )
    };
    let mut args = [
        mesh.clone(),
        velocity(1.0),
        q(1.0, "Time"),
        Value::int(1),
        q(6.0, "Temperature"),
        q(12.0, "Temperature"),
    ];
    let advection = "num_grid_legacy_fvm_advection1";
    let check_field = |result: &Value, values: [f64; 3]| {
        let Value::Tensor(tensor) = result else {
            panic!("expected a measured field")
        };
        assert_eq!(tensor.shape, [3]);
        assert_eq!(tensor.elem, TypeExpr::Float.with_dim(field.clone()));
        assert_eq!(tensor.data.iter().map(number).collect::<Vec<_>>(), values);
    };
    check_field(&call(&session, advection, &args), [4.0, 3.0, 6.0]);
    args[1] = velocity(-1.0);
    check_field(&call(&session, advection, &args), [3.5, 4.5, 9.0]);
    args[2] = q(3.0, "Time");
    assert!(
        session
            .call_value(advection, &args, &CallCtx::new())
            .is_err()
    );
    args[2] = q(1.0, "Time");
    args[4] = q(6.0, "Length");
    assert!(
        session
            .call_value(advection, &args, &CallCtx::new())
            .is_err()
    );
    check_field(
        &call(&session, "num_grid_legacy_fvm_mesh_data", &[mesh]),
        [3.5, 4.5, 9.0],
    );
}

#[test]
fn fvm_mesh1_keeps_sample_order_when_normalized_domain_extent_overflows() {
    let session = Session::new();
    let frame = call(
        &session,
        "num_grid_frame_new",
        &[q(0.5, "Length"), q(1.0, "Time"), q(1.0, "Temperature")],
    );
    let mesh = call(
        &session,
        "num_grid_legacy_fvm_mesh1_new",
        &[
            frame,
            v(vec![Value::int(3)], TypeExpr::Int),
            measured(vec![f64::MAX], vec![1], DimExpr::base("Length")),
            measured(vec![2.0, 7.0, -5.0], vec![3], DimExpr::base("Temperature")),
        ],
    );
    assert_eq!(
        call(
            &session,
            "num_grid_legacy_fvm_num_cells",
            std::slice::from_ref(&mesh)
        ),
        Value::int(3),
    );
    let data = call(&session, "num_grid_legacy_fvm_mesh_data", &[mesh]);
    let Value::Tensor(data) = data else {
        panic!("expected a measured field")
    };
    assert_eq!(
        data.data.iter().map(number).collect::<Vec<_>>(),
        [2.0, 7.0, -5.0]
    );
}

#[test]
fn fvm_advection1_does_not_commit_an_unrepresentable_measured_result() {
    let session = Session::new();
    let frame = call(
        &session,
        "num_grid_frame_new",
        &[q(1.0, "Length"), q(1.0, "Time"), q(3.0, "Temperature")],
    );
    let mesh = call(
        &session,
        "num_grid_legacy_fvm_mesh1_new",
        &[
            frame,
            v(vec![Value::int(3)], TypeExpr::Int),
            measured(vec![3.0], vec![1], DimExpr::base("Length")),
            measured(vec![0.0; 3], vec![3], DimExpr::base("Temperature")),
        ],
    );
    let error = session
        .call_value(
            "num_grid_legacy_fvm_advection1",
            &[
                mesh.clone(),
                measured(
                    vec![1.0],
                    vec![1],
                    DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time").pow(-1)]),
                ),
                q(1.0, "Time"),
                Value::int(1),
                q(f64::MAX, "Temperature"),
                q(0.0, "Temperature"),
            ],
            &CallCtx::new(),
        )
        .unwrap_err();
    // MAX / 3 is finite, but rounding that normalized number back by 3 overflows.
    assert!(error.to_string().contains("overflows Float"), "{error}");
    let result = call(&session, "num_grid_legacy_fvm_mesh_data", &[mesh]);
    let Value::Tensor(field) = result else {
        panic!("expected a measured field");
    };
    assert_eq!(field.shape, [3]);
    assert_eq!(
        field.elem,
        TypeExpr::Float.with_dim(DimExpr::base("Temperature"))
    );
    assert_eq!(field.data.iter().map(number).collect::<Vec<_>>(), [0.0; 3]);
}

#[test]
fn fvm_burgers1_matches_one_measured_lax_friedrichs_step() {
    let session = Session::new();
    let velocity = DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time").pow(-1)]);
    let reference = Value::Dyn(DynValue {
        value: Box::new(Value::Float(2.0)),
        dim: velocity.clone(),
    });
    let frame = call(
        &session,
        "num_grid_frame_new",
        &[q(4.0, "Length"), q(8.0, "Time"), reference],
    );
    let mesh = call(
        &session,
        "num_grid_legacy_fvm_mesh1_new",
        &[
            frame,
            v(vec![Value::int(3)], TypeExpr::Int),
            measured(vec![3.0], vec![1], DimExpr::base("Length")),
            measured(vec![1.0, 2.0, 3.0], vec![3], velocity.clone()),
        ],
    );
    let name = "num_grid_legacy_fvm_burgers1";
    let result = call(
        &session,
        name,
        &[mesh.clone(), q(0.1, "Time"), Value::int(1)],
    );
    let Value::Tensor(field) = &result else {
        panic!("expected a measured field")
    };
    assert_eq!(field.shape, [3]);
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(velocity), &field.elem)
            .accepted()
    );
    // One step using fluxes [0.5, -3.75, -1.75, 4.5] with zero-gradient ghosts.
    for (cell, expected) in field.data.iter().zip([1.425, 1.8, 2.375]) {
        assert!((number(cell) - expected).abs() < 1e-14);
    }
    assert!(
        session
            .call_value(
                name,
                &[mesh.clone(), q(1.0, "Time"), Value::int(1)],
                &CallCtx::new()
            )
            .is_err()
    );
    assert_eq!(
        call(&session, "num_grid_legacy_fvm_mesh_data", &[mesh]),
        result
    );
}

#[test]
fn fvm_shallow_water1_keeps_a_lake_at_rest_and_rejects_dry_momentum() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let momentum = DimExpr::Mul(vec![length.clone().pow(2), DimExpr::base("Time").pow(-1)]);
    let gravity = DimExpr::Mul(vec![length.clone(), DimExpr::base("Time").pow(-2)]);
    let name = "num_grid_legacy_fvm_shallow_water1";
    let mut args = [
        measured(vec![2.0; 3], vec![3], length.clone()),
        measured(vec![0.0; 3], vec![3], momentum.clone()),
        v(vec![Value::int(3)], TypeExpr::Int),
        measured(vec![1.0], vec![1], length.clone()),
        q(0.125, "Time"),
        Value::int(1),
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(1.0)),
            dim: gravity,
        }),
    ];
    let result = call(&session, name, &args);
    let Value::Record(result) = result else {
        panic!("expected a shallow-water record")
    };
    for (key, expected, dim) in [
        ("depth", 2.0, length.clone()),
        ("momentum_x", 0.0, momentum.clone()),
    ] {
        let Value::Tensor(field) = result.get(key).unwrap() else {
            panic!("expected a measured field")
        };
        assert_eq!(field.shape, [3]);
        assert!(
            rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(dim), &field.elem)
                .accepted()
        );
        assert_eq!(
            field.data.iter().map(number).collect::<Vec<_>>(),
            [expected; 3]
        );
    }
    args[0] = measured(vec![2.0, 0.0, 2.0], vec![3], length);
    args[1] = measured(vec![0.0, 1.0, 0.0], vec![3], momentum);
    let error = session
        .call_value(name, &args, &CallCtx::new())
        .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("dry cells require zero momentum"),
        "{error}"
    );
}

#[test]
fn fvm_boundary_masks_reject_unprotected_edges_before_parallel_native_access() {
    let s = Session::new();
    let frame = call(
        &s,
        "num_grid_frame_new",
        &[q(1., "Length"), q(1., "Time"), q(1., "Temperature")],
    );
    let shape = v(vec![Value::int(3), Value::int(3)], TypeExpr::Int);
    let lengths = measured(vec![3., 3.], vec![2], DimExpr::base("Length"));
    let init = measured(vec![1.; 9], vec![9], DimExpr::base("Temperature"));
    let m = call(
        &s,
        "num_grid_legacy_fvm_mesh2_new",
        &[frame, shape, lengths, init],
    );
    let velocity = measured(
        vec![0.1, 0.],
        vec![2],
        DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time").pow(-1)]),
    );
    let mask = v(vec![Value::Bool(false); 9], TypeExpr::Bool);
    assert!(
        s.call_value(
            "num_grid_legacy_fvm_advection2",
            &[m, velocity, q(0.1, "Time"), Value::int(1), mask],
            &CallCtx::new()
        )
        .is_err()
    );
}

#[test]
fn fvm_advection3_uses_unmarked_faces_and_validates_mask_layout() {
    let session = Session::new();
    let frame = call(
        &session,
        "num_grid_frame_new",
        &[q(2.0, "Length"), q(4.0, "Time"), q(3.0, "Temperature")],
    );
    let [width, height, depth] = [5, 4, 3];
    let cells = width * height * depth;
    let dimensions = v(
        [width, height, depth]
            .map(|n| Value::int(n as i64))
            .to_vec(),
        TypeExpr::Int,
    );
    let lengths = measured(vec![5.0, 4.0, 3.0], vec![3], DimExpr::base("Length"));
    let field = DimExpr::base("Temperature");
    let mesh = call(
        &session,
        "num_grid_legacy_fvm_mesh3_new",
        &[
            frame,
            dimensions,
            lengths,
            measured(vec![6.0; cells], vec![cells], field.clone()),
        ],
    );
    let velocity = measured(
        vec![1.0, 0.0, 0.0],
        vec![3],
        DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time").pow(-1)]),
    );
    let mask = |shape, mark: Option<usize>| {
        Value::Tensor(TensorValue {
            shape,
            elem: TypeExpr::Bool,
            data: (0..cells).map(|i| Value::Bool(mark == Some(i))).collect(),
        })
    };
    let advance = |mask| {
        session.call_value(
            "num_grid_legacy_fvm_advection3",
            &[
                mesh.clone(),
                velocity.clone(),
                q(1.0, "Time"),
                Value::int(1),
                mask,
            ],
            &CallCtx::new(),
        )
    };
    let malformed = advance(mask(vec![width, height, depth], None)).unwrap_err();
    assert!(format!("{malformed:?}").contains("3D boundary mask shape"));
    let short = v(vec![Value::Bool(false); cells - 1], TypeExpr::Bool);
    assert!(advance(short).is_err());
    let numbers = v(vec![Value::int(0); cells], TypeExpr::Int);
    assert!(advance(numbers).is_err());
    let Value::Tensor(initial) = call(
        &session,
        "num_grid_legacy_fvm_mesh_data",
        std::slice::from_ref(&mesh),
    ) else {
        panic!("a flat measured field is expected")
    };
    assert!(initial.data.iter().all(|value| number(value) == 6.0));
    let Value::Tensor(after) = advance(mask(vec![depth, height, width], Some(2))).unwrap() else {
        panic!("a flat measured field is expected")
    };
    assert_eq!(after.shape, [cells]);
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(field), &after.elem)
            .accepted()
    );
    assert_eq!(
        after.data[..width].iter().map(number).collect::<Vec<_>>(),
        [0.0, 6.0, 0.0, 6.0, 6.0]
    );
    assert_eq!(
        after.data[cells - width..]
            .iter()
            .map(number)
            .collect::<Vec<_>>(),
        [0.0, 6.0, 6.0, 6.0, 6.0]
    );
    let Value::Tensor(after) = advance(mask(vec![cells], None)).unwrap() else {
        panic!("a flat measured field is expected")
    };
    assert_eq!(
        after.data[..width].iter().map(number).collect::<Vec<_>>(),
        [0.0, 0.0, 6.0, 0.0, 6.0]
    );
    assert_eq!(
        Value::Tensor(after),
        call(&session, "num_grid_legacy_fvm_mesh_data", &[mesh])
    );
}

#[test]
fn fvm_advection3_scenario_returns_unscaled_samples_and_fixed_outer_cells() {
    let Value::Tensor(result) = call(
        &Session::new(),
        "num_grid_legacy_fvm_scenario_advection3",
        &[],
    ) else {
        panic!("a flat dimensionless field is expected")
    };
    assert_eq!(result.shape, [27_000]);
    assert_eq!(result.elem, TypeExpr::Float);
    for (i, value) in result.data.iter().enumerate() {
        let (x, y, z) = (i % 30, (i / 30) % 30, i / 900);
        if [x, y, z]
            .into_iter()
            .any(|coordinate| matches!(coordinate, 0 | 29))
        {
            assert_eq!(number(value), 0.0);
        }
    }
    assert!(result.data.iter().any(|value| number(value) > 0.0));
}

#[test]
fn tiny_fdtd_model_restores_field_units_post_update_times_and_rejects_bad_contracts() {
    let s = Session::new();
    let electric = DimExpr::Mul(vec![
        DimExpr::base("Mass"),
        DimExpr::base("Length"),
        DimExpr::base("Current").pow(-1),
        DimExpr::base("Time").pow(-3),
    ]);
    let quantity = |value, dimension| {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(value)),
            dim: dimension,
        })
    };
    let mut args = vec![
        Value::int(4),
        Value::int(4),
        Value::int(6),
        v(vec![Value::int(2), Value::int(2)], TypeExpr::Int),
        quantity(0.125, DimExpr::base("Time").pow(-1)),
        q(1.0, "Length"),
        quantity(
            1.0,
            DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time").pow(-1)]),
        ),
        quantity(2.0, electric.clone()),
    ];
    let Value::Record(result) = call(&s, "num_native_fdtd_model_run", &args) else {
        panic!("FDTD model returns its measured fields and times")
    };
    let Value::Tensor(times) = result.get("times").unwrap() else {
        panic!("times tensor")
    };
    assert_eq!(times.shape, vec![2]);
    assert_eq!(number(&times.data[0]), 0.5);
    assert_eq!(number(&times.data[1]), 3.0);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(DimExpr::base("Time")),
            &times.elem,
        )
        .accepted()
    );
    let Value::Tensor(field) = result.get("electric_field").unwrap() else {
        panic!("electric field tensor")
    };
    assert_eq!(field.shape, vec![2, 4, 4]);
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(electric), &field.elem,)
            .accepted()
    );
    let native = rssn::sim::models_fdtd_electrodynamics::run_fdtd_simulation(
        &rssn::sim::models_fdtd_electrodynamics::FdtdParameters {
            width: 4,
            height: 4,
            time_steps: 6,
            source_pos: (2, 2),
            source_freq: 0.125 * 0.5,
        },
    );
    assert!(native.iter().flat_map(|a| a.iter()).any(|x| *x != 0.0));
    for (actual, normalized) in field.data.iter().zip(native.iter().flat_map(|a| a.iter())) {
        assert_eq!(number(actual), 2.0 * normalized);
        assert!(matches!(actual, Value::Dyn(_)));
    }
    args[4] = q(0.125, "Length");
    assert!(
        s.call_value("num_native_fdtd_model_run", &args, &CallCtx::new())
            .is_err()
    );
    let filename = std::env::temp_dir().join(format!(
        "rssn-fdtd-rejected-{}-{}.npy",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos(),
    ));
    assert!(!filename.exists());
    assert!(
        s.call_value(
            "num_native_fdtd_model_save",
            &[
                Value::int(4),
                Value::int(6),
                Value::String(filename.to_string_lossy().into_owned()),
                Value::String("measured-electric-field".into()),
            ],
            &CallCtx::new(),
        )
        .is_err()
    );
    assert!(
        !filename.exists(),
        "representation rejection must precede IO"
    );
}

#[test]
fn fdtd_final_field_uses_the_last_step_and_retains_physical_units() {
    use rssn::sim::models_fdtd_electrodynamics::FdtdParameters;
    use rssn::sim::models_fdtd_electrodynamics::run_fdtd_final_field;

    let s = Session::new();
    let electric = DimExpr::Mul(vec![
        DimExpr::base("Mass"),
        DimExpr::base("Length"),
        DimExpr::base("Current").pow(-1),
        DimExpr::base("Time").pow(-3),
    ]);
    let quantity = |value, dim| {
        Value::Dyn(DynValue {
            value: Box::new(Value::Float(value)),
            dim,
        })
    };
    let source = |row, column| v(vec![Value::int(row), Value::int(column)], TypeExpr::Int);
    let mut args = vec![
        Value::int(4),
        Value::int(4),
        Value::int(7),
        source(2, 2),
        quantity(0.25, DimExpr::base("Time").pow(-1)),
        q(2.0, "Length"),
        quantity(
            1.0,
            DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time").pow(-1)]),
        ),
        quantity(2.0, electric.clone()),
    ];
    let name = "num_native_fdtd_model_run_final";
    let time_dimension = DimExpr::base("Time");
    let field_dimension = TypeExpr::Float.with_dim(electric);
    let validate = |result: Value, time, shape: &[usize]| {
        let Value::Record(result) = result else {
            panic!("FDTD result record")
        };
        let t = result.get("time").unwrap();
        let Value::Tensor(field) = result.get("electric_field").unwrap() else {
            panic!("FDTD field tensor")
        };
        assert_eq!(number(t), time);
        assert_eq!(t.dimension(), Some(&time_dimension));
        assert!(rssn::kernel::overload::matches_type(&field_dimension, &field.elem).accepted());
        assert_eq!(field.shape, shape);
        assert!(
            field
                .data
                .iter()
                .all(|value| matches!(value, Value::Dyn(_)))
        );
        field.data.iter().map(number).collect::<Vec<_>>()
    };

    let field = validate(call(&s, name, &args), 7.0, &[4, 4]);
    let params = FdtdParameters {
        width: 4,
        height: 4,
        time_steps: 7,
        source_pos: (2, 2),
        source_freq: 0.25,
    };
    let native = run_fdtd_final_field(&params).unwrap();
    assert_eq!(
        field,
        native.iter().map(|&value| 2.0 * value).collect::<Vec<_>>()
    );
    let Value::Record(sampled) = call(&s, "num_native_fdtd_model_run", &args) else {
        panic!()
    };
    let Value::Tensor(retained) = sampled.get("electric_field").unwrap() else {
        panic!()
    };
    assert!(
        field
            .iter()
            .zip(&retained.data[16..])
            .any(|(last, sampled)| *last != number(sampled))
    );

    args[0] = Value::int(1);
    args[1] = Value::int(1);
    args[2] = Value::int(2);
    args[3] = source(0, 0);
    let one_cell = validate(call(&s, name, &args), 2.0, &[1, 1]);
    let pulse = (-8.41_f64).exp() * 0.95 * 2.0;
    assert!((one_cell[0] - pulse).abs() < 1e-15);
    args[2] = Value::int(0);
    assert_eq!(validate(call(&s, name, &args), 0.0, &[1, 1]), [0.0]);

    args[2] = Value::int(5);
    args[4] = quantity(1e308, DimExpr::base("Time").pow(-1));
    args[5] = q(4.0, "Length");
    assert_eq!(validate(call(&s, name, &args), 10.0, &[1, 1]), [0.0]);
    args[3] = source(0, 1);
    assert!(s.call_value(name, &args, &CallCtx::new()).is_err());
    args[3] = v(vec![Value::int(0)], TypeExpr::Int);
    assert!(s.call_value(name, &args, &CallCtx::new()).is_err());
    args[3] = source(0, 0);
    args[4] = q(1.0, "Length");
    assert!(s.call_value(name, &args, &CallCtx::new()).is_err());

    args[0] = Value::int(128);
    args[1] = Value::int(128);
    args[2] = Value::int(21);
    args[4] = quantity(0.0, DimExpr::base("Time").pow(-1));
    assert!(
        s.call_value("num_native_fdtd_model_run", &args, &CallCtx::new())
            .is_err()
    );
    let no_source = validate(call(&s, name, &args), 42.0, &[128, 128]);
    assert!(no_source.iter().all(|&value| value == 0.0));
}
