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
