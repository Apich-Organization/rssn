//! Native grid callbacks retain coherent reference scales and their compiled parents.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::ComplexValue;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
fn call(
    s: &Session,
    n: &str,
    a: &[Value],
) -> Result<Value, OpError> {
    s.call_value(n, a, &CallCtx::new())
}
fn vec(
    v: Vec<Value>,
    elem: TypeExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![v.len()],
        data: v,
        elem,
    })
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
fn frame(s: &Session) -> Value {
    call(
        s,
        "num_grid_frame_new",
        &[
            measured(2.0, DimExpr::Base("Length".into())),
            measured(2.0, DimExpr::Base("Time".into())),
            measured(3.0, DimExpr::Base("Temperature".into())),
        ],
    )
    .unwrap()
}
fn compile(
    s: &Session,
    expr: &str,
    names: &[&str],
    refs: Vec<Value>,
    output: Value,
) -> Value {
    let backend = call(s, "num_backend_interpreter", &[]).unwrap();
    let refs = refs
        .into_iter()
        .map(|v| call(s, "num_backend_reference", &[v]).unwrap())
        .collect();
    let output = call(s, "num_backend_reference", &[output]).unwrap();
    call(
        s,
        "num_backend_compile",
        &[
            backend,
            Value::String(expr.into()),
            vec(
                names.iter().map(|s| Value::String((*s).into())).collect(),
                TypeExpr::String,
            ),
            vec(
                refs,
                TypeExpr::Opaque("rssn/native/backend-scalar/v1".into()),
            ),
            output,
            Value::int(10000),
        ],
    )
    .unwrap()
}
fn field<'a>(
    v: &'a Value,
    n: &str,
) -> &'a Value {
    let Value::Record(r) = v else {
        panic!("record")
    };
    &r.fields.iter().find(|(k, _)| k == n).unwrap().1
}
fn number(v: &Value) -> f64 {
    if let Value::Dyn(d) = v {
        number(&d.value)
    } else {
        v.as_f64().unwrap()
    }
}

#[test]
fn native_spectral_modes_preserve_the_reference_and_advance_signed_frequencies() {
    let session = Session::new();
    let operator = "num_grid_sm_advance_advection_diffusion_3d_modes";
    let modes = vec(
        vec![Value::Complex(ComplexValue::Float { re: 1.0, im: 0.0 }); 3],
        TypeExpr::Complex,
    );
    let mut args = [
        frame(&session),
        vec(
            [3, 1, 1].into_iter().map(Value::int).collect(),
            TypeExpr::Int,
        ),
        vec(
            [std::f64::consts::TAU / 3.0, 1.0, 1.0]
                .into_iter()
                .map(Value::Float)
                .collect(),
            TypeExpr::Float,
        ),
        modes.clone(),
        vec(
            [1.0, 0.0, 0.0].into_iter().map(Value::Float).collect(),
            TypeExpr::Float,
        ),
        Value::Float(0.25),
        Value::Float(0.5),
        Value::int(2),
    ];
    let result = call(&session, operator, &args).unwrap();
    assert_eq!(field(&result, "shape"), &args[1]);
    let reference = field(&result, "reference");
    assert_eq!(number(reference), 3.0);
    assert_eq!(reference.dimension(), Some(&DimExpr::base("Temperature")));
    let Value::Tensor(data) = field(&result, "data") else {
        panic!("coefficient vector");
    };
    assert_eq!(data.shape, [3]);
    for (mode, (want_re, want_im)) in
        data.data
            .iter()
            .zip([(1.0, 0.0), (0.515625, -0.875), (0.515625, 0.875)])
    {
        let Value::Complex(ComplexValue::Float { re, im }) = mode else {
            panic!("complex coefficient");
        };
        assert!((re - want_re).abs() < 1e-14 && (im - want_im).abs() < 1e-14);
    }
    args[7] = Value::int(0);
    assert_eq!(
        field(&call(&session, operator, &args).unwrap(), "data"),
        &modes
    );
    args[7] = Value::int(8193);
    assert!(call(&session, operator, &args).is_err());
}

#[test]
fn timed_method_of_lines_restores_physical_time_and_outward_gradient() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let temperature = DimExpr::base("Temperature");
    let dim =
        |base: &DimExpr, power| DimExpr::Mul(vec![temperature.clone(), base.clone().pow(power)]);
    let at = |n, dimension: &DimExpr| measured(n, dimension.clone());
    let bound = |x| {
        vec(
            vec![Value::Float(x)],
            TypeExpr::Float.with_dim(length.clone()),
        )
    };
    let grid = call(
        &session,
        "num_analysis_pde_grid",
        &[
            bound(0.0),
            bound(2.0),
            vec(vec![Value::int(5)], TypeExpr::Int),
            at(2.0, &length),
            at(4.0, &time),
            at(3.0, &temperature),
        ],
    )
    .unwrap();
    let gradient = dim(&length, -1);
    let rhs = compile(
        &session,
        "0.5*x + 1 + 0.2*(ux - 0.5*(1+t)) + 0.1*uxx",
        &["t", "x", "u", "ux", "uxx"],
        vec![
            at(1.0, &time),
            at(1.0, &length),
            at(1.0, &temperature),
            at(1.0, &gradient),
            at(1.0, &dim(&length, -2)),
        ],
        at(1.0, &dim(&time, -1)),
    );
    let initial = compile(
        &session,
        "0.5*x",
        &["x"],
        vec![at(1.0, &length)],
        at(1.0, &temperature),
    );
    let left = compile(
        &session,
        "t",
        &["t"],
        vec![at(1.0, &time)],
        at(1.0, &temperature),
    );
    let right = compile(
        &session,
        "0.5*(1+t)",
        &["t"],
        vec![at(1.0, &time)],
        at(1.0, &gradient),
    );
    let solve = |method: &str, left: Value, budget| {
        call(
            &session,
            "num_analysis_pde_method_of_lines_1d_timed",
            &[
                grid.clone(),
                vec(
                    vec![
                        Value::String("Dirichlet".into()),
                        Value::String("Neumann".into()),
                    ],
                    TypeExpr::String,
                ),
                left,
                right.clone(),
                at(1.0, &time),
                Value::String(method.into()),
                Value::Float(1e-8),
                at(1e-10, &temperature),
                at(0.125, &time),
                at(0.0, &time),
                at(0.25, &time),
                Value::int(200),
                Value::int(budget),
                rhs.clone(),
                initial.clone(),
            ],
        )
    };
    for method in ["Radau5", "Bdf3"] {
        let field_handle = solve(method, left.clone(), 10_000).unwrap();
        let result = call(&session, "num_analysis_pde_field_data", &[field_handle]).unwrap();
        assert_eq!(number(field(&result, "time")), 1.0);
        let Value::Tensor(values) = field(&result, "values") else {
            panic!("expected field vector")
        };
        assert_eq!(values.shape, [5]);
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Float.with_dim(temperature.clone()),
                &values.elem
            )
            .accepted()
        );
        for (actual, expected) in values.data.iter().zip([1.0, 1.5, 2.0, 2.5, 3.0]) {
            assert!((number(actual) - expected).abs() < 1e-7);
        }
    }
    assert!(solve("Radau5", initial.clone(), 10_000).is_err());
    assert!(matches!(
        solve("Radau5", left, 1),
        Err(OpError::Domain(message)) if message.contains("compiled function evaluation budget")
    ));
}

#[test]
fn prepared_cartesian_poisson_uses_one_measured_source_evaluation_per_point() {
    let s = Session::new();
    let length = DimExpr::base("Length");
    let temperature = DimExpr::base("Temperature");
    let plan = call(
        &s,
        "num_grid_fem_cartesian_prepare",
        &[frame(&s), vec(vec![Value::int(2); 2], TypeExpr::Int)],
    )
    .unwrap();
    let shape = call(
        &s,
        "num_grid_fem_cartesian_nodal_shape",
        std::slice::from_ref(&plan),
    )
    .unwrap();
    assert_eq!(shape, vec(vec![Value::int(3); 2], TypeExpr::Int));
    let source = compile(
        &s,
        "1",
        &["x", "y"],
        vec![measured(2.0, length.clone()); 2],
        measured(
            0.75,
            DimExpr::Mul(vec![temperature.clone(), length.pow(-2)]),
        ),
    );
    let name = "num_grid_fem_cartesian_solve_prepared";
    assert!(call(&s, name, &[plan.clone(), source.clone(), Value::int(15)]).is_err());
    let result = call(&s, name, &[plan, source, Value::int(16)]).unwrap();
    let reference = field(&result, "value_reference");
    assert_eq!(reference.dimension(), Some(&temperature));
    let Value::Tensor(tensor) = field(&result, "data") else {
        panic!("Poisson values")
    };
    assert_eq!(tensor.shape, [9]);
    assert!((number(&tensor.data[4]) * number(reference) - 0.28125).abs() < 1e-13);
}

#[test]
fn prepared_triangles_reuse_geometry_for_independent_measured_problems() {
    let s = Session::new();
    let grid_frame = frame(&s);
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let temperature = DimExpr::base("Temperature");
    let mesh = call(
        &s,
        "num_triangle_rectangle",
        &[
            grid_frame.clone(),
            vec(
                [0.0, 1.0, 0.0, 1.0].into_iter().map(Value::Float).collect(),
                TypeExpr::Float,
            ),
            Value::int(2),
            Value::int(2),
        ],
    )
    .unwrap();
    let plan = call(&s, "num_triangle_prepare", &[mesh.clone(), Value::int(2)]).unwrap();
    assert_eq!(
        call(
            &s,
            "num_triangle_prepared_dofs",
            std::slice::from_ref(&plan)
        )
        .unwrap(),
        Value::int(25)
    );
    assert!(call(&s, "num_triangle_prepare", &[mesh, Value::int(3)]).is_err());
    let compile_on_grid = |expression, output| {
        compile(
            &s,
            expression,
            &["x", "y"],
            vec![measured(2.0, length.clone()); 2],
            output,
        )
    };
    let k = compile_on_grid(
        "1",
        measured(
            2.0,
            DimExpr::Mul(vec![length.clone().pow(2), time.clone().pow(-1)]),
        ),
    );
    let c = compile_on_grid("0", measured(0.5, time.clone().pow(-1)));
    let f = compile_on_grid(
        "0",
        measured(1.5, DimExpr::Mul(vec![temperature.clone(), time.pow(-1)])),
    );
    for expression in ["1+x+2*y", "2-2*x+y"] {
        let exact = compile_on_grid(expression, measured(3.0, temperature.clone()));
        let problem = call(
            &s,
            "num_triangle_problem",
            &[
                grid_frame.clone(),
                k.clone(),
                c.clone(),
                f.clone(),
                exact.clone(),
                Value::int(1000),
            ],
        )
        .unwrap();
        let solution = call(
            &s,
            "num_triangle_solve_prepared",
            &[plan.clone(), problem.clone()],
        )
        .unwrap();
        let err = call(
            &s,
            "num_triangle_max_nodal_error",
            &[solution, exact, Value::int(100)],
        )
        .unwrap();
        assert_eq!(err.dimension(), Some(&temperature));
        assert!(number(&err) < 1e-10);
        assert!(
            call(
                &Session::new(),
                "num_triangle_solve_prepared",
                &[plan.clone(), problem]
            )
            .is_err()
        );
    }
}

#[test]
fn cutcell_native_simplex_keeps_representable_extreme_fractions() {
    let session = Session::new();
    let op = "num_grid_simplex_positive_fraction";
    let values = |samples: &[f64]| {
        vec(
            samples.iter().copied().map(Value::Float).collect(),
            TypeExpr::Float,
        )
    };
    let big = f64::MAX;
    let tiny = f64::from_bits(1);
    assert_eq!(
        call(&session, op, &[values(&[big, big, -big, -big])]).unwrap(),
        Value::Float(0.5)
    );
    assert_eq!(
        call(&session, op, &[values(&[tiny, tiny, -1.0])]).unwrap(),
        Value::Float(f64::from_bits(2))
    );
    assert!(call(&session, op, &[values(&[1.0, f64::INFINITY])]).is_err());
    assert!(call(&session, op, &[values(&[1.0; 5])]).is_err());
}

#[test]
fn cutcell_compiled_level_sets_use_coherent_coordinates_and_preserve_budgets() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let coordinate_reference = measured(8.0, length.clone());
    let plane = compile(
        &session,
        "x - y",
        &["x", "y", "z"],
        vec![coordinate_reference.clone(); 3],
        coordinate_reference.clone(),
    );
    let points = |x: [f64; 3]| vec(x.into_iter().map(Value::Float).collect(), TypeExpr::Float);
    let axes = |x: &[usize]| vec(x.iter().copied().map(Value::int).collect(), TypeExpr::Int);
    let segment = "num_grid_open_length_fraction";
    let mut segment_args = [
        frame(&session),
        plane.clone(),
        points([0.0, 0.25, 0.0]),
        Value::int(0),
        Value::Float(1.0),
        Value::int(256),
    ];
    assert_eq!(
        call(&session, segment, &segment_args).unwrap(),
        Value::Float(0.75)
    );
    segment_args[5] = Value::int(16);
    let error = call(&session, segment, &segment_args).unwrap_err();
    assert!(error.to_string().contains("budget exhausted"), "{error}");

    let box_clip = "num_grid_open_fraction";
    let mut box_args = [
        frame(&session),
        plane,
        points([0.0; 3]),
        points([1.0, 1.0, 0.0]),
        axes(&[1, 0]),
        Value::int(2),
        Value::int(512),
    ];
    assert_eq!(
        call(&session, box_clip, &box_args).unwrap(),
        Value::Float(0.5)
    );
    box_args[4] = axes(&[0, 0]);
    assert!(call(&session, box_clip, &box_args).is_err());
    box_args[4] = axes(&[0, 1]);
    box_args[5] = Value::int(10);
    assert!(call(&session, box_clip, &box_args).is_err());

    segment_args[1] = compile(
        &session,
        "x - y",
        &["x", "y", "z"],
        vec![
            coordinate_reference.clone(),
            coordinate_reference.clone(),
            measured(1.0, DimExpr::base("Time")),
        ],
        coordinate_reference,
    );
    segment_args[5] = Value::int(256);
    assert!(call(&session, segment, &segment_args).is_err());
}

#[test]
fn cutcell_compiled_distances_keep_nonzero_subnormal_signs() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let field = compile(
        &session,
        "1",
        &["x", "y", "z"],
        vec![measured(1.0, length.clone()); 3],
        measured(f64::from_bits(1), length),
    );
    let zero = vec(vec![Value::Float(0.0); 3], TypeExpr::Float);
    assert_eq!(
        call(
            &session,
            "num_grid_open_length_fraction",
            &[
                frame(&session),
                field.clone(),
                zero.clone(),
                Value::int(0),
                Value::Float(1.0),
                Value::int(17),
            ]
        )
        .unwrap(),
        Value::Float(1.0)
    );
    assert_eq!(
        call(
            &session,
            "num_grid_open_fraction",
            &[
                frame(&session),
                field,
                zero,
                vec(vec![Value::Float(1.0); 3], TypeExpr::Float),
                vec(
                    vec![Value::int(0), Value::int(1), Value::int(2)],
                    TypeExpr::Int
                ),
                Value::int(0),
                Value::int(9),
            ]
        )
        .unwrap(),
        Value::Float(1.0)
    );
}

#[test]
fn legacy_fvm_lax_friedrichs_flux_checks_compiled_physical_domains() {
    let session = Session::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let state = DimExpr::Mul(vec![length.clone(), time.clone().pow(-1)]);
    let flux = DimExpr::Mul(vec![length.clone().pow(2), time.clone().pow(-2)]);
    let input_reference = measured(1.0, state.clone());
    let callback = compile(
        &session,
        "u * u",
        &["u"],
        vec![input_reference.clone()],
        measured(1.0, flux.clone()),
    );
    let name = "num_grid_legacy_fvm_lax_friedrichs_flux";
    let mut arguments = [
        measured(2.0, state.clone()),
        measured(3.0, state),
        measured(0.5, time),
        measured(1.0, length.clone()),
        callback,
    ];
    let result = call(&session, name, &arguments).unwrap();
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(flux), &result.type_expr())
            .accepted()
    );
    assert_eq!(number(&result), 5.5);
    arguments[4] = compile(
        &session,
        "u * u",
        &["u"],
        vec![input_reference],
        measured(1.0, length),
    );
    assert!(call(&session, name, &arguments).is_err());
}

#[test]
fn prepared_simplex_reuses_geometry_and_keeps_independent_measures() {
    let s = Session::new();
    let length = DimExpr::Base("Length".into());
    let time = DimExpr::Base("Time".into());
    let reference = measured(2.0, length.clone());
    let vertices = Value::Tensor(TensorValue {
        shape: vec![3, 2],
        elem: TypeExpr::Float,
        data: [0.0, 0.0, 1.0, 0.0, 0.0, 1.0].map(Value::Float).to_vec(),
    });
    let plan = call(
        &s,
        "num_linear_simplex_map_new",
        &[vertices.clone(), reference.clone()],
    )
    .unwrap();
    let area = call(
        &s,
        "num_linear_simplex_map_volume",
        std::slice::from_ref(&plan),
    )
    .unwrap();
    assert_eq!(
        area,
        measured(2.0, length.clone().pow(Exp::Rat(Rational::int(2))))
    );
    let ordinate = measured(3.0, time.clone());
    let program = compile(
        &s,
        "1",
        &["x", "y"],
        vec![reference.clone(); 2],
        ordinate.clone(),
    );
    let controls = [Value::Float(1e-11), Value::Float(0.0), Value::int(1_000)];
    let retained = [
        vec![plan, ordinate.clone(), program.clone()],
        controls.to_vec(),
    ]
    .concat();
    let one_shot = [
        vec![vertices, reference, ordinate, program],
        controls.to_vec(),
    ]
    .concat();
    let out = call(&s, "num_linear_integrate_simplex_map", &retained).unwrap();
    assert_eq!(
        out,
        call(&s, "num_linear_integrate_simplex", &one_shot).unwrap()
    );
    let value = field(&out, "value");
    assert!((number(value) - 6.0).abs() < 1e-12);
    let Value::Dyn(value) = value else {
        panic!("measured integral");
    };
    assert_eq!(
        value.dim,
        DimExpr::Mul(vec![length.pow(Exp::Rat(Rational::int(2))), time])
    );
}

#[test]
fn periodic_helmholtz_batch_checks_physical_forcing_and_preserves_row_order() {
    let s = Session::new();
    let grid = call(
        &s,
        "num_grid_fourier_new",
        &[
            frame(&s),
            vec(vec![Value::int(4), Value::int(8)], TypeExpr::Int),
            vec(vec![Value::Float(1.25), Value::Float(2.0)], TypeExpr::Float),
        ],
    )
    .unwrap();
    let inverse_square = DimExpr::Base("Length".into()).pow(Exp::Rat(Rational::int(-2)));
    let temperature = DimExpr::Base("Temperature".into());
    let forcing_unit = DimExpr::Mul(vec![inverse_square.clone(), temperature.clone()]);
    let source = |rows, columns, data, dim| {
        Value::Tensor(TensorValue {
            shape: vec![rows, columns],
            elem: TypeExpr::Float.with_dim(dim),
            data,
        })
    };
    let data = [1.5, 3.0]
        .into_iter()
        .flat_map(|v| vec![Value::Float(v); 32])
        .collect::<Vec<_>>();
    let run = |sources, sigma| {
        call(
            &s,
            "num_grid_fourier_helmholtz_batch",
            &[grid.clone(), sources, sigma],
        )
    };
    let sigma = measured(0.5, inverse_square.clone());
    let sources = source(2, 32, data.clone(), forcing_unit.clone());
    let out = run(sources.clone(), sigma.clone()).unwrap();
    let Value::Tensor(shape) = field(&out, "shape") else {
        panic!("batch shape");
    };
    assert_eq!(
        shape.data.iter().map(number).collect::<Vec<_>>(),
        [2.0, 4.0, 8.0]
    );
    let Value::Tensor(values) = field(&out, "data") else {
        panic!("flattened batch data");
    };
    assert_eq!(values.shape, [64]);
    assert!(
        values.data[..32]
            .iter()
            .all(|v| (number(v) - 1.0).abs() < 1e-12)
    );
    assert!(
        values.data[32..]
            .iter()
            .all(|v| (number(v) - 2.0).abs() < 1e-12)
    );
    assert_eq!(
        field(&out, "value_reference"),
        &measured(3.0, temperature.clone())
    );
    assert!(run(sources.clone(), measured(0.0, inverse_square)).is_err());
    assert!(run(sources, Value::Float(0.5)).is_err());
    assert!(run(source(2, 32, data.clone(), temperature), sigma.clone()).is_err());
    assert!(run(source(4, 16, data, forcing_unit.clone()), sigma.clone()).is_err());
    assert!(run(source(0, 32, Vec::new(), forcing_unit), sigma).is_err());
}

#[test]
fn periodic_helmholtz_rejects_a_reference_that_erases_nonzero_inputs() {
    let s = Session::new();
    let length = DimExpr::Base("Length".into());
    let field_dimension = DimExpr::Base("Temperature".into());
    let shift_dimension = length.clone().pow(Exp::Rat(Rational::int(-2)));
    let forcing_dimension = DimExpr::Mul(vec![shift_dimension.clone(), field_dimension.clone()]);
    let tiny_frame = call(
        &s,
        "num_grid_frame_new",
        &[
            measured(1e-200, length),
            measured(1.0, DimExpr::Base("Time".into())),
            measured(1.0, field_dimension),
        ],
    )
    .unwrap();
    let grid = call(
        &s,
        "num_grid_fourier_new",
        &[
            tiny_frame,
            vec(vec![Value::int(4)], TypeExpr::Int),
            vec(vec![Value::Float(1.0)], TypeExpr::Float),
        ],
    )
    .unwrap();
    for (forcing, shift) in [(0.0, 1.0), (1.0, 1e100)] {
        let sources = Value::Tensor(TensorValue {
            shape: vec![1, 4],
            data: vec![Value::Float(forcing); 4],
            elem: TypeExpr::Float.with_dim(forcing_dimension.clone()),
        });
        let error = call(
            &s,
            "num_grid_fourier_helmholtz_batch",
            &[
                grid.clone(),
                sources,
                measured(shift, shift_dimension.clone()),
            ],
        )
        .unwrap_err();
        assert!(error.to_string().contains("underflows"));
    }
}

#[test]
fn periodic_burgers_uses_the_velocity_reference_and_conserves_cell_averages() {
    let s = Session::new();
    let length = DimExpr::Base("Length".into());
    let time = DimExpr::Base("Time".into());
    let velocity = DimExpr::Mul(vec![
        length.clone(),
        time.clone().pow(Exp::Rat(Rational::int(-1))),
    ]);
    let velocity_frame = call(
        &s,
        "num_grid_frame_new",
        &[
            measured(2.0, length),
            measured(2.0, time),
            measured(2.0, velocity.clone()),
        ],
    )
    .unwrap();
    let args = |grid_frame, dt, boundary: &str| {
        [
            grid_frame,
            vec(
                [1.0, 2.0, 3.0, 4.0].map(Value::Float).to_vec(),
                TypeExpr::Float,
            ),
            vec(vec![Value::int(4)], TypeExpr::Int),
            vec(vec![Value::Float(1.0)], TypeExpr::Float),
            Value::Float(dt),
            Value::int(1),
            vec(vec![Value::String(boundary.into())], TypeExpr::String),
        ]
    };
    let name = "num_grid_physics_fvm_burgers_nd_with_boundaries";
    let out = call(&s, name, &args(velocity_frame.clone(), 0.0625, "Periodic")).unwrap();
    let Value::Tensor(data) = field(&out, "data") else {
        panic!("Burgers field");
    };
    assert_eq!(
        data.data.iter().map(number).collect::<Vec<_>>(),
        [3.375, 1.75, 2.625, 2.25]
    );
    assert_eq!(field(&out, "value_reference"), &measured(2.0, velocity));
    assert!(call(&s, name, &args(velocity_frame.clone(), 0.25, "Periodic")).is_err());
    assert!(call(&s, name, &args(velocity_frame, 0.0625, "Reflecting")).is_err());
    assert!(call(&s, name, &args(frame(&s), 0.0625, "Periodic")).is_err());
}

#[test]
fn sampled_dirichlet_wave_restores_displacements_at_signed_steps() {
    let s = Session::new();
    let length = DimExpr::Base("Length".into());
    let time = DimExpr::Base("Time".into());
    let temperature = DimExpr::Base("Temperature".into());
    let inv_time = time.clone().pow(Exp::Rat(Rational::int(-1)));
    let temperature_rate = DimExpr::Mul(vec![temperature.clone(), inv_time.clone()]);
    let speed = DimExpr::Mul(vec![length.clone(), inv_time]);
    let initial = call(
        &s,
        "num_grid_legacy_fdm_with_value",
        &[
            frame(&s),
            vec(vec![Value::int(3), Value::int(3)], TypeExpr::Int),
            measured(9.0, temperature.clone()),
        ],
    )
    .unwrap();
    let run = |rate_dim, cadence| {
        call(
            &s,
            "num_grid_legacy_fdm_solve_wave_dirichlet_sampled",
            &[
                initial.clone(),
                vec(
                    vec![Value::Float(2.0); 9],
                    TypeExpr::Float.with_dim(rate_dim),
                ),
                vec(
                    vec![Value::Float(1.0); 2],
                    TypeExpr::Float.with_dim(length.clone()),
                ),
                measured(0.0, speed.clone()),
                measured(-0.25, time.clone()),
                Value::int(5),
                Value::int(cadence),
            ],
        )
    };
    let result = run(temperature_rate.clone(), 2).unwrap();
    let Value::Tensor(steps) = field(&result, "steps") else {
        panic!("wave steps")
    };
    assert_eq!(steps.data, [0, 2, 4, 5].map(Value::int));
    let Value::Tensor(shape) = field(&result, "grid_shape") else {
        panic!("grid shape")
    };
    assert_eq!(shape.data, [3, 3].map(Value::int));
    let Value::Tensor(values) = field(&result, "data") else {
        panic!("sampled displacement")
    };
    assert_eq!(values.shape, [4, 9]);
    assert_eq!(values.elem, TypeExpr::Float.with_dim(temperature.clone()));
    for (row, center) in values
        .data
        .as_chunks::<9>()
        .0
        .iter()
        .zip([9.0, 8.0, 7.0, 6.5])
    {
        assert!((number(&row[4]) - center).abs() < 1e-12);
        assert!(
            row.iter()
                .enumerate()
                .all(|(index, value)| index == 4 || number(value) == 0.0)
        );
    }
    assert!(run(temperature_rate, 0).is_err());
    assert!(run(temperature, 2).is_err());
}

#[test]
fn crank_nicolson_adi_preserves_fixed_temperature_edges() {
    let s = Session::new();
    let name = "num_grid_cnm_solve_heat_equation_2d_cn_adi_fixed_boundary";
    let mut args = [
        frame(&s),
        vec(vec![Value::int(3); 2], TypeExpr::Int),
        vec([1.0, 2.0].map(Value::Float).to_vec(), TypeExpr::Float),
        vec(
            [1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0]
                .map(Value::Float)
                .to_vec(),
            TypeExpr::Float,
        ),
        Value::Float(1.0),
        Value::Float(0.5),
        Value::int(1),
    ];
    let step = call(&s, name, &args).unwrap();
    let Value::Tensor(data) = field(&step, "data") else {
        panic!("normalized temperature grid");
    };
    assert!((number(&data.data[4]) - 20.0 / 27.0).abs() < 1e-14);
    assert!(
        data.data
            .iter()
            .enumerate()
            .all(|(i, row)| i == 4 || number(row) == 1.0)
    );
    assert_eq!(
        field(&step, "value_reference"),
        &measured(3.0, DimExpr::Base("Temperature".into()))
    );
    args[5] = Value::Float(0.0);
    assert_eq!(field(&call(&s, name, &args).unwrap(), "data"), &args[3]);
    args[5] = Value::Float(-0.5);
    assert!(call(&s, name, &args).is_err());
}

#[test]
fn sampled_transport_keeps_its_frame_and_budgets_only_retained_snapshots() {
    let s = Session::new();
    let eq = call(
        &s,
        "num_grid_transport_equation",
        &[
            frame(&s),
            Value::String("Diffusion".into()),
            vec(Vec::new(), TypeExpr::Float),
            Value::Float(0.25),
        ],
    )
    .unwrap();
    let args = [
        eq,
        vec(vec![Value::int(31); 2], TypeExpr::Int),
        vec(vec![Value::Float(1.0); 2], TypeExpr::Float),
        vec(vec![Value::Float(2.0); 31 * 31], TypeExpr::Float),
        Value::Float(0.125),
        Value::int(100),
    ];
    let solve = "num_grid_dimensional_solve_transport_sampled_nd";
    let sampled = |cadence| {
        let mut sampled_args = args.to_vec();
        sampled_args.push(Value::int(cadence));
        call(&s, solve, &sampled_args)
    };
    assert!(sampled(0).is_err());
    assert!(sampled(1).is_err());
    assert!(call(&s, "num_grid_dimensional_solve_transport_nd", &args).is_err());
    let history = sampled(64).unwrap();
    let Value::Tensor(shape) = field(&history, "shape") else {
        panic!("expected snapshot shape")
    };
    assert_eq!(
        shape.data.iter().map(number).collect::<Vec<_>>(),
        [3.0, 961.0]
    );
    let Value::Tensor(data) = field(&history, "data") else {
        panic!("expected normalized snapshot data")
    };
    assert_eq!(data.shape, [3 * 961]);
    assert!(data.data.iter().all(|value| number(value) == 2.0));
    assert_eq!(
        field(&history, "length_reference"),
        &measured(2.0, DimExpr::Base("Length".into()))
    );
    assert_eq!(
        field(&history, "time_reference"),
        &measured(2.0, DimExpr::Base("Time".into()))
    );
    assert_eq!(
        field(&history, "value_reference"),
        &measured(3.0, DimExpr::Base("Temperature".into()))
    );
}

#[test]
fn retained_time_boundary_evaluates_coherent_seconds_and_holds_parent() {
    let s = Session::new();
    let f = frame(&s);
    let program = compile(
        &s,
        "1+t",
        &["t"],
        vec![measured(1.0, DimExpr::Base("Time".into()))],
        measured(3.0, DimExpr::Base("Temperature".into())),
    );
    let datum = call(
        &s,
        "num_grid_datum_time",
        &[
            f,
            program.clone(),
            Value::int(2),
            measured(3.0, DimExpr::Base("Temperature".into())),
        ],
    )
    .unwrap();
    assert!(call(&s, "num_backend_release", std::slice::from_ref(&program)).is_err());
    let point = vec(vec![Value::Float(0.0); 3], TypeExpr::Float);
    let value = call(
        &s,
        "num_grid_datum_eval",
        &[datum.clone(), Value::Float(0.5), point.clone()],
    )
    .unwrap();
    assert_eq!(number(&value), 6.0);
    assert!(matches!(value, Value::Dyn(_)));
    call(
        &s,
        "num_grid_datum_eval",
        &[datum.clone(), Value::Float(0.0), point.clone()],
    )
    .unwrap();
    assert!(
        call(
            &s,
            "num_grid_datum_eval",
            &[datum.clone(), Value::Float(0.0), point]
        )
        .is_err()
    );
    call(&s, "num_backend_release", &[datum]).unwrap();
    assert!(call(&s, "num_backend_release", &[program]).is_ok());
}

#[test]
fn continuum_multigrid_keeps_extreme_forcing_and_callback_diagnostics() {
    use rssn::kernel::value::RecordValue;

    let session = Session::new();
    let faces = call(
        &session,
        "num_grid_faces_dirichlet_zero",
        &[frame(&session)],
    )
    .unwrap();
    let grid = call(
        &session,
        "num_grid_structured_unit",
        &[faces, Value::int(1), Value::int(9)],
    )
    .unwrap();
    let length = DimExpr::base("Length");
    let rhs = DimExpr::Mul(vec![DimExpr::base("Temperature"), length.clone().pow(-2)]);
    let opts = Value::Record(RecordValue {
        fields: vec![
            ("cycle".into(), Value::String("Fmg".into())),
            ("pre".into(), Value::int(2)),
            ("post".into(), Value::int(2)),
            ("max_cycles".into(), Value::int(8)),
            ("tol".into(), Value::Float(1e-10)),
            ("coarsest_nodes".into(), Value::int(9)),
        ],
    });
    let solve = |expression: &str, budget| {
        let source = compile(
            &session,
            expression,
            &["x", "y", "z"],
            vec![measured(1.0, length.clone()); 3],
            measured(1.0, rhs.clone()),
        );
        call(
            &session,
            "num_grid_solve_helmholtz_multigrid",
            &[
                grid.clone(),
                Value::Float(0.0),
                source,
                Value::int(budget),
                opts.clone(),
            ],
        )
    };
    for (expression, expected) in [("1e-200", 5e-201), ("1e200", 5e199)] {
        let solution = solve(expression, 100).unwrap();
        let result = call(&session, "num_grid_multigrid_result", &[solution]).unwrap();
        assert_eq!(field(&result, "converged"), &Value::Bool(true));
        let physical = field(&result, "field");
        let Value::Tensor(values) = field(physical, "data") else {
            panic!("field tensor")
        };
        assert_eq!(values.shape, [9]);
        let center = number(&values.data[4]) * number(field(physical, "value_reference"));
        assert!((center / expected - 1.0).abs() < 1e-10);
    }
    let zero = solve("0", 100).unwrap();
    let factor = call(&session, "num_grid_multigrid_convergence_factor", &[zero]).unwrap();
    assert_eq!(number(&factor), 0.0);
    assert!(
        matches!(solve("1", 1), Err(OpError::Domain(message)) if message.contains("compiled function evaluation budget"))
    );
}

#[test]
fn continuum_multigrid_keeps_both_fmg_samples_at_the_cycle_limit() {
    use rssn::kernel::value::RecordValue;

    let session = Session::new();
    let faces = call(
        &session,
        "num_grid_faces_dirichlet_zero",
        &[frame(&session)],
    )
    .unwrap();
    let grid = call(
        &session,
        "num_grid_structured_unit",
        &[faces, Value::int(1), Value::int(9)],
    )
    .unwrap();
    let options = Value::Record(RecordValue {
        fields: vec![
            ("cycle".into(), Value::String("Fmg".into())),
            ("pre".into(), Value::int(0)),
            ("post".into(), Value::int(0)),
            ("max_cycles".into(), Value::int(512)),
            ("tol".into(), Value::Float(1e-10)),
            ("coarsest_nodes".into(), Value::int(3)),
        ],
    });
    let hierarchy = call(
        &session,
        "num_grid_multigrid_new",
        &[grid, Value::Float(0.0), options],
    )
    .unwrap();
    let source = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0];
    let solution = call(
        &session,
        "num_grid_multigrid_solve",
        &[
            hierarchy,
            vec(
                source.into_iter().map(Value::Float).collect(),
                TypeExpr::Float,
            ),
            Value::Float(0.0),
        ],
    )
    .unwrap();
    let result = call(&session, "num_grid_multigrid_result", &[solution]).unwrap();
    assert_eq!(field(&result, "converged"), &Value::Bool(false));
    let Value::Tensor(history) = field(&result, "residuals") else {
        panic!("residual tensor")
    };
    assert_eq!(history.shape, [514]);
    assert!(history.data.iter().all(|value| number(value) == 1.0));
}

#[test]
fn native_helmholtz_source_keeps_physical_laplacian_units_and_exact_grid_profile() {
    let s = Session::new();
    let f = frame(&s);
    let b = call(&s, "num_grid_faces_dirichlet_zero", &[f]).unwrap();
    let g = call(
        &s,
        "num_grid_structured_unit",
        &[b, Value::int(1), Value::int(9)],
    )
    .unwrap();
    let forcing = DimExpr::Mul(vec![
        DimExpr::Base("Temperature".into()),
        DimExpr::Inverse(Box::new(DimExpr::Mul(vec![
            DimExpr::Base("Length".into()),
            DimExpr::Base("Length".into()),
        ]))),
    ]);
    let program = compile(
        &s,
        "1",
        &["x", "y", "z"],
        vec![measured(1.0, DimExpr::Base("Length".into())); 3],
        measured(0.75, forcing),
    );
    let out = call(
        &s,
        "num_grid_solve_helmholtz",
        &[
            g,
            Value::Float(0.0),
            program,
            Value::int(100),
            Value::String("Direct".into()),
            Value::Float(1e-12),
        ],
    )
    .unwrap();
    let Value::Tensor(data) = field(&out, "data") else {
        panic!("data")
    };
    assert!((number(&data.data[4]) - 0.125).abs() < 1e-12);
    assert_eq!(number(field(&out, "value_reference")), 3.0);
}
