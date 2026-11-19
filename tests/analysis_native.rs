//! Regression checks for new analysis transport and genuine native method selection.
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
fn tensor(
    values: &[f64],
    dimension: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), dimension),
        data: values.iter().map(|v| Value::Float(*v)).collect(),
    })
}
fn values(value: &Value) -> Vec<f64> {
    let Value::Tensor(t) = value else {
        panic!("expected numeric vector: {value:?}");
    };
    t.data
        .iter()
        .map(|v| {
            let v = if let Value::Dyn(d) = v {
                d.value.as_ref()
            } else {
                v
            };
            v.as_f64().unwrap()
        })
        .collect()
}
#[test]
fn physical_wavelet_bands_roundtrip_without_rectangularizing_ragged_levels() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let dimension = DimExpr::Base("Length".into());
    let original = [1.0, 2.0, -3.0, 4.0, 5.0, -6.0, 7.0, 8.0];
    let bands = session
        .call_value(
            "num_analysis_wavelets_wavedec",
            &[
                tensor(&original, dimension.clone()),
                Value::String("Haar".into()),
                Value::int(2),
            ],
            &ctx,
        )
        .unwrap();
    let Value::Record(record) = &bands else {
        panic!("wavelet bands require explicit data/offsets");
    };
    assert_eq!(
        values(record.get("offsets").unwrap()),
        vec![0.0, 2.0, 4.0, 8.0]
    );
    let restored = session
        .call_value(
            "num_analysis_wavelets_waverec",
            &[bands, Value::String("Haar".into())],
            &ctx,
        )
        .unwrap();
    for (a, b) in values(&restored).iter().zip(original) {
        assert!((a - b).abs() < 1e-12);
    }
    let TypeExpr::Tensor { elem, .. } = ctx.value_type(&restored).unwrap() else {
        panic!()
    };
    assert_eq!(*elem, TypeExpr::Dim(Box::new(TypeExpr::Float), dimension));
}
#[test]
fn nonlinear_conjugate_gradient_entry_uses_its_named_native_algorithm() {
    use argmin::core::OptimizationResult;
    use argmin::core::Solver;
    use argmin::core::State;
    fn solver_name<O, S, I: State>(result: &OptimizationResult<O, S, I>) -> &str
    where
        S: Solver<O, I>,
    {
        result.solver.name()
    }
    let options = rssn::kernels::optimize::OptimizationConfig {
        max_iters: 50,
        tolerance: 1e-12,
        problem_type: rssn::kernels::optimize::ProblemType::Sphere,
        dimension: 2,
    };
    let result = rssn::kernels::optimize::EquationOptimizer::auto_solve_conjugate_gradient(
        rssn::kernels::optimize::Sphere,
        ndarray::Array1::from(vec![4.0, -3.0]),
        &options,
    )
    .unwrap();
    assert_eq!(solver_name(&result), "Nonlinear Conjugate Gradient");
    assert!(result.state.get_best_cost() < 1e-10);
}

fn integers(v: &[usize]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![v.len()],
        elem: TypeExpr::Int,
        data: v.iter().map(|n| Value::int(*n)).collect(),
    })
}

fn scalar(
    value: f64,
    dimension: DimExpr,
) -> Value {
    Value::Dyn(rssn::kernel::value::DynValue {
        dim: dimension,
        value: Box::new(Value::Float(value)),
    })
}
fn dynamic_vector(data: Vec<Value>) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        elem: TypeExpr::Dyn(rssn::kernel::signature::RepExpr {
            repr: "term".into(),
            dim: None,
        }),
        data,
    })
}
fn strings(values: &[&str]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        elem: TypeExpr::String,
        data: values.iter().map(|s| Value::String((*s).into())).collect(),
    })
}
fn ode_options(time: &DimExpr) -> Value {
    Value::Record(rssn::kernel::value::RecordValue {
        fields: vec![
            ("rtol".into(), Value::Float(1e-7)),
            ("atol".into(), Value::Float(1e-9)),
            ("h0".into(), scalar(0.01, time.clone())),
            ("hmin".into(), scalar(0.0, time.clone())),
            ("hmax".into(), scalar(0.1, time.clone())),
            ("max_steps".into(), Value::int(2000)),
        ],
    })
}
fn coherent(value: &Value) -> f64 {
    match value {
        | Value::Dyn(x) => coherent(&x.value),
        | Value::Quantity(x) => coherent(&x.value),
        | _ => value.as_f64().unwrap(),
    }
}
#[test]
fn native_ode_mixed_physical_state_retains_true_methods_and_budgets() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let time = DimExpr::Base("Time".into());
    let length = DimExpr::Base("Length".into());
    let velocity = DimExpr::Mul(vec![
        length.clone(),
        DimExpr::Inverse(Box::new(time.clone())),
    ]);
    let acceleration = DimExpr::Mul(vec![
        velocity.clone(),
        DimExpr::Inverse(Box::new(time.clone())),
    ]);
    let tape = session
        .call_value(
            "num_native_tape_lower",
            &[
                strings(&["v", "-q"]),
                strings(&["t", "q", "v"]),
                strings(&[]),
                dynamic_vector(vec![
                    scalar(1.0, time.clone()),
                    scalar(1.0, length.clone()),
                    scalar(1.0, velocity.clone()),
                ]),
                dynamic_vector(vec![]),
                dynamic_vector(vec![
                    scalar(1.0, velocity.clone()),
                    scalar(1.0, acceleration),
                ]),
            ],
            &ctx,
        )
        .unwrap();
    let rhs = session
        .call_value(
            "num_analysis_ode_rhs_new",
            &[tape, dynamic_vector(vec![])],
            &ctx,
        )
        .unwrap();
    let initial = dynamic_vector(vec![
        scalar(1.0, length.clone()),
        scalar(0.0, velocity.clone()),
    ]);
    let empty = dynamic_vector(vec![]);
    let args = |controls: Value, budget: usize| {
        vec![
            rhs.clone(),
            initial.clone(),
            scalar(0.0, time.clone()),
            scalar(1.0, time.clone()),
            ode_options(&time),
            controls,
            empty.clone(),
            Value::int(budget),
        ]
    };
    let fixed = session
        .call_value(
            "num_analysis_ode_fixed",
            &args(integers(&[32, 2]), 500),
            &ctx,
        )
        .unwrap();
    let last = session
        .call_value("num_analysis_ode_last", &[fixed], &ctx)
        .unwrap();
    let Value::Tensor(last) = last else {
        panic!("physical state vector expected")
    };
    assert!((coherent(&last.data[0]) - 1.0f64.cos()).abs() < 2e-7);
    assert!((coherent(&last.data[1]) + 1.0f64.sin()).abs() < 2e-7);
    // The legacy adaptive function cannot silently ignore an explicit h0.
    let unsupported =
        session.call_value("num_analysis_ode_adaptive", &args(integers(&[]), 500), &ctx);
    assert!(
        matches!(unsupported,Err(rssn::kernel::OpError::Domain(ref message)) if message.contains("does not expose h0/hmin/hmax"))
    );
    // The stiff facade calls native Radau with its exact native tape Jacobian.
    let radau = session
        .call_value(
            "num_analysis_ode_radau_jac",
            &args(integers(&[1]), 20000),
            &ctx,
        )
        .unwrap();
    let at_half = session
        .call_value(
            "num_analysis_ode_eval",
            &[radau, scalar(0.5, time.clone())],
            &ctx,
        )
        .unwrap();
    let Value::Variant(result) = at_half else {
        panic!("native dense evaluation returns Option")
    };
    assert_eq!(result.tag, "Some");
    let Value::Tensor(state) = &result.payload[0] else {
        panic!()
    };
    assert!((coherent(&state.data[0]) - 0.5f64.cos()).abs() < 2e-6);
    let exhausted =
        session.call_value("num_analysis_ode_fixed", &args(integers(&[32, 2]), 1), &ctx);
    assert!(matches!(
        exhausted,
        Err(rssn::kernel::OpError::BudgetExceeded)
    ));
    // An initial velocity carrying length units must be rejected before native work.
    let mut wrong = args(integers(&[2, 2]), 100);
    wrong[1] = dynamic_vector(vec![scalar(1.0, length.clone()), scalar(0.0, length)]);
    assert!(
        session
            .call_value("num_analysis_ode_fixed", &wrong, &ctx)
            .is_err()
    );
}

#[test]
fn native_embedded_runge_kutta_stops_on_nonfinite_oracle_without_nan_step_loop() {
    use rssn::sim::physics_rkm::BogackiShampine23;
    use rssn::sim::physics_rkm::CashKarp45;
    use rssn::sim::physics_rkm::DormandPrince54;
    use rssn::sim::physics_rkm::OdeSystem;
    struct Invalid;
    impl OdeSystem for Invalid {
        fn dim(&self) -> usize {
            1
        }

        fn eval(
            &self,
            _t: f64,
            _y: &[f64],
            out: &mut [f64],
        ) {
            out.fill(f64::NAN);
        }
    }
    assert_eq!(
        DormandPrince54::new()
            .solve(&Invalid, &[1.0], (0.0, 1.0), 0.1, (1e-6, 1e-9))
            .len(),
        1
    );
    assert_eq!(
        CashKarp45::default()
            .solve(&Invalid, &[1.0], (0.0, 1.0), 0.1, (1e-6, 1e-9))
            .len(),
        1
    );
    assert_eq!(
        BogackiShampine23::default()
            .solve(&Invalid, &[1.0], (0.0, 1.0), 0.1, (1e-6, 1e-9))
            .len(),
        1
    );
}

#[test]
fn native_fixed_delay_exact_cell_lookup_and_short_final_step_are_valid() {
    let (times, states) = rssn::kernels::ode_adaptive::dde_rk4(
        |_t, _state, lag| vec![lag[0]],
        |_t| vec![1.0],
        0.2,
        0.0,
        0.25,
        0.2,
    )
    .unwrap();
    assert_eq!(times.len(), 3);
    assert!((times.last().unwrap() - 0.25).abs() < 1e-15);
    assert!(states.iter().flatten().all(|v| v.is_finite()));
}

#[test]
fn yee_source_owns_its_charged_program_after_original_handle_release() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name: &str, args: &[Value]| {
        session
            .call_value(name, args, &ctx)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"))
    };
    let l = DimExpr::base("Length");
    let t = DimExpr::base("Time");
    let m = DimExpr::base("Mass");
    let charge = DimExpr::Mul(vec![DimExpr::base("Current"), t.clone()]);
    let inverse = |d| DimExpr::Inverse(Box::new(d));
    let electric = DimExpr::Mul(vec![
        m.clone(),
        l.clone(),
        inverse(charge.clone()),
        inverse(t.clone()),
        inverse(t.clone()),
    ]);
    let impedance = DimExpr::Mul(vec![
        m,
        l.clone(),
        l.clone(),
        inverse(t.clone()),
        inverse(charge.clone()),
        inverse(charge.clone()),
    ]);
    let current = DimExpr::Mul(vec![
        charge,
        inverse(t.clone()),
        inverse(l.clone()),
        inverse(l.clone()),
    ]);
    let options = Value::Record(rssn::kernel::value::RecordValue {
        fields: vec![
            ("cfl".into(), Value::Float(0.5)),
            ("pml_cells".into(), Value::int(0)),
            ("pml_order".into(), Value::Float(3.0)),
            ("pml_kappa_max".into(), Value::Float(1.0)),
            ("dm_min_area".into(), Value::Float(0.01)),
            ("pml_alpha_max".into(), scalar(0.0, inverse(t.clone()))),
            ("boundary".into(), Value::String("PEC".into())),
        ],
    });
    let mut config_args = vec![
        integers(&[8]),
        tensor(&[1.0], l.clone()),
        scalar(1.0, l.clone()),
        scalar(1.0, DimExpr::Mul(vec![l.clone(), inverse(t.clone())])),
        scalar(1.0, impedance),
        scalar(1.0, electric.clone()),
        options,
        dynamic_vector(vec![]),
        dynamic_vector(vec![]),
        Value::int(1000),
    ];
    let config = call("num_analysis_yee_config_new", &config_args);
    let solver = call("num_analysis_yee_new", &[config]);
    let input_ref = call("num_backend_reference", &[scalar(1.0, t)]);
    let output_ref = call("num_backend_reference", &[scalar(1.0, current)]);
    let program = call(
        "num_backend_compile_current",
        &[
            Value::String("1".into()),
            strings(&["time"]),
            dynamic_vector(vec![input_ref]),
            output_ref,
            Value::int(10000),
        ],
    );
    call(
        "num_analysis_yee_add_source",
        &[
            solver.clone(),
            integers(&[3, 0, 0]),
            Value::String("Ey".into()),
            program.clone(),
            Value::int(1000),
        ],
    );
    call("num_backend_release", &[program]);
    call("num_analysis_yee_step", std::slice::from_ref(&solver));
    let field = call(
        "num_analysis_yee_field",
        &[solver, Value::String("Ey".into())],
    );
    let Value::Tensor(field) = field else {
        panic!("physical electric field array expected")
    };
    assert!(field.data.iter().any(|v| coherent(v).abs() > 1e-8));
    let expected = TypeExpr::Float.with_dim(electric);
    assert!(rssn::kernel::overload::matches_type(&expected, &field.elem).accepted());
    // A sphere phi=x²+y²+z²-r² is a valid Length² zero level set.
    // Only build the lightweight config: no second field-grid allocation.
    let coordinate_ref = call("num_backend_reference", &[scalar(1.0, l.clone())]);
    let level_ref = call(
        "num_backend_reference",
        &[scalar(1.0, DimExpr::Mul(vec![l.clone(), l]))],
    );
    let sphere = call(
        "num_backend_compile_current",
        &[
            Value::String("x*x+y*y+z*z-1".into()),
            strings(&["x", "y", "z"]),
            dynamic_vector(vec![
                coordinate_ref.clone(),
                coordinate_ref.clone(),
                coordinate_ref,
            ]),
            level_ref,
            Value::int(10000),
        ],
    );
    config_args[7] = dynamic_vector(vec![Value::Record(rssn::kernel::value::RecordValue {
        fields: vec![
            ("kind".into(), Value::String("pec_body".into())),
            ("program".into(), sphere),
        ],
    })]);
    let _sphere_config = call("num_analysis_yee_config_new", &config_args);
}

#[test]
fn native_dual_min_integer_power_preserves_negative_base_derivative_parity() {
    let x = rssn::kernels::autodiff::Dual::var(-1.0).powi(i32::MIN);
    assert_eq!(x.val, 1.0);
    assert_eq!(x.deriv, 2147483648.0);
    let y = rssn::kernels::autodiff::Dual::var(1.0).powi(i32::MIN);
    assert_eq!(y.deriv, -2147483648.0);
    let ordinary = rssn::kernels::autodiff::Dual::var(-2.0).powi(3);
    assert_eq!(ordinary.val, -8.0);
    assert_eq!(ordinary.deriv, 12.0);
}

#[test]
fn native_ttl_accepts_negative_binding_for_valid_unbound_constant_potential() {
    let s = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name: &str, args: &[Value]| {
        s.call_value(name, args, &ctx)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"))
    };
    let l = DimExpr::base("Length");
    let m = DimExpr::base("Mass");
    let t = DimExpr::base("Time");
    let momentum = DimExpr::Mul(vec![
        m.clone(),
        l.clone(),
        DimExpr::Inverse(Box::new(t.clone())),
    ]);
    let energy = DimExpr::Mul(vec![
        m.clone(),
        l.clone(),
        l.clone(),
        DimExpr::Inverse(Box::new(t.clone())),
        DimExpr::Inverse(Box::new(t.clone())),
    ]);
    let potential = call(
        "num_native_tape_lower",
        &[
            strings(&["-1"]),
            strings(&["x"]),
            strings(&[]),
            dynamic_vector(vec![scalar(1.0, l.clone())]),
            dynamic_vector(vec![]),
            dynamic_vector(vec![scalar(1.0, energy.clone())]),
        ],
    );
    let state = call(
        "num_analysis_symplectic_ttl_new",
        &[
            tensor(&[0.0], l.clone()),
            tensor(&[2.0], momentum.clone()),
            tensor(&[1.0], m.clone()),
            scalar(0.0, t.clone()),
            scalar(1.0, l),
            scalar(1.0, m),
            scalar(1.0, t.clone()),
        ],
    );
    let ds = scalar(0.1, DimExpr::Mul(vec![energy.clone(), t]));
    call(
        "num_analysis_symplectic_ttl_step",
        &[
            state.clone(),
            scalar(-1.0, energy.clone()),
            ds.clone(),
            Value::int(2),
            potential.clone(),
            dynamic_vector(vec![]),
            Value::int(100),
        ],
    );
    let result = call(
        "num_analysis_symplectic_ttl_data",
        std::slice::from_ref(&state),
    );
    let Value::Record(result) = result else {
        panic!("physical TTL state expected")
    };
    assert!((values(result.get("positions").unwrap())[0] - 0.2).abs() < 1e-14);
    assert!((values(result.get("momenta").unwrap())[0] - 2.0).abs() < 1e-14);
    assert!((coherent(result.get("time").unwrap()) - 0.1).abs() < 1e-14);
    // T=2, so B=-3 is outside the native positive logarithmic kinetic domain.
    let invalid = s.call_value(
        "num_analysis_symplectic_ttl_step",
        &[
            state,
            scalar(-3.0, energy),
            ds,
            Value::int(2),
            potential,
            dynamic_vector(vec![]),
            Value::int(100),
        ],
        &ctx,
    );
    assert!(
        matches!(invalid,Err(rssn::kernel::OpError::Domain(ref message)) if message.contains("T+B>0"))
    );
}
