//! Regression checks for new analysis transport and genuine native method selection.
use rssn::kernel::CallCtx;
use rssn::kernel::NativeCallback;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::manifest::SignatureSource;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::Signature;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use std::sync::Arc;
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
fn measured_chain_acceleration_follows_the_input_body_order() {
    use rssn::kernel::value::DynValue;
    let s = rssn::Session::new();
    let dims = [
        DimExpr::base("Length"),
        DimExpr::base("Mass"),
        DimExpr::base("Time"),
    ];
    let [l, m, t] = dims.clone();
    let delta = DimExpr::Mul(vec![l.clone(), t.clone().pow(-2)]);
    let measured = |v, dim| {
        Value::Dyn(DynValue {
            dim,
            value: Box::new(Value::Float(v)),
        })
    };
    let rows = |data: [f64; 6], dimension| {
        Value::Tensor(TensorValue {
            shape: vec![2, 3],
            elem: TypeExpr::Float.with_dim(dimension),
            data: data.into_iter().map(Value::Float).collect(),
        })
    };
    let arguments = [
        measured(
            8.0,
            DimExpr::Mul(vec![l.clone().pow(3), m.clone().pow(-1), t.clone().pow(-2)]),
        ),
        tensor(&[1.0, 3.0], m),
        rows([0.0, 0.0, 0.0, 2.0, 0.0, 0.0], l.clone()),
        rows([0.0; 6], DimExpr::Mul(vec![l, t.pow(-1)])),
        measured(4.0, dims[0].clone()),
        measured(2.0, dims[1].clone()),
        measured(2.0, dims[2].clone()),
    ];
    let ctx = CallCtx::new();
    let chain = s
        .call_value("num_analysis_nbody_chain_new", &arguments, &ctx)
        .unwrap();
    let result = s
        .call_value("num_analysis_nbody_chain_accelerations", &[chain], &ctx)
        .unwrap();
    let Value::Tensor(field) = &result else {
        panic!("measured acceleration tensor")
    };
    assert_eq!(field.shape, [2, 3]);
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(delta), &field.elem)
            .accepted()
    );
    assert_eq!(values(&result), [6.0, 0.0, 0.0, -2.0, 0.0, 0.0]);
}

#[test]
fn fixed_spherical_coordinates_preserve_shape_and_near_axis_angles() {
    let session = rssn::Session::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &CallCtx::new());
    let point = |values: &[f64]| {
        Value::Tensor(TensorValue {
            shape: vec![values.len()],
            elem: TypeExpr::Float,
            data: values.iter().copied().map(Value::Float).collect(),
        })
    };
    let forward = "num_analysis_coordinates_cartesian_to_spherical";
    let inverse = "num_analysis_coordinates_spherical_to_cartesian";
    let sphere = call(forward, &[point(&[1e-10, 0.0, 1.0])]).unwrap();
    assert!((values(&sphere)[1] - 1e-10).abs() < 1e-24);
    let restored = call(inverse, &[sphere]).unwrap();
    let Value::Tensor(result) = &restored else {
        panic!("Cartesian output requires three coordinates");
    };
    assert_eq!(result.shape, [3]);
    let components = values(&restored);
    assert!((components[0] - 1e-10).abs() < 1e-24);
    assert_eq!(components[2], 1.0);
    assert!(call(forward, &[point(&[1.0, 2.0])]).is_err());
    assert!(
        call(
            inverse,
            &[tensor(&[1.0, 0.0, 0.0], DimExpr::Base("Length".into()))]
        )
        .is_err()
    );
}

#[test]
fn complex_reciprocal_restores_the_inverse_dimension_at_extreme_magnitudes() {
    use rssn::kernel::value::ComplexValue;
    use rssn::kernel::value::DynValue;

    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let maximum = f64::MAX;
    let length = DimExpr::Base("Length".into());
    let input = Value::Dyn(DynValue {
        dim: length.clone(),
        value: Box::new(Value::Complex(ComplexValue::Float {
            re: maximum,
            im: maximum,
        })),
    });
    let name = "num_analysis_complex_scalar_method_recip";
    let result = session.call_value(name, &[input], &ctx).unwrap();
    let expected = TypeExpr::Complex.with_dim(length.pow(-1));
    assert!(rssn::kernel::overload::matches_type(&expected, &result.type_expr()).accepted());
    let Value::Dyn(measured) = result else {
        panic!("expected measured reciprocal");
    };
    let Value::Complex(ComplexValue::Float { re, im }) = *measured.value else {
        panic!("expected complex payload");
    };
    assert_eq!(re, maximum.recip() * 0.5);
    assert_eq!(im, -re);
    assert!(
        session
            .call_value(
                name,
                &[Value::Complex(ComplexValue::Float { re: 0.0, im: 0.0 })],
                &ctx
            )
            .is_err()
    );
}

#[test]
fn prepared_wavelet_keeps_measured_samples_and_requires_compatible_bands() {
    let session = rssn::Session::new();
    let context = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &context);
    let plan = call(
        "num_analysis_wavelets_dwt_plan_new",
        &[
            Value::String("Db4".into()),
            Value::String("Reflect".into()),
            Value::int(5),
        ],
    )
    .unwrap();
    assert_eq!(
        call(
            "num_analysis_wavelets_dwt_plan_signal_len",
            std::slice::from_ref(&plan)
        )
        .unwrap(),
        Value::int(5)
    );
    let length = DimExpr::Base("Length".into());
    let input = [2.0, -3.0, 4.0, 1.0, -1.0];
    let bands = call(
        "num_analysis_wavelets_dwt_plan_analyze",
        &[plan.clone(), tensor(&input, length.clone())],
    )
    .unwrap();
    let Value::Record(bands) = bands else {
        panic!("expected approximation and detail");
    };
    let approximation = bands.get("approximation").unwrap();
    let detail = bands.get("detail").unwrap();
    assert_eq!(
        call(
            "num_analysis_wavelets_dwt_plan_coeff_len",
            std::slice::from_ref(&plan)
        )
        .unwrap(),
        Value::int(values(detail).len())
    );
    let name = "num_analysis_wavelets_dwt_plan_synthesize";
    let result = call(name, &[plan.clone(), approximation.clone(), detail.clone()]).unwrap();
    let TypeExpr::Tensor { elem, .. } = context.value_type(&result).unwrap() else {
        panic!("expected a measured signal");
    };
    assert_eq!(*elem, TypeExpr::Float.with_dim(length));
    for (actual, expected) in values(&result).into_iter().zip(input) {
        assert!((actual - expected).abs() < 1e-12);
    }
    let wrong = tensor(&values(detail), DimExpr::Base("Time".into()));
    assert!(call(name, &[plan, approximation.clone(), wrong]).is_err());
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

#[test]
fn sampled_chebyshev_descriptor_preserves_coordinate_and_ordinate_dimensions() {
    let s = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name, args: &[Value]| s.call_value(name, args, &ctx);
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let name = "num_native_cheb_from_values";
    let samples = tensor(&[5.0, 3.0, 1.0], length.clone());
    let model = call(
        name,
        &[
            samples.clone(),
            scalar(0.0, time.clone()),
            scalar(2.0, time.clone()),
        ],
    )
    .unwrap();
    let result = call("num_native_cheb_eval", &[model, scalar(0.5, time.clone())]).unwrap();
    assert!((coherent(&result) - 2.0).abs() < 1e-13);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(length.clone()),
            &result.type_expr()
        )
        .accepted()
    );
    assert!(call(name, &[samples, scalar(0.0, time), scalar(2.0, length)]).is_err());
}

#[test]
fn prepared_pic_assignment_restores_charge_density_and_field_units() {
    let s = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name, args: &[Value]| s.call_value(name, args, &ctx);
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let charge = DimExpr::Mul(vec![DimExpr::base("Current"), time.clone()]);
    let area_permittivity = DimExpr::Mul(vec![
        charge.clone().pow(2),
        time.clone().pow(2),
        DimExpr::base("Mass").pow(-1),
        length.clone().pow(-1),
    ]);
    let grid = call(
        "num_analysis_pic_grid_new",
        &[
            Value::int(8),
            scalar(8.0, length.clone()),
            scalar(2.0, length.clone()),
            scalar(1.0, time.clone()),
            scalar(2.0, charge.clone()),
            scalar(1.0, area_permittivity),
        ],
    )
    .unwrap();
    let plan = call(
        "num_analysis_pic_assignment_new",
        &[
            grid,
            tensor(&[0.0, 1.5], length.clone()),
            Value::String("CIC".into()),
        ],
    )
    .unwrap();
    let density = call(
        "num_analysis_pic_assignment_deposit",
        &[plan.clone(), tensor(&[2.0, 4.0], charge.clone())],
    )
    .unwrap();
    assert_eq!(values(&density), [2.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    let Value::Tensor(rho) = density else {
        panic!("expected a mesh tensor");
    };
    let expected = TypeExpr::Float.with_dim(DimExpr::Mul(vec![charge.clone(), length.pow(-1)]));
    assert!(rssn::kernel::overload::matches_type(&expected, &rho.elem).accepted());
    let field = tensor(&[0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0], time.clone());
    let sampled = call(
        "num_analysis_pic_assignment_gather",
        &[plan.clone(), field, scalar(2.0, time.clone())],
    )
    .unwrap();
    assert_eq!(values(&sampled), [0.0, 1.5]);
    let Value::Tensor(sampled) = sampled else {
        panic!("expected a particle tensor");
    };
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(time), &sampled.elem)
            .accepted()
    );
    assert!(
        call(
            "num_analysis_pic_assignment_deposit",
            &[plan, tensor(&[2.0], charge)]
        )
        .is_err()
    );
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

#[test]
fn prepared_newton_cg_reuses_a_tape_after_the_original_is_released() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name: &str, values: &[Value]| session.call_value(name, values, &ctx);
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let mass = DimExpr::base("Mass");
    let tape = call(
        "num_native_tape_lower",
        &[
            strings(&["(x-p)*(x-p)+(y-2)*(y-2)"]),
            strings(&["x", "y"]),
            strings(&["p"]),
            dynamic_vector(vec![scalar(1.0, length.clone()), scalar(1.0, time.clone())]),
            dynamic_vector(vec![scalar(1.0, length.clone())]),
            dynamic_vector(vec![scalar(5.0, mass.clone())]),
        ],
    )
    .unwrap();
    let make = "num_analysis_newton_cg_tape_plan_new";
    let run = "num_analysis_newton_cg_tape_plan_minimize";
    assert!(call(make, &[tape.clone(), Value::int(1)]).is_err());
    let plan = call(make, &[tape.clone(), Value::int(0)]).unwrap();
    let start = dynamic_vector(vec![scalar(0.0, length.clone()), scalar(0.0, time.clone())]);
    let params = |target| dynamic_vector(vec![scalar(target, length.clone())]);
    let options = Value::Record(rssn::kernel::value::RecordValue {
        fields: vec![
            ("gradient_tolerance".into(), Value::Float(1e-10)),
            ("max_iterations".into(), Value::int(24)),
            ("max_cg_iterations".into(), Value::int(8)),
            ("cg_relative_tolerance".into(), Value::Float(0.1)),
            ("initial_radius".into(), Value::Float(1.0)),
            ("max_radius".into(), Value::Float(10.0)),
            ("min_radius".into(), Value::Float(1e-10)),
            ("acceptance_threshold".into(), Value::Float(0.1)),
        ],
    });
    let source_result = call(
        "num_analysis_newton_cg_tape",
        &[
            tape.clone(),
            Value::int(0),
            start.clone(),
            params(3.0),
            options.clone(),
        ],
    )
    .unwrap();
    call("num_native_release", &[tape]).unwrap();
    let first = call(
        run,
        &[plan.clone(), start.clone(), params(3.0), options.clone()],
    )
    .unwrap();
    assert_eq!(first, source_result);
    let wrong = dynamic_vector(vec![scalar(4.0, time.clone())]);
    assert!(call(run, &[plan.clone(), start.clone(), wrong, options.clone()]).is_err());
    let args = [plan.clone(), start, params(4.0), options];
    let foreign = rssn::Session::new();
    assert!(foreign.call_value(run, &args, &ctx).is_err());
    let second = call(run, &args).unwrap();
    let Value::Record(second) = second else {
        panic!("expected a solution record")
    };
    assert_eq!(
        second.get("status"),
        Some(&Value::String("FirstOrderConverged".into()))
    );
    let Some(Value::Tensor(point)) = second.get("x") else {
        panic!("expected the measured optimum")
    };
    let [Value::Dyn(x), Value::Dyn(y)] = point.data.as_slice() else {
        panic!("expected two physical coordinates");
    };
    assert_eq!(&x.dim, &length);
    assert_eq!(&y.dim, &time);
    assert!((coherent(&point.data[0]) - 4.0).abs() < 1e-9);
    assert!((coherent(&point.data[1]) - 2.0).abs() < 1e-9);
    let value = second.get("value").unwrap();
    let Value::Dyn(measured_value) = value else {
        panic!("expected a physical objective value");
    };
    assert_eq!(&measured_value.dim, &mass);
    assert!(coherent(value) < 1e-18);
    call("num_native_release", &[plan]).unwrap();
    assert!(call(run, &args).is_err());
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
fn native_fused_regression_restores_loss_and_each_gradient_dimension() {
    let s = rssn::Session::new();
    let ctx = CallCtx::new();
    let time = DimExpr::base("Time");
    let length = DimExpr::base("Length");
    let x = Value::Tensor(TensorValue {
        shape: vec![3, 1],
        elem: TypeExpr::Float.with_dim(time.clone()),
        data: vec![Value::Float(2.0), Value::Float(4.0), Value::Float(6.0)],
    });
    let model = s
        .call_value(
            "num_analysis_linear_regression_new",
            &[
                x,
                tensor(&[8.0, 12.0, 16.0], length.clone()),
                scalar(2.0, time.clone()),
                scalar(4.0, length.clone()),
            ],
            &ctx,
        )
        .unwrap();
    let coefficients = dynamic_vector(vec![
        scalar(0.0, length.clone()),
        scalar(
            2.0,
            DimExpr::Mul(vec![length.clone(), time.clone().pow(-1)]),
        ),
    ]);
    let value_args = [model.clone(), coefficients];
    let Value::Record(result) = s
        .call_value(
            "num_analysis_regression_value_and_gradient",
            &value_args,
            &ctx,
        )
        .unwrap()
    else {
        panic!("expected fused regression record")
    };
    let value = result.get("value").unwrap();
    assert_eq!(coherent(value), 8.0);
    assert_eq!(
        *value,
        s.call_value("num_analysis_regression_cost", &value_args, &ctx)
            .unwrap()
    );
    let Value::Tensor(gradient) = result.get("gradient").unwrap() else {
        panic!("expected measured gradient")
    };
    assert_eq!(gradient.shape, vec![2]);
    for (entry, expected, dimension) in [
        (&gradient.data[0], -4.0, length.clone()),
        (
            &gradient.data[1],
            -16.0,
            DimExpr::Mul(vec![length.clone(), time]),
        ),
    ] {
        assert_eq!(coherent(entry), expected);
        let Value::Dyn(measured) = entry else {
            panic!("gradient unit erased")
        };
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Float.with_dim(dimension),
                &TypeExpr::Float.with_dim(measured.dim.clone())
            )
            .accepted()
        );
    }
    assert!(
        s.call_value(
            "num_analysis_regression_value_and_gradient",
            &[
                model,
                dynamic_vector(vec![scalar(0.0, length.clone()), scalar(2.0, length)])
            ],
            &ctx,
        )
        .is_err()
    );
}

#[test]
fn sampled_implicit_sde_preserves_states_units_and_callback_budget() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name: &str, args: &[Value]| {
        session
            .call_value(name, args, &ctx)
            .unwrap_or_else(|error| panic!("{name}: {error}"))
    };
    let time = DimExpr::base("Time");
    let length = DimExpr::base("Length");
    let reference = |dim| scalar(1.0, dim);
    let rate = DimExpr::Mul(vec![
        length.clone(),
        DimExpr::Inverse(Box::new(time.clone())),
    ]);
    let diffusion = DimExpr::Mul(vec![
        length.clone(),
        DimExpr::Inverse(Box::new(time.clone().pow(Exp::Rat(Rational::new(1, 2))))),
    ]);
    let tape = |expression, output| {
        call(
            "num_native_tape_lower",
            &[
                strings(&[expression]),
                strings(&["t", "x"]),
                strings(&[]),
                dynamic_vector(vec![reference(time.clone()), reference(length.clone())]),
                dynamic_vector(vec![]),
                dynamic_vector(vec![reference(output)]),
            ],
        )
    };
    let drift = call(
        "num_analysis_ode_rhs_new",
        &[tape("1", rate), dynamic_vector(vec![])],
    );
    let model = call(
        "num_analysis_sde_model_new",
        &[
            drift,
            tape("0", diffusion),
            dynamic_vector(vec![]),
            Value::int(1),
        ],
    );
    let rng = call("num_native_rng_new", &[Value::int(42)]);
    let noise = call(
        "num_analysis_sde_sample_noise",
        &[
            rng,
            Value::int(1),
            Value::int(5),
            scalar(0.25, time.clone()),
            scalar(2.0, time.clone()),
            Value::int(2),
        ],
    );
    let args = [
        model,
        noise,
        dynamic_vector(vec![scalar(3.0, length.clone())]),
        scalar(0.0, time),
        Value::String("BackwardEuler".into()),
        Value::Float(1.0),
    ];
    let full = call(
        "num_analysis_sde_drift_implicit_core",
        &[args.as_slice(), &[Value::int(200)]].concat(),
    );
    let sampled =
        |stride, budget| [args.as_slice(), &[Value::int(stride), Value::int(budget)]].concat();
    let op = "num_analysis_sde_drift_implicit_sampled";
    let saved = call(op, &sampled(2, 200));
    let (Value::Tensor(full), Value::Tensor(saved)) = (full, saved) else {
        panic!("expected state sequences")
    };
    assert_eq!(saved.shape, [4]);
    for (state, index) in saved.data.iter().zip([0, 2, 4, 5]) {
        assert_eq!(state, &full.data[index]);
        let Value::Tensor(state) = state else {
            panic!("expected a physical state")
        };
        let Value::Dyn(component) = &state.data[0] else {
            panic!("expected a dimensioned state")
        };
        assert_eq!(component.dim, length);
    }
    assert!(session.call_value(op, &sampled(0, 200), &ctx).is_err());
    assert!(session.call_value(op, &sampled(2, 1), &ctx).is_err());
}

#[test]
fn gauss_kronrod_breakpoints_expose_declared_physical_callback_contracts() {
    let op = "num_analysis_quadrature_gk_breakpoints_scaled";
    let session = rssn::Session::new();
    let manifest = session.manifest();
    let descriptor = manifest.ops.iter().find(|entry| entry.name == op).unwrap();
    assert_eq!(descriptor.signature_source, SignatureSource::Declared);
    assert_eq!(
        descriptor.rust_path.as_deref(),
        Some("crate::kernels::quadrature::integrate_gk_breakpoints"),
    );
    assert!(descriptor.doc.contains("known interior breakpoints"));
    assert!(descriptor.higher_order);

    let mut ctx = CallCtx::new();
    let callback = ctx
        .add_callback(Arc::new(NativeCallback::new(
            |args: &[Value]| {
                let x = coherent(&args[0]);
                Ok(Value::Float(if 0.25 < x && x < 0.5 {
                    1.0
                } else {
                    0.0
                }))
            },
            Signature::new(vec![TypeExpr::Float], TypeExpr::Float),
        )))
        .unwrap();
    let length = DimExpr::Base("Length".into());
    let time = DimExpr::Base("Time".into());
    let args = |breakpoints: [f64; 2], budget| {
        [
            tensor(
                &[0.0, 1.0, 1e-13, 1e-13, breakpoints[0], breakpoints[1]],
                DimExpr::Dimless,
            ),
            integers(&[0, 3]),
            scalar(2.0, time.clone()),
            scalar(3.0, length.clone()),
            Value::int(budget),
        ]
    };
    let result = session
        .call_ho(op, &args([0.25, 0.5], 45), &[callback], &ctx)
        .unwrap();
    let Value::Record(result) = result else {
        panic!("expected the integral report")
    };
    let Value::Dyn(value) = result.get("value").unwrap() else {
        panic!("expected physical value")
    };
    assert!((coherent(&value.value) - 1.5).abs() < 1e-13);
    assert_eq!(value.dim, DimExpr::Mul(vec![time.clone(), length.clone()]));
    assert_eq!(result.get("evaluations"), Some(&Value::int(45)));
    let Value::Variant(converged) = result.get("converged").unwrap() else {
        panic!("expected a convergence flag")
    };
    assert_eq!(converged.tag, "Some");
    assert_eq!(converged.payload, [Value::Bool(true)]);
    assert!(matches!(
        session.call_ho(op, &args([0.25, 0.5], 44), &[callback], &ctx),
        Err(OpError::BudgetExceeded),
    ));
    assert!(matches!(
        session.call_ho(op, &args([0.5, 0.25], 45), &[callback], &ctx),
        Err(OpError::Domain(message)) if message.contains("strictly increasing"),
    ));
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
    let adjacent = 1.0f64.next_up();
    let short_interval = |step, maxeval| {
        session.call_value(
            "num_analysis_physics_dormand",
            &[
                rhs.clone(),
                initial.clone(),
                scalar(1.0, time.clone()),
                scalar(adjacent, time.clone()),
                scalar(step, time.clone()),
                Value::Float(1e-6),
                Value::Float(1e-9),
                Value::int(maxeval),
            ],
            &ctx,
        )
    };
    assert!(matches!(
        short_interval(2e-20, 7),
        Err(OpError::Domain(message)) if message.contains("stopped before requested endpoint")
    ));
    assert!(matches!(
        short_interval(adjacent - 1.0, 6),
        Err(OpError::BudgetExceeded)
    ));
    let result = short_interval(adjacent - 1.0, 7).unwrap();
    let Value::Tensor(rows) = result else {
        panic!("expected two time–state rows, got {result:?}");
    };
    assert_eq!(rows.shape, [2]);
    let Value::Record(last) = &rows.data[1] else {
        panic!("expected a time–state record");
    };
    assert_eq!(
        coherent(last.get("time").unwrap()).to_bits(),
        adjacent.to_bits()
    );
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
    let final_args = |steps, method, budget| {
        [
            rhs.clone(),
            initial.clone(),
            scalar(0.0, time.clone()),
            scalar(1.0, time.clone()),
            Value::int(steps),
            Value::int(method),
            Value::int(budget),
        ]
    };
    let final_state = |steps, method, budget| {
        session.call_value(
            "num_analysis_ode_fixed_final",
            &final_args(steps, method, budget),
            &ctx,
        )
    };
    assert_eq!(final_state(32, 2, 128).unwrap(), Value::Tensor(last));
    assert!(matches!(
        final_state(32, 2, 127),
        Err(rssn::kernel::OpError::BudgetExceeded)
    ));
    assert!(final_state(0, 2, 128).is_err());
    assert!(final_state(32, 3, 128).is_err());
    let Value::Tensor(long_run) = final_state(24_000, 0, 24_000).unwrap() else {
        panic!("expected the final mixed physical state after a larger run");
    };
    assert!((coherent(&long_run.data[0]) - 1.0f64.cos()).abs() < 2e-5);
    assert!((coherent(&long_run.data[1]) + 1.0f64.sin()).abs() < 2e-5);
    let mut wrong_endpoint = final_args(32, 2, 128);
    wrong_endpoint[3] = scalar(1.0, length.clone());
    assert!(
        session
            .call_value("num_analysis_ode_fixed_final", &wrong_endpoint, &ctx)
            .is_err()
    );
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
fn checked_explicit_euler_restores_signed_times_and_exact_callback_budgets() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let time = DimExpr::base("Time");
    let length = DimExpr::base("Length");
    let velocity = DimExpr::Mul(vec![length.clone(), time.clone().pow(-1)]);
    let empty = dynamic_vector(vec![]);
    let tape = session
        .call_value(
            "num_native_tape_lower",
            &[
                strings(&["1.0"]),
                strings(&["t", "y"]),
                strings(&[]),
                dynamic_vector(vec![scalar(1.0, time.clone()), scalar(1.0, length.clone())]),
                empty.clone(),
                dynamic_vector(vec![scalar(1.0, velocity)]),
            ],
            &ctx,
        )
        .unwrap();
    let rhs = session
        .call_value("num_analysis_ode_rhs_new", &[tape, empty], &ctx)
        .unwrap();
    let solve = |method: &str, step, end, budget| {
        session.call_value(
            "num_analysis_physics_explicit_euler",
            &[
                rhs.clone(),
                dynamic_vector(vec![scalar(4.0, length.clone())]),
                scalar(2.0, time.clone()),
                scalar(end, time.clone()),
                scalar(step, time.clone()),
                Value::String(method.into()),
                Value::int(budget),
            ],
            &ctx,
        )
    };
    for (method, budget) in [("Forward", 3), ("Midpoint", 6), ("Heun", 6)] {
        assert!(matches!(
            solve(method, -0.75, 0.0, budget - 1),
            Err(OpError::BudgetExceeded)
        ));
        let Value::Tensor(rows) = solve(method, -0.75, 0.0, budget).unwrap() else {
            panic!("expected a measured history");
        };
        assert_eq!(rows.shape, [4]);
        for (row, (t, y)) in
            rows.data
                .iter()
                .zip([(2.0, 4.0), (1.25, 3.25), (0.5, 2.5), (0.0, 2.0)])
        {
            let Value::Record(row) = row else {
                panic!("expected time and state")
            };
            assert_eq!(coherent(row.get("time").unwrap()), t);
            assert_eq!(values(row.get("state").unwrap()), [y]);
        }
    }
    let Value::Tensor(zero_span) = solve("Heun", 0.125, 2.0, 1).unwrap() else {
        panic!("expected initial state only");
    };
    assert_eq!(zero_span.shape, [1]);
    assert!(solve("Midpoint", 0.75, 0.0, 6).is_err());
    assert!(solve("Forward", -1e-7, 0.0, 1).is_err());
    assert!(solve("Invalid", -0.75, 0.0, 6).is_err());
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
        charge.clone(),
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
    let input_ref = call("num_backend_reference", &[scalar(1.0, t.clone())]);
    let output_ref = call("num_backend_reference", &[scalar(1.0, current)]);
    let program = call(
        "num_backend_compile_current",
        &[
            Value::String("1".into()),
            strings(&["time"]),
            dynamic_vector(vec![input_ref.clone()]),
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
        &[solver.clone(), Value::String("Ey".into())],
    );
    let Value::Tensor(field) = field else {
        panic!("physical electric field array expected")
    };
    assert!(field.data.iter().any(|v| coherent(v).abs() > 1e-8));
    let expected = TypeExpr::Float.with_dim(electric.clone());
    assert!(rssn::kernel::overload::matches_type(&expected, &field.elem).accepted());
    let sample_name = "num_analysis_yee_plane_wave_sample";
    let absent = call(sample_name, &[solver.clone(), Value::int(0)]);
    assert!(matches!(absent, Value::Variant(v) if v.tag == "None"));
    let electric_ref = call("num_backend_reference", &[scalar(1.0, electric.clone())]);
    let wave = call(
        "num_backend_compile_current",
        &[
            Value::String("2".into()),
            strings(&["time"]),
            dynamic_vector(vec![input_ref]),
            electric_ref,
            Value::int(10000),
        ],
    );
    call(
        "num_analysis_yee_set_plane_wave",
        &[
            solver.clone(),
            Value::int(2),
            Value::int(5),
            Value::Float(0.0),
            wave,
            Value::int(10),
            Value::int(1000),
        ],
    );
    call("num_analysis_yee_step", std::slice::from_ref(&solver));
    let sample = call(sample_name, &[solver.clone(), Value::int(0)]);
    let Value::Variant(sample) = sample else {
        panic!("expected optional wave sample");
    };
    assert_eq!(sample.tag, "Some");
    let [Value::Record(sample)] = sample.payload.as_slice() else {
        panic!("expected wave fields");
    };
    let state = call("num_analysis_yee_state", std::slice::from_ref(&solver));
    let Value::Record(state) = state else {
        panic!("expected Yee state");
    };
    let time = coherent(state.get("time").unwrap());
    let dt = coherent(state.get("dt").unwrap());
    let amplitude = sample.get("electric").unwrap();
    assert_eq!(coherent(amplitude), 2.0);
    assert_eq!(coherent(sample.get("electric_time").unwrap()), time);
    assert_eq!(
        coherent(sample.get("magnetic_time").unwrap()),
        time - 0.5 * dt
    );
    let magnetic = DimExpr::Mul(vec![charge, inverse(t.clone()), inverse(l.clone())]);
    for (field, dim) in [
        ("electric", electric),
        ("magnetic", magnetic),
        ("electric_time", t.clone()),
        ("magnetic_time", t),
    ] {
        let expected = TypeExpr::Float.with_dim(dim);
        assert!(
            rssn::kernel::overload::matches_type(
                &expected,
                &sample.get(field).unwrap().type_expr()
            )
            .accepted()
        );
    }
    let outside = call(sample_name, &[solver, Value::int(10_000)]);
    assert!(matches!(outside, Value::Variant(v) if v.tag == "None"));
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
fn sampled_symplectic_retains_step_indices_and_mixed_state_dimensions() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |name: &str, args: &[Value]| {
        session
            .call_value(name, args, &ctx)
            .unwrap_or_else(|e| panic!("{name}: {e:?}"))
    };
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let speed = DimExpr::Mul(vec![
        length.clone(),
        DimExpr::Inverse(Box::new(time.clone())),
    ]);
    let rate = DimExpr::Mul(vec![
        speed.clone(),
        DimExpr::Inverse(Box::new(time.clone())),
    ]);
    let tape = call(
        "num_native_tape_lower",
        &[
            strings(&["-x"]),
            strings(&["x"]),
            strings(&[]),
            dynamic_vector(vec![scalar(1.0, length.clone())]),
            dynamic_vector(vec![]),
            dynamic_vector(vec![scalar(1.0, rate)]),
        ],
    );
    let acceleration = call(
        "num_analysis_mechanics_acceleration_new",
        &[tape, dynamic_vector(vec![]), scalar(1.0, time.clone())],
    );
    let args = [
        acceleration,
        dynamic_vector(vec![
            scalar(1.0, length.clone()),
            scalar(0.0, speed.clone()),
        ]),
        scalar(-0.1, time),
        Value::int(5),
        Value::String("Verlet".into()),
        Value::int(2),
        Value::int(30),
    ];
    let op = "num_analysis_symplectic_integrate_sampled";
    let Value::Record(result) = call(op, &args) else {
        panic!("expected sampled positions and velocities");
    };
    assert_eq!(result.get("steps"), Some(&integers(&[0, 2, 4, 5])));
    let Some(Value::Tensor(states)) = result.get("states") else {
        panic!("expected the ordered state vector");
    };
    assert_eq!(states.shape, [4]);
    assert_eq!(values(&states.data[0]), [1.0, 0.0]);
    let Value::Tensor(final_state) = &states.data[3] else {
        panic!("expected the final state");
    };
    let [Value::Dyn(q), Value::Dyn(p)] = final_state.data.as_slice() else {
        panic!("expected physical position and velocity");
    };
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(length),
            &TypeExpr::Float.with_dim(q.dim.clone()),
        )
        .accepted()
    );
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(speed),
            &TypeExpr::Float.with_dim(p.dim.clone()),
        )
        .accepted()
    );
    assert!((coherent(&final_state.data[0]) - 0.5_f64.cos()).abs() < 1e-3);
    assert!((coherent(&final_state.data[1]) - 0.5_f64.sin()).abs() < 1e-3);
    let mut invalid = args.clone();
    invalid[5] = Value::int(0);
    assert!(session.call_value(op, &invalid, &ctx).is_err());
    let mut limited = args;
    limited[6] = Value::int(1);
    assert!(matches!(
        session.call_value(op, &limited, &ctx),
        Err(OpError::BudgetExceeded)
    ));
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

#[test]
fn prepared_l2_measures_independently_dimensioned_fields_and_keeps_extreme_normalization() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |op, arguments: &[Value]| session.call_value(op, arguments, &ctx).unwrap();
    let new = "num_functional_l2_grid_new";
    let norm = "num_functional_l2_grid_norm";
    let normalize = "num_functional_l2_grid_norm_and_normalize";
    let time = DimExpr::base("Time");
    let length = DimExpr::base("Length");
    let half_time = time.clone().pow(Exp::Rat(Rational::new(1, 2)));
    let grid = call(new, &[tensor(&[0.0, 1.0, 4.0], time.clone())]);
    for (scale, dim) in [(3.0, length.clone()), (4.0, DimExpr::base("Mass"))] {
        let result = call(norm, &[grid.clone(), tensor(&[scale; 3], dim.clone())]);
        let Value::Dyn(result) = result else {
            panic!("expected measured norm")
        };
        assert!((coherent(&result.value) - 2.0 * scale).abs() < 1e-14);
        assert!(
            rssn::kernel::overload::matches_type(
                &TypeExpr::Float.with_dim(DimExpr::Mul(vec![dim, half_time.clone()])),
                &TypeExpr::Float.with_dim(result.dim),
            )
            .accepted()
        );
    }
    for (width, scale, norm_tag, expected) in [
        (1e308, 1e308, "None", 1e-154),
        (1e-308, 1e-308, "Some", 1e154),
        (0.0, 1.0, "Some", 0.0),
    ] {
        let grid = call(new, &[tensor(&[0.0, width], time.clone())]);
        let out = call(normalize, &[grid, tensor(&[scale; 2], length.clone())]);
        let Value::Record(out) = out else {
            panic!("expected a norm record")
        };
        let Some(Value::Variant(norm)) = out.get("norm") else {
            panic!("expected a norm option")
        };
        let Some(Value::Variant(normalized)) = out.get("normalized") else {
            panic!("expected a normalized option")
        };
        assert_eq!(norm.tag, norm_tag);
        if width == 0.0 {
            assert_eq!(normalized.tag, "None");
        } else {
            assert_eq!(normalized.tag, "Some");
            let Value::Tensor(result) = &normalized.payload[0] else {
                panic!("expected normalized ordinates")
            };
            assert!(
                rssn::kernel::overload::matches_type(
                    &TypeExpr::Float.with_dim(half_time.clone().pow(-1)),
                    &result.elem
                )
                .accepted()
            );
            assert!(
                result
                    .data
                    .iter()
                    .all(|x| (coherent(x) / expected - 1.0).abs() < 2e-15)
            );
        }
    }
    assert!(
        session
            .call_value(norm, &[grid, tensor(&[1.0], length)], &ctx)
            .is_err()
    );
}

#[test]
fn native_jump_batch_preserves_stream_order_and_rolls_back_on_error() {
    let session = rssn::Session::new();
    let ctx = CallCtx::new();
    let call = |op, arguments: &[Value]| session.call_value(op, arguments, &ctx);
    let length = DimExpr::base("Length");
    let make_dist = |low: &[f64], high: &[f64], scale: f64| {
        call(
            "num_analysis_stochastic_jump_distribution",
            &[
                Value::String("Uniform".into()),
                tensor(low, DimExpr::Dimless),
                tensor(high, DimExpr::Dimless),
                tensor(&[], DimExpr::Dimless),
                scalar(scale, length.clone()),
                Value::int(low.len()),
            ],
        )
        .unwrap()
    };
    let distribution = make_dist(&[0.0, 1.0], &[1.0, 2.0], 2.0);
    let overflows = make_dist(&[2.0, 2.0], &[3.0, 3.0], f64::MAX);
    let make_stream = || call("num_analysis_stochastic_stream_new", &[Value::int(1729)]).unwrap();
    let batched_stream = make_stream();
    let scalar_stream = make_stream();
    let batch = "num_analysis_stochastic_jump_sample_batch";
    let next = "num_analysis_stochastic_next_unit";
    let empty = call(
        batch,
        &[distribution.clone(), batched_stream.clone(), Value::int(0)],
    )
    .unwrap();
    let Value::Tensor(empty) = empty else {
        panic!("expected a rank-two empty batch")
    };
    assert_eq!(empty.shape, [0, 2]);
    assert!(empty.data.is_empty());
    assert!(
        call(
            batch,
            &[
                distribution.clone(),
                batched_stream.clone(),
                Value::int(262_144)
            ]
        )
        .is_err()
    );
    assert!(call(batch, &[overflows, batched_stream.clone(), Value::int(2)]).is_err());
    assert_eq!(
        call(next, &[batched_stream.clone()]).unwrap(),
        call(next, &[scalar_stream.clone()]).unwrap()
    );
    let out = call(
        batch,
        &[distribution.clone(), batched_stream.clone(), Value::int(5)],
    )
    .unwrap();
    let Value::Tensor(out) = out else {
        panic!("expected a rank-two sample batch")
    };
    assert_eq!(out.shape, [5, 2]);
    assert_eq!(out.elem, TypeExpr::Float.with_dim(length));
    let expected = (0..5)
        .flat_map(|_| {
            let row = call(
                "num_analysis_stochastic_jump_sample",
                &[distribution.clone(), scalar_stream.clone()],
            )
            .unwrap();
            let Value::Tensor(row) = row else {
                panic!("expected one sample row")
            };
            row.data
        })
        .collect::<Vec<_>>();
    assert_eq!(out.data, expected);
    assert_eq!(
        call(next, &[batched_stream]).unwrap(),
        call(next, &[scalar_stream]).unwrap()
    );
}
