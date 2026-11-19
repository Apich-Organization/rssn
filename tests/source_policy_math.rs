//! Source policies preserve measured function data and fixed SI physical laws.
#[path = "support/dimensions.rs"]
mod test_dimensions;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
fn dim(name: &str) -> DimExpr {
    DimExpr::Base(name.into())
}
fn scalar(
    value: f64,
    dimension: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(value)),
        dim: dimension,
    })
}
fn vector(
    data: &[f64],
    dimension: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        data: data.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), dimension),
    })
}
fn sample(
    x: &[f64],
    y: &[f64],
    ordinate: DimExpr,
) -> Value {
    Value::Record(RecordValue {
        fields: vec![
            ("item0".into(), vector(x, dim("Time"))),
            ("item1".into(), vector(y, ordinate)),
        ],
    })
}
fn call(
    name: &str,
    args: &[Value],
) -> Value {
    Session::new()
        .call_value(name, args, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{name}: {e:?}"))
}
fn units(value: &Value) -> DimExpr {
    match value {
        | Value::Dyn(v) => v.dim.clone(),
        | Value::Tensor(t) => {
            match &t.elem {
                | TypeExpr::Dim(_, d) => d.clone(),
                | _ => DimExpr::Dimless,
            }
        },
        | _ => DimExpr::Dimless,
    }
}
fn assert_units(
    value: &Value,
    expected: DimExpr,
) {
    assert_eq!(
        test_dimensions::factors(&units(value)),
        test_dimensions::factors(&expected)
    );
}
fn number(value: &Value) -> f64 {
    match value {
        | Value::Float(v) => *v,
        | Value::Dyn(v) => number(&v.value),
        | _ => panic!("number: {value:?}"),
    }
}
fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(r) = value else {
        panic!("record")
    };
    &r.fields.iter().find(|(n, _)| n == name).unwrap().1
}
#[test]
fn sampled_norms_projection_and_normalization_keep_independent_units() {
    let f = sample(&[0., 1., 2.], &[2., 2., 2.], dim("Length"));
    let l1 = call("num_functional_l1_norm", std::slice::from_ref(&f));
    assert_eq!(number(&l1), 4.);
    assert_units(&l1, DimExpr::Mul(vec![dim("Time"), dim("Length")]));
    let l2 = call("num_functional_l2_norm", std::slice::from_ref(&f));
    assert!((number(&l2) - 8f64.sqrt()).abs() < 1e-14);
    assert_units(
        &l2,
        DimExpr::Mul(vec![
            dim("Time").pow(Exp::Rat(Rational::new(1, 2))),
            dim("Length"),
        ]),
    );
    let g = sample(&[0., 1., 2.], &[3., 3., 3.], dim("Mass"));
    let projected = call("num_functional_project", &[f.clone(), g]);
    assert_units(field(&projected, "item0"), dim("Time"));
    assert_units(field(&projected, "item1"), dim("Length"));
    let normalized = call("num_functional_normalize", &[f]);
    assert_units(
        field(&normalized, "item1"),
        dim("Time").pow(Exp::Rat(Rational::new(-1, 2))),
    );
    assert!(
        Session::new()
            .call_value(
                "num_functional_inner_product",
                &[
                    sample(&[0., 1., 2.], &[1., 1., 1.], dim("Length")),
                    sample(&[0., 1.1, 2.], &[1., 1., 1.], dim("Mass"))
                ],
                &CallCtx::new()
            )
            .is_err()
    );
}
#[test]
fn fixed_si_physics_checks_inputs_and_restores_exact_output_dimensions() {
    let wavelength = scalar(500e-9, dim("Length"));
    let energy = call("num_classical_photon_energy", &[wavelength]);
    assert_units(
        &energy,
        DimExpr::Mul(vec![
            dim("Mass"),
            dim("Length").pow(Exp::Rat(Rational::new(2, 1))),
            dim("Time").pow(Exp::Rat(Rational::new(-2, 1))),
        ]),
    );
    assert!((number(&energy) / 3.972891714297857e-19 - 1.).abs() < 1e-12);
    assert!(
        Session::new()
            .call_value(
                "num_classical_photon_energy",
                &[scalar(1., dim("Time"))],
                &CallCtx::new()
            )
            .is_err()
    );
    let tiny_photon = call(
        "num_classical_photon_energy",
        &[scalar(1e-25, dim("Length"))],
    );
    assert!((number(&tiny_photon) / 1.9864458571489286 - 1.0).abs() < 1e-13);
    let tiny_temperature = call(
        "num_classical_wien_displacement_wavelength",
        &[scalar(1e-15, dim("Temperature"))],
    );
    assert!((number(&tiny_temperature) / 2.897771955e12 - 1.0).abs() < 1e-14);
    assert!(rssn::sim::classical::heisenberg_position_uncertainty(0.0).is_infinite());
    assert!(rssn::sim::classical::heisenberg_position_uncertainty(-1.0).is_nan());
    let uncertainty = call(
        "num_classical_heisenberg_position_uncertainty",
        &[scalar(
            1e-45,
            DimExpr::Mul(vec![
                dim("Mass"),
                dim("Length"),
                dim("Time").pow(Exp::Rat(Rational::new(-1, 1))),
            ]),
        )],
    );
    assert!((number(&uncertainty) / 5.272859085e10 - 1.0).abs() < 1e-14);
    assert_units(&uncertainty, dim("Length"));
    let gas = call(
        "num_classical_ideal_gas_pressure",
        &[
            scalar(1., dim("Amount")),
            scalar(300., dim("Temperature")),
            scalar(1., dim("Length").pow(Exp::Rat(Rational::new(3, 1)))),
        ],
    );
    assert!((number(&gas) - 2494.3387854).abs() < 1e-5);
}
#[test]
fn source_integer_and_quadrature_policies_are_executable() {
    assert_eq!(
        call("num_number_theory_gcd", &[Value::int(84), Value::int(30)]),
        Value::int(6)
    );
    let power = call(
        "num_number_theory_mod_pow",
        &[Value::int(1u128 << 100), Value::int(2), Value::int(17)],
    );
    assert_eq!(power, Value::int(1));
    let quantile = call(
        "num_robust_quantile",
        &[
            vector(&[1., 2., 3., 4.], dim("Length")),
            Value::Float(0.5),
            Value::int(7),
        ],
    );
    assert_eq!(number(&quantile), 2.5);
    assert_units(&quantile, dim("Length"));
    let rule = call("num_quadrature_gauss_legendre_rule", &[Value::int(3)]);
    let Value::Tensor(weights) = field(&rule, "item1") else {
        panic!("weights")
    };
    assert!((weights.data.iter().map(number).sum::<f64>() - 2.).abs() < 1e-13);
}

#[test]
fn weighted_statistics_preserve_measurement_units() {
    let samples = vector(&[1., 2., 3.], dim("Length"));
    let weights = vector(&[1., 1., 2.], DimExpr::Dimless);
    let quantile = call(
        "num_robust_weighted_quantile",
        &[samples.clone(), weights.clone(), Value::Float(0.5)],
    );
    assert_eq!(number(&quantile), 2.0);
    assert_units(&quantile, dim("Length"));
    let ecdf = call(
        "num_robust_weighted_ecdf",
        &[samples.clone(), weights.clone(), scalar(2., dim("Length"))],
    );
    assert_eq!(number(&ecdf), 0.5);
    assert_units(&ecdf, DimExpr::Dimless);

    let session = Session::new();
    assert!(
        session
            .call_value(
                "num_robust_weighted_quantile",
                &[samples.clone(), samples.clone(), Value::Float(0.5)],
                &CallCtx::new(),
            )
            .is_err()
    );
    assert!(
        session
            .call_value(
                "num_robust_weighted_ecdf",
                &[samples, weights, scalar(2., dim("Time"))],
                &CallCtx::new(),
            )
            .is_err()
    );
}

#[test]
fn native_integer_polynomial_product_reports_overflow() {
    let integers = |data: &[i64]| {
        Value::Tensor(TensorValue {
            shape: vec![data.len()],
            data: data.iter().copied().map(Value::int).collect(),
            elem: TypeExpr::Int,
        })
    };
    let product = call(
        "num_fftconv_poly_mul_i64_checked",
        &[integers(&[i64::MAX, i64::MAX]), integers(&[1, -1])],
    );
    let Value::Variant(product) = product else {
        panic!("expected the checked product variant")
    };
    let [Value::Tensor(coefficients)] = product.payload.as_slice() else {
        panic!("expected a coefficient vector")
    };
    assert_eq!(product.tag, "Some");
    assert_eq!(
        coefficients.data,
        [Value::int(i64::MAX), Value::int(0), Value::int(-i64::MAX)]
    );
    let overflow = call(
        "num_fftconv_poly_mul_i64_checked",
        &[integers(&[i64::MAX]), integers(&[2])],
    );
    assert!(matches!(overflow, Value::Variant(v) if v.tag == "None" && v.payload.is_empty()));
}

#[test]
fn native_rng_jumps_mutate_the_session_owned_stream() {
    let session = Session::new();
    let context = CallCtx::new();
    for (jump, next) in [
        ("num_native_rng_jump", 0x2107_D23F_5380_538B_u64),
        ("num_native_rng_long_jump", 0x7089_19B1_47F7_8AF3),
    ] {
        let handle = session
            .call_value("num_native_rng_new", &[Value::int(0)], &context)
            .unwrap();
        assert_eq!(
            session
                .call_value(jump, std::slice::from_ref(&handle), &context)
                .unwrap(),
            Value::Bool(true)
        );
        assert_eq!(
            session
                .call_value("num_native_rng_next_u64", &[handle], &context)
                .unwrap(),
            Value::int(next)
        );
    }
}

#[test]
fn explicit_correlation_radii_share_coordinate_units() {
    let name = "num_fractal_geometry_and_chaos_correlation_dimension_with_radii";
    let points = |rows, data| {
        Value::Tensor(TensorValue {
            shape: vec![rows, 2],
            data,
            elem: TypeExpr::Dim(Box::new(TypeExpr::Float), dim("Length")),
        })
    };
    let samples = points(
        3,
        [0., 0., 1., 0., 3., 0.]
            .into_iter()
            .map(Value::Float)
            .collect(),
    );
    let length = vector(&[2., 4.], dim("Length"));
    let result = call(name, &[samples.clone(), length.clone()]);
    let Value::Variant(result) = result else {
        panic!("expected a correlation estimate")
    };
    let [estimate] = result.payload.as_slice() else {
        panic!("expected two radii to determine a slope")
    };
    assert_eq!(result.tag, "Some");
    assert!((number(estimate) - 3_f64.ln() / 2_f64.ln()).abs() < 1e-14);
    assert_units(estimate, DimExpr::Dimless);

    let session = Session::new();
    let ctx = CallCtx::new();
    assert!(
        session
            .call_value(name, &[samples, vector(&[2., 4.], dim("Time"))], &ctx)
            .is_err()
    );
    assert!(
        session
            .call_value(
                name,
                &[points(1025, vec![Value::Float(0.); 2050]), length],
                &ctx
            )
            .is_err()
    );
}
