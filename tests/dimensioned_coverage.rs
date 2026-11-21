//! Regressions for source-generated physical contracts and portable model APIs.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::Exp;
use rssn::kernel::signature::Rational;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;
fn dim(name: &str) -> DimExpr {
    DimExpr::Base(name.into())
}
fn tensor(
    shape: &[usize],
    data: &[f64],
    d: DimExpr,
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        data: data.iter().copied().map(Value::Float).collect(),
        elem: if d == DimExpr::Dimless {
            TypeExpr::Float
        } else {
            TypeExpr::Dim(Box::new(TypeExpr::Float), d)
        },
    })
}
fn scalar(
    x: f64,
    d: DimExpr,
) -> Value {
    Value::Dyn(DynValue {
        value: Box::new(Value::Float(x)),
        dim: d,
    })
}
fn field<'a>(
    v: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(r) = v else {
        panic!("record expected: {v:?}")
    };
    &r.fields.iter().find(|(n, _)| n == name).unwrap().1
}
fn call(
    s: &Session,
    name: &str,
    args: &[Value],
) -> Value {
    s.call_value(name, args, &CallCtx::new())
        .unwrap_or_else(|e| panic!("{name}: {e:?}"))
}
#[path = "support/dimensions.rs"]
mod test_dimensions;
fn assert_dimension(
    value: &Value,
    expected: DimExpr,
) {
    let actual = match value {
        | Value::Dyn(v) => v.dim.clone(),
        | Value::Tensor(t) => {
            match &t.elem {
                | TypeExpr::Dim(_, d) => d.clone(),
                | _ => DimExpr::Dimless,
            }
        },
        | _ => DimExpr::Dimless,
    };
    assert_eq!(
        test_dimensions::factors(&actual),
        test_dimensions::factors(&expected)
    );
}
#[test]
fn measured_interpolation_fft_statistics_and_matrix_functions() {
    let s = Session::new();
    let x = tensor(&[3], &[0.0, 1.0, 2.0], dim("Time"));
    let y = tensor(&[3], &[0.0, 2.0, 4.0], dim("Length"));
    let out = call(
        &s,
        "num_interp_ext_pchip",
        &[x.clone(), y.clone(), scalar(0.5, dim("Time"))],
    );
    assert_dimension(&out, dim("Length"));
    assert!(
        s.call_value(
            "num_interp_ext_pchip",
            &[x, y.clone(), scalar(0.5, dim("Length"))],
            &CallCtx::new()
        )
        .is_err()
    );
    assert_dimension(
        &call(&s, "num_signal_dct", std::slice::from_ref(&y)),
        dim("Length"),
    );
    let none = Value::Variant(VariantValue {
        tag: "None".into(),
        payload: vec![],
    });
    let density = call(&s, "num_stats_ext_kde", &[y.clone(), y, none]);
    assert_dimension(&density, DimExpr::Inverse(Box::new(dim("Length"))));
    let a = tensor(
        &[2, 2],
        &[4.0, 0.0, 0.0, 9.0],
        dim("Length").pow(Exp::Rat(Rational::new(2, 1))),
    );
    assert_dimension(
        &call(&s, "num_matfun_sqrtm", std::slice::from_ref(&a)),
        dim("Length"),
    );
    let schur = call(&s, "num_schur_schur", &[a]);
    assert_dimension(field(&schur, "q"), DimExpr::Dimless);
    assert_dimension(
        field(&schur, "t"),
        dim("Length").pow(Exp::Rat(Rational::new(2, 1))),
    );
    let fit = call(
        &s,
        "num_stats_ext_ols",
        &[
            Value::int(3),
            Value::int(1),
            tensor(&[3], &[1.0, 2.0, 3.0], dim("Time")),
            tensor(&[3], &[3.0, 5.0, 7.0], dim("Length")),
            Value::Bool(true),
        ],
    );
    assert_dimension(
        field(&fit, "coefficients"),
        DimExpr::Mul(vec![dim("Length"), DimExpr::Inverse(Box::new(dim("Time")))]),
    );
    assert_dimension(field(&fit, "intercept"), dim("Length"));
}
#[test]
fn opaque_models_roundtrip_with_distinct_axis_and_value_units() {
    let s = Session::new();
    let points = tensor(&[3, 2], &[0.0, 0.0, 2.0, 0.0, 0.0, 2.0], dim("Length"));
    let model = call(&s, "num_kdtree_build", &[points]);
    let out = call(
        &s,
        "num_kdtree_nearest",
        &[
            model.clone(),
            tensor(&[2], &[1.8, 0.0], dim("Length")),
            Value::int(2),
        ],
    );
    assert_dimension(field(&out, "distances"), dim("Length"));
    assert!(
        s.call_value(
            "num_kdtree_nearest",
            &[model, tensor(&[2], &[1.0, 0.0], dim("Time")), Value::int(1)],
            &CallCtx::new()
        )
        .is_err()
    );
    let model = call(
        &s,
        "num_bspline_build",
        &[
            tensor(&[3], &[0.0, 1.0, 2.0], dim("Time")),
            tensor(&[3], &[0.0, 2.0, 4.0], dim("Length")),
            Value::int(2),
        ],
    );
    let out = call(
        &s,
        "num_bspline_query",
        &[model.clone(), scalar(0.5, dim("Time"))],
    );
    assert_dimension(&out, dim("Length"));
    let mut malformed = model;
    if let Value::Record(r) = &mut malformed
        && let Value::Opaque(p) = &mut r.fields.iter_mut().find(|(n, _)| n == "payload").unwrap().1
    {
        p.bytes = vec![255; 8];
    }
    assert!(matches!(
        s.call_value(
            "num_bspline_query",
            &[malformed, scalar(0.5, dim("Time"))],
            &CallCtx::new()
        ),
        Err(OpError::Domain(_))
    ));
}
#[test]
fn normalized_ad_and_polynomial_adapters_restore_units() {
    let s = Session::new();
    let args = vec![
        Value::String("x^3".into()),
        Value::Tensor(TensorValue {
            shape: vec![1],
            data: vec![Value::String("x".into())],
            elem: TypeExpr::String,
        }),
        tensor(&[1], &[4.0], dim("Length")),
        scalar(2.0, dim("Length")),
        scalar(3.0, dim("Time")),
        Value::String("reverse".into()),
    ];
    let gradient = call(&s, "num_ad_gradient_scaled", &args);
    assert_dimension(
        &gradient,
        DimExpr::Mul(vec![dim("Time"), DimExpr::Inverse(Box::new(dim("Length")))]),
    );
    let Value::Tensor(g) = gradient else { panic!() };
    let Value::Dyn(v) = &g.data[0] else { panic!() };
    assert_eq!(*v.value, Value::Float(18.0));
    let coefficients = call(
        &s,
        "num_polyfit_normalized",
        &[
            tensor(&[3], &[0.0, 2.0, 4.0], dim("Time")),
            tensor(&[3], &[1.0, 3.0, 5.0], dim("Length")),
            Value::int(1),
            scalar(2.0, dim("Time")),
        ],
    );
    let value = call(
        &s,
        "num_polyval_normalized",
        &[
            coefficients,
            scalar(3.0, dim("Time")),
            scalar(2.0, dim("Time")),
        ],
    );
    assert_dimension(&value, dim("Length"));
}

#[test]
fn shape_dependent_determinant_and_density_units_are_exact() {
    let session = Session::new();
    let a = tensor(&[2, 2], &[2.0, 0.0, 0.0, 3.0], dim("Length"));
    let determinant = call(&session, "num_dense_det_units", &[a]);
    assert_dimension(
        &determinant,
        dim("Length").pow(Exp::Rat(Rational::new(2, 1))),
    );
    let density = call(
        &session,
        "num_distributions_multivariate_normal_pdf",
        &[
            tensor(&[2], &[0.0, 0.0], dim("Length")),
            tensor(
                &[2, 2],
                &[1.0, 0.0, 0.0, 1.0],
                dim("Length").pow(Exp::Rat(Rational::new(2, 1))),
            ),
            tensor(&[2], &[0.0, 0.0], dim("Length")),
        ],
    );
    assert_dimension(&density, dim("Length").pow(Exp::Rat(Rational::new(-2, 1))));
}
