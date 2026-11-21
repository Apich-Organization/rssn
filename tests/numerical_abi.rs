use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::rules::numerical_abi::EXPORTS;
use rssn::rules::numerical_abi::GENERATED_EXPORTS;
use std::collections::BTreeSet;

fn vector(values: &[f64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}
fn matrix(
    rows: usize,
    cols: usize,
    values: &[f64],
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![rows, cols],
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float,
    })
}
fn call(
    session: &Session,
    name: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(name, args, &CallCtx::new())
}

#[test]
fn carlson_and_jacobi_cross_the_typed_value_abi() {
    let session = Session::new();
    let rg = call(
        &session,
        "num_special_ext_carlson_rg",
        &[Value::Float(4.0), Value::Float(4.0), Value::Float(4.0)],
    )
    .unwrap();
    assert_eq!(rg, Value::Float(2.0));

    let u = 0.75_f64;
    let amplitude = call(
        &session,
        "num_special_ext_jacobi_amplitude",
        &[Value::Float(u), Value::Float(0.0)],
    )
    .unwrap();
    assert_eq!(amplitude, Value::Float(u));

    let Value::Tensor(samples) = call(
        &session,
        "num_special_ext_jacobi_elliptic_batch",
        &[vector(&[0.0, u]), Value::Float(0.0)],
    )
    .unwrap() else {
        panic!("expected a matrix of Jacobi samples")
    };
    assert_eq!(samples.shape, [2, 3]);
    assert_eq!(samples.elem, TypeExpr::Float);
    for (got, expected) in samples
        .data
        .iter()
        .zip([0.0, 1.0, 1.0, u.sin(), u.cos(), 1.0])
    {
        assert!((got.as_f64().unwrap() - expected).abs() < 1e-15);
    }
}

#[test]
fn levin_condition_profile_preserves_term_dimensions_and_emits_unitless_orders() {
    let session = Session::new();
    let op = "num_analysis_convergence_levin_condition_numbers";
    let args = [
        measured_geometry(vec![4], &[1.0, -0.5, 1.0 / 3.0, -0.25], "Length"),
        Value::String("U".into()),
        Value::Float(1.0),
    ];
    let Value::Tensor(profile) = call(&session, op, &args).unwrap() else {
        panic!("expected the Levin condition number by order");
    };
    assert_eq!(profile.shape, [4]);
    assert_eq!(profile.elem, TypeExpr::Float);
    assert_eq!(profile.data[0], Value::Float(1.0));
    assert_eq!(
        profile.data.last(),
        Some(
            &call(
                &session,
                "num_analysis_convergence_levin_condition_number",
                &args
            )
            .unwrap()
        ),
    );
    let mut invalid = args;
    invalid[1] = Value::String("invalid".into());
    assert!(call(&session, op, &invalid).is_err());
}

#[test]
fn custom_histogram_checks_shared_dimensions_and_returns_exact_counts() {
    let session = Session::new();
    let op = "num_descriptive_histogram_edges";
    let samples = measured_geometry(vec![4], &[0.0, 2.0, 5.0, 5.0], "Length");
    let edges = measured_geometry(vec![3], &[-1.0, 2.0, 5.0], "Length");
    let Value::Tensor(counts) = call(&session, op, &[samples.clone(), edges]).unwrap() else {
        panic!("expected bin counts");
    };
    assert_eq!(counts.shape, [2]);
    assert_eq!(counts.elem, TypeExpr::Int);
    assert_eq!(counts.data, [Value::int(1), Value::int(3)]);
    for invalid in [
        measured_geometry(vec![3], &[-1.0, 2.0, 5.0], "Time"),
        measured_geometry(vec![3], &[-1.0, 2.0, 2.0], "Length"),
        measured_geometry(vec![3], &[-1.0, 2.0, 4.0], "Length"),
    ] {
        assert!(call(&session, op, &[samples.clone(), invalid]).is_err());
    }
}

#[test]
fn symmetric_matrix_powers_check_domain_and_preserve_zero_eigenvalues() {
    let session = Session::new();
    let op = "num_matfun_powm_psd";
    let positive_semidefinite = matrix(2, 2, &[0.0, 0.0, 0.0, 16.0]);
    let Value::Tensor(half_power) = call(
        &session,
        op,
        &[positive_semidefinite.clone(), Value::Float(0.5)],
    )
    .unwrap() else {
        panic!("expected a matrix power");
    };
    assert_eq!(half_power.shape, [2, 2]);
    assert_eq!(half_power.data, [0.0, 0.0, 0.0, 4.0].map(Value::Float));
    assert!(call(&session, op, &[positive_semidefinite, Value::Float(-1.0)]).is_err());
    let measured = measured_geometry(vec![2, 2], &[4.0, 0.0, 0.0, 16.0], "Length");
    assert!(call(&session, op, &[measured, Value::Float(0.5)]).is_err());
}

#[test]
fn shape_preserving_interpolators_match_query_and_value_dimensions() {
    use rssn::kernel::signature::DimExpr;

    let s = Session::new();
    let points = measured_geometry(vec![3], &[0.0, 1.0, 2.0], "Time");
    let values = measured_geometry(vec![3], &[1.0, 3.0, 5.0], "Length");
    let queries = measured_geometry(vec![3], &[1.5, 0.25, 2.0], "Time");
    let wrong = measured_geometry(vec![3], &[1.5, 0.25, 2.0], "Length");
    let expected = TypeExpr::Float.with_dim(DimExpr::base("Length"));
    for name in ["num_interp_ext_akima_many", "num_interp_ext_pchip_many"] {
        let Value::Tensor(out) =
            call(&s, name, &[points.clone(), values.clone(), queries.clone()]).unwrap()
        else {
            panic!("expected interpolated samples");
        };
        assert_eq!(out.shape, [3]);
        assert!(rssn::kernel::overload::matches_type(&expected, &out.elem).accepted());
        for (actual, expected) in out.data.iter().map(geometry_scalar).zip([4.0, 1.5, 5.0]) {
            assert!((actual - expected).abs() < 1e-14);
        }
        assert!(call(&s, name, &[points.clone(), values.clone(), wrong.clone()]).is_err());
    }
}

#[test]
fn consecutive_mobius_values_and_congruence_classes_are_exact() {
    let s = Session::new();
    let result = call(&s, "num_number_ext_mobius_sieve", &[Value::int(10)]).unwrap();
    let Value::Tensor(result) = result else {
        panic!("expected Möbius values")
    };
    assert_eq!(result.shape, [11]);
    assert_eq!(result.elem, TypeExpr::Int);
    assert_eq!(
        result.data,
        [0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1].map(Value::int)
    );

    let op = "num_number_ext_solve_linear_congruence";
    let Value::Variant(solution) =
        call(&s, op, &[Value::int(14), Value::int(30), Value::int(100)]).unwrap()
    else {
        panic!("expected an optional congruence class")
    };
    assert_eq!(solution.tag, "Some");
    assert_eq!(
        geometry_field(&solution.payload[0], "item0"),
        &Value::int(45)
    );
    assert_eq!(
        geometry_field(&solution.payload[0], "item1"),
        &Value::int(50)
    );
    let Value::Variant(no_solution) =
        call(&s, op, &[Value::int(6), Value::int(5), Value::int(15)]).unwrap()
    else {
        panic!("expected an optional congruence class")
    };
    assert_eq!(no_solution.tag, "None");
    assert!(no_solution.payload.is_empty());
    assert!(call(&s, op, &[Value::int(1), Value::int(0), Value::int(0)]).is_err());
}
#[test]
fn descriptors_register_exact_manifest_overloads_and_unique_keys() {
    let session = Session::new();
    let manifest = session.manifest();
    let mut names = BTreeSet::new();
    assert!(!GENERATED_EXPORTS.is_empty());
    for export in EXPORTS {
        assert!(names.insert(export.name));
        let info = manifest.op(export.name).expect(export.name);
        assert_eq!(info.rule_set, export.namespace);
        assert_eq!(info.rust_path.as_deref(), Some(export.rust_path));
        // The manifest keeps the declared text first and may complete it
        // from the rustdoc of `rust_path` (`doc_source` says which).
        assert!(
            info.doc.starts_with(export.doc.trim()),
            "{}: {:?} does not start with the declared {:?}",
            export.name,
            info.doc,
            export.doc
        );
        assert_eq!(info.stability, export.stability);
        assert_eq!(info.overloads.len(), 1);
        assert_eq!(info.keys.len(), 1);
        assert_eq!(
            info.signature_source,
            rssn::kernel::manifest::SignatureSource::Declared
        );
    }
}

#[test]
fn analytic_signal_and_zoom_fft_cross_the_typed_value_abi() {
    use rssn::kernel::value::ComplexValue;

    let session = Session::new();
    let Value::Tensor(analytic) = call(
        &session,
        "num_signal_analytic_signal",
        &[vector(&[1.0, 0.0, -1.0, 0.0])],
    )
    .unwrap() else {
        panic!("analytic signal must return a complex tensor")
    };
    assert_eq!(analytic.shape, [4]);
    assert_eq!(analytic.elem, TypeExpr::Complex);
    let samples: Vec<_> = analytic
        .data
        .into_iter()
        .map(|value| {
            let Value::Complex(value) = value else {
                panic!("expected a complex sample")
            };
            value.to_f64()
        })
        .collect();
    for ((re, im), (expected_re, expected_im)) in
        samples
            .into_iter()
            .zip([(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)])
    {
        assert!((re - expected_re).abs() <= 1e-12);
        assert!((im - expected_im).abs() <= 1e-12);
    }

    let impulse = Value::Tensor(TensorValue {
        shape: vec![2],
        data: vec![
            Value::Complex(ComplexValue::Float { re: 1.0, im: -0.25 }),
            Value::Complex(ComplexValue::Float { re: 0.0, im: 0.0 }),
        ],
        elem: TypeExpr::Complex,
    });
    let Value::Tensor(spectrum) = call(
        &session,
        "num_signal_zoom_fft",
        &[
            impulse,
            Value::Float(0.07),
            Value::Float(0.21),
            Value::int(5),
        ],
    )
    .unwrap() else {
        panic!("zoom FFT must return a complex tensor")
    };
    assert_eq!(spectrum.shape, [5]);
    assert_eq!(spectrum.elem, TypeExpr::Complex);
    for sample in spectrum.data {
        let Value::Complex(sample) = sample else {
            panic!("expected a complex frequency bin")
        };
        let (re, im) = sample.to_f64();
        assert!((re - 1.0).abs() <= 1e-12);
        assert!((im + 0.25).abs() <= 1e-12);
    }
}

#[test]
fn typed_signal_resampling_and_inverse_stft_reconstruct_samples() {
    let session = Session::new();
    let samples = vector(&[1.0, 2.0, 3.0, 4.0]);
    let window = vector(&[1.0, 1.0]);
    let hop = Value::int(2);
    let spectrum = call(
        &session,
        "num_analysis_signal_stft",
        &[samples.clone(), window.clone(), hop.clone()],
    )
    .unwrap();
    assert_eq!(
        call(
            &session,
            "num_signal_inverse_stft",
            &[spectrum, window, hop]
        )
        .unwrap(),
        samples
    );
    assert_eq!(
        call(
            &session,
            "num_signal_upfirdn",
            &[
                vector(&[1.0, 2.0]),
                vector(&[1.0, 1.0]),
                Value::int(2),
                Value::int(1)
            ]
        )
        .unwrap(),
        vector(&[1.0, 1.0, 2.0, 2.0])
    );
    assert_eq!(
        call(&session, "num_special_bessel_i0e", &[Value::Float(0.0)]).unwrap(),
        Value::Float(1.0)
    );
}

#[test]
fn checked_fft_convolution_preserves_both_operand_dimensions() {
    use rssn::kernel::signature::DimExpr;

    let session = Session::new();
    let args = [
        measured_geometry(vec![2], &[2.0, 3.0], "Length"),
        measured_geometry(vec![2], &[1.0, -1.0], "Time"),
    ];
    let op = "num_fftconv_convolve_overlap_save";
    let Value::Tensor(result) = call(&session, op, &args).unwrap() else {
        panic!("expected a vector of quantities");
    };
    assert_eq!(result.shape, [3]);
    let product = DimExpr::Mul(vec![DimExpr::base("Length"), DimExpr::base("Time")]);
    assert!(
        rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(product), &result.elem)
            .accepted()
    );
    for (sample, expected) in result.data.iter().zip([2.0, 1.0, -3.0]) {
        assert!((geometry_scalar(sample) - expected).abs() < 1e-12);
    }
}

#[test]
fn welch_record_fields_preserve_their_independent_dimensions() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::QuantityValue;

    let session = Session::new();
    let hertz = session
        .manifest()
        .units
        .iter()
        .find(|unit| unit.name == "Hz")
        .unwrap()
        .clone();
    let frequency = hertz.dimension.clone();
    let sample_rate = Value::Quantity(QuantityValue {
        value: Box::new(Value::Float(4.0)),
        unit: hertz,
    });
    let spectrum = call(
        &session,
        "num_signal_welch_psd",
        &[
            measured_geometry(vec![4], &[1.0; 4], "Length"),
            vector(&[1.0; 4]),
            Value::int(4),
            sample_rate,
        ],
    )
    .unwrap();
    let Value::Record(record) = spectrum else {
        panic!("expected a Welch power record")
    };
    let Value::Tensor(frequencies) = record.get("frequencies").unwrap() else {
        panic!("expected a frequency vector")
    };
    let Value::Tensor(density) = record.get("density").unwrap() else {
        panic!("expected a density vector")
    };
    let density_dimension =
        DimExpr::Mul(vec![DimExpr::base("Length").pow(2), DimExpr::base("Time")]);
    for (tensor, expected) in [(frequencies, frequency), (density, density_dimension)] {
        assert!(
            rssn::kernel::overload::matches_type(&TypeExpr::Float.with_dim(expected), &tensor.elem)
                .accepted()
        );
    }
    assert!((geometry_scalar(&frequencies.data[1]) - 1.0).abs() < 1e-12);
    assert!((geometry_scalar(&density.data[0]) - 1.0).abs() < 1e-12);
}

#[test]
fn time_series_lag_operators_keep_covariance_units_and_normalize_correlations() {
    use rssn::kernel::signature::DimExpr;

    let session = Session::new();
    let samples = measured_geometry(vec![2], &[1.0, -1.0], "Length");
    let args = [samples.clone(), Value::int(1)];
    let Value::Tensor(covariance) = call(&session, "num_signal_autocovariance", &args).unwrap()
    else {
        panic!("expected a vector of covariances");
    };
    assert_eq!(covariance.shape, [2]);
    let squared_length = TypeExpr::Float.with_dim(DimExpr::base("Length").pow(2));
    assert!(rssn::kernel::overload::matches_type(&squared_length, &covariance.elem).accepted());
    for (value, expected) in covariance.data.iter().zip([1.0, -0.5]) {
        assert!((geometry_scalar(value) - expected).abs() < 1e-12);
    }

    for name in [
        "num_signal_autocorrelation",
        "num_signal_partial_autocorrelation",
    ] {
        assert_eq!(call(&session, name, &args).unwrap(), vector(&[1.0, -0.5]));
        assert!(matches!(
            call(&session, name, &[samples.clone(), Value::int(2)]),
            Err(OpError::Domain(_))
        ));
    }
}

#[test]
fn numeric_families_dispatch_through_name_and_key() {
    let session = Session::new();
    let manifest = session.manifest();
    let args = [vector(&[1.0, 2.0]), vector(&[3.0, 4.0])];
    let key = &manifest.op("num_dot").unwrap().keys[0];
    assert_eq!(
        call(&session, "num_dot", &args).unwrap(),
        Value::Float(11.0)
    );
    assert_eq!(
        session.call_key(key, &args, &CallCtx::new()).unwrap(),
        Value::Float(11.0)
    );
    assert_eq!(
        call(&session, "num_mean", &[vector(&[1.0, 2.0, 3.0])]).unwrap(),
        Value::Float(2.0)
    );
    assert_eq!(
        call(&session, "num_sigmoid", &[Value::Float(0.0)]).unwrap(),
        Value::Float(0.5)
    );
    assert_eq!(
        call(
            &session,
            "num_matmul",
            &[matrix(1, 2, &[2.0, 3.0]), matrix(2, 1, &[4.0, 5.0])]
        )
        .unwrap(),
        matrix(1, 1, &[23.0])
    );
}
#[test]
fn numeric_adapters_reject_bad_shapes_and_non_finite_values() {
    let session = Session::new();
    assert!(matches!(
        call(&session, "num_dot", &[vector(&[1.0]), vector(&[1.0, 2.0])]),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(&session, "num_norm", &[vector(&[f64::NAN])]),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(
            &session,
            "num_matmul",
            &[matrix(1, 2, &[1.0, 2.0]), matrix(1, 2, &[1.0, 2.0])]
        ),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(
            &session,
            "num_cholesky",
            &[matrix(2, 2, &[1.0, 2.0, 2.0, 1.0])]
        ),
        Err(OpError::Domain(_))
    ));
    let malformed = matrix(2, 2, &[1.0]);
    assert!(matches!(
        call(&session, "num_transpose", &[malformed]),
        Err(OpError::Domain(_))
    ));
    let transposed = call(&session, "num_transpose", &[matrix(2, 1, &[4.0, 5.0])]).unwrap();
    assert_eq!(
        call(
            &session,
            "num_matmul",
            &[transposed, matrix(2, 1, &[2.0, 3.0])]
        )
        .unwrap(),
        matrix(1, 1, &[23.0])
    );
}

#[test]
fn structured_codecs_preserve_schema_and_decompositions() {
    let session = Session::new();
    let manifest = session.manifest();
    for (name, fields) in [
        ("num_svd", vec!["u", "s", "v"]),
        ("num_qr", vec!["q", "r"]),
        ("num_eigen_symmetric", vec!["values", "vectors"]),
    ] {
        let info = manifest.op(name).expect(name);
        let TypeExpr::Record(schema) = &info.overloads[0].ret else {
            panic!("{name} has no record signature");
        };
        assert_eq!(
            schema
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            fields
        );
        let result = call(&session, name, &[matrix(2, 2, &[3.0, 0.0, 0.0, 2.0])]).unwrap();
        let Value::Record(record) = result else {
            panic!("{name} has no record result");
        };
        assert_eq!(
            record
                .fields
                .iter()
                .map(|(name, _)| name.as_str())
                .collect::<Vec<_>>(),
            fields
        );
        for ((_, field_ty), (_, field)) in schema.iter().zip(&record.fields) {
            assert!(rssn::kernel::overload::matches_type(field_ty, &field.type_expr()).accepted());
        }
    }
    let Value::Record(svd) =
        call(&session, "num_svd", &[matrix(2, 2, &[3.0, 0.0, 0.0, 2.0])]).unwrap()
    else {
        panic!("svd record");
    };
    let (Value::Tensor(u), Value::Tensor(s), Value::Tensor(v)) = (
        svd.get("u").unwrap(),
        svd.get("s").unwrap(),
        svd.get("v").unwrap(),
    ) else {
        panic!("svd tensors");
    };
    for i in 0..2 {
        for j in 0..2 {
            let mut actual = 0.0;
            for k in 0..2 {
                actual += u.data[i * 2 + k].as_f64().unwrap()
                    * s.data[k].as_f64().unwrap()
                    * v.data[j * 2 + k].as_f64().unwrap();
            }
            let expected = [3.0, 0.0, 0.0, 2.0][i * 2 + j];
            assert!((actual - expected).abs() < 1e-12);
        }
    }
}
#[test]
fn integer_and_multiple_argument_codecs_check_domains() {
    let session = Session::new();
    let Value::Float(fact) = call(&session, "num_factorial", &[Value::int(5)]).unwrap() else {
        panic!("factorial");
    };
    assert!((fact - 120.0).abs() < 1e-9);
    let Value::Float(comb) =
        call(&session, "num_binomial", &[Value::int(6), Value::int(2)]).unwrap()
    else {
        panic!("binomial");
    };
    assert!((comb - 15.0).abs() < 1e-9);
    assert!(matches!(
        call(&session, "num_factorial", &[Value::int(-1)]),
        Err(OpError::Domain(_))
    ));
    assert!(call(&session, "num_factorial", &[Value::Float(2.5)]).is_err());
    assert!(matches!(
        call(
            &session,
            "num_rank_tol",
            &[matrix(1, 1, &[1.0]), Value::Float(-1.0)]
        ),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(
            &session,
            "num_lstsq",
            &[matrix(2, 1, &[1.0, 2.0]), vector(&[1.0])]
        ),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn documented_infinity_is_not_confused_with_invalid_nan() {
    let session = Session::new();
    let condition = call(&session, "num_cond", &[matrix(2, 2, &[1.0, 0.0, 0.0, 0.0])]).unwrap();
    assert!(matches!(condition, Value::Float(value) if value == f64::INFINITY));
    assert!(matches!(
        call(&session, "num_geometric_mean", &[vector(&[-1.0, 2.0])]),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn tensor_train_export_keeps_shapes_cores_and_error_report() {
    let session = Session::new();
    let shape = Value::Tensor(TensorValue {
        shape: vec![2],
        data: vec![Value::int(2), Value::int(2)],
        elem: TypeExpr::Int,
    });
    let args = [
        vector(&[2.0, 0.0, 0.0, 2.0]),
        shape,
        Value::Float(0.0),
        Value::int(2),
    ];
    let Value::Record(result) = call(&session, "num_tt_svd", &args).unwrap() else {
        panic!("tensor train record");
    };
    assert_eq!(result.get("shape"), Some(&args[1]));
    assert_eq!(result.get("tolerance_met"), Some(&Value::Bool(true)));
    assert!(matches!(result.get("error"),Some(Value::Float(x)) if *x < 1e-12));
    let Some(Value::Tensor(cores)) = result.get("cores") else {
        panic!("core list");
    };
    assert_eq!(cores.shape, vec![2]);
    assert_eq!(result.get("scale"), Some(&Value::Float(2.0)));
    let (Value::Tensor(left), Value::Tensor(right)) = (&cores.data[0], &cores.data[1]) else {
        panic!("core tensors");
    };
    let rank = left.shape[2];
    for i in 0..2 {
        for j in 0..2 {
            let value = (0..rank)
                .map(|r| {
                    left.data[i * rank + r].as_f64().unwrap()
                        * right.data[r * 2 + j].as_f64().unwrap()
                })
                .sum::<f64>()
                * 2.0;
            assert!((value - if i == j { 2.0 } else { 0.0 }).abs() < 1e-12);
        }
    }
    assert!(
        cores
            .data
            .iter()
            .all(|v| matches!(v,Value::Tensor(t) if t.shape.len() == 3))
    );
}

#[test]
fn dimensional_policies_preserve_and_combine_witnesses() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    let session = Session::new();
    let tensor = |unit: &str, xs: &[f64]| {
        Value::Tensor(TensorValue {
            shape: vec![xs.len()],
            data: xs.iter().copied().map(Value::Float).collect(),
            elem: TypeExpr::Dim(Box::new(TypeExpr::Float), DimExpr::Base(unit.into())),
        })
    };
    let a = tensor("Length", &[1.0, 2.0]);
    let b = tensor("Time", &[3.0, 4.0]);
    let Value::Dyn(DynValue { value, dim }) =
        call(&session, "num_dot", &[a.clone(), b.clone()]).unwrap()
    else {
        panic!("dot witness");
    };
    assert_eq!(*value, Value::Float(11.0));
    assert_eq!(
        dim,
        DimExpr::Mul(vec![
            DimExpr::Base("Length".into()),
            DimExpr::Base("Time".into())
        ])
    );
    let Value::Dyn(DynValue { value, dim }) =
        call(&session, "num_mean", std::slice::from_ref(&a)).unwrap()
    else {
        panic!("mean witness");
    };
    assert_eq!(*value, Value::Float(1.5));
    assert_eq!(dim, DimExpr::Base("Length".into()));
    assert!(matches!(
        call(&session, "num_vector_add", &[a, b]),
        Err(OpError::Domain(_))
    ));
}

#[test]
fn quantities_are_converted_to_coherent_units_before_numerical_work() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    use rssn::kernel::value::QuantityValue;
    let session = Session::new();
    let mut unit = session
        .manifest()
        .units
        .iter()
        .find(|u| u.name == "m")
        .expect("metre")
        .clone();
    unit.scale.numer = 1000.into();
    let argument = Value::Tensor(TensorValue {
        shape: vec![1],
        elem: TypeExpr::Dim(Box::new(TypeExpr::Float), DimExpr::Base("Length".into())),
        data: vec![Value::Quantity(QuantityValue {
            value: Box::new(Value::Float(2.0)),
            unit,
        })],
    });
    let Value::Dyn(DynValue { value, dim }) = call(&session, "num_mean", &[argument]).unwrap()
    else {
        panic!("quantity result witness");
    };
    assert_eq!(*value, Value::Float(2000.0));
    assert_eq!(dim, DimExpr::Base("Length".into()));
}

#[test]
fn generated_source_bindings_dispatch_and_obey_output_policies() {
    let session = Session::new();
    let product = call(
        &session,
        "num_cross_product",
        &[vector(&[1.0, 0.0, 0.0]), vector(&[0.0, 1.0, 0.0])],
    )
    .expect("3D cross product");
    let Value::Tensor(product) = product else {
        panic!("expected tensor")
    };
    assert_eq!(
        product.data,
        vec![Value::Float(0.0), Value::Float(0.0), Value::Float(1.0)]
    );
    assert!(
        call(
            &session,
            "num_cross_product",
            &[vector(&[1.0, 0.0]), vector(&[0.0, 1.0])]
        )
        .is_err()
    );
    assert_eq!(
        call(&session, "num_inverse_erfc", &[Value::Float(0.0)]).unwrap(),
        Value::Float(f64::INFINITY)
    );
    assert!(call(&session, "num_inverse_erfc", &[Value::Float(-1.0)]).is_err());
    let Value::Float(result) = call(&session, "num_erfi", &[Value::Float(1.0)]).unwrap() else {
        panic!("expected scalar")
    };
    assert!((result - 1.6504257588).abs() < 1e-8);
    let sample = vector(&[9.0, 1.0, 3.0]);
    assert_eq!(
        call(&session, "num_median", std::slice::from_ref(&sample)).unwrap(),
        Value::Float(3.0)
    );
    assert_eq!(
        sample,
        vector(&[9.0, 1.0, 3.0]),
        "mutable Rust kernels must not mutate caller tensors"
    );
    assert!(call(&session, "num_skewness", &[vector(&[1.0, 2.0])]).is_err());
    assert_eq!(
        call(
            &session,
            "num_lp_norm",
            &[vector(&[3.0, 4.0]), Value::Float(2.0)]
        )
        .unwrap(),
        Value::Float(5.0)
    );
    assert!(
        call(
            &session,
            "num_lp_norm",
            &[vector(&[3.0]), Value::Float(0.0)]
        )
        .is_err()
    );
    let Value::Tensor(scaled) = call(
        &session,
        "num_scalar_mul",
        &[vector(&[2.0, 3.0]), Value::Float(4.0)],
    )
    .unwrap() else {
        panic!("expected tensor")
    };
    assert_eq!(scaled.data, vec![Value::Float(8.0), Value::Float(12.0)]);
}

#[test]
fn generated_structured_policies_are_exhaustive_and_honest_about_unbound_types() {
    use rssn::rules::numerical_abi::codecs::STRUCTURED_POLICIES;
    use rssn::rules::numerical_abi::dimensions::FieldDimension;
    let svd = STRUCTURED_POLICIES
        .iter()
        .find(|policy| policy.rust_path == "crate::kernels::dense::Svd")
        .unwrap();
    assert!(svd.bound);
    assert_eq!(
        svd.fields.iter().map(|(name, _)| *name).collect::<Vec<_>>(),
        ["u", "s", "v"]
    );
    assert!(
        svd.fields
            .iter()
            .all(|(_, dimension)| matches!(dimension, FieldDimension::Formula(_)))
    );
    let integral = STRUCTURED_POLICIES
        .iter()
        .find(|policy| policy.rust_path == "crate::kernels::quadrature::Integral")
        .unwrap();
    assert!(integral.bound);
    assert!(matches!(integral.fields[0].1, FieldDimension::Formula(_)));
    assert!(matches!(integral.fields[1].1, FieldDimension::Formula(_)));
    assert!(matches!(integral.fields[2].1, FieldDimension::NonNumeric));
    let iterative = STRUCTURED_POLICIES
        .iter()
        .find(|policy| {
            policy.rust_path == "crate::kernels::krylov::IterResult"
                && policy.contract == "ITERATIVE_LINEAR_SOLVE"
        })
        .unwrap();
    assert!(iterative.bound);
    assert!(matches!(iterative.fields[0].1, FieldDimension::Formula(_)));
    assert!(matches!(iterative.fields[1].1, FieldDimension::Formula(_)));
    assert!(matches!(iterative.fields[2].1, FieldDimension::NonNumeric));
}

fn csr(
    indptr: &[i64],
    indices: &[i64],
    values: &[f64],
) -> Value {
    use rssn::kernel::value::RecordValue;
    let ints = |items: &[i64]| {
        Value::Tensor(TensorValue {
            shape: vec![items.len()],
            data: items.iter().copied().map(Value::int).collect(),
            elem: TypeExpr::Int,
        })
    };
    Value::Record(RecordValue {
        fields: vec![
            ("rows".into(), Value::int(2)),
            ("cols".into(), Value::int(2)),
            ("indptr".into(), ints(indptr)),
            ("indices".into(), ints(indices)),
            ("values".into(), vector(values)),
        ],
    })
}
fn result_field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected record")
    };
    &record
        .fields
        .iter()
        .find(|(field, _)| field == name)
        .unwrap()
        .1
}
#[test]
fn structural_csr_solvers_check_residuals_and_manifest_origins() {
    let session = Session::new();
    for name in ["num_csr_cg", "num_csr_gmres", "num_csr_bicgstab"] {
        let mut args = vec![
            csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]),
            vector(&[2.0, 8.0]),
            Value::Float(1e-10),
            Value::int(40),
        ];
        if name == "num_csr_gmres" {
            args.push(Value::int(2));
        }
        let result = call(&session, name, &args).unwrap();
        assert_eq!(result_field(&result, "converged"), &Value::Bool(true));
        assert!(result_field(&result, "residual").as_f64().unwrap() < 1e-10);
        let key = session.manifest().op(name).unwrap().keys[0].clone();
        assert_eq!(
            session.call_key(&key, &args, &CallCtx::new()).unwrap(),
            result
        );
    }
}
#[test]
fn structural_csr_rejects_bad_rows_columns_and_controls() {
    let session = Session::new();
    let good_rhs = vector(&[2.0, 8.0]);
    for matrix in [
        csr(&[0, 2, 1], &[0, 1], &[2.0, 4.0]),
        csr(&[0, 1, 2], &[2, 1], &[2.0, 4.0]),
        csr(&[0, 2, 2], &[0, 0], &[2.0, 4.0]),
        csr(&[0, 1, 3], &[0, 1], &[2.0, 4.0]),
    ] {
        assert!(
            call(
                &session,
                "num_csr_cg",
                &[matrix, good_rhs.clone(), Value::Float(1e-8), Value::int(10)]
            )
            .is_err()
        );
    }
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[
                csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]),
                vector(&[1.0]),
                Value::Float(1e-8),
                Value::int(10)
            ]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "num_csr_gmres",
            &[
                csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]),
                good_rhs.clone(),
                Value::Float(1e-8),
                Value::int(10),
                Value::int(0)
            ]
        )
        .is_err()
    );
    let mut oversized = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let Value::Record(ref mut record) = oversized else {
        unreachable!()
    };
    record.fields[0].1 = Value::int(-1);
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[
                oversized,
                good_rhs.clone(),
                Value::Float(1e-8),
                Value::int(10)
            ]
        )
        .is_err()
    );
    let mut duplicate = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let Value::Record(ref mut record) = duplicate else {
        unreachable!()
    };
    record.fields[1].0 = "rows".into();
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[duplicate, good_rhs, Value::Float(1e-8), Value::int(10)]
        )
        .is_err()
    );
}

#[test]
fn structural_csr_zero_rhs_and_iteration_limit_are_honest() {
    let session = Session::new();
    let matrix = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let zero = call(
        &session,
        "num_csr_cg",
        &[
            matrix.clone(),
            vector(&[0.0, 0.0]),
            Value::Float(1e-10),
            Value::int(5),
        ],
    )
    .unwrap();
    assert_eq!(result_field(&zero, "converged"), &Value::Bool(true));
    assert_eq!(result_field(&zero, "residual").as_f64(), Some(0.0));
    let limited = call(
        &session,
        "num_csr_cg",
        &[
            matrix,
            vector(&[2.0, 8.0]),
            Value::Float(1e-14),
            Value::int(1),
        ],
    )
    .unwrap();
    assert_eq!(result_field(&limited, "converged"), &Value::Bool(false));
    assert!(result_field(&limited, "residual").as_f64().unwrap() > 1e-14);
}
#[test]
fn structural_csr_preserves_units_and_rejects_conflicting_element_witnesses() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    let length = DimExpr::Base("Length".into());
    let time = DimExpr::Base("Time".into());
    let mut matrix = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let Value::Record(ref mut fields) = matrix else {
        unreachable!()
    };
    let Value::Tensor(values) = &mut fields.fields[4].1 else {
        unreachable!()
    };
    values.elem = TypeExpr::Dim(Box::new(TypeExpr::Float), length.clone());
    let mut rhs = vector(&[2.0, 8.0]);
    let Value::Tensor(rhs_values) = &mut rhs else {
        unreachable!()
    };
    rhs_values.elem = TypeExpr::Dim(Box::new(TypeExpr::Float), time.clone());
    let session = Session::new();
    let result = call(
        &session,
        "num_csr_cg",
        &[
            matrix.clone(),
            rhs.clone(),
            Value::Float(1e-10),
            Value::int(10),
        ],
    )
    .unwrap();
    let Value::Tensor(x) = result_field(&result, "x") else {
        panic!("solution tensor")
    };
    let TypeExpr::Dim(_, dimension) = &x.elem else {
        panic!("solution unit")
    };
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(dimension.clone()),
            &TypeExpr::Float.with_dim(DimExpr::Mul(vec![
                time.clone(),
                DimExpr::Inverse(Box::new(length.clone()))
            ]))
        )
        .accepted()
    );
    let Value::Record(ref mut fields) = matrix else {
        unreachable!()
    };
    let Value::Tensor(values) = &mut fields.fields[4].1 else {
        unreachable!()
    };
    values.data[0] = Value::Dyn(DynValue {
        value: Box::new(Value::Float(2.0)),
        dim: time,
    });
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[
                matrix.clone(),
                rhs.clone(),
                Value::Float(1e-10),
                Value::int(10)
            ]
        )
        .is_err()
    );
    let matrix = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let bad_tol = Value::Dyn(DynValue {
        value: Box::new(Value::Float(1e-10)),
        dim: length,
    });
    assert!(
        call(
            &session,
            "num_csr_cg",
            &[matrix, rhs, bad_tol, Value::int(10)]
        )
        .is_err()
    );
}

#[test]
fn structural_csr_bindings_solve_uniformly_tiny_and_huge_systems() {
    let session = Session::new();
    for scale in [f64::from_bits(1), 1e-320, 1e-300, 1e300] {
        for name in ["num_csr_cg", "num_csr_gmres", "num_csr_bicgstab"] {
            let mut args = vec![
                csr(&[0, 1, 2], &[0, 1], &[scale, scale]),
                vector(&[scale, -scale]),
                Value::Float(1e-10),
                Value::int(12),
            ];
            if name == "num_csr_gmres" {
                args.push(Value::int(2));
            }
            let result = call(&session, name, &args).unwrap();
            assert_eq!(
                result_field(&result, "converged"),
                &Value::Bool(true),
                "{name}, scale {scale:e}"
            );
            let Value::Tensor(x) = result_field(&result, "x") else {
                panic!("solution tensor")
            };
            assert!(
                (x.data[0].as_f64().unwrap() - 1.0).abs() < 1e-10,
                "{name}, scale {scale:e}"
            );
            assert!(
                (x.data[1].as_f64().unwrap() + 1.0).abs() < 1e-10,
                "{name}, scale {scale:e}"
            );
        }
    }
}

#[test]
fn generated_contract_table_preserves_dimension_rules_across_families() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    let session = Session::new();
    let dimensioned = |name: &str, xs: &[f64]| {
        Value::Tensor(TensorValue {
            shape: vec![xs.len()],
            data: xs.iter().copied().map(Value::Float).collect(),
            elem: TypeExpr::Dim(Box::new(TypeExpr::Float), DimExpr::Base(name.into())),
        })
    };
    let length = dimensioned("Length", &[1.0, 3.0]);
    let time = dimensioned("Time", &[1.0, 0.0]);
    let Value::Dyn(variance) =
        call(&session, "num_variance", std::slice::from_ref(&length)).unwrap()
    else {
        panic!("variance quantity")
    };
    let expected = DimExpr::Mul(vec![
        DimExpr::Base("Length".into()),
        DimExpr::Base("Length".into()),
    ]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(variance.dim),
            &TypeExpr::Float.with_dim(expected)
        )
        .accepted()
    );
    assert!(matches!(
        call(
            &session,
            "num_cosine_similarity",
            &[length.clone(), time.clone()]
        )
        .unwrap(),
        Value::Float(_)
    ));
    let Value::Tensor(projected) = call(&session, "num_project", &[length.clone(), time]).unwrap()
    else {
        panic!("projection")
    };
    assert!(
        rssn::kernel::overload::matches_type(
            &projected.elem,
            &TypeExpr::Float.with_dim(DimExpr::Base("Length".into()))
        )
        .accepted()
    );
    let empty = dimensioned("Length", &[]);
    let Value::Dyn(mean) = call(&session, "num_mean", &[empty]).unwrap() else {
        panic!("empty mean witness")
    };
    assert_eq!(mean.dim, DimExpr::Base("Length".into()));
    let bad_scale = Value::Dyn(DynValue {
        value: Box::new(Value::Float(2.0)),
        dim: DimExpr::Base("Time".into()),
    });
    assert!(call(&session, "num_lp_norm", &[length, bad_scale]).is_err());
}

#[test]
fn source_codecs_reject_dimensions_hidden_inside_optional_arguments() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    use rssn::kernel::value::VariantValue;
    let session = Session::new();
    let optional = |value| {
        Value::Variant(VariantValue {
            tag: "Some".into(),
            payload: vec![value],
        })
    };
    let args = [
        vector(&[0.0, 1.0, 2.0]),
        vector(&[1.0]),
        optional(Value::Float(0.5)),
    ];
    let valid = call(&session, "num_stats_ext_kde", &args);
    assert!(
        valid.is_ok(),
        "valid KDE optional bandwidth rejected: {valid:?}"
    );
    let none = Value::Variant(VariantValue {
        tag: "None".into(),
        payload: vec![],
    });
    let inferred = call(
        &session,
        "num_stats_ext_kde",
        &[args[0].clone(), args[1].clone(), none],
    );
    assert!(
        inferred.is_ok(),
        "valid KDE None bandwidth rejected: {inferred:?}"
    );
    for malformed in [
        Value::Variant(VariantValue {
            tag: "None".into(),
            payload: vec![Value::Float(0.5)],
        }),
        Value::Variant(VariantValue {
            tag: "Some".into(),
            payload: vec![],
        }),
        Value::Variant(VariantValue {
            tag: "Unknown".into(),
            payload: vec![],
        }),
        optional(Value::Bool(true)),
    ] {
        assert!(
            call(
                &session,
                "num_stats_ext_kde",
                &[args[0].clone(), args[1].clone(), malformed]
            )
            .is_err()
        );
    }
    let test_result = call(
        &session,
        "num_stats_ext_t_test_one_sample",
        &[vector(&[1.0, 2.0, 3.0]), Value::Float(2.0)],
    )
    .unwrap();
    let Value::Record(test_result) = test_result else {
        panic!("expected test-result record")
    };
    assert_eq!(test_result.get("df"), Some(&optional(Value::Float(2.0))));
    let ks_result = call(
        &session,
        "num_stats_ext_ks_test_two_sample",
        &[vector(&[1.0, 2.0]), vector(&[1.0, 2.0])],
    )
    .unwrap();
    let Value::Record(ks_result) = ks_result else {
        panic!("expected test-result record")
    };
    assert_eq!(
        ks_result.get("df"),
        Some(&Value::Variant(VariantValue {
            tag: "None".into(),
            payload: vec![]
        }))
    );
    let dimensioned = Value::Dyn(DynValue {
        value: Box::new(Value::Float(0.5)),
        dim: DimExpr::Base("Length".into()),
    });
    let bad_args = [args[0].clone(), args[1].clone(), optional(dimensioned)];
    let rejected = call(&session, "num_stats_ext_kde", &bad_args);
    assert!(
        rejected.is_err(),
        "dimensioned optional bandwidth was accepted: {rejected:?}"
    );
}

#[test]
fn sparse_eigen_and_amg_contracts_preserve_physical_dimensions_and_verified_residuals() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    let length = DimExpr::Base("Length".into());
    let time = DimExpr::Base("Time".into());
    let mut a = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let Value::Record(ref mut record) = a else {
        unreachable!()
    };
    let Value::Tensor(values) = &mut record.fields[4].1 else {
        unreachable!()
    };
    values.elem = TypeExpr::Dim(Box::new(TypeExpr::Float), length.clone());
    let session = Session::new();
    for name in ["num_csr_eigs_sym", "num_csr_eigs"] {
        let result = call(
            &session,
            name,
            &[
                a.clone(),
                Value::int(1),
                Value::String("largest_magnitude".into()),
                Value::Float(1e-8),
                Value::int(30),
                Value::int(0),
            ],
        )
        .unwrap();
        assert_eq!(result_field(&result, "converged"), &Value::Bool(true));
        for field in ["values_re", "residuals"] {
            let Value::Tensor(values) = result_field(&result, field) else {
                panic!("tensor")
            };
            assert_eq!(
                values.elem,
                TypeExpr::Dim(Box::new(TypeExpr::Float), length.clone())
            );
        }
        let Value::Tensor(vectors) = result_field(&result, "vectors_re") else {
            panic!("vectors")
        };
        assert_eq!(vectors.elem, TypeExpr::Float);
        assert_eq!(vectors.shape, vec![2, 1]);
    }
    let sigma = Value::Dyn(DynValue {
        value: Box::new(Value::Float(3.8)),
        dim: length.clone(),
    });
    let result = call(
        &session,
        "num_csr_eigs_shift_invert",
        &[
            a.clone(),
            Value::int(1),
            sigma,
            Value::Bool(true),
            Value::Float(1e-8),
            Value::int(30),
            Value::int(0),
        ],
    )
    .unwrap();
    let Value::Tensor(values) = result_field(&result, "values_re") else {
        panic!("values")
    };
    assert!((sparse_number(&values.data[0]) - 4.0).abs() < 1e-8);
    let mut rhs = vector(&[2.0, 8.0]);
    let Value::Tensor(ref mut v) = rhs else {
        unreachable!()
    };
    v.elem = TypeExpr::Dim(Box::new(TypeExpr::Float), time);
    for name in [
        "num_csr_amg_solve",
        "num_csr_amg_general_solve",
        "num_csr_cg_amg",
    ] {
        let result = call(
            &session,
            name,
            &[a.clone(), rhs.clone(), Value::Float(1e-8), Value::int(30)],
        )
        .unwrap();
        assert_eq!(result_field(&result, "converged"), &Value::Bool(true));
        let Value::Tensor(x) = result_field(&result, "x") else {
            panic!("x")
        };
        assert!(matches!(&x.elem, TypeExpr::Dim(..)));
        assert!((sparse_number(&x.data[1]) - 2.0).abs() < 1e-8);
        assert!(result_field(&result, "residual").as_f64().unwrap() < 1e-8);
    }
}

#[test]
fn sparse_eigen_adapters_reject_wrong_symmetry_controls_and_shift_units() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;
    let session = Session::new();
    let nonsymmetric = csr(&[0, 2, 3], &[0, 1, 1], &[2.0, 1.0, 4.0]);
    assert!(
        call(
            &session,
            "num_csr_eigs_sym",
            &[
                nonsymmetric,
                Value::int(1),
                Value::String("largest_real".into()),
                Value::Float(1e-8),
                Value::int(20),
                Value::int(0)
            ]
        )
        .is_err()
    );
    let good = csr(&[0, 1, 2], &[0, 1], &[2.0, 4.0]);
    let wrong = Value::Dyn(DynValue {
        value: Box::new(Value::Float(3.8)),
        dim: DimExpr::Base("Time".into()),
    });
    assert!(
        call(
            &session,
            "num_csr_eigs_shift_invert",
            &[
                good.clone(),
                Value::int(1),
                wrong,
                Value::Bool(true),
                Value::Float(1e-8),
                Value::int(20),
                Value::int(0)
            ]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "num_csr_eigs",
            &[
                good,
                Value::int(0),
                Value::String("largest_real".into()),
                Value::Float(1e-8),
                Value::int(20),
                Value::int(0)
            ]
        )
        .is_err()
    );
}

fn sparse_number(value: &Value) -> f64 {
    match value {
        | Value::Dyn(v) => sparse_number(&v.value),
        | other => other.as_f64().unwrap(),
    }
}

fn geometry_field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected geometry result record")
    };
    record.get(name).unwrap()
}
fn geometry_scalar(value: &Value) -> f64 {
    match value {
        | Value::Dyn(v) => geometry_scalar(&v.value),
        | _ => value.as_f64().unwrap(),
    }
}
fn measured_geometry(
    shape: Vec<usize>,
    values: &[f64],
    unit: &str,
) -> Value {
    Value::Tensor(TensorValue {
        shape,
        data: values.iter().copied().map(Value::Float).collect(),
        elem: TypeExpr::Float.with_dim(rssn::kernel::signature::DimExpr::Base(unit.into())),
    })
}
#[test]
fn scaled_euclidean_operators_preserve_physical_results() {
    use rssn::kernel::signature::DimExpr;
    let s = Session::new();
    let input = |values: &[f64], unit| measured_geometry(vec![values.len()], values, unit);
    let tiny = input(&[3e-200, 4e-200], "Length");
    let huge = input(&[3e200, 4e200], "Length");
    let zero = input(&[0.0, 0.0], "Length");
    let length = TypeExpr::Float.with_dim(DimExpr::Base("Length".into()));
    for (name, args, expected) in [
        ("num_norm", vec![tiny], 5e-200),
        ("num_distance", vec![huge.clone(), zero], 5e200),
    ] {
        let result = call(&s, name, &args).unwrap();
        assert!(rssn::kernel::overload::matches_type(&length, &result.type_expr()).accepted());
        assert!((geometry_scalar(&result) / expected - 1.0).abs() <= 2e-15);
    }
    let right_angle = call(
        &s,
        "num_angle",
        &[huge.clone(), input(&[-4e-200, 3e-200], "Time")],
    )
    .unwrap();
    assert!((right_angle.as_f64().unwrap() - std::f64::consts::FRAC_PI_2).abs() <= 2e-15);
    let Value::Tensor(direction) = call(&s, "num_normalize", &[huge]).unwrap() else {
        panic!("expected a tensor of normalized coordinates")
    };
    assert_eq!(direction.elem, TypeExpr::Float);
    for (value, expected) in direction.data.iter().zip([0.6, 0.8]) {
        assert!((geometry_scalar(value) - expected).abs() <= 2e-15);
    }
}

#[test]
fn reflection_restores_dimensions_after_a_compensated_overflow() {
    let s = Session::new();
    let tiny = f64::from_bits(1);
    let incoming = measured_geometry(vec![2], &[f64::MAX, tiny], "Length");
    let Value::Tensor(output) = call(&s, "num_reflect", &[incoming, vector(&[1.0, 0.0])]).unwrap()
    else {
        panic!("expected a reflected vector");
    };
    assert_eq!(
        output.elem,
        TypeExpr::Float.with_dim(rssn::kernel::signature::DimExpr::Base("Length".into()))
    );
    assert_eq!(
        geometry_scalar(&output.data[0]).to_bits(),
        (-f64::MAX).to_bits()
    );
    assert_eq!(geometry_scalar(&output.data[1]).to_bits(), 1);
}

#[test]
fn finite_interval_roots_share_measured_bounds_and_include_endpoints() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;

    let session = Session::new();
    let polynomial = call(
        &session,
        "num_native_polynomial_normalized_new",
        &[vector(&[1.0, -0.25, 0.0])],
    )
    .unwrap();
    let measured = |coordinate, dim: &str| {
        Value::Dyn(DynValue {
            dim: DimExpr::base(dim),
            value: Box::new(Value::Float(coordinate)),
        })
    };
    let length = |coordinate| measured(coordinate, "Length");
    let endpoint_roots = [
        polynomial,
        length(0.0),
        length(1.0),
        length(1e-10),
        length(4.0),
    ];
    let op = "num_native_real_roots_find_roots_in";
    let Value::Tensor(result) = call(&session, op, &endpoint_roots).unwrap() else {
        panic!("expected a measured root vector");
    };
    assert_eq!(result.shape, [2]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(DimExpr::base("Length")),
            &result.elem,
        )
        .accepted()
    );
    assert_eq!(
        result.data.iter().map(geometry_scalar).collect::<Vec<_>>(),
        [0.0, 1.0]
    );
    let mut invalid = endpoint_roots.clone();
    invalid[2] = measured(1.0, "Time");
    assert!(call(&session, op, &invalid).is_err());
    let mut invalid = endpoint_roots;
    invalid[1] = length(2.0);
    assert!(call(&session, op, &invalid).is_err());
}

#[test]
fn portable_kdtree_batch_preserves_query_order_and_coordinate_units() {
    let session = Session::new();
    let points = measured_geometry(vec![3, 2], &[0.0, 0.0, 2.0, 0.0, 0.0, 2.0], "Length");
    let model = call(&session, "num_kdtree_build", &[points]).unwrap();
    let queries = [1.75, 0.0, 0.0, 1.5, 0.25, 0.0];
    let op = "num_kdtree_nearest_batch";
    let result = call(
        &session,
        op,
        &[
            model.clone(),
            measured_geometry(vec![3, 2], &queries, "Length"),
        ],
    )
    .unwrap();
    let Value::Tensor(indices) = geometry_field(&result, "indices") else {
        panic!("expected indices");
    };
    assert_eq!(indices.elem, TypeExpr::Int);
    assert_eq!(indices.data, [Value::int(1), Value::int(2), Value::int(0)]);
    let Value::Tensor(distances) = geometry_field(&result, "distances") else {
        panic!("expected distances");
    };
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(rssn::kernel::signature::DimExpr::base("Length")),
            &distances.elem,
        )
        .accepted()
    );
    assert_eq!(
        distances
            .data
            .iter()
            .map(geometry_scalar)
            .collect::<Vec<_>>(),
        [0.25, 0.5, 0.25],
    );
    for invalid in [
        measured_geometry(vec![3, 2], &queries, "Time"),
        measured_geometry(vec![2, 3], &queries, "Length"),
    ] {
        assert!(call(&session, op, &[model.clone(), invalid]).is_err());
    }
    let empty = call(
        &session,
        op,
        &[model, measured_geometry(vec![0, 2], &[], "Length")],
    )
    .unwrap();
    assert!(matches!(
        geometry_field(&empty, "indices"),
        Value::Tensor(t) if t.shape == [0] && t.data.is_empty()
    ));
}

#[test]
fn portable_bvh_roundtrips_dimensions_and_rejects_corrupt_or_inconsistent_definitions() {
    let session = Session::new();
    let lower = measured_geometry(vec![2, 3], &[0.0, 0.0, 0.0, 10.0, 0.0, 0.0], "Length");
    let upper = measured_geometry(vec![2, 3], &[1.0, 1.0, 1.0, 11.0, 1.0, 1.0], "Length");
    let model = call(&session, "num_bvh_build", &[lower, upper]).unwrap();
    let serialized = serde_json::to_vec(&model).unwrap();
    let model: Value = serde_json::from_slice(&serialized).unwrap();
    let query = measured_geometry(vec![3], &[2.0, 0.5, 0.5], "Length");
    let nearest = call(&session, "num_bvh_nearest", &[model.clone(), query]).unwrap();
    let Value::Tensor(indices) = geometry_field(&nearest, "indices") else {
        panic!("indices")
    };
    assert_eq!(indices.data, [Value::int(0)]);
    let Value::Tensor(distances) = geometry_field(&nearest, "distances") else {
        panic!("distances")
    };
    assert_eq!(geometry_scalar(&distances.data[0]), 1.0);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(rssn::kernel::signature::DimExpr::Base("Length".into())),
            &distances.elem
        )
        .accepted()
    );
    let origin = measured_geometry(vec![3], &[-1.0, 0.5, 0.5], "Length");
    let hits = call(
        &session,
        "num_bvh_ray",
        &[model.clone(), origin.clone(), vector(&[1e-300, 0.0, 0.0])],
    )
    .unwrap();
    let Value::Tensor(hits) = hits else {
        panic!("hits")
    };
    let mut hit_indices: Vec<_> = hits.data.iter().map(geometry_scalar).collect();
    hit_indices.sort_by(f64::total_cmp);
    assert_eq!(hit_indices, [0.0, 1.0]);
    assert!(
        call(
            &session,
            "num_bvh_ray",
            &[
                model.clone(),
                origin.clone(),
                measured_geometry(vec![3], &[1.0, 0.0, 0.0], "Length")
            ]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "num_bvh_ray",
            &[model.clone(), origin.clone(), vector(&[0.0, 0.0, 0.0])]
        )
        .is_err()
    );
    for corrupt_version in [true, false] {
        let mut corrupt = model.clone();
        let Value::Record(record) = &mut corrupt else {
            unreachable!()
        };
        let Value::Opaque(payload) = &mut record
            .fields
            .iter_mut()
            .find(|(name, _)| name == "payload")
            .unwrap()
            .1
        else {
            unreachable!()
        };
        if corrupt_version {
            payload.kind = "rssn/bvh-defining-data/v99".into();
        } else {
            payload.bytes.push(0);
        }
        assert!(call(&session, "num_bvh_nearest", &[corrupt, origin.clone()]).is_err());
    }
    assert!(
        call(
            &session,
            "num_bvh_build",
            &[
                measured_geometry(vec![1, 3], &[2.0, 0.0, 0.0], "Length"),
                measured_geometry(vec![1, 3], &[1.0, 0.0, 0.0], "Length")
            ]
        )
        .is_err()
    );
}
#[test]
fn voronoi_topology_is_flat_and_only_vertices_carry_coordinate_units() {
    let session = Session::new();
    let sites = measured_geometry(
        vec![5, 2],
        &[0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0, 2.0, 1.0, 1.0],
        "Length",
    );
    let diagram = call(&session, "num_voronoi_2d", &[sites]).unwrap();
    let Value::Tensor(vertices) = geometry_field(&diagram, "vertices") else {
        panic!("vertices")
    };
    assert_eq!(vertices.shape, [4, 2]);
    assert!(
        rssn::kernel::overload::matches_type(
            &TypeExpr::Float.with_dim(rssn::kernel::signature::DimExpr::Base("Length".into())),
            &vertices.elem
        )
        .accepted()
    );
    let Value::Tensor(directions) = geometry_field(&diagram, "ray_directions") else {
        panic!("ray directions")
    };
    assert_eq!(directions.elem, TypeExpr::Float);
    let Value::Tensor(offsets) = geometry_field(&diagram, "cell_offsets") else {
        panic!("offsets")
    };
    assert_eq!(offsets.shape, [6]);
    let Value::Tensor(indices) = geometry_field(&diagram, "cell_indices") else {
        panic!("cell indices")
    };
    assert_eq!(
        offsets.data.last().unwrap(),
        &Value::int(indices.data.len())
    );
    assert!(indices.data.iter().all(|v| geometry_scalar(v) < 4.0));
    assert!(
        call(
            &session,
            "num_voronoi_2d",
            &[matrix(2, 2, &[0.0, 0.0, 1.0, 1.0])]
        )
        .is_err()
    );
}

#[test]
fn kruskal_wallis_accepts_ragged_groups_with_common_units() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::RecordValue;
    let session = Session::new();
    let input = |offsets: &[usize], values: &[f64]| {
        let offsets = Value::Tensor(TensorValue {
            shape: vec![offsets.len()],
            elem: TypeExpr::Int,
            data: offsets.iter().copied().map(Value::int).collect(),
        });
        let values = Value::Tensor(TensorValue {
            shape: vec![values.len()],
            elem: TypeExpr::Float.with_dim(DimExpr::base("Length")),
            data: values.iter().copied().map(Value::Float).collect(),
        });
        Value::Record(RecordValue {
            fields: vec![("offsets".into(), offsets), ("values".into(), values)],
        })
    };
    let call =
        |argument| session.call_value("num_stats_kruskal_wallis", &[argument], &CallCtx::new());
    let result = call(input(
        &[0, 2, 5, 9],
        &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
    ))
    .unwrap();
    let Value::Record(result) = result else {
        panic!("expected a rank test record")
    };
    assert!((result.get("statistic").unwrap().as_f64().unwrap() - 7.0).abs() < 2e-14);
    assert!((result.get("p_value").unwrap().as_f64().unwrap() - (-3.5_f64).exp()).abs() < 2e-14);
    let Some(Value::Variant(df)) = result.get("df") else {
        panic!("expected degrees of freedom")
    };
    assert_eq!(df.tag, "Some");
    assert_eq!(df.payload, [Value::Float(2.0)]);
    assert!(call(input(&[0, 2, 2, 3], &[1.0, 2.0, 3.0])).is_err());
    assert!(call(input(&[0, 2, 5, 9], &[1.0; 9])).is_err());
}

#[test]
fn bell_batch_keeps_input_order_and_localizes_overflow() {
    let session = Session::new();
    let indices = Value::Tensor(TensorValue {
        shape: vec![6],
        data: [218, 219, 0, 5, 218, u64::MAX]
            .into_iter()
            .map(Value::int)
            .collect(),
        elem: TypeExpr::Int,
    });
    let value = call(&session, "num_combinatorics_bell_batch", &[indices]).unwrap();
    let Value::Tensor(value) = value else {
        panic!("expected a vector of optional counts")
    };
    assert_eq!(value.shape, [6]);
    let count = |index| {
        let Value::Variant(value) = &value.data[index] else {
            panic!("expected an optional Bell number")
        };
        match (value.tag.as_str(), value.payload.as_slice()) {
            | ("None", []) => None,
            | ("Some", [Value::Float(count)]) => Some(*count),
            | _ => panic!("invalid count option"),
        }
    };
    assert!(count(0).is_some_and(|x| x.is_finite() && x > 1e300));
    assert_eq!(count(1), None);
    assert_eq!(count(2), Some(1.0));
    assert_eq!(count(3), Some(52.0));
    assert_eq!(count(4), count(0));
    assert_eq!(count(5), None);
}

#[test]
fn native_prime_sqrt_plan_checks_modulus_and_reuses_roots() {
    let session = Session::new();
    let new = "num_number_theory_prime_sqrt_plan_new";
    let sqrt = "num_number_theory_prime_sqrt_plan_sqrt";
    let modulus = "num_number_theory_prime_sqrt_plan_modulus";
    let invalid = call(&session, new, &[Value::int(49)]).unwrap();
    assert!(matches!(invalid, Value::Variant(value) if value.tag == "None"));
    let prepared = call(&session, new, &[Value::int(17)]).unwrap();
    let Value::Variant(prepared) = prepared else {
        panic!("expected an optional plan")
    };
    assert_eq!(prepared.tag, "Some");
    let plan = &prepared.payload[0];
    assert_eq!(
        call(&session, modulus, std::slice::from_ref(plan)).unwrap(),
        Value::int(17)
    );
    for value in [0, 1, 2, 3, 4, 5, 16, 18, u64::MAX] {
        let expected = call(
            &session,
            "num_number_theory_sqrt_mod_prime",
            &[Value::int(value), Value::int(17)],
        )
        .unwrap();
        let actual = call(&session, sqrt, &[plan.clone(), Value::int(value)]).unwrap();
        assert_eq!(actual, expected);
    }
    let other_session = Session::new();
    assert!(call(&other_session, modulus, std::slice::from_ref(plan)).is_err());
}
