//! Portable generic models invoke native methods, reject forged data and preserve exact arithmetic.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
fn matrix(
    values: &[f64],
    rows: usize,
    cols: usize,
) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![rows, cols],
        elem: TypeExpr::Float,
        data: values.iter().copied().map(Value::Float).collect(),
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
fn portable_matrix_methods_validate_kind_indices_and_allocation_before_native_calls() {
    let session = Session::new();
    let model = call(
        &session,
        "num_matrix_float_new",
        &[matrix(&[1.0, 2.0, 3.0, 4.0], 2, 2), Value::int(0)],
    )
    .unwrap();
    assert!(
        call(
            &session,
            "num_matrix_float_get",
            &[model.clone(), Value::int(2), Value::int(0)]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "num_matrix_float_zeros",
            &[Value::int(usize::MAX), Value::int(2), Value::Float(1.0)]
        )
        .is_err()
    );
    let Value::Record(mut forged) = model.clone() else {
        panic!("record")
    };
    forged
        .fields
        .iter_mut()
        .find(|(name, _)| name == "format")
        .unwrap()
        .1 = Value::String("rssn/matrix-float/v999".into());
    assert!(call(&session, "num_matrix_float_trace", &[Value::Record(forged)]).is_err());
    let mut malformed = matrix(&[1.0], 1, 1);
    if let Value::Tensor(t) = &mut malformed {
        t.shape = vec![2, 2];
    }
    assert!(
        call(
            &session,
            "num_matrix_float_new",
            &[malformed, Value::int(0)]
        )
        .is_err()
    );
    assert!(
        call(
            &session,
            "num_matrix_real_new",
            &[matrix(&[1.0], 1, 1), Value::int(0)]
        )
        .is_err()
    );
}
#[test]
fn generic_exact_matrix_uses_updated_pivots_and_never_discards_small_nonzero_rationals() {
    use num_rational::BigRational;
    use rssn::kernels::matrix::Matrix;
    let values = [1, 1, 0, 1, 1, 1, 0, 1, 0]
        .into_iter()
        .map(|n| BigRational::from_integer(n.into()))
        .collect();
    let a = Matrix::new(3, 3, values);
    assert_eq!(
        a.determinant_lu().unwrap(),
        BigRational::from_integer((-1).into())
    );
    let tiny = BigRational::new(1.into(), num_bigint::BigInt::from(10).pow(100));
    let mut a = Matrix::new(
        2,
        2,
        vec![
            tiny.clone(),
            BigRational::from_integer(0.into()),
            BigRational::from_integer(0.into()),
            tiny,
        ],
    );
    assert_eq!(a.rref().unwrap(), 2);
}
#[test]
fn linear_chebyshev_roots_and_constant_polynomial_remainders_are_valid_models() {
    let roots = rssn::kernels::approx::ChebSeries {
        coeffs: vec![1.0, 2.0],
        a: 0.0,
        b: 2.0,
    }
    .roots()
    .unwrap();
    assert_eq!(roots, vec![0.5]);
    let session = Session::new();
    let vector = |data: &[f64]| matrix(data, 1, data.len());
    let values = |data: &[f64]| {
        let Value::Tensor(mut t) = vector(data) else {
            unreachable!()
        };
        t.shape = vec![data.len()];
        Value::Tensor(t)
    };
    let a = call(
        &session,
        "num_polynomial_model_new",
        &[values(&[1.0, 2.0]), Value::Float(1.0)],
    )
    .unwrap();
    let b = call(
        &session,
        "num_polynomial_model_new",
        &[values(&[1.0]), Value::Float(1.0)],
    )
    .unwrap();
    let division = call(&session, "num_polynomial_model_long_division", &[a, b]).unwrap();
    let Value::Record(record) = division else {
        panic!("record")
    };
    let remainder = &record
        .fields
        .iter()
        .find(|(name, _)| name == "remainder")
        .unwrap()
        .1;
    assert!(
        call(
            &session,
            "num_polynomial_model_eval",
            &[remainder.clone(), Value::Float(1.0)]
        )
        .is_ok()
    );
}

#[test]
fn measured_polynomial_fit_composition_and_batch_preserve_axis_and_ordinate_units() {
    use rssn::kernel::signature::DimExpr;
    use rssn::kernel::value::DynValue;

    let session = Session::new();
    let length = DimExpr::base("Length");
    let time = DimExpr::base("Time");
    let mass = DimExpr::base("Mass");
    let scalar = |value, dim: &DimExpr| {
        Value::Dyn(DynValue {
            dim: dim.clone(),
            value: Box::new(Value::Float(value)),
        })
    };
    let samples = |xs: &[f64], dim: &DimExpr| {
        Value::Tensor(TensorValue {
            shape: vec![xs.len()],
            elem: TypeExpr::Float.with_dim(dim.clone()),
            data: xs.iter().copied().map(Value::Float).collect(),
        })
    };
    // outer(length) = 3 * length + 1 in mass units.
    let outer = call(
        &session,
        "num_polynomial_model_interpolate_lagrange",
        &[
            samples(&[0.0, 2.0], &length),
            samples(&[1.0, 7.0], &mass),
            scalar(2.0, &length),
        ],
    )
    .unwrap();
    // inner(time) = 2 * time in length units, with a different coordinate scale.
    let inner = call(
        &session,
        "num_polynomial_model_new",
        &[samples(&[8.0, 0.0], &length), scalar(4.0, &time)],
    )
    .unwrap();
    let composed = call(
        &session,
        "num_polynomial_model_compose",
        &[outer.clone(), inner.clone()],
    )
    .unwrap();
    let name = "num_polynomial_model_eval_batch";
    for points in [
        vec![],
        vec![3.0, -1.0, 0.0],
        (0..257).map(f64::from).collect(),
    ] {
        let Value::Tensor(result) =
            call(&session, name, &[composed.clone(), samples(&points, &time)]).unwrap()
        else {
            panic!("expected an ordinate tensor");
        };
        assert_eq!(result.elem, TypeExpr::Float.with_dim(mass.clone()));
        assert_eq!(result.shape, [points.len()]);
        for (value, x) in result.data.iter().zip(points) {
            let Value::Dyn(measured) = value else {
                panic!("expected a measured ordinate")
            };
            assert_eq!(measured.dim, mass);
            assert_eq!(measured.value.as_f64(), Some(6.0 * x + 1.0));
        }
    }
    assert!(call(&session, name, &[composed, samples(&[1.0], &length)]).is_err());
    assert!(call(&session, "num_polynomial_model_compose", &[inner, outer]).is_err());
}

#[test]
fn monic_polynomial_at_the_degree_limit_can_be_reopened_and_evaluated() {
    let session = Session::new();
    let samples = |x: &[f64]| {
        Value::Tensor(TensorValue::of(
            vec![x.len()],
            x.iter().copied().map(Value::Float).collect(),
            TypeExpr::Float,
        ))
    };
    let make = |roots: &[f64]| {
        call(
            &session,
            "num_polynomial_model_from_roots",
            &[samples(roots), Value::Float(2.0)],
        )
    };
    let polynomial = make(&[0.0; 128]).unwrap();
    let at = |model, x| {
        call(
            &session,
            "num_polynomial_model_eval",
            &[model, Value::Float(x)],
        )
        .unwrap()
    };
    assert!((at(polynomial, 2.0).as_f64().unwrap() - 1.0).abs() < 1e-12);
    assert!(make(&[0.0; 129]).is_err());
    let repeated = make(&[2.0, 2.0]).unwrap();
    assert_eq!(at(repeated.clone(), 2.0).as_f64(), Some(0.0));
    assert_eq!(at(repeated.clone(), 4.0).as_f64(), Some(1.0));
    for name in ["mul_schoolbook", "mul_karatsuba", "mul_fft"] {
        let product = call(
            &session,
            &format!("num_polynomial_model_{name}"),
            &[repeated.clone(), repeated.clone()],
        )
        .unwrap();
        assert!((at(product, 6.0).as_f64().unwrap() - 16.0).abs() < 1e-11);
    }
}

#[test]
fn checked_polynomial_division_binds_the_new_method_and_keeps_finite_remainders() {
    let session = Session::new();
    let polynomial = |coefficients: &[f64]| {
        let coefficients = Value::Tensor(TensorValue {
            shape: vec![coefficients.len()],
            elem: TypeExpr::Float,
            data: coefficients.iter().copied().map(Value::Float).collect(),
        });
        call(
            &session,
            "num_polynomial_model_new",
            &[coefficients, Value::Float(1.0)],
        )
        .unwrap()
    };
    let dividend = polynomial(&[2.0, f64::MAX]);
    let divisor = polynomial(&[0.0, 1.0, f64::MAX]);
    let name = "num_polynomial_model_try_long_division";
    let Value::Record(result) = call(&session, name, &[dividend.clone(), divisor]).unwrap() else {
        panic!("expected the quotient and remainder");
    };
    for (field, coefficient) in [("quotient", 2.0), ("remainder", -f64::MAX)] {
        let model = result.get(field).unwrap();
        let value = call(
            &session,
            "num_polynomial_model_eval",
            &[model.clone(), Value::Float(0.0)],
        )
        .unwrap();
        assert_eq!(value, Value::Float(coefficient));
    }
    assert!(call(&session, name, &[dividend, polynomial(&[0.0])]).is_err());
    let manifest = session.manifest();
    assert_eq!(
        manifest.op(name).unwrap().rust_path.as_deref(),
        Some("crate::kernels::polynomial::Polynomial::try_long_division")
    );
}

#[test]
fn portable_factors_reconstruct_native_solvers_and_reject_kind_permutation_forgery() {
    use rssn::kernel::value::RecordValue;
    let session = Session::new();
    let factor = call(
        &session,
        "num_dense_lu_model_build",
        &[matrix(&[2.0, 0.0, 0.0, 4.0], 2, 2)],
    )
    .unwrap();
    let mut rhs = matrix(&[6.0, 8.0], 1, 2);
    if let Value::Tensor(t) = &mut rhs {
        t.shape = vec![2];
    }
    let x = call(&session, "num_dense_lu_model_solve", &[factor, rhs]).unwrap();
    let Value::Tensor(t) = x else {
        panic!("vector")
    };
    assert_eq!(t.data, vec![Value::Float(3.0), Value::Float(2.0)]);
    let ints = |data: &[usize]| {
        Value::Tensor(TensorValue {
            shape: vec![data.len()],
            elem: TypeExpr::Int,
            data: data.iter().copied().map(Value::int).collect(),
        })
    };
    let mut data = matrix(&[2.0, 4.0], 1, 2);
    if let Value::Tensor(t) = &mut data {
        t.shape = vec![2];
    }
    let csr = Value::Record(RecordValue {
        fields: vec![
            ("rows".into(), Value::int(2)),
            ("cols".into(), Value::int(2)),
            ("indptr".into(), ints(&[0, 1, 2])),
            ("indices".into(), ints(&[0, 1])),
            ("values".into(), data),
        ],
    });
    assert!(
        call(
            &session,
            "num_csr_model_permute",
            &[csr.clone(), ints(&[0, 0])]
        )
        .is_err()
    );
    let ilu = call(&session, "num_sparse_ilu_model_build", &[csr]).unwrap();
    assert!(call(&session, "num_sparse_lu_model_dim", &[ilu]).is_err());
}
