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
