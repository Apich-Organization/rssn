//! Prepared prime fields and exact sparse native objects.
use num_bigint::BigInt;
use num_rational::BigRational;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernels::sparse_exact::Gfp;
use rssn::kernels::sparse_exact::word_primes;

const PRIME: u64 = 18_446_744_073_709_551_557;
const PREFIX: &str = "num_linear_sparse_exact_";

fn call(
    session: &Session,
    operation: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    session.call_value(operation, args, &CallCtx::new())
}

fn exact(
    session: &Session,
    operation: &str,
    args: &[Value],
) -> Result<Value, OpError> {
    call(session, &format!("{PREFIX}{operation}"), args)
}

fn tensor<const N: usize>(
    shape: [usize; N],
    elem: TypeExpr,
    data: Vec<Value>,
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        elem,
        data,
    })
}

fn rational(
    values: impl IntoIterator<Item = BigRational>,
    shape: &[usize],
) -> Value {
    Value::Tensor(TensorValue {
        shape: shape.to_vec(),
        elem: TypeExpr::Real,
        data: values.into_iter().map(Value::rational).collect(),
    })
}

fn r(
    n: impl Into<BigInt>,
    d: impl Into<BigInt>,
) -> BigRational {
    BigRational::new(n.into(), d.into())
}

fn rational_matrix(
    session: &Session,
    rows: usize,
    columns: usize,
    values: impl IntoIterator<Item = BigRational>,
) -> Value {
    call(
        session,
        "num_linear_rational_sparse_from_dense",
        &[
            Value::int(rows),
            Value::int(columns),
            rational(values, &[rows * columns]),
        ],
    )
    .unwrap()
}

fn variant(
    value: Value,
    tag: &str,
) -> Vec<Value> {
    let Value::Variant(value) = value else {
        panic!("expected an Option")
    };
    assert_eq!(value.tag, tag);
    value.payload
}

fn field<'a>(
    value: &'a Value,
    name: &str,
) -> &'a Value {
    let Value::Record(record) = value else {
        panic!("expected a record")
    };
    &record.fields.iter().find(|(key, _)| key == name).unwrap().1
}

#[test]
fn prepared_fields_check_primality_and_keep_the_full_word() {
    let s = Session::new();
    for invalid_modulus in [0_u64, 1, 15, 341, 3_825_123_056_546_413_051] {
        assert!(exact(&s, "field_new", &[Value::int(invalid_modulus)]).is_err());
    }
    let f = exact(&s, "field_new", &[Value::int(PRIME)]).unwrap();
    let fargs = std::slice::from_ref(&f);
    assert_eq!(exact(&s, "field_prime", fargs).unwrap(), Value::int(PRIME));
    let inputs = [f.clone(), Value::int(u64::MAX), Value::int(u64::MAX - 7)];
    assert_eq!(exact(&s, "field_add", &inputs).unwrap(), Value::int(109));
    assert_eq!(exact(&s, "field_sub", &inputs).unwrap(), Value::int(7));
    assert_eq!(exact(&s, "field_mul", &inputs).unwrap(), Value::int(2_958));
    assert_eq!(
        exact(
            &s,
            "field_reduce",
            &[f.clone(), Value::int(-(BigInt::from(1) << 256usize))]
        )
        .unwrap(),
        Value::int(Gfp::new(PRIME).reduce(&(-(BigInt::from(1) << 256usize))))
    );
    let operands = [PRIME - 1, 17, u64::MAX, PRIME - 7];
    let data = tensor(
        [4],
        TypeExpr::Int,
        operands.into_iter().map(Value::int).collect(),
    );
    let result = exact(&s, "field_batch_inv", &[f.clone(), data]).unwrap();
    let inner = variant(result, "Some");
    let [Value::Tensor(inverses)] = inner.as_slice() else {
        panic!("expected a tensor of inverse residues")
    };
    for (inverse, value) in inverses.data.iter().zip(operands) {
        assert_eq!(
            exact(
                &s,
                "field_mul",
                &[f.clone(), inverse.clone(), Value::int(value)]
            )
            .unwrap(),
            Value::int(1)
        );
    }
    let bad = tensor([2], TypeExpr::Int, vec![Value::int(1), Value::int(PRIME)]);
    assert!(
        variant(
            exact(&s, "field_batch_inv", &[f.clone(), bad]).unwrap(),
            "None"
        )
        .is_empty()
    );
    assert!(exact(&Session::new(), "field_prime", fargs).is_err());
    let plan = exact(&s, "field_reducer", fargs).unwrap();
    assert_eq!(
        call(&s, "num_native_release", fargs).unwrap(),
        Value::Bool(true)
    );
    assert_eq!(
        call(&s, "num_modular_barrett_modulus", &[plan]).unwrap(),
        Value::int(PRIME)
    );
}

#[test]
fn rational_reduction_and_prime_rref_keep_exact_residues() {
    let s = Session::new();
    let f = exact(&s, "field_new", &[Value::int(7)]).unwrap();
    let a = rational_matrix(
        &s,
        2,
        3,
        [r(1, 2), r(0, 1), r(-3, 4), r(0, 1), r(2, 3), r(1, 1)],
    );
    let matrix = variant(
        exact(&s, "rational_to_prime", &[a, f.clone()]).unwrap(),
        "Some",
    )
    .pop()
    .unwrap();
    let dense = |m| exact(&s, "prime_to_dense", &[m]).unwrap();
    assert_eq!(
        dense(matrix.clone()),
        tensor(
            [2, 3],
            TypeExpr::Int,
            [4, 0, 1, 0, 3, 1].map(Value::int).to_vec()
        )
    );
    let reduced = exact(&s, "prime_rref", &[matrix]).unwrap();
    assert_eq!(
        field(&reduced, "pivots"),
        &tensor([2], TypeExpr::Int, [0, 1].map(Value::int).to_vec())
    );
    assert_eq!(
        dense(field(&reduced, "matrix").clone()),
        tensor(
            [2, 3],
            TypeExpr::Int,
            [1, 0, 2, 0, 1, 5].map(Value::int).to_vec()
        )
    );
    let bad = rational_matrix(&s, 1, 1, [r(1, 7)]);
    assert!(variant(exact(&s, "rational_to_prime", &[bad, f]).unwrap(), "None").is_empty());
}

#[test]
fn rational_rref_and_cached_order_validate_shapes() {
    let s = Session::new();
    let coeff = BigInt::from(1) << 384usize;
    let a = rational_matrix(
        &s,
        2,
        3,
        [
            r(3, 1),
            r(0, 1),
            r(coeff.clone(), 1),
            r(0, 1),
            r(0, 1),
            r(0, 1),
        ],
    );
    assert!(exact(&s, "rcm_plan", std::slice::from_ref(&a)).is_err());
    let result = exact(&s, "rational_rref", &[a]).unwrap();
    let dense = call(
        &s,
        "num_linear_rational_sparse_to_dense",
        std::slice::from_ref(field(&result, "matrix")),
    )
    .unwrap();
    assert_eq!(
        field(&result, "pivots"),
        &tensor([1], TypeExpr::Int, vec![Value::int(0)])
    );
    assert_eq!(
        field(&dense, "values"),
        &rational(
            [r(1, 1), r(0, 1), r(coeff, 3), r(0, 1), r(0, 1), r(0, 1)],
            &[6]
        )
    );
    let a = rational_matrix(&s, 2, 2, [r(1, 1), r(1, 1), r(1, 1), r(0, 1)]);
    let args = std::slice::from_ref(&a);
    let plan = exact(&s, "rcm_plan", args).unwrap();
    assert_eq!(plan, exact(&s, "rcm_plan", args).unwrap());
    let Value::Tensor(perm) = field(&plan, "permutation") else {
        panic!("permutation vector")
    };
    let Value::Tensor(inverse) = field(&plan, "inverse") else {
        panic!("inverse permutation")
    };
    for (index, old) in perm.data.iter().enumerate() {
        let Value::Int(old) = old else {
            panic!("exact row index")
        };
        assert_eq!(
            inverse.data[usize::try_from(old).unwrap()],
            Value::int(index)
        );
    }
}

#[test]
fn dixon_finds_a_compatible_prime_and_replays_large_exact_right_hand_sides() {
    let s = Session::new();
    let primes = word_primes().take(4).map(BigInt::from).product::<BigInt>();
    let a = rational_matrix(&s, 1, 1, [r(primes.clone(), 1)]);
    let rhs = rational([r(1, 1), r(&primes * 3, 1)], &[2, 1]);
    let got = variant(exact(&s, "dixon_multi", &[a, rhs.clone()]).unwrap(), "Some")
        .pop()
        .unwrap();
    assert_eq!(got, rational([r(1, primes), r(3, 1)], &[2, 1]));
    let singular = rational_matrix(&s, 1, 1, [r(0, 1)]);
    assert!(
        variant(
            exact(&s, "dixon_multi", &[singular.clone(), rhs]).unwrap(),
            "None"
        )
        .is_empty()
    );
    assert!(
        exact(
            &s,
            "dixon_multi",
            &[singular, rational([r(1, 1), r(2, 1)], &[1, 2])]
        )
        .is_err()
    );
    let costly = rational_matrix(&s, 1, 1, [r(BigInt::from(1) << 9_000usize, 1)]);
    assert!(matches!(
        exact(&s, "dixon_multi", &[costly, rational([r(1, 1)], &[1, 1])]),
        Err(OpError::BudgetExceeded)
    ));
}

#[test]
fn sparse_exact_manifest_declares_the_ten_missing_sources() {
    let manifest = Session::new().manifest();
    let entries = manifest
        .ops
        .iter()
        .filter(|entry| entry.name.starts_with(PREFIX))
        .collect::<Vec<_>>();
    assert_eq!(entries.len(), 16);
    assert_eq!(
        entries
            .iter()
            .filter(|entry| entry.rust_path.is_some())
            .count(),
        10
    );
    assert!(entries.iter().all(|entry| {
        entry.signature_source == rssn::kernel::manifest::SignatureSource::Declared
    }));
}

#[test]
fn sparse_rref_preserves_zero_row_rectangular_shapes() {
    let session = Session::new();
    let matrix = rational_matrix(&session, 0, 3, []);
    assert!(exact(&session, "rcm_plan", std::slice::from_ref(&matrix)).is_err());
    let rref = exact(&session, "rational_rref", &[matrix]).unwrap();
    assert_eq!(field(&rref, "pivots"), &tensor([0], TypeExpr::Int, vec![]));
    let dense = call(
        &session,
        "num_linear_rational_sparse_to_dense",
        std::slice::from_ref(field(&rref, "matrix")),
    )
    .unwrap();
    assert_eq!(field(&dense, "rows"), &Value::int(0));
    assert_eq!(field(&dense, "cols"), &Value::int(3));
}
