//! Dimensionless exact number contracts through the public numerical ABI.
use num_bigint::BigInt;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::RecordValue;
use rssn::kernel::value::TensorValue;
use rssn::kernel::value::VariantValue;

fn ints(values: impl IntoIterator<Item = impl Into<BigInt>>) -> Value {
    let data: Vec<_> = values.into_iter().map(Value::int).collect();
    Value::Tensor(TensorValue::of(vec![data.len()], data, TypeExpr::Int))
}

fn option(value: Option<Value>) -> Value {
    Value::Variant(VariantValue {
        tag: if value.is_some() {
            "Some"
        } else {
            "None"
        }
        .into(),
        payload: value.into_iter().collect(),
    })
}

#[test]
fn exact_number_exports_keep_signed_bigints_and_option_contracts() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let call = |name, values: &[Value]| session.call_value(name, values, &ctx).unwrap();
    let scale = BigInt::from(1) << 128_u32;
    let a: BigInt = (-6_i32) * &scale;
    let b: BigInt = 15_i32 * &scale;
    let args = [Value::int(a.clone()), Value::int(b.clone())];
    let gcd = call("num_number_exact_bigint_gcd", &args);
    assert_eq!(gcd, Value::int(3_i32 * &scale));
    let Value::Record(RecordValue { fields }) = call("num_number_exact_bigint_xgcd", &args) else {
        panic!("expected the extended GCD record");
    };
    let [
        (_, Value::Int(divisor)),
        (_, Value::Int(x)),
        (_, Value::Int(y)),
    ] = fields.as_slice()
    else {
        panic!("expected three exact integer fields");
    };
    assert_eq!(Value::int(divisor.clone()), gcd);
    assert_eq!(a * x + b * y, *divisor);

    let root = (BigInt::from(1) << 64_u32) + 1_u8;
    let square = &root * &root;
    let cube = &square * &root;
    assert_eq!(
        call(
            "num_number_exact_bigint_isqrt",
            &[Value::int(&square + &root)]
        ),
        option(Some(Value::int(root.clone())))
    );
    assert_eq!(
        call(
            "num_number_exact_bigint_nth_root",
            &[Value::int(-cube - 1_u8), Value::int(3)]
        ),
        option(Some(Value::int(-root)))
    );
    assert_eq!(
        call("num_number_exact_bigint_isqrt", &[Value::int(-1)]),
        option(None)
    );
    assert_eq!(
        call(
            "num_number_exact_bigint_nth_root",
            &[Value::int(-1), Value::int(2)]
        ),
        option(None)
    );
    assert_eq!(
        call(
            "num_number_exact_jacobi_bigint",
            &[Value::int(7), Value::int(21)]
        ),
        option(Some(Value::int(0)))
    );
    assert_eq!(
        call(
            "num_number_exact_jacobi_bigint",
            &[Value::int(5), Value::int(20)]
        ),
        option(None)
    );
    assert_eq!(
        call("num_number_exact_strong_lucas_bigint", &[Value::int(97)]),
        Value::Bool(true)
    );
    assert_eq!(
        call("num_number_exact_strong_lucas_bigint", &[Value::int(91)]),
        Value::Bool(false)
    );

    let modulus: BigInt = 4_i32 * &scale + 1_u8;
    assert_eq!(
        call(
            "num_number_exact_bigint_mod_inverse",
            &[Value::int(2), Value::int(modulus)]
        ),
        option(Some(Value::int(2_i32 * &scale + 1_u8)))
    );
    assert_eq!(
        call(
            "num_number_exact_bigint_mod_inverse",
            &[Value::int(6), Value::int(15)]
        ),
        option(None)
    );
}

#[test]
fn exact_vectors_validate_width_shapes_and_kernel_work() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let call = |name, args: &[Value]| session.call_value(name, args, &ctx);
    let batch = "num_number_ext_batch_mod_pow";
    for (modulus, expected) in [(23, [10, 1, 21]), (24, [21, 1, 0])] {
        let args = [ints([-3, 7, 42]), ints([5, 0, 17]), Value::int(modulus)];
        assert_eq!(call(batch, &args).unwrap(), ints(expected));
    }
    let empty = ints([] as [u64; 0]);
    assert_eq!(
        call(batch, &[empty.clone(), empty.clone(), Value::int(7)]).unwrap(),
        empty
    );
    let bad_batch = [
        [ints([5, 7]), ints([1]), Value::int(9)],
        [ints([5]), ints([-1]), Value::int(9)],
        [ints([5]), ints([1]), Value::int(0)],
        [ints([5]), ints([u128::MAX]), Value::int(9)],
    ];
    for args in bad_batch {
        assert!(matches!(call(batch, &args), Err(OpError::Domain(_))));
    }
    let malformed = Value::Tensor(TensorValue::of(vec![2], vec![Value::int(7)], TypeExpr::Int));
    assert!(call(batch, &[malformed, ints([1]), Value::int(9)]).is_err());

    let garner = "num_number_exact_bigint_crt_garner";
    assert_eq!(
        call(garner, &[ints([2, 3, 4]), ints([5, 7, 11])]).unwrap(),
        option(Some(Value::int(367)))
    );
    for pair in [
        [empty.clone(), empty.clone()],
        [ints([1]), ints([2, 3])],
        [ints([1]), ints([0])],
        [ints([0, 0]), ints([4, 6])],
    ] {
        assert_eq!(call(garner, &pair).unwrap(), option(None));
    }
    let too_many = ints(std::iter::repeat_n(1, 1025));
    assert!(matches!(
        call(garner, &[too_many.clone(), too_many]),
        Err(OpError::Domain(message)) if message.contains("work/allocation")
    ));

    let crt = "num_number_ext_crt_bigint";
    let large: BigInt = BigInt::from(1) << 130_u32;
    let first = &large - 1_u8;
    assert_eq!(
        call(
            crt,
            &[
                ints([first.clone(), BigInt::from(3)]),
                ints([large.clone(), BigInt::from(6)])
            ]
        )
        .unwrap(),
        Value::Record(RecordValue {
            fields: vec![
                ("item0".into(), Value::int(first)),
                ("item1".into(), Value::int(3_i32 * large))
            ]
        })
    );
    assert!(matches!(
        call(crt, &[empty.clone(), empty]),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(crt, &[ints([0, 1]), ints([4, 6])]),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        call(crt, &[ints([0]), ints([4, 6])]),
        Err(OpError::Domain(_))
    ));
    let too_wide: BigInt = BigInt::from(1) << 65_536_u32;
    assert!(matches!(
        call(crt, &[ints([BigInt::from(0)]), ints([too_wide])]),
        Err(OpError::Domain(message)) if message.contains("coefficient bit budget")
    ));

    let totient = "num_number_ext_totient_sieve";
    assert_eq!(call(totient, &[Value::int(0)]).unwrap(), ints([0]));
    assert_eq!(
        call(totient, &[Value::int(10)]).unwrap(),
        ints([0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4])
    );
    assert!(matches!(
        call(totient, &[Value::int(usize::MAX)]),
        Err(OpError::Domain(_))
    ));
    let measured = Value::Dyn(DynValue {
        value: Box::new(Value::int(2)),
        dim: DimExpr::base("Length"),
    });
    assert!(call("num_number_exact_bigint_gcd", &[measured, Value::int(3)]).is_err());
}
