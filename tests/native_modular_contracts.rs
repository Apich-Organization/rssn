//! ABI regressions for typed exact modular arithmetic.
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;

fn residues(values: &[u64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })
}

#[test]
fn vector_arithmetic_keeps_full_width_residues() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let m = u64::MAX;
    let a = residues(&[m - 1, 2, m, m - 2, 3]);
    let b = residues(&[m - 2, 3, 7, m - 3, 4]);
    let pair = [a.clone(), b.clone(), Value::int(m)];
    for (op, expected) in [
        (
            "num_finite_field_slice_add_mod",
            residues(&[m - 3, 5, 7, m - 5, 7]),
        ),
        (
            "num_finite_field_slice_sub_mod",
            residues(&[1, m - 1, m - 7, 1, m - 1]),
        ),
        ("num_finite_field_dot_mod", Value::int(26)),
    ] {
        assert_eq!(
            session.call_value(op, &pair, &ctx).unwrap(),
            expected,
            "{op}"
        );
    }
    assert_eq!(
        session
            .call_value(
                "num_finite_field_slice_scale_mod",
                &[a.clone(), Value::int(m - 1), Value::int(m)],
                &ctx
            )
            .unwrap(),
        residues(&[1, m - 2, 0, 2, m - 3])
    );
    assert_eq!(
        session
            .call_value(
                "num_finite_field_slice_axpy_mod",
                &[Value::int(m - 1), a, b, Value::int(m)],
                &ctx
            )
            .unwrap(),
        residues(&[m - 1, 1, 7, m - 1, 1])
    );
}

#[test]
fn modular_vectors_check_shape_modulus_and_invertibility() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let strict = "num_finite_field_batch_inverse_mod";
    let skip = "num_finite_field_batch_inverse_skip_zeros";
    assert_eq!(
        session
            .call_value(strict, &[residues(&[2, 3, 4]), Value::int(7)], &ctx)
            .unwrap(),
        residues(&[4, 5, 2])
    );
    assert_eq!(
        session
            .call_value(strict, &[residues(&[]), Value::int(7)], &ctx)
            .unwrap(),
        residues(&[])
    );
    let composite = [residues(&[7, 4, 15, 18, 8]), Value::int(15)];
    assert!(matches!(
        session.call_value(strict, &composite, &ctx),
        Err(OpError::Domain(_))
    ));
    assert_eq!(
        session.call_value(skip, &composite, &ctx).unwrap(),
        residues(&[13, 4, 0, 0, 2])
    );
    assert!(matches!(
        session.call_value(skip, &[residues(&[1]), Value::int(1)], &ctx),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        session.call_value(
            "num_finite_field_dot_mod",
            &[residues(&[1]), residues(&[]), Value::int(7)],
            &ctx
        ),
        Err(OpError::Domain(_))
    ));
    let dimensioned = Value::Dyn(DynValue {
        value: Box::new(Value::int(2)),
        dim: DimExpr::base("Length"),
    });
    assert!(
        session
            .call_value(
                "num_finite_field_slice_scale_mod",
                &[residues(&[1]), dimensioned, Value::int(7)],
                &ctx
            )
            .is_err()
    );
}

#[test]
fn wide_modular_arithmetic_preserves_bigint_values_and_checks_range() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let m = u128::MAX;
    for (op, args, expected) in [
        (
            "num_number_theory_mul_mod_u128",
            [Value::int(m - 1), Value::int(m - 2), Value::int(m)],
            Value::int(2),
        ),
        (
            "num_number_theory_mod_pow_u128",
            [Value::int(m - 2), Value::int(3), Value::int(m)],
            Value::int(m - 8),
        ),
        (
            "num_number_theory_mod_pow_u128",
            [Value::int(m), Value::int(0), Value::int(1)],
            Value::int(0),
        ),
    ] {
        assert_eq!(
            session.call_value(op, &args, &ctx).unwrap(),
            expected,
            "{op}"
        );
    }
    for name in [
        "num_number_theory_mul_mod",
        "num_number_theory_mul_mod_u128",
        "num_number_theory_mod_pow_u128",
    ] {
        assert!(matches!(
            session.call_value(name, &[Value::int(1), Value::int(1), Value::int(0)], &ctx),
            Err(OpError::Domain(_))
        ));
    }
    let two_to_96 = 1_u128 << 96;
    let gcd = "num_number_theory_gcd_u128";
    assert_eq!(
        session
            .call_value(
                gcd,
                &[Value::int(two_to_96), Value::int(3 * two_to_96)],
                &ctx
            )
            .unwrap(),
        Value::int(two_to_96)
    );
    assert!(matches!(
        session.call_value(gcd, &[Value::int(-1), Value::int(1)], &ctx),
        Err(OpError::Domain(_))
    ));
    let outside = num_bigint::BigInt::from(u128::MAX) + 1_u32;
    assert!(matches!(
        session.call_value(gcd, &[Value::int(outside), Value::int(1)], &ctx),
        Err(OpError::Domain(_))
    ));
}
