//! Exactness and domain checks for prime groups and NTT convolution.

use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::OpError;
use rssn::kernel::Value;
use rssn::kernel::signature::DimExpr;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::DynValue;
use rssn::kernel::value::TensorValue;
use rssn::kernels::finite_field;
use rssn::kernels::modular;
use rssn::kernels::number_theory;

const GOLDILOCKS: u64 = 18_446_744_069_414_584_321;

fn integer_vector(values: &[u64]) -> Value {
    Value::Tensor(TensorValue {
        shape: vec![values.len()],
        data: values.iter().copied().map(Value::int).collect(),
        elem: TypeExpr::Int,
    })
}

fn schoolbook(
    a: &[u64],
    b: &[u64],
    modulus: u64,
) -> Vec<u64> {
    let mut product = vec![0; a.len() + b.len() - 1];
    for (i, &left) in a.iter().enumerate() {
        for (j, &right) in b.iter().enumerate() {
            let m = u128::from(modulus);
            let term = u128::from(left) * u128::from(right) % m;
            product[i + j] = ((u128::from(product[i + j]) + term) % m) as u64;
        }
    }
    product
}

#[test]
fn prime_generators_agree_with_general_roots_and_ntt_contexts() {
    // p - 1 = 2^32 * 3 * 5 * 17 * 257 * 65537 for the Goldilocks prime.
    for candidate in 2..=7 {
        let has_full_order = [2, 3, 5, 17, 257, 65_537].into_iter().all(|factor| {
            modular::Barrett::new(GOLDILOCKS)
                .unwrap()
                .pow(candidate, (GOLDILOCKS - 1) / factor)
                != 1
        });
        assert_eq!(has_full_order, candidate == 7);
    }
    assert_eq!(finite_field::primitive_root_prime(GOLDILOCKS), Some(7));
    assert_eq!(number_theory::primitive_root(GOLDILOCKS), Some(7));
    assert_eq!(finite_field::primitive_root_prime(2), Some(1));
    assert_eq!(finite_field::primitive_root_prime(9), None);
    assert_eq!(number_theory::primitive_root(9), Some(2));
    assert_eq!(number_theory::primitive_root(18), Some(5));
    assert_eq!(number_theory::primitive_root(15), None);

    for len in [1, 8, 64, 256] {
        let plan = finite_field::ntt_plan(GOLDILOCKS, len).unwrap();
        if len > 2 {
            let omega = plan.mont.from_mont(plan.roots_mont[1]);
            assert_eq!(plan.mont.pow(omega, len as u64), 1);
            assert_ne!(plan.mont.pow(omega, len as u64 / 2), 1);
        }
        let again = finite_field::ntt_plan(GOLDILOCKS, len).unwrap();
        assert!(std::sync::Arc::ptr_eq(&plan, &again));
    }
}

#[test]
fn ntt_matches_schoolbook_across_transform_lengths_and_full_word_coefficients() {
    for (left_len, right_len) in [(1, 1), (5, 9), (257, 65), (8_193, 3)] {
        let a: Vec<_> = (0..left_len).map(|i| u64::MAX - 37 * i).collect();
        let b: Vec<_> = (0..right_len).map(|i| GOLDILOCKS - 1 - 17 * i).collect();
        let expected = schoolbook(&a, &b, GOLDILOCKS);
        assert_eq!(
            finite_field::poly_mul_ntt_mod(&a, &b, GOLDILOCKS),
            Some(expected)
        );
    }
    assert_eq!(
        finite_field::poly_mul_ntt_mod(&[1, 2, 3], &[4, 5, 6], 7),
        None
    );
}

#[test]
fn euclidean_coefficients_and_one_coefficient_inverses_remain_exact() {
    for (a, b) in [
        (0, 0),
        (7, 0),
        (0, 7),
        (u64::MAX, u64::MAX - 2),
        (GOLDILOCKS, 2),
    ] {
        let (gcd, left, right) = modular::extended_gcd_u64(a, b);
        assert_eq!(gcd, number_theory::gcd(a, b));
        assert_eq!(
            i128::from(a) * left + i128::from(b) * right,
            i128::from(gcd)
        );
    }
    for (value, modulus) in [
        (3, 25),
        (4, 15),
        (u64::MAX, 2),
        (u64::MAX - 1, u64::MAX),
        (2, GOLDILOCKS),
    ] {
        let inverse = modular::mod_inv_u64(value, modulus).unwrap();
        assert!(inverse < modulus);
        assert_eq!(
            u128::from(value) * u128::from(inverse) % u128::from(modulus),
            1
        );
    }
    assert_eq!(modular::mod_inv_u64(6, 15), None);
    assert_eq!(modular::mod_inv_u64(1, 1), None);
}

#[test]
fn exact_value_operations_validate_domains_without_erasing_zero() {
    let session = Session::new();
    let ctx = CallCtx::new();
    let jacobi = "num_number_theory_jacobi_symbol";
    let legendre = "num_number_theory_legendre_symbol";
    let mul = "num_finite_field_poly_mul_ntt_mod";

    assert_eq!(
        session
            .call_value(jacobi, &[Value::int(-2), Value::int(7)], &ctx)
            .unwrap(),
        Value::int(-1)
    );
    assert_eq!(
        session
            .call_value(legendre, &[Value::int(u64::MAX), Value::int(3)], &ctx)
            .unwrap(),
        Value::int(0)
    );
    assert!(matches!(
        session.call_value(legendre, &[Value::int(3), Value::int(9)], &ctx),
        Err(OpError::Domain(_))
    ));
    assert!(matches!(
        session.call_value(jacobi, &[Value::int(2), Value::int(0)], &ctx),
        Err(OpError::Domain(_))
    ));
    assert_eq!(
        session
            .call_value(
                mul,
                &[
                    integer_vector(&[GOLDILOCKS - 1, 2, 3]),
                    integer_vector(&[4, 5]),
                    Value::int(GOLDILOCKS)
                ],
                &ctx
            )
            .unwrap(),
        integer_vector(&[GOLDILOCKS - 4, 3, 22, 15])
    );
    for args in [
        [
            integer_vector(&[1, 2, 3]),
            integer_vector(&[1, 2, 3]),
            Value::int(7),
        ],
        [integer_vector(&[]), integer_vector(&[1]), Value::int(9)],
        [
            integer_vector(&[1]),
            integer_vector(&[1]),
            Value::Dyn(DynValue {
                value: Box::new(Value::int(17)),
                dim: DimExpr::base("Length"),
            }),
        ],
    ] {
        assert!(session.call_value(mul, &args, &ctx).is_err());
    }
}
