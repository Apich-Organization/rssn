//! Exact generalized congruences and arithmetic across the word-size boundary.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::One;
use rssn::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::TypeExpr;
use rssn::kernel::Value;
use rssn::kernel::signature::NatExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernels::number_theory::primes_sieve;
use rssn::rules::number_theory::number_theory;

fn integers(values: impl IntoIterator<Item = impl Into<BigInt>>) -> Value {
    let data: Vec<_> = values
        .into_iter()
        .map(|value| Value::Int(value.into()))
        .collect();
    Value::Tensor(TensorValue {
        shape: vec![data.len()],
        elem: TypeExpr::Int,
        data,
    })
}

fn crt_arguments(
    remainders: impl IntoIterator<Item = impl Into<BigInt>>,
    moduli: impl IntoIterator<Item = impl Into<BigInt>>,
) -> [Value; 2] {
    [integers(remainders), integers(moduli)]
}

#[test]
fn generalized_crt_matches_brute_congruences() {
    let session = Session::with_rules(&[number_theory()]).unwrap();
    let ctx = CallCtx::new();
    for (left, right) in [(1_i64, 1_i64), (2, 4), (3, 5), (4, 6), (6, 8), (7, 7)] {
        let period = (1..=left * right)
            .find(|value| value % left == 0 && value % right == 0)
            .unwrap();
        for first in [-9_i64, -1, 0, 4] {
            for second in [-4_i64, 0, 5] {
                let expected = (0..period).find(|value| {
                    (value - first).rem_euclid(left) == 0 && (value - second).rem_euclid(right) == 0
                });
                let args = crt_arguments([first, second], [left, right]);
                let result = session.apply_diag("crt_general", &args, &ctx);
                match expected {
                    | Some(expected) => assert_eq!(result.unwrap(), integers([expected, period])),
                    | None => assert!(result.is_err(), "{first} {second} {left} {right}"),
                }
            }
        }
    }
}

#[test]
fn generalized_crt_exports_exact_integer_vectors_and_handles_large_shared_factors() {
    let session = Session::with_rules(&[number_theory()]).unwrap();
    let ctx = CallCtx::new();
    let overloads = session.overloads("crt_general");
    assert_eq!(overloads.len(), 1);
    let (key, signature) = &overloads[0];
    assert_eq!(
        signature.ret,
        TypeExpr::Tensor {
            rank: NatExpr::Lit(1),
            shape: Some(vec![NatExpr::Lit(2)]),
            elem: Box::new(TypeExpr::Int),
        }
    );

    let scale = BigInt::one() << 200_u32;
    let period = &scale * 15_u32;
    let moduli = [scale.clone(), &scale * 3, &scale * 5];
    let remainders = [BigInt::from(-7), &moduli[1] - 7, -&moduli[2] - 7];
    for order in [[0, 1, 2], [2, 0, 1], [1, 2, 0]] {
        let args = crt_arguments(
            order.map(|index| remainders[index].clone()),
            order.map(|index| moduli[index].clone()),
        );
        assert_eq!(
            session.resolve_overload("crt_general", &args).unwrap(),
            *key
        );
        assert_eq!(
            session.call_key(key, &args, &ctx).unwrap(),
            integers([&period - 7, period.clone()])
        );
    }
    for args in [
        crt_arguments([0_i64], [0_i64]),
        crt_arguments([0_i64], [-7_i64]),
        crt_arguments([1_i64], [6_i64, 8]),
        crt_arguments([0_i64; 0], [0_i64; 0]),
    ] {
        assert!(session.apply_diag("crt_general", &args, &ctx).is_err());
    }

    assert_eq!(
        session
            .call_value("crt", &crt_arguments([-1_i64, 2], [5_i64, 7]), &ctx)
            .unwrap(),
        Value::Int(BigInt::from(9))
    );
    assert!(matches!(
        session.call_value("crt", &crt_arguments([5_i64, 5], [6_i64, 8]), &ctx),
        Ok(Value::Term(_))
    ));
}

#[test]
fn exact_signed_bezout_coefficients_and_long_two_adic_factors() {
    let session = Session::with_rules(&[number_theory()]).unwrap();
    let ctx = CallCtx::new();
    let mut fibonacci = (BigInt::from(0), BigInt::one());
    for _ in 0..650 {
        fibonacci = (fibonacci.1.clone(), fibonacci.0 + fibonacci.1);
    }
    for (a, b) in [
        (BigInt::from(i64::MIN), BigInt::from(0)),
        (BigInt::from(0), BigInt::from(0)),
        (BigInt::from(-33), BigInt::from(99)),
        (-fibonacci.0, fibonacci.1),
    ] {
        let result = session
            .call_value(
                "egcd",
                &[Value::Int(a.clone()), Value::Int(b.clone())],
                &ctx,
            )
            .unwrap();
        let Value::Tensor(result) = result else {
            panic!("expected three exact Bezout integers: {result:?}");
        };
        let [Value::Int(g), Value::Int(s), Value::Int(t)] = result.data.as_slice() else {
            panic!("expected three exact Bezout integers: {result:?}");
        };
        assert_eq!(*g, a.gcd(&b));
        assert_eq!(a * s + b * t, *g);
    }

    let n = (BigInt::one() << 65_539_u32) + BigInt::from(5);
    for (exponent, expected) in [(65_536_u32, 1), (65_537, -1)] {
        let args = [Value::Int(BigInt::one() << exponent), Value::Int(n.clone())];
        assert_eq!(
            session.call_value("jacobi", &args, &ctx).unwrap(),
            Value::Int(BigInt::from(expected))
        );
    }
}

#[test]
fn coprime_product_tree_keeps_the_canonical_result_across_worker_counts() {
    let primes = primes_sieve(5_000);
    for (count, power) in [(513, 1), (73, 12)] {
        let mut moduli: Vec<BigInt> = primes
            .iter()
            .take(count)
            .map(|&prime| BigInt::from(prime).pow(power))
            .collect();
        moduli.splice(7..7, [BigInt::one(), BigInt::one()]);
        let period: BigInt = moduli.iter().product();
        let first = &period / 7_u32 + 5_u32;
        let remainders: Vec<BigInt> = moduli.iter().map(|m| (&first % m) - m).collect();
        let args = crt_arguments(remainders, moduli.clone());
        for workers in [1, 4] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(workers)
                .build()
                .unwrap();
            pool.install(|| {
                let session = Session::with_rules(&[number_theory()]).unwrap();
                assert_eq!(
                    session.call_value("crt", &args, &CallCtx::new()).unwrap(),
                    Value::Int(first.clone())
                );
            });
        }

        let last = moduli.len() - 1;
        moduli[last] = moduli[0].clone();
        let args = crt_arguments(vec![0_i64; moduli.len()], moduli);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        pool.install(|| {
            let session = Session::with_rules(&[number_theory()]).unwrap();
            assert!(matches!(
                session.call_value("crt", &args, &CallCtx::new()),
                Ok(Value::Term(_))
            ));
        });
    }
}
