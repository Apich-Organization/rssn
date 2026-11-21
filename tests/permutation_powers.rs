//! Permutation powering through exact signed and fixed-width inputs.

use num_bigint::BigInt;
use rssn::api::Session;
use rssn::kernel::CallCtx;
use rssn::kernel::Value;
use rssn::kernel::signature::TypeExpr;
use rssn::kernel::value::TensorValue;
use rssn::kernels::permgroup::compose;
use rssn::kernels::permgroup::cycle_type;
use rssn::kernels::permgroup::element_order;
use rssn::kernels::permgroup::identity;
use rssn::kernels::permgroup::power;
use rssn::kernels::permgroup::power_exact;

const CYCLES: [u32; 9] = [1, 2, 0, 4, 5, 6, 7, 3, 8];
const ROTATED: [u32; 9] = [1, 2, 0, 5, 6, 7, 3, 4, 8];

fn wide_exponent() -> BigInt {
    -((BigInt::from(1) << 512_usize) + BigInt::from(7))
}

#[test]
fn signed_powers_reduce_each_cycle_and_fix_singletons() {
    assert_eq!(power_exact(&CYCLES, &wide_exponent()).unwrap(), ROTATED);
    assert_eq!(
        power_exact(&CYCLES, &BigInt::from(-1)).unwrap(),
        [2, 0, 1, 7, 3, 4, 5, 6, 8]
    );
    assert_eq!(power_exact(&CYCLES, &BigInt::from(0)).unwrap(), identity(9));
    assert_eq!(cycle_type(&CYCLES), [1, 3, 5]);
    assert_eq!(element_order(&CYCLES), 15_u32.into());
    for invalid in [[0, 0], [0, 2]] {
        assert!(power_exact(&invalid, &BigInt::from(0)).is_none());
        assert!(power_exact(&invalid, &wide_exponent()).is_none());
    }
    assert!(power_exact(&[], &wide_exponent()).unwrap().is_empty());
}

fn squaring(
    p: &[u32],
    mut exponent: u64,
) -> Vec<u32> {
    let mut base = p.to_vec();
    let mut out = identity(p.len());
    while exponent != 0 {
        if exponent % 2 == 1 {
            out = compose(&out, &base);
        }
        base = compose(&base, &base);
        exponent /= 2;
    }
    out
}

#[test]
fn cycle_rotations_match_composition_on_each_small_permutation() {
    fn permutations(
        p: &mut [u32],
        fixed: usize,
    ) {
        if fixed == p.len() {
            for exponent in [0, 1, 2, 3, 4, 5, 17, 255, 1 << 42, u64::MAX] {
                let expected = squaring(p, exponent);
                assert_eq!(power(p, exponent), expected);
                assert_eq!(power_exact(p, &BigInt::from(exponent)).unwrap(), expected);
            }
            return;
        }
        for next in fixed..p.len() {
            p.swap(fixed, next);
            permutations(p, fixed + 1);
            p.swap(fixed, next);
        }
    }
    for n in 0..=5 {
        permutations(&mut identity(n), 0);
    }
}

#[test]
fn native_permutation_power_accepts_exact_signed_input_and_rejects_invalid_domains() {
    let session = Session::new();
    let call = |images: &[u32], exponent| {
        session.call_value(
            "num_permgroup_power_exact",
            &[
                Value::Tensor(TensorValue {
                    shape: vec![images.len()],
                    elem: TypeExpr::Int,
                    data: images.iter().copied().map(Value::int).collect(),
                }),
                Value::Int(exponent),
            ],
            &CallCtx::new(),
        )
    };
    let Value::Tensor(result) = call(&CYCLES, wide_exponent()).unwrap() else {
        panic!("expected an exact permutation tensor");
    };
    assert_eq!(result.elem, TypeExpr::Int);
    assert_eq!(result.data, ROTATED.map(Value::int));
    assert!(call(&[0, 0], 2.into()).is_err());
    assert!(call(&CYCLES, BigInt::from(1) << 65_536_usize).is_err());
}
