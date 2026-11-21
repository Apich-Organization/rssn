//! Exact rational and modular Chebyshev values against independent references.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use rssn::kernels::special_exact::chebyshev_t_exact;
use rssn::kernels::special_exact::chebyshev_tu_mod;
use rssn::kernels::special_exact::chebyshev_u_exact;

fn q(
    n: i64,
    d: i64,
) -> BigRational {
    BigRational::new(n.into(), d.into())
}

fn ordinary(
    n: u32,
    x: &BigRational,
    second_kind: bool,
) -> BigRational {
    if n == 0 {
        return q(1, 1);
    }
    let two_x = x * q(2, 1);
    let mut previous = q(1, 1);
    let mut current = if second_kind {
        two_x.clone()
    } else {
        x.clone()
    };
    for _ in 1..n {
        (previous, current) = (current.clone(), &two_x * &current - previous);
    }
    current
}

#[test]
fn exact_values_match_the_scalar_recurrences() {
    for x in [
        q(3, 7),
        q(-10, 3),
        q(4, 5),
        q(-1, 2),
        q(0, 1),
        q(1, 2),
        q(1, 1),
        q(-1, 1),
        BigRational::new_raw(10.into(), (-14).into()),
    ] {
        for n in [0, 1, 2, 3, 4, 6, 7, 11, 16, 33, 65] {
            assert_eq!(chebyshev_t_exact(n, &x), ordinary(n, &x, false), "{n}, {x}");
            assert_eq!(chebyshev_u_exact(n, &x), ordinary(n, &x, true), "{n}, {x}");
        }
    }
}

#[test]
fn bounded_periodic_cases_support_the_maximum_degree() {
    let n = u32::MAX;
    assert_eq!(chebyshev_t_exact(n, &q(0, 1)), q(0, 1));
    assert_eq!(chebyshev_u_exact(n, &q(0, 1)), q(0, 1));
    assert_eq!(chebyshev_t_exact(n, &q(1, 1)), q(1, 1));
    assert_eq!(chebyshev_t_exact(n, &q(-1, 1)), q(-1, 1));
    assert_eq!(chebyshev_u_exact(n, &q(1, 1)), q(4_294_967_296, 1));
    assert_eq!(chebyshev_u_exact(n, &q(-1, 1)), q(-4_294_967_296, 1));
    assert_eq!(chebyshev_t_exact(n, &q(1, 2)), q(-1, 1));
    assert_eq!(chebyshev_t_exact(n, &q(-1, 2)), q(1, 1));
    assert_eq!(chebyshev_u_exact(n, &q(1, 2)), q(-1, 1));
    assert_eq!(chebyshev_u_exact(n, &q(-1, 2)), q(1, 1));
    assert_eq!(chebyshev_t_exact(n - 1, &q(1, 2)), q(-1, 2));
    assert_eq!(chebyshev_u_exact(n - 1, &q(1, 2)), q(0, 1));
}

type Matrix = [[BigInt; 2]; 2];

fn matrix_product(
    a: &Matrix,
    b: &Matrix,
    modulus: &BigInt,
) -> Matrix {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| (&a[i][0] * &b[0][j] + &a[i][1] * &b[1][j]).mod_floor(modulus))
    })
}

fn modular_oracle(
    mut n: u64,
    x: i64,
    modulus: u64,
) -> (u64, u64) {
    let modulus = BigInt::from(modulus);
    let x = BigInt::from(x);
    let mut matrix = [[2_i64 * &x, (-1).into()], [1.into(), 0.into()]];
    let mut power: Matrix = [[1.into(), 0.into()], [0.into(), 1.into()]];
    while n != 0 {
        if n % 2 == 1 {
            power = matrix_product(&power, &matrix, &modulus);
        }
        matrix = matrix_product(&matrix, &matrix, &modulus);
        n >>= 1;
    }
    // Transition the respective initial states (T0,T-1) = (1,x)
    // and (U0,U-1) = (1,0) rather than forming either full polynomial.
    let t = (&power[0][0] + &power[0][1] * x).mod_floor(&modulus);
    let u = power[0][0].mod_floor(&modulus);
    (t.to_u64().unwrap(), u.to_u64().unwrap())
}

#[test]
fn modular_evaluation_handles_full_width_even_and_composite_moduli() {
    assert_eq!(chebyshev_tu_mod(30, -5, 0), None);
    for modulus in [1, 2, 8, 35, 97, 1_u64 << 63, u64::MAX - 1, u64::MAX] {
        for n in [0, 1, 2, 17, 63, 1_u64 << 63, u64::MAX] {
            for x in [0, -1, 6, i64::MIN, i64::MAX] {
                assert_eq!(
                    chebyshev_tu_mod(n, x, modulus),
                    Some(modular_oracle(n, x, modulus)),
                    "degree={n}, x={x}, modulus={modulus}"
                );
            }
        }
    }
}
