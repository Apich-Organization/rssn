//! Coverage tests for `rssn::kernels::special_exact`.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use num_traits::Zero;
use rssn::kernels::special_exact::chebyshev_t_exact;
use rssn::kernels::special_exact::chebyshev_tu_mod;
use rssn::kernels::special_exact::chebyshev_u_exact;
use rssn::kernels::special_exact::hermite_h_exact;
use rssn::kernels::special_exact::laguerre_l_exact;
use rssn::kernels::special_exact::legendre_p_exact;

fn rat(
    n: i64,
    d: i64,
) -> BigRational {
    BigRational::new(BigInt::from(n), BigInt::from(d))
}

#[test]
fn legendre_hermite_laguerre_closed_forms_and_recurrences() {
    let x = rat(1, 3);
    // P_0(x)=1, P_1(x)=x, P_2(x)=(3x^2-1)/2, P_3(x)=(5x^3-3x)/2
    assert_eq!(legendre_p_exact(0, &x), BigRational::one());
    assert_eq!(legendre_p_exact(1, &x), x);
    assert_eq!(
        legendre_p_exact(2, &x),
        (rat(3, 1) * &x * &x - BigRational::one()) / rat(2, 1)
    );
    assert_eq!(
        legendre_p_exact(3, &x),
        (rat(5, 1) * &x * &x * &x - rat(3, 1) * &x) / rat(2, 1)
    );
    // Endpoint identities: P_n(1)=1, P_n(-1)=(-1)^n.
    for n in 0..=6 {
        assert_eq!(legendre_p_exact(n, &rat(1, 1)), rat(1, 1));
        let sign = if n % 2 == 0 { 1 } else { -1 };
        assert_eq!(legendre_p_exact(n, &rat(-1, 1)), rat(sign, 1));
    }

    // Hermite H_n(x): H_0=1, H_1=2x, H_2=4x^2-2, H_3=8x^3-12x
    let x2 = rat(-2, 5);
    assert_eq!(hermite_h_exact(0, &x2), BigRational::one());
    assert_eq!(hermite_h_exact(1, &x2), rat(2, 1) * &x2);
    assert_eq!(hermite_h_exact(2, &x2), rat(4, 1) * &x2 * &x2 - rat(2, 1));
    assert_eq!(
        hermite_h_exact(3, &x2),
        rat(8, 1) * &x2 * &x2 * &x2 - rat(12, 1) * &x2
    );

    // Laguerre L_n(0) = 1 and L_2(x) = (x^2 - 4x + 2)/2
    for n in 0..=6 {
        assert_eq!(
            laguerre_l_exact(n, &BigRational::zero()),
            BigRational::one()
        );
    }
    assert_eq!(
        laguerre_l_exact(2, &x),
        (&x * &x - rat(4, 1) * &x + rat(2, 1)) / rat(2, 1)
    );
}

#[test]
fn chebyshev_exact_and_modular_identities() {
    let points = [
        rat(0, 1),
        rat(1, 1),
        rat(-1, 1),
        rat(1, 2),
        rat(-1, 2),
        rat(2, 5),
        rat(-3, 7),
    ];
    for x in &points {
        let mut t_prev = BigRational::one();
        let mut t_curr = x.clone();
        let mut u_prev = BigRational::one();
        let mut u_curr = rat(2, 1) * x;
        assert_eq!(chebyshev_t_exact(0, x), t_prev);
        assert_eq!(chebyshev_t_exact(1, x), t_curr);
        assert_eq!(chebyshev_u_exact(0, x), u_prev);
        assert_eq!(chebyshev_u_exact(1, x), u_curr);
        for n in 2..=12 {
            let t_next = rat(2, 1) * x * &t_curr - &t_prev;
            let u_next = rat(2, 1) * x * &u_curr - &u_prev;
            assert_eq!(chebyshev_t_exact(n, x), t_next, "T_{n}({x})");
            assert_eq!(chebyshev_u_exact(n, x), u_next, "U_{n}({x})");
            t_prev = t_curr;
            t_curr = t_next;
            u_prev = u_curr;
            u_curr = u_next;
        }
    }

    assert_eq!(chebyshev_tu_mod(10, 3, 0), None);
    assert_eq!(chebyshev_tu_mod(10, 3, 1), Some((0, 0)));
    let modulus = 1_000_000_007_u64;
    for (n, x) in [(0_u64, 5_i64), (1, -3), (7, 4), (11, -7)] {
        let (tm, um) = chebyshev_tu_mod(n, x, modulus).unwrap();
        let te = chebyshev_t_exact(n as u32, &rat(x, 1));
        let ue = chebyshev_u_exact(n as u32, &rat(x, 1));
        let m = BigInt::from(modulus);
        let t_rem = ((te.numer() % &m) + &m) % &m;
        let u_rem = ((ue.numer() % &m) + &m) % &m;
        assert_eq!(BigInt::from(tm), t_rem);
        assert_eq!(BigInt::from(um), u_rem);
    }
}
