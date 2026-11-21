//! Exact Pauli-matrix checks for three-dimensional multivector inverses.

use num_complex::Complex;
use num_rational::BigRational;
use num_traits::ToPrimitive;
use num_traits::Zero;
use rssn::kernels::geometric_algebra::Multivector3D;

type Exact = Complex<BigRational>;

fn from_coefficients([s, v1, v2, v3, b12, b23, b31, pss]: [f64; 8]) -> Multivector3D {
    Multivector3D::new(s, v1, v2, v3, b12, b23, b31, pss)
}

fn coefficients(a: Multivector3D) -> [f64; 8] {
    [a.s, a.v1, a.v2, a.v3, a.b12, a.b23, a.b31, a.pss]
}

fn exact(
    re: f64,
    im: f64,
) -> Exact {
    Exact::new(
        BigRational::from_float(re).unwrap(),
        BigRational::from_float(im).unwrap(),
    )
}

fn pauli(a: Multivector3D) -> [Exact; 4] {
    let z = exact(a.s, a.pss);
    let w = exact(a.v3, a.b12);
    let x = exact(a.v1, a.b23);
    let iy = exact(-a.b31, a.v2);
    [&z + &w, &x - &iy, &x + &iy, z - w]
}

fn from_pauli([a, b, c, d]: [Exact; 4]) -> Option<[f64; 8]> {
    let half = exact(0.5, 0.0);
    let z = (&a + &d) * &half;
    let w = (&a - &d) * &half;
    let x = (&b + &c) * &half;
    let y = (c - b) * exact(0.0, -0.5);
    let value = |x: BigRational| x.to_f64().filter(|y| y.is_finite());
    Some([
        value(z.re)?,
        value(x.re)?,
        value(y.re)?,
        value(w.re)?,
        value(w.im)?,
        value(x.im)?,
        value(y.im)?,
        value(z.im)?,
    ])
}

fn pauli_inverse(a: Multivector3D) -> Option<[f64; 8]> {
    let [a, b, c, d] = pauli(a);
    let determinant = &a * &d - &b * &c;
    if determinant.is_zero() {
        return None;
    }
    from_pauli([d, -b, -c, a].map(|entry| entry / &determinant))
}

fn assert_inverse_ulps(a: Multivector3D) {
    let expected = pauli_inverse(a);
    let actual = a.inv().map(coefficients);
    match (actual, expected) {
        | (None, None) => {},
        | (Some(actual), Some(expected)) => {
            for (value, reference) in actual.into_iter().zip(expected) {
                assert!(
                    value == reference
                        || (value.is_sign_negative() == reference.is_sign_negative()
                            && value.to_bits().abs_diff(reference.to_bits()) <= 2),
                    "{a:?}: {value:e} != {reference:e}",
                );
            }
        },
        | pair => panic!("{a:?}: {pair:?}"),
    }
}

#[test]
fn central_product_matches_exact_pauli_scaling() {
    let a = from_coefficients([0.5, -2.0, 3.0, -4.0, 5.0, 6.0, -7.0, 8.0]);
    for (re, im) in [(1.0, 0.0), (0.0, 1.0), (-0.5, 0.75), (0.0, 0.0)] {
        let factor = exact(re, im);
        let reference = from_pauli(pauli(a).map(|value| value * &factor)).unwrap();
        assert_eq!(coefficients(a.mul_central(re, im)), reference);
        let centre = from_coefficients([re, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, im]);
        assert_eq!(a.mul_central(re, im), a * centre);
        assert_eq!(a.mul_central(re, im), centre * a);
    }
}

#[test]
fn inversion_keeps_small_terms_in_a_large_near_singular_element() {
    for a in [
        [1e200, 1e200, 0.0, 0.0, 0.0, 0.0, 0.0, 1e-200],
        [1e300, 1e300, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        [1.0, 1.0_f64.next_down(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0_f64.next_up(), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [
            2.0_f64.powi(-200),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            f64::MIN_POSITIVE,
        ],
        [
            2.0_f64.powi(200),
            2.0_f64.powi(-200),
            0.0,
            0.0,
            0.0,
            2.0_f64.powi(-200),
            0.0,
            0.0,
        ],
    ] {
        let a = from_coefficients(a);
        assert!(a.inv().is_some());
        assert_inverse_ulps(a);
    }
}

#[test]
fn exact_singularity_and_inverse_overflow_are_rejected() {
    let rounded_overflow = 2.356_029_373_680_492_3e-293_f64;
    for a in [
        [3.0, 1.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0],
        [3e200, 3e200, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0; 8],
        [f64::MIN_POSITIVE * 0.125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [
            rounded_overflow,
            rounded_overflow.next_down(),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    ] {
        let a = from_coefficients(a);
        assert!(a.inv().is_none());
        assert_inverse_ulps(a);
    }
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            from_coefficients([1.0, 0.0, 0.0, 0.0, value, 0.0, 0.0, 0.0])
                .inv()
                .is_none()
        );
    }
    for value in [f64::MIN_POSITIVE * 0.5, f64::MAX] {
        let a = from_coefficients([value, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(a.inv().is_some());
        assert_inverse_ulps(a);
    }
}

#[test]
fn normwise_accuracy_across_three_common_scales() {
    let mut state = 1009_u32;
    for _ in 0..32 {
        let base = [(); 8].map(|()| {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            f64::from(state % 19) - 9.0
        });
        for exponent in [-800, 0, 800] {
            let scale = 2.0_f64.powi(exponent);
            let a = from_coefficients(base.map(|x| x * scale));
            let expected = pauli_inverse(a);
            let actual = a.inv().map(coefficients);
            match (actual, expected) {
                | (None, None) => {},
                | (Some(actual), Some(expected)) => {
                    let maximum = expected.iter().map(|x| x.abs()).fold(0.0, f64::max);
                    for (value, reference) in actual.into_iter().zip(expected) {
                        assert!(((value - reference) / maximum).abs() < 3e-14, "{a:?}");
                    }
                },
                | pair => panic!("{a:?}: {pair:?}"),
            }
        }
    }
}
