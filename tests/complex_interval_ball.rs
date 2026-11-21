//! Exact rational checks for complex boxes and shared polynomial evaluation.

use num_complex::Complex64;
use num_rational::BigRational;
use rssn::backend::tape::Interval;
use rssn::kernels::ball::CTerm;
use rssn::kernels::ball::Cball;
use rssn::kernels::ball::eval_jacobian_terms;
use rssn::kernels::ball::eval_terms;
use rssn::kernels::ball::eval_value_and_jacobian_terms;
use rssn::kernels::ball::jacobian_at_terms;
use rssn::kernels::ball::values_and_jacobian_at_terms;
use rssn::kernels::ball::values_at_terms;

type Rational = BigRational;
type Gaussian = (Rational, Rational);

fn rational(x: f64) -> Rational {
    Rational::from_float(x).unwrap_or_else(|| panic!("nonfinite reference input {x}"))
}

fn gaussian(z: Complex64) -> Gaussian {
    (rational(z.re), rational(z.im))
}

fn multiply(
    a: &Gaussian,
    b: &Gaussian,
) -> Gaussian {
    (&a.0 * &b.0 - &a.1 * &b.1, &a.0 * &b.1 + &a.1 * &b.0)
}

fn power(
    z: Complex64,
    mut exponent: u32,
) -> Gaussian {
    let mut result = gaussian(Complex64::new(1.0, 0.0));
    let mut factor = gaussian(z);
    while exponent != 0 {
        if exponent & 1 != 0 {
            result = multiply(&result, &factor);
        }
        exponent >>= 1;
        if exponent != 0 {
            factor = multiply(&factor, &factor);
        }
    }
    result
}

fn assert_encloses_real(
    interval: Interval,
    value: &Rational,
) {
    assert!(
        interval.lo == f64::NEG_INFINITY || rational(interval.lo) <= *value,
        "lower bound misses the reference: {interval:?}, {value}"
    );
    assert!(
        interval.hi == f64::INFINITY || rational(interval.hi) >= *value,
        "upper bound misses the reference: {interval:?}, {value}"
    );
}

fn assert_encloses(
    ball: Cball,
    value: &Gaussian,
) {
    assert_encloses_real(ball.re, &value.0);
    assert_encloses_real(ball.im, &value.1);
}

#[test]
fn modulus_bounds_round_outward_at_normal_and_extreme_scales() {
    let least = f64::from_bits(1);
    for (re, im) in [(1.0, 1.0), (1.0, 5.0), (least, least), (f64::MAX, f64::MAX)] {
        let ball = Cball::point(re, im);
        let lower = rational(ball.modulus_lower());
        let expected = rational(re).pow(2) + rational(im).pow(2);
        assert!(&lower * &lower <= expected, "{re}+{im}i");
        let upper = ball.modulus_upper();
        if upper.is_finite() {
            let upper = rational(upper);
            assert!(&upper * &upper >= expected, "{re}+{im}i");
        }
    }

    let rectangle = Cball::rectangle((-3.0, 5.0), (1.0, 4.0));
    assert_eq!(rectangle.modulus_lower(), 1.0);
    assert!(rational(rectangle.modulus_upper()).pow(2) >= rational(41.0));
    assert_eq!(Cball::point(least, 0.0).centre().re, least);
    assert_eq!(
        Cball::rectangle((-f64::MAX, f64::MAX), (0.0, 0.0)).half_width(),
        f64::MAX
    );
}

fn quotient(
    a: Complex64,
    b: Complex64,
) -> Gaussian {
    let a = gaussian(a);
    let b = gaussian(b);
    let denominator = &b.0 * &b.0 + &b.1 * &b.1;
    let numerator = multiply(&a, &(b.0, -b.1));
    (numerator.0 / &denominator, numerator.1 / &denominator)
}

#[test]
fn division_encloses_nonzero_axis_crossings_and_extreme_rectangles() {
    let numerator = Cball::rectangle((1.0, 2.0), (-0.5, 1.0));
    let divisor = Cball::rectangle((-10.0, 10.0), (1.0, 1.0));
    let q = numerator
        .div(&divisor)
        .unwrap_or_else(|| panic!("the divisor excludes zero"));
    for a in [Complex64::new(1.0, -0.5), Complex64::new(2.0, 1.0)] {
        for b in [
            Complex64::new(-10.0, 1.0),
            Complex64::new(0.0, 1.0),
            Complex64::new(10.0, 1.0),
        ] {
            assert_encloses(q, &quotient(a, b));
        }
    }
    assert!(
        numerator
            .div(&Cball::rectangle((-1.0, 1.0), (0.0, 2.0)))
            .is_none()
    );

    for scale in [1e-200, 1e200] {
        let a = Complex64::new(scale, -2.0 * scale);
        let b = Complex64::new(-2.0 * scale, -3.0 * scale);
        let q = Cball::of(a)
            .div(&Cball::of(b))
            .unwrap_or_else(|| panic!("nonzero point divisor"));
        assert!(q.is_finite() && q.half_width() < 2e-14, "{q:?}");
        assert_encloses(q, &quotient(a, b));
    }
}

fn exact_term(
    term: &CTerm,
    point: &[Complex64],
    derivative: Option<usize>,
) -> Gaussian {
    let mut result = gaussian(term.coef);
    for (column, (&z, &degree)) in point.iter().zip(&term.exps).enumerate() {
        let factor = if derivative == Some(column) {
            if degree == 0 {
                return gaussian(Complex64::new(0.0, 0.0));
            }
            let (re, im) = power(z, degree - 1);
            let multiplier = rational(f64::from(degree));
            (re * &multiplier, im * &multiplier)
        } else {
            power(z, degree)
        };
        result = multiply(&result, &factor);
    }
    if derivative.is_some_and(|column| column >= term.exps.len()) {
        return gaussian(Complex64::new(0.0, 0.0));
    }
    result
}

fn exact_row(
    row: &[CTerm],
    point: &[Complex64],
    derivative: Option<usize>,
) -> Gaussian {
    row.iter()
        .map(|term| exact_term(term, point, derivative))
        .fold(gaussian(Complex64::new(0.0, 0.0)), |sum, term| {
            (sum.0 + term.0, sum.1 + term.1)
        })
}

#[test]
fn fused_evaluation_encloses_nondyadic_products_and_irregular_rows() {
    let system = vec![
        vec![
            CTerm {
                coef: Complex64::new(0.7, -1.2),
                exps: vec![2, 4, 1],
            },
            CTerm::real(-0.9, vec![2, 4, 0, 1]),
            CTerm::real(1.7, vec![0, 1]),
        ],
        vec![],
        vec![
            CTerm::real(0.3, vec![3, 0, 2]),
            CTerm::real(-0.2, vec![0, 0, 1]),
            CTerm::real(1.0, vec![]),
        ],
    ];
    let box_point = [
        Cball::rectangle((-0.6, 0.3), (-0.2, 0.4)),
        Cball::rectangle((0.2, 0.7), (-0.3, 0.1)),
        Cball::rectangle((-0.1, 0.2), (-0.4, 0.3)),
    ];
    let (values, jacobian) = eval_value_and_jacobian_terms(&system, &box_point);
    assert_eq!(values, eval_terms(&system, &box_point));
    assert_eq!(jacobian, eval_jacobian_terms(&system, &box_point));
    for z in [
        Complex64::new(-0.6, -0.2),
        Complex64::new(0.0, 0.0),
        Complex64::new(0.3, 0.4),
    ] {
        for w in [Complex64::new(0.2, -0.3), Complex64::new(0.7, 0.1)] {
            let point = [z, w, Complex64::new(0.0, 0.0)];
            let (point_values, point_jacobian) = values_and_jacobian_at_terms(&system, &point);
            assert_eq!(point_values, values_at_terms(&system, &point));
            assert_eq!(point_jacobian, jacobian_at_terms(&system, &point));
            for (row, equation) in system.iter().enumerate() {
                assert_encloses(values[row], &exact_row(equation, &point, None));
                for column in 0..point.len() {
                    assert_encloses(
                        jacobian[row][column],
                        &exact_row(equation, &point, Some(column)),
                    );
                }
            }
        }
    }
}

#[test]
fn full_unsigned_exponents_and_roots_do_not_require_coordinate_division() {
    let exponent = u32::MAX;
    let system = vec![
        vec![
            CTerm::real(1.0, vec![exponent, 0]),
            CTerm::real(3.0, vec![0, 1]),
        ],
        vec![CTerm::real(1.0, vec![0, exponent])],
    ];
    let point = [Complex64::new(0.0, 1.0), Complex64::new(0.0, 0.0)];
    let (value, jacobian) = values_and_jacobian_at_terms(&system, &point);
    assert_eq!(
        value,
        vec![Complex64::new(0.0, -1.0), Complex64::new(0.0, 0.0)]
    );
    assert_eq!(jacobian[0][0], Complex64::new(-f64::from(exponent), 0.0));
    assert_eq!(jacobian[0][1], Complex64::new(3.0, 0.0));
    assert_eq!(jacobian[1], vec![Complex64::new(0.0, 0.0); 2]);

    let (boxes, slopes) = eval_value_and_jacobian_terms(&system, &point.map(Cball::of));
    for (row, equation) in system.iter().enumerate() {
        assert_encloses(boxes[row], &exact_row(equation, &point, None));
        for column in 0..2 {
            assert_encloses(
                slopes[row][column],
                &exact_row(equation, &point, Some(column)),
            );
        }
    }
    assert_eq!(Cball::zero().powi(exponent), Cball::zero());
    assert_eq!(Cball::zero().powi(0), Cball::one());
}
