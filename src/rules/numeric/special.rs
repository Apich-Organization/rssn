//! Special functions of a complex argument, as operators.
//!
//! A complex number has no term encoding, so a complex argument is passed as
//! its real and imaginary parts and a complex result comes back as
//! `list(re, im)`.
//!
//! | operator | value |
//! |---|---|
//! | `lgamma_complex(re, im)` | `ln Gamma(z)` on the branch continuous off the non-positive real axis |
//! | `gamma_complex(re, im)` | `Gamma(z)` |
//! | `digamma_complex(re, im)` | `psi(z) = Gamma'(z) / Gamma(z)` |
//! | `lambert_w_complex(re, im, k)` | branch `k` of the Lambert `W`, the solution of `w exp(w) = z` |
//! | `bessel_j_complex(nu, re, im)` | `J_nu(z)` for real order, principal branch |
//! | `bessel_i_complex(nu, re, im)` | `I_nu(z)` for real order, principal branch |
//! | `erf_complex(re, im)` | the error function `erf(z)` |
//! | `fresnel_s_complex(re, im)`, `fresnel_c_complex(re, im)` | the Fresnel integrals `S(z)` and `C(z)` in the `sin(pi t^2 / 2)` normalisation |
//!
//! The arguments must be finite numbers of moderate size; a pole, a branch
//! point whose value is infinite, or any non-finite result leaves the call as
//! a request instead of returning an infinity.

use num_complex::Complex64;

use super::bridge::Data;
use super::bridge::Entry;
use super::bridge::require;
use super::bridge::scalar;
use crate::kernel::diag::DiagKind;
use crate::kernels::fresnel;
use crate::kernels::special_complex;

/// Largest real or imaginary part of an argument.
const MAX_PART: f64 = 1.0e6;
/// Largest modulus the Bessel functions accept.
const MAX_BESSEL: f64 = 1.0e3;
/// Largest order magnitude of a Bessel function.
const MAX_ORDER: f64 = 100.0;

/// The complex number whose parts are the arguments `at` and `at + 1`.
fn point(
    a: &[Data],
    at: usize,
) -> Option<Complex64> {
    let (re, im) = (scalar(a, at)?, scalar(a, at + 1)?);
    require(
        re.abs() <= MAX_PART && im.abs() <= MAX_PART,
        Some(at),
        DiagKind::Domain,
        format_args!("the real and imaginary parts must not exceed {MAX_PART:e} in magnitude"),
    )?;
    Some(Complex64::new(re, im))
}

/// A complex answer as `list(re, im)`; a non-finite one is no answer.
fn pair(value: Complex64) -> Option<Data> {
    require(
        value.re.is_finite() && value.im.is_finite(),
        None,
        DiagKind::Domain,
        format_args!("the value is not finite (a pole, or an overflow)"),
    )?;
    Some(Data::vector(vec![value.re, value.im]))
}

fn lgamma_complex(a: &[Data]) -> Option<Data> {
    pair(special_complex::ln_gamma_complex(point(a, 0)?))
}

fn gamma_complex(a: &[Data]) -> Option<Data> {
    pair(special_complex::gamma_complex(point(a, 0)?))
}

fn digamma_complex(a: &[Data]) -> Option<Data> {
    pair(special_complex::digamma_complex(point(a, 0)?))
}

fn lambert_w_complex(a: &[Data]) -> Option<Data> {
    let z = point(a, 0)?;
    let k = scalar(a, 2)?;
    (k.fract() == 0.0 && k.abs() <= 50.0).then_some(())?;
    pair(special_complex::lambert_w_complex(z, k as i32))
}

/// Order and argument of a Bessel function.
fn bessel_arguments(a: &[Data]) -> Option<(f64, Complex64)> {
    let nu = scalar(a, 0)?;
    let z = point(a, 1)?;
    (nu.abs() <= MAX_ORDER && z.norm() <= MAX_BESSEL).then_some((nu, z))
}

fn bessel_j_complex(a: &[Data]) -> Option<Data> {
    let (nu, z) = bessel_arguments(a)?;
    pair(special_complex::bessel_j_complex(nu, z))
}

fn bessel_i_complex(a: &[Data]) -> Option<Data> {
    let (nu, z) = bessel_arguments(a)?;
    pair(special_complex::bessel_i_complex(nu, z))
}

/// Largest modulus the error function and Fresnel integrals accept; beyond
/// it `exp(-z^2)` over- or underflows and the series loses all digits.
const MAX_ERF: f64 = 1.0e3;

fn erf_argument(a: &[Data]) -> Option<Complex64> {
    point(a, 0).filter(|z| z.norm() <= MAX_ERF)
}

fn erf_complex(a: &[Data]) -> Option<Data> {
    pair(fresnel::erf_complex(erf_argument(a)?))
}

fn fresnel_s_complex(a: &[Data]) -> Option<Data> {
    let z = erf_argument(a)?;
    let ((re, im), _) = fresnel::fresnel_sc_complex(z.re, z.im);
    pair(Complex64::new(re, im))
}

fn fresnel_c_complex(a: &[Data]) -> Option<Data> {
    let z = erf_argument(a)?;
    let (_, (re, im)) = fresnel::fresnel_sc_complex(z.re, z.im);
    pair(Complex64::new(re, im))
}

pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "lgamma_complex",
        params: "ss",
        eval: lgamma_complex,
    },
    Entry {
        name: "gamma_complex",
        params: "ss",
        eval: gamma_complex,
    },
    Entry {
        name: "digamma_complex",
        params: "ss",
        eval: digamma_complex,
    },
    Entry {
        name: "lambert_w_complex",
        params: "sss",
        eval: lambert_w_complex,
    },
    Entry {
        name: "bessel_j_complex",
        params: "sss",
        eval: bessel_j_complex,
    },
    Entry {
        name: "bessel_i_complex",
        params: "sss",
        eval: bessel_i_complex,
    },
    Entry {
        name: "erf_complex",
        params: "ss",
        eval: erf_complex,
    },
    Entry {
        name: "fresnel_s_complex",
        params: "ss",
        eval: fresnel_s_complex,
    },
    Entry {
        name: "fresnel_c_complex",
        params: "ss",
        eval: fresnel_c_complex,
    },
];

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use num_complex::Complex64;

    use crate::rules::numeric::numeric;
    use crate::rules::testing::reduce_with;
    use crate::rules::testing::simplify;

    fn run(src: &str) -> Vec<f64> {
        let text = simplify(&[numeric()], src);
        text.replace(['(', ')', ','], " ")
            .split_whitespace()
            .filter(|word| *word != "list")
            .map(|word| word.parse().unwrap_or(f64::NAN))
            .collect()
    }

    fn complex(src: &str) -> Complex64 {
        let v = run(src);
        assert_eq!(v.len(), 2, "{src}: {v:?}");
        Complex64::new(v[0], v[1])
    }

    fn near(
        got: Complex64,
        want: Complex64,
        tol: f64,
    ) {
        assert!(
            (got - want).norm() <= tol * (1.0 + want.norm()),
            "{got} vs {want}"
        );
    }

    #[test]
    fn gamma_matches_closed_forms_and_the_reflection_formula() {
        near(
            complex("gamma_complex(5, 0)"),
            Complex64::new(24.0, 0.0),
            1e-12,
        );
        near(
            complex("lgamma_complex(0.5, 0)"),
            Complex64::new(PI.sqrt().ln(), 0.0),
            1e-12,
        );
        // Gamma(1 + i) = -i Gamma(i) and Gamma(i) is the tabulated value.
        near(
            complex("gamma_complex(1, 1)"),
            Complex64::new(0.498_015_668_118_356, -0.154_949_828_301_811),
            1e-12,
        );
        // Gamma(z) Gamma(1 - z) = pi / sin(pi z), a second implementation.
        let z = Complex64::new(0.3, 0.4);
        let product = complex("gamma_complex(0.3, 0.4)") * complex("gamma_complex(0.7, -0.4)");
        near(product, PI / (z * PI).sin(), 1e-11);
        // exp(lgamma) is gamma.
        near(
            complex("lgamma_complex(2.5, 1.5)").exp(),
            complex("gamma_complex(2.5, 1.5)"),
            1e-11,
        );
    }

    #[test]
    fn digamma_satisfies_its_recurrence() {
        near(
            complex("digamma_complex(1, 0)"),
            Complex64::new(-std::f64::consts::EULER_GAMMA, 0.0),
            1e-12,
        );
        let z = Complex64::new(1.3, 0.7);
        near(
            complex("digamma_complex(2.3, 0.7)"),
            complex("digamma_complex(1.3, 0.7)") + 1.0 / z,
            1e-11,
        );
    }

    #[test]
    fn lambert_w_inverts_w_exp_w() {
        near(
            complex("lambert_w_complex(1, 0, 0)"),
            Complex64::new(0.567_143_290_409_783_9, 0.0),
            1e-12,
        );
        for (re, im, k) in [
            (0.5, 2.0, 0),
            (-1.0, 0.3, 0),
            (2.0, -1.0, -1),
            (1.0, 1.0, 1),
        ] {
            let w = complex(&format!("lambert_w_complex({re}, {im}, {k})"));
            near(w * w.exp(), Complex64::new(re, im), 1e-9);
        }
    }

    #[test]
    fn bessel_functions_match_tabulated_values_and_half_orders() {
        near(
            complex("bessel_j_complex(0, 2, 0)"),
            Complex64::new(0.223_890_779_141_235_67, 0.0),
            1e-12,
        );
        near(
            complex("bessel_i_complex(0, 1, 0)"),
            Complex64::new(1.266_065_877_752_008_2, 0.0),
            1e-12,
        );
        // J_{1/2}(x) = sqrt(2 / (pi x)) sin x.
        let x = 3.7_f64;
        near(
            complex("bessel_j_complex(0.5, 3.7, 0)"),
            Complex64::new((2.0 / (PI * x)).sqrt() * x.sin(), 0.0),
            1e-11,
        );
        // I_nu(z) = i^-nu J_nu(i z); for nu = 1 that is -i J_1(i z).
        let jz = complex("bessel_j_complex(1, 0, 1.5)");
        near(
            complex("bessel_i_complex(1, 1.5, 0)"),
            Complex64::new(0.0, -1.0) * jz,
            1e-11,
        );
    }

    #[test]
    fn error_function_and_fresnel_integrals_match_known_values() {
        // erf(1) and the odd symmetry erf(-z) = -erf(z).
        near(
            complex("erf_complex(1, 0)"),
            Complex64::new(0.842_700_792_949_714_9, 0.0),
            1e-12,
        );
        near(
            complex("erf_complex(-0.5, -0.7)"),
            -complex("erf_complex(0.5, 0.7)"),
            1e-12,
        );
        // erf(i y) = i erfi(y); erfi(1) = 1.650425758797543.
        near(
            complex("erf_complex(0, 1)"),
            Complex64::new(0.0, 1.650_425_758_797_543),
            1e-12,
        );
        // On the real axis S(1) and C(1) are the tabulated values.
        near(
            complex("fresnel_s_complex(1, 0)"),
            Complex64::new(0.438_259_147_390_354_8, 0.0),
            1e-10,
        );
        near(
            complex("fresnel_c_complex(1, 0)"),
            Complex64::new(0.779_893_400_376_822_8, 0.0),
            1e-10,
        );
        // C and S are odd; both are built from erf, so a complex point agrees
        // with the closed form C(z) + i S(z) = (1 + i)/2 erf((1 - i) sqrt(pi) z / 2).
        let z = Complex64::new(0.6, 0.4);
        let c = complex("fresnel_c_complex(0.6, 0.4)");
        let s = complex("fresnel_s_complex(0.6, 0.4)");
        let a = PI.sqrt() / 2.0;
        let e = complex(&format!(
            "erf_complex({}, {})",
            (Complex64::new(a, -a) * z).re,
            (Complex64::new(a, -a) * z).im
        ));
        near(
            c + Complex64::new(0.0, 1.0) * s,
            Complex64::new(0.5, 0.5) * e,
            1e-11,
        );
        for src in [
            "erf_complex(1e4, 0)",
            "erf_complex(x, 0)",
            "fresnel_s_complex(0, 5e3)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }

    #[test]
    fn poles_symbols_and_huge_arguments_are_left_as_requests() {
        for src in [
            "gamma_complex(0, 0)",
            "gamma_complex(-3, 0)",
            "gamma_complex(a, 0)",
            "gamma_complex(1e300, 0)",
            "lambert_w_complex(1, 0, 0.5)",
            "bessel_j_complex(1e6, 1, 0)",
            "bessel_j_complex(0, 1e9, 0)",
        ] {
            let (text, reduced) = reduce_with(&[numeric()], src, &[]);
            assert!(!reduced, "{text}");
        }
    }
}
