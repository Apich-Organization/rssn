//! # Fresnel integrals
//!
//! The normalised Fresnel integrals (DLMF 7.2.7)
//!
//! ```text
//! S(x) = ∫_0^x sin(π t² / 2) dt,      C(x) = ∫_0^x cos(π t² / 2) dt,
//! ```
//!
//! odd, entire, with `S(±∞) = C(±∞) = ±1/2`, plus the auxiliary functions
//! `f` and `g` of DLMF 7.2.10, and the complex-argument versions through the
//! error function of a complex argument ([`erf_complex`], [`faddeeva_w`]).
//!
//! ## Real argument
//!
//! * `|x| <= 1.5`: the Maclaurin series in `t = π x² / 2`, which has no
//!   serious cancellation on that range;
//! * `|x| > 1.5`: the auxiliary functions from the Laplace-type continued
//!   fraction (modified Lentz, no cancellation), then
//!   `C = 1/2 + f sin θ - g cos θ`, `S = 1/2 - f cos θ - g sin θ` with
//!   `θ = π x² / 2` reduced *exactly*: `x²` is split into a double with its
//!   `fma` rounding residual, reduced modulo 4, and the remainder rotated
//!   by whole quadrants, so the phase is accurate for any `|x| <= 1e17`.
//!
//! Both are accurate to a few units in the last place away from the zeros of
//! the function (where only absolute accuracy of about `1e-16` is possible).
//!
//! ## Complex argument
//!
//! `C(z) + i S(z) = (1+i)/2 · erf(√π/2 (1-i) z)` holds for real `z` only; for
//! complex `z` the sine and cosine integrals are the combinations of
//! `erf(√π/2 (1 ± i) z)` given at [`fresnel_sc_complex`]. The complex error
//! function uses its Maclaurin series for `|z| < 1.5` and Weideman's rational
//! approximation of the Faddeeva function otherwise (about `1e-14` relative).

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;
use std::sync::OnceLock;

use num_complex::Complex64;

/// Below this the Maclaurin series is used, above it the continued fraction.
const SERIES_LIMIT: f64 = 1.5;
/// Beyond this `|x|` the correction to `±1/2` is below half an ulp.
const SATURATION: f64 = 1.0e17;

/// `sin(π u / 2)` and `cos(π u / 2)` for a reduced `u` in about `[0, 4]`.
fn quarter_turns(u: f64) -> (f64, f64) {
    let k = u.round();
    let phi = FRAC_PI_2 * (u - k);
    let (sp, cp) = phi.sin_cos();
    // k mod 4 in 0..=3 (u is non-negative here, but stay safe for tiny negatives).
    let q = k.rem_euclid(4.0) as i32;
    match q {
        | 0 => (sp, cp),
        | 1 => (cp, -sp),
        | 2 => (-sp, -cp),
        | _ => (-cp, sp),
    }
}

/// `(sin(π x² / 2), cos(π x² / 2))` with the phase reduced exactly.
#[must_use]
pub fn sin_cos_half_pi_square(x: f64) -> (f64, f64) {
    let ax = x.abs();
    if !ax.is_finite() {
        return (f64::NAN, f64::NAN);
    }
    if ax >= 1.0e154 {
        // x is a multiple of a huge power of two: x² ≡ 0 (mod 4).
        return (0.0, 1.0);
    }
    let p = ax * ax;
    let e = ax.mul_add(ax, -p);
    let u = (p % 4.0) + e;
    quarter_turns(u)
}

/// Maclaurin series; returns `(S, C)` for `x >= 0` small.
fn series(x: f64) -> (f64, f64) {
    let t = FRAC_PI_2 * x * x;
    let t2 = t * t;
    // S = x Σ (-1)^n t^(2n+1) / ((2n+1)! (4n+3))
    // C = x Σ (-1)^n t^(2n)   / ((2n)!   (4n+1))
    let mut sine_term = t;
    let mut cosine_term = 1.0;
    let mut s = t / 3.0;
    let mut c = 1.0;
    let mut n = 1.0_f64;
    for _ in 0..60 {
        sine_term *= -t2 / ((2.0 * n) * (2.0 * n + 1.0));
        cosine_term *= -t2 / ((2.0 * n - 1.0) * (2.0 * n));
        let ds = sine_term / (4.0 * n + 3.0);
        let dc = cosine_term / (4.0 * n + 1.0);
        s += ds;
        c += dc;
        if ds.abs() < 1e-18 * s.abs() && dc.abs() < 1e-18 * c.abs() {
            break;
        }
        n += 1.0;
    }
    (x * s, x * c)
}

/// `g(x) + i f(x) = x · h` with `h` the Lentz-evaluated continued fraction, `x > 0`.
fn aux_continued_fraction(x: f64) -> (f64, f64) {
    let pix2 = PI * x * x;
    let tiny = 1.0e-300;
    let mut b = Complex64::new(1.0, -pix2);
    let mut c = Complex64::new(1.0 / tiny, 0.0);
    let mut d = Complex64::new(1.0, 0.0) / b;
    let mut h = d;
    let mut n = -1.0_f64;
    for _ in 2..5000 {
        n += 2.0;
        let a = -n * (n + 1.0);
        b += 4.0;
        d = Complex64::new(1.0, 0.0) / (a * d + b);
        c = b + a / c;
        let delta = c * d;
        h *= delta;
        if (delta.re - 1.0).abs() + delta.im.abs() < 1.0e-16 {
            break;
        }
    }
    (x * h.im, x * h.re)
}

/// `(f, g)` for `x > 0` large enough for the continued fraction.
fn aux_positive(x: f64) -> (f64, f64) {
    aux_continued_fraction(x)
}

/// Fresnel integrals `(S(x), C(x))` together.
#[must_use]
pub fn fresnel_sc(x: f64) -> (f64, f64) {
    if x.is_nan() {
        return (f64::NAN, f64::NAN);
    }
    let ax = x.abs();
    let (s, c) = if ax <= SERIES_LIMIT {
        series(ax)
    } else if ax >= SATURATION {
        (0.5, 0.5)
    } else {
        let (f, g) = aux_positive(ax);
        let (sin_t, cos_t) = sin_cos_half_pi_square(ax);
        (
            0.5 - f.mul_add(cos_t, g * sin_t),
            0.5 + f.mul_add(sin_t, -(g * cos_t)),
        )
    };
    if x.is_sign_negative() {
        (-s, -c)
    } else {
        (s, c)
    }
}

/// Fresnel sine integral `S(x) = ∫_0^x sin(π t² / 2) dt`.
#[must_use]
pub fn fresnel_s(x: f64) -> f64 {
    fresnel_sc(x).0
}

/// Fresnel cosine integral `C(x) = ∫_0^x cos(π t² / 2) dt`.
#[must_use]
pub fn fresnel_c(x: f64) -> f64 {
    fresnel_sc(x).1
}

/// Auxiliary functions `(f(x), g(x))` of DLMF 7.2.10:
/// `f = (1/2 - S) cos θ - (1/2 - C) sin θ`, `g = (1/2 - C) cos θ + (1/2 - S) sin θ`,
/// `θ = π x² / 2`.
#[must_use]
pub fn fresnel_fg(x: f64) -> (f64, f64) {
    if x.is_nan() {
        return (f64::NAN, f64::NAN);
    }
    if x.is_infinite() {
        return if x > 0.0 {
            (0.0, 0.0)
        } else {
            (f64::INFINITY, f64::INFINITY)
        };
    }
    if x > SERIES_LIMIT {
        if x >= SATURATION {
            return (1.0 / (PI * x), 1.0 / (PI * PI * x * x * x));
        }
        return aux_positive(x);
    }
    let (s, c) = fresnel_sc(x);
    let (sin_t, cos_t) = sin_cos_half_pi_square(x);
    let (hs, hc) = (0.5 - s, 0.5 - c);
    (
        hs.mul_add(cos_t, -(hc * sin_t)),
        hc.mul_add(cos_t, hs * sin_t),
    )
}

/// Auxiliary function `f(x)` of DLMF 7.2.10 (`≈ 1/(π x)` for large `x`).
#[must_use]
pub fn fresnel_f(x: f64) -> f64 {
    fresnel_fg(x).0
}

/// Auxiliary function `g(x)` of DLMF 7.2.10 (`≈ 1/(π² x³)` for large `x`).
#[must_use]
pub fn fresnel_g(x: f64) -> f64 {
    fresnel_fg(x).1
}

// ----------------------------------------------------------------------
// Complex error function
// ----------------------------------------------------------------------

const WEIDEMAN_N: usize = 48;

/// Weideman's expansion coefficients (highest degree first) and scale `L`.
fn weideman() -> &'static (Vec<f64>, f64) {
    static TABLE: OnceLock<(Vec<f64>, f64)> = OnceLock::new();
    TABLE.get_or_init(|| {
        let n = WEIDEMAN_N;
        let m = 2 * n;
        let two_m = 2 * m;
        let l = (n as f64 / 2.0_f64.sqrt()).sqrt();
        // shifted[j] = F(k) with k = j for j < M, k = j - 2M for j >= M, F(-M) = 0.
        let mut shifted = vec![0.0_f64; two_m];
        for (j, slot) in shifted.iter_mut().enumerate() {
            let k = if j < m {
                j as i64
            } else {
                j as i64 - two_m as i64
            };
            if k == -(m as i64) {
                continue;
            }
            let theta = k as f64 * PI / m as f64;
            let t = l * (theta / 2.0).tan();
            *slot = (-t * t).exp() * (l * l + t * t);
        }
        let mut coefficients = Vec::with_capacity(n);
        for index in (1..=n).rev() {
            let mut sum = 0.0;
            for (j, value) in shifted.iter().enumerate() {
                let angle = PI * ((index * j) % two_m) as f64 / m as f64;
                sum += value * angle.cos();
            }
            coefficients.push(sum / two_m as f64);
        }
        (coefficients, l)
    })
}

/// Faddeeva function `w(z) = exp(-z²) erfc(-i z)` for `Im z >= 0`.
fn faddeeva_upper(z: Complex64) -> Complex64 {
    let (coefficients, l) = weideman();
    let i = Complex64::new(0.0, 1.0);
    let denominator = Complex64::new(*l, 0.0) - i * z;
    let zz = (Complex64::new(*l, 0.0) + i * z) / denominator;
    let mut p = Complex64::new(0.0, 0.0);
    for &a in coefficients {
        p = p * zz + a;
    }
    2.0 * p / (denominator * denominator) + (1.0 / PI.sqrt()) / denominator
}

/// The Faddeeva function `w(z) = exp(-z²) erfc(-i z)`, any `z` (about `1e-14`).
#[must_use]
pub fn faddeeva_w(z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if z.im >= 0.0 {
        faddeeva_upper(z)
    } else {
        // w(z) = 2 exp(-z²) - w(-z)
        2.0 * (-(z * z)).exp() - faddeeva_upper(-z)
    }
}

/// Error function of a complex argument.
///
/// Maclaurin series for `|z| < 1.5`, otherwise `1 - exp(-z²) w(i z)` after
/// folding to `Re z >= 0` with `erf(-z) = -erf(z)`. Accurate to roughly
/// `1e-14` relative to `max(1, |erf z|)`.
#[must_use]
pub fn erf_complex(z: Complex64) -> Complex64 {
    if z.re.is_nan() || z.im.is_nan() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    if z.im == 0.0 {
        return Complex64::new(crate::kernels::special::erf_numerical(z.re), z.im);
    }
    if z.re == 0.0 {
        // erf(i y) = i erfi(y)
        return Complex64::new(0.0, crate::kernels::special::erfi(z.im));
    }
    if z.norm() < 1.5 {
        let z2 = z * z;
        let mut term = z;
        let mut sum = z;
        for n in 1..80 {
            let nf = f64::from(n);
            term = -term * z2 / nf;
            let delta = term / (2.0 * nf + 1.0);
            sum += delta;
            if delta.norm() < 1e-18 * sum.norm() {
                break;
            }
        }
        return sum * (2.0 / PI.sqrt());
    }
    if z.re < 0.0 {
        return -erf_complex(-z);
    }
    let i = Complex64::new(0.0, 1.0);
    Complex64::new(1.0, 0.0) - (-(z * z)).exp() * faddeeva_upper(i * z)
}

/// `(S(z), C(z))` for a complex argument `z = re + i·im`, returned as
/// `((S.re, S.im), (C.re, C.im))`.
///
/// With `a = √π/2`, `E₁ = erf(a (1-i) z)`, `E₂ = erf(a (1+i) z)`:
/// `C(z) = (1+i)/4 · E₁ + (1-i)/4 · E₂`,
/// `S(z) = (1-i)/4 · E₁ + (1+i)/4 · E₂`.
/// Real arguments take the real kernels, so they keep full accuracy.
#[must_use]
pub fn fresnel_sc_complex(
    re: f64,
    im: f64,
) -> ((f64, f64), (f64, f64)) {
    if im == 0.0 {
        let (s, c) = fresnel_sc(re);
        return ((s, 0.0), (c, 0.0));
    }
    if re == 0.0 {
        // S(i y) = -i S(y) and C(i y) = i C(y) by oddness and S(iy) = i^3 S(y)... (see tests)
        let (s, c) = fresnel_sc(im);
        return ((0.0, -s), (0.0, c));
    }
    let z = Complex64::new(re, im);
    let a = PI.sqrt() / 2.0;
    let e1 = erf_complex(Complex64::new(a, -a) * z);
    let e2 = erf_complex(Complex64::new(a, a) * z);
    let quarter = Complex64::new(0.25, 0.25);
    let quarter_conj = Complex64::new(0.25, -0.25);
    let c = quarter * e1 + quarter_conj * e2;
    let s = quarter_conj * e1 + quarter * e2;
    ((s.re, s.im), (c.re, c.im))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[track_caller]
    fn rel(
        got: f64,
        want: f64,
        tolerance: f64,
    ) {
        assert!(
            (got - want).abs() <= tolerance * want.abs().max(1e-300),
            "got {got:e}, want {want:e}"
        );
    }

    const REAL: [(f64, f64, f64); 16] = [
        (0.1, 0.00052358954761221068668, 0.099997532627085073601),
        (0.5, 0.064732432859999277611, 0.49234422587144639288),
        (1.0, 0.43825914739035476608, 0.77989340037682282947),
        (1.5, 0.69750496008209301308, 0.44526117603982153506),
        (1.6, 0.63888768350938083462, 0.36546168344048765296),
        (2.0, 0.3434156783636982422, 0.4882534060753407545),
        (2.5, 0.61918175581959293611, 0.45741300964177704525),
        (3.5, 0.4152480119724375239, 0.53257243502800084533),
        (5.0, 0.49919138191711688675, 0.5636311887040122311),
        (7.25, 0.47194472715389316354, 0.53376656924787114232),
        (10.0, 0.4681699785848822404, 0.49989869420551572361),
        (37.3, 0.49624741109995594099, 0.4923355750873330416),
        (100.0, 0.49681690114783755327, 0.49999989867881789756),
        (1e5, 0.49999681690113816209, 0.49999999999999989868),
        (123456.789, 0.50000246370825017869, 0.49999923985222636448),
        (1e8, 0.49999999681690113816, 0.5),
    ];
    const CPLX: [(f64, f64, f64, f64, f64, f64); 8] = [
        (
            1.0,
            1.0,
            -2.06188821919484047,
            2.06188821919484047,
            2.55579377810243902,
            2.55579377810243902,
        ),
        (
            0.3,
            0.4,
            -0.0619690054654396111,
            0.0228870671898499883,
            0.300558613614952577,
            0.407737046485807309,
        ),
        (
            2.0,
            -1.0,
            -15.5877511044045873,
            36.7254648839914384,
            -36.2256879928816502,
            -16.0878713741254804,
        ),
        (
            -1.5,
            0.5,
            -1.59589959661936143,
            -0.216460933721620534,
            -0.290344529899153447,
            -1.07958259611780916,
        ),
        (
            3.0,
            3.0,
            -51417357578.0630208,
            51417357578.0630208,
            51417357578.5630208,
            51417357578.5630208,
        ),
        (
            0.5,
            2.0,
            0.400910040561988269,
            1.28573215843460146,
            1.78079635294636956,
            0.0949915944867347364,
        ),
        (
            5.0,
            0.1,
            0.498389031514127687,
            0.146500224666963458,
            0.659704145460982847,
            0.000971856588047086373,
        ),
        (
            0.001,
            0.002,
            -5.7595865315839641e-9,
            -1.04719755122225181e-9,
            0.000999999999989883676,
            0.00200000000000937617,
        ),
    ];

    const AUX: [(f64, f64, f64); 6] = [
        (0.5, 0.39920505852570224, 0.173642699613237748),
        (1.0, 0.279893400376822829, 0.0617408526096452339),
        (2.0, 0.156584321636301758, 0.0117465939246592455),
        (3.5, 0.090765558315310835, 0.00234017563172892039),
        (10.0, 0.0318300214151177596, 0.000101305794484276386),
        (100.0, 0.00318309885216244673, 1.01321182102440532e-7),
    ];

    #[test]
    fn fresnel_real_reference_values() {
        for &(x, s, c) in &REAL {
            rel(fresnel_s(x), s, 2e-15);
            rel(fresnel_c(x), c, 2e-15);
            rel(fresnel_s(-x), -s, 2e-15);
            rel(fresnel_c(-x), -c, 2e-15);
        }
    }

    #[test]
    fn fresnel_limits_and_specials() {
        assert_eq!(fresnel_sc(0.0), (0.0, 0.0));
        assert_eq!(fresnel_s(f64::INFINITY), 0.5);
        assert_eq!(fresnel_c(f64::NEG_INFINITY), -0.5);
        assert!(fresnel_s(f64::NAN).is_nan() && fresnel_c(f64::NAN).is_nan());
        assert_eq!(fresnel_c(1e300), 0.5);
        assert_eq!(fresnel_s(-1e20), -0.5);
        assert!(fresnel_s(-0.0).is_sign_negative());
        // Continuity across the series / continued-fraction switch.
        for x in [1.4999999, 1.5000001] {
            let (s, c) = fresnel_sc(x);
            let (s0, c0) = fresnel_sc(1.5);
            assert!((s - s0).abs() < 1e-6 && (c - c0).abs() < 1e-6);
        }
    }

    #[test]
    fn fresnel_auxiliary_functions() {
        for &(x, f, g) in &AUX {
            rel(fresnel_f(x), f, 1e-13);
            rel(fresnel_g(x), g, 1e-13);
        }
        // Large x: f ~ 1/(pi x).
        rel(fresnel_f(1e6), 1.0 / (PI * 1e6), 1e-10);
    }

    #[test]
    fn fresnel_matches_derivative_by_quadrature() {
        // S(b) - S(a) by composite Gauss-Legendre of sin(pi t^2 / 2).
        let nodes = [
            -0.906_179_845_938_664,
            -0.538_469_310_105_683,
            0.0,
            0.538_469_310_105_683,
            0.906_179_845_938_664,
        ];
        let weights = [
            0.236_926_885_056_189,
            0.478_628_670_499_366,
            0.568_888_888_888_889,
            0.478_628_670_499_366,
            0.236_926_885_056_189,
        ];
        let (a, b, steps) = (0.3_f64, 4.1_f64, 400);
        let h = (b - a) / f64::from(steps);
        let mut sum = 0.0;
        for k in 0..steps {
            let mid = a + h * (f64::from(k) + 0.5);
            for (n, w) in nodes.iter().zip(weights) {
                let t = mid + 0.5 * h * n;
                sum += 0.5 * h * w * (FRAC_PI_2 * t * t).sin();
            }
        }
        assert!((fresnel_s(b) - fresnel_s(a) - sum).abs() < 1e-11);
    }

    #[test]
    fn complex_erf_reference() {
        let v = erf_complex(Complex64::new(1.0, 1.0));
        assert!((v.re - 1.316_151_281_697_947_6).abs() < 1e-13);
        assert!((v.im - 0.190_453_469_237_834_7).abs() < 1e-13);
        let w = faddeeva_w(Complex64::new(0.0, 0.0));
        assert!((w.re - 1.0).abs() < 1e-14 && w.im.abs() < 1e-14);
        // erf(-z) = -erf(z), erf(conj z) = conj erf(z)
        let z = Complex64::new(0.7, -2.2);
        let (p, q) = (erf_complex(z), erf_complex(-z));
        assert!((p + q).norm() < 1e-13);
        assert!((erf_complex(z.conj()) - p.conj()).norm() < 1e-13);
    }

    #[test]
    fn fresnel_complex_reference_values() {
        for &(re, im, sr, si, cr, ci) in &CPLX {
            let ((a, b), (c, d)) = fresnel_sc_complex(re, im);
            let scale = |x: f64, y: f64| x.hypot(y);
            let (s_want, c_want) = (scale(sr, si), scale(cr, ci));
            assert!(
                (a - sr).hypot(b - si) <= 1e-12 * s_want.max(1e-3),
                "S({re},{im}) = {a}+{b}i vs {sr}+{si}i"
            );
            assert!(
                (c - cr).hypot(d - ci) <= 1e-12 * c_want.max(1e-3),
                "C({re},{im}) = {c}+{d}i vs {cr}+{ci}i"
            );
        }
    }

    #[test]
    fn fresnel_imaginary_axis_symmetry() {
        for y in [0.2, 1.0, 2.5] {
            let ((sr, si), (cr, ci)) = fresnel_sc_complex(0.0, y);
            assert_eq!((sr, si), (0.0, -fresnel_s(y)));
            assert_eq!((cr, ci), (0.0, fresnel_c(y)));
            // The generic complex path agrees up to a tiny real part.
            let ((sr2, si2), (cr2, ci2)) = fresnel_sc_complex(1e-300, y);
            assert!(sr2.abs() < 1e-12 && (si2 + fresnel_s(y)).abs() < 1e-12);
            assert!(cr2.abs() < 1e-12 && (ci2 - fresnel_c(y)).abs() < 1e-12);
        }
    }
}
