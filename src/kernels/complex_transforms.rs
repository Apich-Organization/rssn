//! Numerical methods for complex integral transforms.
//!
//! Includes Mellin and inverse Mellin transforms, Talbot, Weideman and Gaver–Stehfest Laplace
//! inversion, Mellin–Barnes integration, the Meijer G-function, and Hilbert, Hankel, and Borel
//! transforms.

use crate::kernels::integrate::gauss_kronrod_any;
use crate::kernels::special_complex;
use crate::kernels::special_complex::gamma_complex;
use num_complex::Complex;
use std::f64::consts::PI;

fn integrate_complex_ab<F>(
    f: F,
    a: f64,
    b: f64,
    tol: f64,
    max_evals: usize,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    let re = gauss_kronrod_any(|t| f(t).re, a, b, tol, max_evals).value;
    let im = gauss_kronrod_any(|t| f(t).im, a, b, tol, max_evals).value;
    Complex::new(re, im)
}

/// Numerical Mellin transform `\mathcal{M}[f](s) = \int_0^\infty x^{s-1} f(x) \, dx`
/// via the substitution `x = e^u` (`u \in [-L, L]`).
pub fn mellin_transform<F>(
    f: F,
    s: Complex<f64>,
    log_bound: f64,
    tol: f64,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    integrate_complex_ab(
        |u| {
            let x = u.exp();
            (s * u).exp() * f(x)
        },
        -log_bound,
        log_bound,
        tol,
        400,
    )
}

/// Numerical inverse Mellin transform
/// $f(x) = \frac{1}{2\pi i} \int_{c - i T}^{c + i T} x^{-s} F(s) \, ds$.
pub fn inverse_mellin_transform<F>(
    big_f: F,
    x: Complex<f64>,
    c: f64,
    t_max: f64,
    tol: f64,
) -> Complex<f64>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let ln_x = x.ln();
    let val = integrate_complex_ab(
        |t| {
            let s = Complex::new(c, t);
            (-s * ln_x).exp() * big_f(s)
        },
        -t_max,
        t_max,
        tol,
        400,
    );
    val / (2.0 * PI)
}

/// Numerical Bromwich Laplace inversion using Weideman's optimal Talbot contour
/// shifted by abscissa `sigma` (so singularities with `Re(s) < sigma` are enclosed).
///
/// # Errors
///
/// Returns an error if time is nonpositive, the shift is non-finite, or there are fewer than four quadrature points.
pub fn bromwich_talbot_shifted<F>(
    big_f: F,
    t: f64,
    sigma: f64,
    n_pts: usize,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    if t <= 0.0 || n_pts < 4 || !sigma.is_finite() {
        return Err("t must be positive, sigma finite, and n_pts >= 4".to_string());
    }
    let nf = n_pts as f64;
    let a = 0.6407_f64;
    let dtheta = 2.0 * PI / nf;
    let mut sum = Complex::new(0.0, 0.0);
    for k in 0..n_pts {
        let theta = -PI + ((k as f64) + 0.5) * dtheta;
        let (re_part, dre_part) = if theta.abs() < 1e-12 {
            (-0.6122 + 0.5017 / a, 0.0)
        } else {
            let at = a * theta;
            let cot = 1.0 / at.tan();
            let csc = 1.0 / at.sin();
            (
                -0.6122 + 0.5017 * theta * cot,
                0.5017 * (cot - at * csc * csc),
            )
        };
        let s = Complex::new(sigma, 0.0) + Complex::new(re_part, 0.2645 * theta) * (nf / t);
        let ds_dtheta = Complex::new(dre_part, 0.2645) * (nf / t);
        sum += (s * t).exp() * big_f(s) * ds_dtheta;
    }
    Ok(sum * dtheta / Complex::new(0.0, 2.0 * PI))
}

/// Numerical Bromwich Laplace inversion using Weideman's optimal Talbot contour:
/// $f(t) = \frac{1}{2\pi i} \int_{\Gamma} e^{s t} F(s) \, ds$ for $t > 0$.
///
/// # Errors
///
/// Returns an error if time is nonpositive or there are fewer than four quadrature points.
pub fn bromwich_talbot<F>(
    big_f: F,
    t: f64,
    n_pts: usize,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    bromwich_talbot_shifted(big_f, t, 0.0, n_pts)
}

/// Numerical Bromwich Laplace inversion using Weideman's parabolic contour
/// shifted by abscissa `sigma`: $s(u) = \sigma + \mu (i u + 1)^2$.
///
/// # Errors
///
/// Returns an error if time is nonpositive, the shift is non-finite, or there are fewer than four quadrature points.
pub fn bromwich_weideman_parabola_shifted<F>(
    big_f: F,
    t: f64,
    sigma: f64,
    n_pts: usize,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    if t <= 0.0 || n_pts < 4 || !sigma.is_finite() {
        return Err("t must be positive, sigma finite, and n_pts >= 4".to_string());
    }
    let half_n = n_pts / 2;
    let half_f = half_n as f64;
    let mu = (PI * half_f) / (12.0 * t);
    let h = 3.0 / half_f;
    let mut sum = Complex::new(0.0, 0.0);
    for k in -(half_n as isize)..=(half_n as isize) {
        let u = (k as f64) * h;
        let w = Complex::new(1.0, u);
        let s = Complex::new(sigma, 0.0) + w * w * mu;
        let ds_du = Complex::new(0.0, 2.0 * mu) * w;
        sum += (s * t).exp() * big_f(s) * ds_du;
    }
    Ok(sum * h / Complex::new(0.0, 2.0 * PI))
}

/// Numerical Bromwich Laplace inversion using Weideman's parabolic contour
/// $s(u) = \mu (i u + 1)^2$ ($u \in \mathbb{R}$) for $t > 0$.
///
/// # Errors
///
/// Returns an error if time is nonpositive or there are fewer than four quadrature points.
pub fn bromwich_weideman_parabola<F>(
    big_f: F,
    t: f64,
    n_pts: usize,
) -> Result<Complex<f64>, String>
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    bromwich_weideman_parabola_shifted(big_f, t, 0.0, n_pts)
}

/// Gaver–Stehfest real-axis inverse Laplace transform for $t > 0$ (`m_even` in `2..=14`).
///
/// # Errors
///
/// Returns an error if time is nonpositive or the order is not even and in `2..=14`.
pub fn gaver_stehfest_inverse_laplace<F>(
    big_f: F,
    t: f64,
    m_even: usize,
) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    if t <= 0.0 || m_even == 0 || !m_even.is_multiple_of(2) || m_even > 14 {
        return Err("t must be positive and m_even must be an even integer in 2..=14".to_string());
    }
    let half = m_even / 2;
    let fact = |n: usize| -> f64 {
        let mut p = 1.0_f64;
        for i in 2..=n {
            p *= i as f64;
        }
        p
    };
    let ln2_t = 2.0_f64.ln() / t;
    let mut sum = 0.0_f64;
    for k in 1..=m_even {
        let j_min = k.div_ceil(2);
        let j_max = k.min(half);
        let mut v_k = 0.0_f64;
        for j in j_min..=j_max {
            let num = (j as f64).powi(half as i32) * fact(2 * j);
            let den = fact(half - j) * fact(j) * fact(j - 1) * fact(k - j) * fact(2 * j - k);
            let sign = if (k + half).is_multiple_of(2) {
                1.0
            } else {
                -1.0
            };
            v_k += sign * num / den;
        }
        sum += v_k * big_f((k as f64) * ln2_t);
    }
    Ok(ln2_t * sum)
}

/// Evaluates a Mellin–Barnes integral.
///
/// Evaluates a Mellin–Barnes integral
/// `\frac{1}{2\pi i} \int_{c - iT}^{c + iT} \frac{\prod \Gamma(\alpha_j + s) \prod \Gamma(\beta_j - s)}{\prod \Gamma(\gamma_j + s) \prod \Gamma(\delta_j - s)} z^{-s} \, ds`
/// along a slightly tilted/vertical contour `\text{Re}(s) = c`.
#[allow(clippy::too_many_arguments)]
#[must_use]
pub fn mellin_barnes_integral(
    num_plus: &[Complex<f64>],
    num_minus: &[Complex<f64>],
    den_plus: &[Complex<f64>],
    den_minus: &[Complex<f64>],
    z: Complex<f64>,
    c: f64,
    t_max: f64,
    tol: f64,
) -> Complex<f64> {
    let ln_z = z.ln();
    let val = integrate_complex_ab(
        |t| {
            let s = Complex::new(c, t);
            let mut num = Complex::new(1.0, 0.0);
            for &a in num_plus {
                num *= gamma_complex(a + s);
            }
            for &b in num_minus {
                num *= gamma_complex(b - s);
            }
            let mut den = Complex::new(1.0, 0.0);
            for &g in den_plus {
                den *= gamma_complex(g + s);
            }
            for &d in den_minus {
                den *= gamma_complex(d - s);
            }
            (num / den) * (-s * ln_z).exp()
        },
        -t_max,
        t_max,
        tol,
        500,
    );
    val / (2.0 * PI)
}

/// Evaluates the Meijer G-function.
///
/// Evaluates the Meijer G-function
/// $G_{p,q}^{m,n}\left(z \middle| \begin{smallmatrix} `a_1,\dots,a_p` \\ `b_1,\dots,b_q` \end{smallmatrix}\right)$
/// via its Mellin–Barnes contour representation along $\text{Re}(s) = c$.
///
/// # Errors
///
/// Returns an error if the numerator indices exceed the lengths of their respective parameter lists.
#[expect(
    clippy::too_many_arguments,
    reason = "The public signature follows the four Meijer G parameter groups and the contour parameters"
)]
pub fn meijer_g(
    a: &[Complex<f64>],
    b: &[Complex<f64>],
    m: usize,
    n: usize,
    z: Complex<f64>,
    c: f64,
    t_max: f64,
    tol: f64,
) -> Result<Complex<f64>, String> {
    if m > b.len() || n > a.len() {
        return Err("Invalid Meijer G indices (m <= q and n <= p required)".to_string());
    }
    let one = Complex::new(1.0, 0.0);
    // Standard convention: int (prod_{j=1}^m Gamma(b_j - s) prod_{j=1}^n Gamma(1 - a_j + s)) /
    // (prod_{j=m+1}^q Gamma(1 - b_j + s) prod_{j=n+1}^p Gamma(a_j - s)) * z^s ds / (2*pi*i)
    // Using s -> -s: z^{-s} with:
    // num_plus: b_j (j < m), num_minus: 1 - a_j (j < n),
    // den_plus: 1 - b_j (j >= m), den_minus: a_j (j >= n).
    let num_plus: Vec<Complex<f64>> = b.iter().take(m).copied().collect();
    let num_minus: Vec<Complex<f64>> = a.iter().take(n).map(|&ak| one - ak).collect();
    let den_plus: Vec<Complex<f64>> = b.iter().skip(m).map(|&bk| one - bk).collect();
    let den_minus: Vec<Complex<f64>> = a.iter().skip(n).copied().collect();
    Ok(mellin_barnes_integral(
        &num_plus, &num_minus, &den_plus, &den_minus, z, c, t_max, tol,
    ))
}

/// Numerical Hilbert transform $\mathcal{H}\[f\](x) = \frac{1}{\pi} \text{PV} \int_{-L}^L \frac{f(t)}{x - t} \, dt$
/// using singularity subtraction on $[-L, L]$.
///
/// # Errors
///
/// Returns an error if the evaluation point lies outside the open integration interval.
pub fn hilbert_transform_real<F>(
    f: F,
    x: f64,
    bound_l: f64,
    tol: f64,
) -> Result<f64, String>
where
    F: Fn(f64) -> f64,
{
    if x.abs() >= bound_l {
        return Err("Evaluation point x must lie strictly inside (-L, L)".to_string());
    }
    let fx = f(x);
    let smooth_int = gauss_kronrod_any(
        |t: f64| {
            let dt = x - t;
            if dt.abs() < 1e-12 {
                0.0
            } else {
                (f(t) - fx) / dt
            }
        },
        -bound_l,
        bound_l,
        tol,
        400,
    )
    .value;
    let pv_const = fx * ((bound_l + x) / (bound_l - x)).ln();
    Ok((smooth_int + pv_const) / PI)
}

/// Numerical Hankel transform of order `nu`:
/// `\mathcal{H}_\nu[f](k) = \int_0^{R_{\max}} f(r) J_\nu(k r) r \, dr`.
pub fn hankel_transform<F>(
    f: F,
    nu: f64,
    k: f64,
    r_max: f64,
    tol: f64,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    integrate_complex_ab(
        |r| {
            let j_val = special_complex::bessel_j_complex(nu, Complex::new(k * r, 0.0));
            f(r) * j_val * r
        },
        0.0,
        r_max,
        tol,
        400,
    )
}

/// Evaluates the truncated Borel transform `B[a](t) = \sum_{n=0}^{N-1} \frac{a_n}{n!} t^n`.
#[must_use]
pub fn borel_transform_series(
    coeffs: &[Complex<f64>],
    t: Complex<f64>,
) -> Complex<f64> {
    let mut sum = Complex::new(0.0, 0.0);
    let mut t_pow_over_fact = Complex::new(1.0, 0.0);
    for (n, &an) in coeffs.iter().enumerate() {
        if n > 0 {
            t_pow_over_fact *= t / (n as f64);
        }
        sum += an * t_pow_over_fact;
    }
    sum
}

/// Numerical Stieltjes transform `S[f](z) = int_0^{r_max} f(t) / (t + z) dt`
/// for `z` off `(-inf, 0]`.
pub fn stieltjes_transform<F>(
    f: F,
    z: Complex<f64>,
    r_max: f64,
    tol: f64,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    integrate_complex_ab(|t| f(t) / (t + z), 0.0, r_max, tol, 400)
}
