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
use super::bridge::require;
use super::bridge::scalar;
use crate::kernel::diag::DiagKind;
use crate::kernels::fresnel;
use crate::kernels::hypergeometric;
use crate::kernels::special_complex;
use crate::kernels::special_complex_ext as ext;
use crate::kernels::special_ext;

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

fn rgamma_complex(a: &[Data]) -> Option<Data> {
    pair(ext::rgamma_complex(point(a, 0)?))
}

fn polygamma_complex(a: &[Data]) -> Option<Data> {
    let n = scalar(a, 0)?;
    if n < 0.0 || n.fract() != 0.0 || n > 20.0 {
        return None;
    }
    pair(ext::polygamma_complex(n as u32, point(a, 1)?))
}

fn beta_complex(a: &[Data]) -> Option<Data> {
    pair(ext::beta_complex(point(a, 0)?, point(a, 2)?))
}

fn ln_beta_complex(a: &[Data]) -> Option<Data> {
    pair(ext::ln_beta_complex(point(a, 0)?, point(a, 2)?))
}

fn lower_gamma_complex(a: &[Data]) -> Option<Data> {
    pair(ext::lower_incomplete_gamma_complex(
        point(a, 0)?,
        point(a, 2)?,
    ))
}

fn upper_gamma_complex(a: &[Data]) -> Option<Data> {
    pair(ext::upper_incomplete_gamma_complex(
        point(a, 0)?,
        point(a, 2)?,
    ))
}

fn incomplete_beta_complex(a: &[Data]) -> Option<Data> {
    pair(ext::incomplete_beta_complex(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn pochhammer_complex(a: &[Data]) -> Option<Data> {
    pair(ext::pochhammer_complex(point(a, 0)?, point(a, 2)?))
}

fn barnes_g_complex(a: &[Data]) -> Option<Data> {
    pair(ext::barnes_g_complex(point(a, 0)?))
}

fn erfc_complex(a: &[Data]) -> Option<Data> {
    pair(ext::erfc_complex(erf_argument(a)?))
}

fn erfcx_complex(a: &[Data]) -> Option<Data> {
    pair(ext::erfcx_complex(erf_argument(a)?))
}

fn erfi_complex(a: &[Data]) -> Option<Data> {
    pair(ext::erfi_complex(erf_argument(a)?))
}

fn faddeeva_w_complex(a: &[Data]) -> Option<Data> {
    pair(ext::voigt_complex(erf_argument(a)?))
}

fn dawson_complex(a: &[Data]) -> Option<Data> {
    pair(ext::dawson_complex(erf_argument(a)?))
}

fn exp_integral_e1_complex(a: &[Data]) -> Option<Data> {
    pair(ext::exp_integral_e1_complex(point(a, 0)?))
}

fn exp_integral_ei_complex(a: &[Data]) -> Option<Data> {
    pair(ext::exp_integral_ei_complex(point(a, 0)?))
}

fn log_integral_complex(a: &[Data]) -> Option<Data> {
    pair(ext::log_integral_complex(point(a, 0)?))
}

fn sin_integral_complex(a: &[Data]) -> Option<Data> {
    pair(ext::sin_cos_integral_complex(point(a, 0)?).0)
}

fn cos_integral_complex(a: &[Data]) -> Option<Data> {
    pair(ext::sin_cos_integral_complex(point(a, 0)?).1)
}

fn sinh_integral_complex(a: &[Data]) -> Option<Data> {
    pair(ext::sinh_cosh_integral_complex(point(a, 0)?).0)
}

fn cosh_integral_complex(a: &[Data]) -> Option<Data> {
    pair(ext::sinh_cosh_integral_complex(point(a, 0)?).1)
}

fn bessel_y_complex(a: &[Data]) -> Option<Data> {
    let (nu, z) = bessel_arguments(a)?;
    pair(ext::bessel_y_complex(nu, z))
}

fn bessel_k_complex(a: &[Data]) -> Option<Data> {
    let (nu, z) = bessel_arguments(a)?;
    pair(ext::bessel_k_complex(nu, z))
}

fn hankel_h1_complex(a: &[Data]) -> Option<Data> {
    let (nu, z) = bessel_arguments(a)?;
    pair(ext::hankel_h1_complex(nu, z))
}

fn hankel_h2_complex(a: &[Data]) -> Option<Data> {
    let (nu, z) = bessel_arguments(a)?;
    pair(ext::hankel_h2_complex(nu, z))
}

fn airy_ai_complex(a: &[Data]) -> Option<Data> {
    pair(ext::airy_ai_complex(point(a, 0)?).0)
}

fn airy_ai_prime_complex(a: &[Data]) -> Option<Data> {
    pair(ext::airy_ai_complex(point(a, 0)?).1)
}

fn airy_bi_complex(a: &[Data]) -> Option<Data> {
    pair(ext::airy_bi_complex(point(a, 0)?).0)
}

fn airy_bi_prime_complex(a: &[Data]) -> Option<Data> {
    pair(ext::airy_bi_complex(point(a, 0)?).1)
}

fn struve_h_complex(a: &[Data]) -> Option<Data> {
    let (nu, z) = bessel_arguments(a)?;
    pair(ext::struve_h_complex(nu, z))
}

fn struve_l_complex(a: &[Data]) -> Option<Data> {
    let (nu, z) = bessel_arguments(a)?;
    pair(ext::struve_l_complex(nu, z))
}

fn riemann_zeta_complex(a: &[Data]) -> Option<Data> {
    pair(ext::riemann_zeta_complex(point(a, 0)?))
}

fn dirichlet_eta_complex(a: &[Data]) -> Option<Data> {
    pair(ext::dirichlet_eta_complex(point(a, 0)?))
}

fn hurwitz_zeta_complex(a: &[Data]) -> Option<Data> {
    pair(ext::hurwitz_zeta_complex(point(a, 0)?, point(a, 2)?))
}

fn polylog_complex_op(a: &[Data]) -> Option<Data> {
    let s = scalar(a, 0)?;
    if s.abs() > 100.0 {
        return None;
    }
    pair(special_ext::polylog_complex(
        Complex64::new(s, 0.0),
        point(a, 1)?,
    ))
}

fn lerch_phi_complex(a: &[Data]) -> Option<Data> {
    pair(ext::lerch_phi_complex(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn elliptic_k_complex(a: &[Data]) -> Option<Data> {
    pair(ext::elliptic_k_complex(point(a, 0)?))
}

fn elliptic_e_complex(a: &[Data]) -> Option<Data> {
    pair(ext::elliptic_e_complex(point(a, 0)?))
}

fn jacobi_sn_complex(a: &[Data]) -> Option<Data> {
    pair(ext::jacobi_elliptic_complex(point(a, 0)?, point(a, 2)?).0)
}

fn jacobi_cn_complex(a: &[Data]) -> Option<Data> {
    pair(ext::jacobi_elliptic_complex(point(a, 0)?, point(a, 2)?).1)
}

fn jacobi_dn_complex(a: &[Data]) -> Option<Data> {
    pair(ext::jacobi_elliptic_complex(point(a, 0)?, point(a, 2)?).2)
}

fn jacobi_theta_complex(a: &[Data]) -> Option<Data> {
    let k = scalar(a, 0)?;
    if k.fract() != 0.0 || !(1.0..=4.0).contains(&k) {
        return None;
    }
    pair(ext::jacobi_theta_complex(
        k as u8,
        point(a, 1)?,
        point(a, 3)?,
    ))
}

fn weierstrass_p_complex(a: &[Data]) -> Option<Data> {
    pair(ext::weierstrass_p_complex(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn dedekind_eta_complex(a: &[Data]) -> Option<Data> {
    pair(ext::dedekind_eta_complex(point(a, 0)?))
}

fn klein_j_complex(a: &[Data]) -> Option<Data> {
    pair(ext::klein_j_complex(point(a, 0)?))
}

fn wright_omega_complex(a: &[Data]) -> Option<Data> {
    pair(ext::wright_omega_complex(point(a, 0)?))
}

fn whittaker_m_complex(a: &[Data]) -> Option<Data> {
    pair(ext::whittaker_m_complex(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn whittaker_w_complex(a: &[Data]) -> Option<Data> {
    pair(ext::whittaker_w_complex(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn parabolic_cylinder_d_complex(a: &[Data]) -> Option<Data> {
    pair(ext::parabolic_cylinder_d_complex(
        point(a, 0)?,
        point(a, 2)?,
    ))
}

fn hyp0f1_complex(a: &[Data]) -> Option<Data> {
    pair(ext::hyp0f1_complex(point(a, 0)?, point(a, 2)?))
}

fn hyp1f1_complex(a: &[Data]) -> Option<Data> {
    pair(hypergeometric::hyp1f1_cplx(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn hyp2f1_complex(a: &[Data]) -> Option<Data> {
    pair(hypergeometric::hyp2f1_cplx(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
        point(a, 6)?,
    ))
}

fn binomial_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::binomial_complex(point(a, 0)?, point(a, 2)?))
}

fn reg_lower_gamma_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::reg_lower_gamma_complex(point(a, 0)?, point(a, 2)?))
}

fn reg_upper_gamma_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::reg_upper_gamma_complex(point(a, 0)?, point(a, 2)?))
}

fn reg_beta_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::reg_beta_complex(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn spence_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::spence_complex(point(a, 0)?))
}

fn bessel_ke_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::bessel_ke_complex(scalar(a, 0)?, point(a, 1)?))
}

fn spherical_jn_complex_op(a: &[Data]) -> Option<Data> {
    let n = scalar(a, 0)?.round() as i32;
    pair(ext::spherical_jn_complex(n, point(a, 1)?))
}

fn spherical_yn_complex_op(a: &[Data]) -> Option<Data> {
    let n = scalar(a, 0)?.round() as i32;
    pair(ext::spherical_yn_complex(n, point(a, 1)?))
}

fn spherical_in_complex_op(a: &[Data]) -> Option<Data> {
    let n = scalar(a, 0)?.round() as i32;
    pair(ext::spherical_in_complex(n, point(a, 1)?))
}

fn spherical_kn_complex_op(a: &[Data]) -> Option<Data> {
    let n = scalar(a, 0)?.round() as i32;
    pair(ext::spherical_kn_complex(n, point(a, 1)?))
}


fn airy_gi_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::airy_gi_complex(point(a, 0)?))
}

fn airy_hi_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::airy_hi_complex(point(a, 0)?))
}

fn elliptic_f_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::elliptic_f_complex(point(a, 0)?, point(a, 2)?))
}

fn elliptic_e_inc_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::elliptic_e_inc_complex(point(a, 0)?, point(a, 2)?))
}

fn elliptic_pi_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::elliptic_pi_complex(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn mittag_leffler_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::mittag_leffler_complex(
        point(a, 0)?,
        point(a, 2)?,
        point(a, 4)?,
    ))
}

fn eisenstein_e4_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::eisenstein_e4_complex(point(a, 0)?))
}

fn eisenstein_e6_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::eisenstein_e6_complex(point(a, 0)?))
}

fn modular_lambda_complex_op(a: &[Data]) -> Option<Data> {
    pair(ext::modular_lambda_complex(point(a, 0)?))
}

// Multi-component kernels share a source path across the operators that
// expose their individual results.
reviewed_entries! {
    (
        "lgamma_complex", "ss", lgamma_complex, special_complex::ln_gamma_complex,
        "`lgamma_complex(re, im)`: `list(re, im)` of `ln Gamma(z)` for `z = re + i*im`, continuous off the non-positive real axis."
    ),
    (
        "gamma_complex", "ss", gamma_complex, special_complex::gamma_complex,
        "`gamma_complex(re, im)`: `list(re, im)` of the complex Gamma function `Gamma(re + i*im)`."
    ),
    (
        "digamma_complex", "ss", digamma_complex, special_complex::digamma_complex,
        "`digamma_complex(re, im)`: `list(re, im)` of the complex digamma function `psi(z) = Gamma'(z) / Gamma(z)`."
    ),
    (
        "lambert_w_complex", "sss", lambert_w_complex, special_complex::lambert_w_complex,
        "`lambert_w_complex(re, im, k)`: `list(re, im)` of branch `k` of the complex Lambert W function solving `w exp(w) = re + i*im`."
    ),
    (
        "bessel_j_complex", "sss", bessel_j_complex, special_complex::bessel_j_complex,
        "`bessel_j_complex(nu, re, im)`: `list(re, im)` of the Bessel function of the first kind `J_nu(re + i*im)` for real order `nu`."
    ),
    (
        "bessel_i_complex", "sss", bessel_i_complex, special_complex::bessel_i_complex,
        "`bessel_i_complex(nu, re, im)`: `list(re, im)` of the modified Bessel function of the first kind `I_nu(re + i*im)` for real order `nu`."
    ),
    (
        "erf_complex", "ss", erf_complex, fresnel::erf_complex,
        "`erf_complex(re, im)`: `list(re, im)` of the complex error function `erf(re + i*im)`."
    ),
    (
        "fresnel_s_complex", "ss", fresnel_s_complex, fresnel::fresnel_sc_complex,
        "`fresnel_s_complex(re, im)`: `list(re, im)` of the complex Fresnel sine integral `S(re + i*im)`."
    ),
    (
        "fresnel_c_complex", "ss", fresnel_c_complex, fresnel::fresnel_sc_complex,
        "`fresnel_c_complex(re, im)`: `list(re, im)` of the complex Fresnel cosine integral `C(re + i*im)`."
    ),
    (
        "rgamma_complex", "ss", rgamma_complex, special_complex_ext::gamma_beta::rgamma_complex,
        "`rgamma_complex(re, im)`: `list(re, im)` of `1 / Gamma(z)`."
    ),
    (
        "polygamma_complex", "sss", polygamma_complex, special_complex_ext::gamma_beta::polygamma_complex,
        "`polygamma_complex(n, re, im)`: `list(re, im)` of `psi^(n)(z)`."
    ),
    (
        "beta_complex", "ssss", beta_complex, special_complex_ext::gamma_beta::beta_complex,
        "`beta_complex(ar, ai, br, bi)`: `list(re, im)` of `B(a, b)`."
    ),
    (
        "ln_beta_complex", "ssss", ln_beta_complex, special_complex_ext::gamma_beta::ln_beta_complex,
        "`ln_beta_complex(ar, ai, br, bi)`: `list(re, im)` of `ln B(a, b)`."
    ),
    (
        "lower_gamma_complex", "ssss", lower_gamma_complex, special_complex_ext::gamma_beta::lower_incomplete_gamma_complex,
        "`lower_gamma_complex(ar, ai, zr, zi)`: `list(re, im)` of lower incomplete gamma `gamma(a, z)`."
    ),
    (
        "upper_gamma_complex", "ssss", upper_gamma_complex, special_complex_ext::gamma_beta::upper_incomplete_gamma_complex,
        "`upper_gamma_complex(ar, ai, zr, zi)`: `list(re, im)` of upper incomplete gamma `Gamma(a, z)`."
    ),
    (
        "incomplete_beta_complex", "ssssss", incomplete_beta_complex, special_complex_ext::gamma_beta::incomplete_beta_complex,
        "`incomplete_beta_complex(ar, ai, br, bi, zr, zi)`: `list(re, im)` of `B_z(a, b)`."
    ),
    (
        "pochhammer_complex", "ssss", pochhammer_complex, special_complex_ext::gamma_beta::pochhammer_complex,
        "`pochhammer_complex(zr, zi, wr, wi)`: `list(re, im)` of `(z)_w`."
    ),
    (
        "barnes_g_complex", "ss", barnes_g_complex, special_complex_ext::gamma_beta::barnes_g_complex,
        "`barnes_g_complex(re, im)`: `list(re, im)` of Barnes G-function `G(z)`."
    ),
    (
        "erfc_complex", "ss", erfc_complex, special_complex_ext::erf_integrals::erfc_complex,
        "`erfc_complex(re, im)`: `list(re, im)` of `erfc(z)`."
    ),
    (
        "erfcx_complex", "ss", erfcx_complex, special_complex_ext::erf_integrals::erfcx_complex,
        "`erfcx_complex(re, im)`: `list(re, im)` of `exp(z^2) * erfc(z)`."
    ),
    (
        "erfi_complex", "ss", erfi_complex, special_complex_ext::erf_integrals::erfi_complex,
        "`erfi_complex(re, im)`: `list(re, im)` of `erfi(z)`."
    ),
    (
        "faddeeva_w_complex", "ss", faddeeva_w_complex, special_complex_ext::erf_integrals::voigt_complex,
        "`faddeeva_w_complex(re, im)`: `list(re, im)` of Faddeeva function `w(z)`."
    ),
    (
        "dawson_complex", "ss", dawson_complex, special_complex_ext::erf_integrals::dawson_complex,
        "`dawson_complex(re, im)`: `list(re, im)` of Dawson integral `D+(z)`."
    ),
    (
        "exp_integral_e1_complex", "ss", exp_integral_e1_complex, special_complex_ext::erf_integrals::exp_integral_e1_complex,
        "`exp_integral_e1_complex(re, im)`: `list(re, im)` of `E_1(z)`."
    ),
    (
        "exp_integral_ei_complex", "ss", exp_integral_ei_complex, special_complex_ext::erf_integrals::exp_integral_ei_complex,
        "`exp_integral_ei_complex(re, im)`: `list(re, im)` of `Ei(z)`."
    ),
    (
        "log_integral_complex", "ss", log_integral_complex, special_complex_ext::erf_integrals::log_integral_complex,
        "`log_integral_complex(re, im)`: `list(re, im)` of `li(z)`."
    ),
    (
        "sin_integral_complex", "ss", sin_integral_complex, special_complex_ext::erf_integrals::sin_cos_integral_complex,
        "`sin_integral_complex(re, im)`: `list(re, im)` of `Si(z)`."
    ),
    (
        "cos_integral_complex", "ss", cos_integral_complex, special_complex_ext::erf_integrals::sin_cos_integral_complex,
        "`cos_integral_complex(re, im)`: `list(re, im)` of `Ci(z)`."
    ),
    (
        "sinh_integral_complex", "ss", sinh_integral_complex, special_complex_ext::erf_integrals::sinh_cosh_integral_complex,
        "`sinh_integral_complex(re, im)`: `list(re, im)` of `Shi(z)`."
    ),
    (
        "cosh_integral_complex", "ss", cosh_integral_complex, special_complex_ext::erf_integrals::sinh_cosh_integral_complex,
        "`cosh_integral_complex(re, im)`: `list(re, im)` of `Chi(z)`."
    ),
    (
        "bessel_y_complex", "sss", bessel_y_complex, special_complex_ext::bessel_airy::bessel_y_complex,
        "`bessel_y_complex(nu, re, im)`: `list(re, im)` of `Y_nu(z)`."
    ),
    (
        "bessel_k_complex", "sss", bessel_k_complex, special_complex_ext::bessel_airy::bessel_k_complex,
        "`bessel_k_complex(nu, re, im)`: `list(re, im)` of `K_nu(z)`."
    ),
    (
        "hankel_h1_complex", "sss", hankel_h1_complex, special_complex_ext::bessel_airy::hankel_h1_complex,
        "`hankel_h1_complex(nu, re, im)`: `list(re, im)` of `H_nu^(1)(z)`."
    ),
    (
        "hankel_h2_complex", "sss", hankel_h2_complex, special_complex_ext::bessel_airy::hankel_h2_complex,
        "`hankel_h2_complex(nu, re, im)`: `list(re, im)` of `H_nu^(2)(z)`."
    ),
    (
        "airy_ai_complex", "ss", airy_ai_complex, special_complex_ext::bessel_airy::airy_ai_complex,
        "`airy_ai_complex(re, im)`: `list(re, im)` of `Ai(z)`."
    ),
    (
        "airy_ai_prime_complex", "ss", airy_ai_prime_complex, special_complex_ext::bessel_airy::airy_ai_complex,
        "`airy_ai_prime_complex(re, im)`: `list(re, im)` of `Ai'(z)`."
    ),
    (
        "airy_bi_complex", "ss", airy_bi_complex, special_complex_ext::bessel_airy::airy_bi_complex,
        "`airy_bi_complex(re, im)`: `list(re, im)` of `Bi(z)`."
    ),
    (
        "airy_bi_prime_complex", "ss", airy_bi_prime_complex, special_complex_ext::bessel_airy::airy_bi_complex,
        "`airy_bi_prime_complex(re, im)`: `list(re, im)` of `Bi'(z)`."
    ),
    (
        "struve_h_complex", "sss", struve_h_complex, special_complex_ext::bessel_airy::struve_h_complex,
        "`struve_h_complex(nu, re, im)`: `list(re, im)` of `H_nu(z)`."
    ),
    (
        "struve_l_complex", "sss", struve_l_complex, special_complex_ext::bessel_airy::struve_l_complex,
        "`struve_l_complex(nu, re, im)`: `list(re, im)` of `L_nu(z)`."
    ),
    (
        "riemann_zeta_complex", "ss", riemann_zeta_complex, special_complex_ext::zeta_polylog::riemann_zeta_complex,
        "`riemann_zeta_complex(re, im)`: `list(re, im)` of `zeta(s)`."
    ),
    (
        "dirichlet_eta_complex", "ss", dirichlet_eta_complex, special_complex_ext::zeta_polylog::dirichlet_eta_complex,
        "`dirichlet_eta_complex(re, im)`: `list(re, im)` of `eta(s)`."
    ),
    (
        "hurwitz_zeta_complex", "ssss", hurwitz_zeta_complex, special_complex_ext::zeta_polylog::hurwitz_zeta_complex,
        "`hurwitz_zeta_complex(sr, si, ar, ai)`: `list(re, im)` of `zeta(s, a)`."
    ),
    (
        "polylog_complex", "sss", polylog_complex_op, special_ext::polylog_complex,
        "`polylog_complex(s, re, im)`: `list(re, im)` of `Li_s(z)`."
    ),
    (
        "lerch_phi_complex", "ssssss", lerch_phi_complex, special_complex_ext::zeta_polylog::lerch_phi_complex,
        "`lerch_phi_complex(zr, zi, sr, si, ar, ai)`: `list(re, im)` of `Phi(z, s, a)`."
    ),
    (
        "elliptic_k_complex", "ss", elliptic_k_complex, special_complex_ext::elliptic_modular::elliptic_k_complex,
        "`elliptic_k_complex(re, im)`: `list(re, im)` of `K(m)`."
    ),
    (
        "elliptic_e_complex", "ss", elliptic_e_complex, special_complex_ext::elliptic_modular::elliptic_e_complex,
        "`elliptic_e_complex(re, im)`: `list(re, im)` of `E(m)`."
    ),
    (
        "jacobi_sn_complex", "ssss", jacobi_sn_complex, special_complex_ext::elliptic_modular::jacobi_elliptic_complex,
        "`jacobi_sn_complex(ur, ui, mr, mi)`: `list(re, im)` of `sn(u | m)`."
    ),
    (
        "jacobi_cn_complex", "ssss", jacobi_cn_complex, special_complex_ext::elliptic_modular::jacobi_elliptic_complex,
        "`jacobi_cn_complex(ur, ui, mr, mi)`: `list(re, im)` of `cn(u | m)`."
    ),
    (
        "jacobi_dn_complex", "ssss", jacobi_dn_complex, special_complex_ext::elliptic_modular::jacobi_elliptic_complex,
        "`jacobi_dn_complex(ur, ui, mr, mi)`: `list(re, im)` of `dn(u | m)`."
    ),
    (
        "jacobi_theta_complex", "sssss", jacobi_theta_complex, special_complex_ext::elliptic_modular::jacobi_theta_complex,
        "`jacobi_theta_complex(k, zr, zi, qr, qi)`: `list(re, im)` of `theta_k(z, q)`."
    ),
    (
        "weierstrass_p_complex", "ssssss", weierstrass_p_complex, special_complex_ext::elliptic_modular::weierstrass_p_complex,
        "`weierstrass_p_complex(zr, zi, w1r, w1i, w2r, w2i)`: `list(re, im)` of `wp(z; w1, w2)`."
    ),
    (
        "dedekind_eta_complex", "ss", dedekind_eta_complex, special_complex_ext::elliptic_modular::dedekind_eta_complex,
        "`dedekind_eta_complex(re, im)`: `list(re, im)` of `eta(tau)`."
    ),
    (
        "klein_j_complex", "ss", klein_j_complex, special_complex_ext::elliptic_modular::klein_j_complex,
        "`klein_j_complex(re, im)`: `list(re, im)` of `j(tau)`."
    ),
    (
        "wright_omega_complex", "ss", wright_omega_complex, special_complex_ext::orthopoly_hyper::wright_omega_complex,
        "`wright_omega_complex(re, im)`: `list(re, im)` of `omega(z)`."
    ),
    (
        "whittaker_m_complex", "ssssss", whittaker_m_complex, special_complex_ext::orthopoly_hyper::whittaker_m_complex,
        "`whittaker_m_complex(kr, ki, mr, mi, zr, zi)`: `list(re, im)` of `M_{kappa, mu}(z)`."
    ),
    (
        "whittaker_w_complex", "ssssss", whittaker_w_complex, special_complex_ext::orthopoly_hyper::whittaker_w_complex,
        "`whittaker_w_complex(kr, ki, mr, mi, zr, zi)`: `list(re, im)` of `W_{kappa, mu}(z)`."
    ),
    (
        "parabolic_cylinder_d_complex", "ssss", parabolic_cylinder_d_complex, special_complex_ext::orthopoly_hyper::parabolic_cylinder_d_complex,
        "`parabolic_cylinder_d_complex(nur, nui, zr, zi)`: `list(re, im)` of `D_nu(z)`."
    ),
    (
        "hyp0f1_complex", "ssss", hyp0f1_complex, special_complex_ext::orthopoly_hyper::hyp0f1_complex,
        "`hyp0f1_complex(br, bi, zr, zi)`: `list(re, im)` of `0F1(; b; z)`."
    ),
    (
        "hyp1f1_complex", "ssssss", hyp1f1_complex, hypergeometric::hyp1f1_cplx,
        "`hyp1f1_complex(ar, ai, br, bi, zr, zi)`: `list(re, im)` of `1F1(a; b; z)`."
    ),
    (
        "hyp2f1_complex", "ssssssss", hyp2f1_complex, hypergeometric::hyp2f1_cplx,
        "`hyp2f1_complex(ar, ai, br, bi, cr, ci, zr, zi)`: `list(re, im)` of `2F1(a, b; c; z)`."
    ),
    (
        "binomial_complex", "ssss", binomial_complex_op, special_complex_ext::gamma_beta::binomial_complex,
        "`binomial_complex(nr, ni, kr, ki)`: `list(re, im)` of `C(n, k)`."
    ),
    (
        "reg_lower_gamma_complex", "ssss", reg_lower_gamma_complex_op, special_complex_ext::gamma_beta::reg_lower_gamma_complex,
        "`reg_lower_gamma_complex(ar, ai, zr, zi)`: `list(re, im)` of `P(a, z)`."
    ),
    (
        "reg_upper_gamma_complex", "ssss", reg_upper_gamma_complex_op, special_complex_ext::gamma_beta::reg_upper_gamma_complex,
        "`reg_upper_gamma_complex(ar, ai, zr, zi)`: `list(re, im)` of `Q(a, z)`."
    ),
    (
        "reg_beta_complex", "ssssss", reg_beta_complex_op, special_complex_ext::gamma_beta::reg_beta_complex,
        "`reg_beta_complex(ar, ai, br, bi, zr, zi)`: `list(re, im)` of `I_z(a, b)`."
    ),
    (
        "spence_complex", "ss", spence_complex_op, special_complex_ext::erf_integrals::spence_complex,
        "`spence_complex(re, im)`: `list(re, im)` of `Li_2(1 - z)`."
    ),
    (
        "bessel_ke_complex", "sss", bessel_ke_complex_op, special_complex_ext::bessel_airy::bessel_ke_complex,
        "`bessel_ke_complex(nu, re, im)`: `list(re, im)` of `exp(z) * K_nu(z)`."
    ),
    (
        "spherical_jn_complex", "sss", spherical_jn_complex_op, special_complex_ext::bessel_airy::spherical_jn_complex,
        "`spherical_jn_complex(n, re, im)`: `list(re, im)` of `j_n(z)`."
    ),
    (
        "spherical_yn_complex", "sss", spherical_yn_complex_op, special_complex_ext::bessel_airy::spherical_yn_complex,
        "`spherical_yn_complex(n, re, im)`: `list(re, im)` of `y_n(z)`."
    ),
    (
        "spherical_in_complex", "sss", spherical_in_complex_op, special_complex_ext::bessel_airy::spherical_in_complex,
        "`spherical_in_complex(n, re, im)`: `list(re, im)` of `i_n(z)`."
    ),
    (
        "spherical_kn_complex", "sss", spherical_kn_complex_op, special_complex_ext::bessel_airy::spherical_kn_complex,
        "`spherical_kn_complex(n, re, im)`: `list(re, im)` of `k_n(z)`."
    ),
    (
        "airy_gi_complex", "ss", airy_gi_complex_op, special_complex_ext::bessel_airy::airy_gi_complex,
        "`airy_gi_complex(re, im)`: `list(re, im)` of Scorer `Gi(z)`."
    ),
    (
        "airy_hi_complex", "ss", airy_hi_complex_op, special_complex_ext::bessel_airy::airy_hi_complex,
        "`airy_hi_complex(re, im)`: `list(re, im)` of Scorer `Hi(z)`."
    ),
    (
        "elliptic_f_complex", "ssss", elliptic_f_complex_op, special_complex_ext::elliptic_modular::elliptic_f_complex,
        "`elliptic_f_complex(pr, pi, mr, mi)`: `list(re, im)` of `F(phi | m)`."
    ),
    (
        "elliptic_e_inc_complex", "ssss", elliptic_e_inc_complex_op, special_complex_ext::elliptic_modular::elliptic_e_inc_complex,
        "`elliptic_e_inc_complex(pr, pi, mr, mi)`: `list(re, im)` of `E(phi | m)`."
    ),
    (
        "elliptic_pi_complex", "ssssss", elliptic_pi_complex_op, special_complex_ext::elliptic_modular::elliptic_pi_complex,
        "`elliptic_pi_complex(nr, ni, pr, pi, mr, mi)`: `list(re, im)` of `Pi(n; phi | m)`."
    ),
    (
        "mittag_leffler_complex", "ssssss", mittag_leffler_complex_op, special_complex_ext::orthopoly_hyper::mittag_leffler_complex,
        "`mittag_leffler_complex(ar, ai, br, bi, zr, zi)`: `list(re, im)` of `E_{alpha, beta}(z)`."
    ),
    (
        "eisenstein_e4_complex", "ss", eisenstein_e4_complex_op, special_complex_ext::elliptic_modular::eisenstein_e4_complex,
        "`eisenstein_e4_complex(re, im)`: `list(re, im)` of `E_4(tau)`."
    ),
    (
        "eisenstein_e6_complex", "ss", eisenstein_e6_complex_op, special_complex_ext::elliptic_modular::eisenstein_e6_complex,
        "`eisenstein_e6_complex(re, im)`: `list(re, im)` of `E_6(tau)`."
    ),
    (
        "modular_lambda_complex", "ss", modular_lambda_complex_op, special_complex_ext::elliptic_modular::modular_lambda_complex,
        "`modular_lambda_complex(re, im)`: `list(re, im)` of `lambda(tau)`."
    ),
}
