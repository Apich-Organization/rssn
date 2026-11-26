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

pub(super) const ENTRIES: &[Entry] = &[
    Entry {
        name: "lgamma_complex",
        params: "ss",
        eval: lgamma_complex,
        doc: "`lgamma_complex(re, im)`: `list(re, im)` of `ln Gamma(z)` for `z = re + i*im`, continuous off the non-positive real axis.",
    },
    Entry {
        name: "gamma_complex",
        params: "ss",
        eval: gamma_complex,
        doc: "`gamma_complex(re, im)`: `list(re, im)` of the complex Gamma function `Gamma(re + i*im)`.",
    },
    Entry {
        name: "digamma_complex",
        params: "ss",
        eval: digamma_complex,
        doc: "`digamma_complex(re, im)`: `list(re, im)` of the complex digamma function `psi(z) = Gamma'(z) / Gamma(z)`.",
    },
    Entry {
        name: "lambert_w_complex",
        params: "sss",
        eval: lambert_w_complex,
        doc: "`lambert_w_complex(re, im, k)`: `list(re, im)` of branch `k` of the complex Lambert W function solving `w exp(w) = re + i*im`.",
    },
    Entry {
        name: "bessel_j_complex",
        params: "sss",
        eval: bessel_j_complex,
        doc: "`bessel_j_complex(nu, re, im)`: `list(re, im)` of the Bessel function of the first kind `J_nu(re + i*im)` for real order `nu`.",
    },
    Entry {
        name: "bessel_i_complex",
        params: "sss",
        eval: bessel_i_complex,
        doc: "`bessel_i_complex(nu, re, im)`: `list(re, im)` of the modified Bessel function of the first kind `I_nu(re + i*im)` for real order `nu`.",
    },
    Entry {
        name: "erf_complex",
        params: "ss",
        eval: erf_complex,
        doc: "`erf_complex(re, im)`: `list(re, im)` of the complex error function `erf(re + i*im)`.",
    },
    Entry {
        name: "fresnel_s_complex",
        params: "ss",
        eval: fresnel_s_complex,
        doc: "`fresnel_s_complex(re, im)`: `list(re, im)` of the complex Fresnel sine integral `S(re + i*im)`.",
    },
    Entry {
        name: "fresnel_c_complex",
        params: "ss",
        eval: fresnel_c_complex,
        doc: "`fresnel_c_complex(re, im)`: `list(re, im)` of the complex Fresnel cosine integral `C(re + i*im)`.",
    },
    Entry {
        name: "rgamma_complex",
        params: "ss",
        eval: rgamma_complex,
        doc: "`rgamma_complex(re, im)`: `list(re, im)` of `1 / Gamma(z)`.",
    },
    Entry {
        name: "polygamma_complex",
        params: "sss",
        eval: polygamma_complex,
        doc: "`polygamma_complex(n, re, im)`: `list(re, im)` of `psi^(n)(z)`.",
    },
    Entry {
        name: "beta_complex",
        params: "ssss",
        eval: beta_complex,
        doc: "`beta_complex(ar, ai, br, bi)`: `list(re, im)` of `B(a, b)`.",
    },
    Entry {
        name: "ln_beta_complex",
        params: "ssss",
        eval: ln_beta_complex,
        doc: "`ln_beta_complex(ar, ai, br, bi)`: `list(re, im)` of `ln B(a, b)`.",
    },
    Entry {
        name: "lower_gamma_complex",
        params: "ssss",
        eval: lower_gamma_complex,
        doc: "`lower_gamma_complex(ar, ai, zr, zi)`: `list(re, im)` of lower incomplete gamma `gamma(a, z)`.",
    },
    Entry {
        name: "upper_gamma_complex",
        params: "ssss",
        eval: upper_gamma_complex,
        doc: "`upper_gamma_complex(ar, ai, zr, zi)`: `list(re, im)` of upper incomplete gamma `Gamma(a, z)`.",
    },
    Entry {
        name: "incomplete_beta_complex",
        params: "ssssss",
        eval: incomplete_beta_complex,
        doc: "`incomplete_beta_complex(ar, ai, br, bi, zr, zi)`: `list(re, im)` of `B_z(a, b)`.",
    },
    Entry {
        name: "pochhammer_complex",
        params: "ssss",
        eval: pochhammer_complex,
        doc: "`pochhammer_complex(zr, zi, wr, wi)`: `list(re, im)` of `(z)_w`.",
    },
    Entry {
        name: "barnes_g_complex",
        params: "ss",
        eval: barnes_g_complex,
        doc: "`barnes_g_complex(re, im)`: `list(re, im)` of Barnes G-function `G(z)`.",
    },
    Entry {
        name: "erfc_complex",
        params: "ss",
        eval: erfc_complex,
        doc: "`erfc_complex(re, im)`: `list(re, im)` of `erfc(z)`.",
    },
    Entry {
        name: "erfcx_complex",
        params: "ss",
        eval: erfcx_complex,
        doc: "`erfcx_complex(re, im)`: `list(re, im)` of `exp(z^2) * erfc(z)`.",
    },
    Entry {
        name: "erfi_complex",
        params: "ss",
        eval: erfi_complex,
        doc: "`erfi_complex(re, im)`: `list(re, im)` of `erfi(z)`.",
    },
    Entry {
        name: "faddeeva_w_complex",
        params: "ss",
        eval: faddeeva_w_complex,
        doc: "`faddeeva_w_complex(re, im)`: `list(re, im)` of Faddeeva function `w(z)`.",
    },
    Entry {
        name: "dawson_complex",
        params: "ss",
        eval: dawson_complex,
        doc: "`dawson_complex(re, im)`: `list(re, im)` of Dawson integral `D+(z)`.",
    },
    Entry {
        name: "exp_integral_e1_complex",
        params: "ss",
        eval: exp_integral_e1_complex,
        doc: "`exp_integral_e1_complex(re, im)`: `list(re, im)` of `E_1(z)`.",
    },
    Entry {
        name: "exp_integral_ei_complex",
        params: "ss",
        eval: exp_integral_ei_complex,
        doc: "`exp_integral_ei_complex(re, im)`: `list(re, im)` of `Ei(z)`.",
    },
    Entry {
        name: "log_integral_complex",
        params: "ss",
        eval: log_integral_complex,
        doc: "`log_integral_complex(re, im)`: `list(re, im)` of `li(z)`.",
    },
    Entry {
        name: "sin_integral_complex",
        params: "ss",
        eval: sin_integral_complex,
        doc: "`sin_integral_complex(re, im)`: `list(re, im)` of `Si(z)`.",
    },
    Entry {
        name: "cos_integral_complex",
        params: "ss",
        eval: cos_integral_complex,
        doc: "`cos_integral_complex(re, im)`: `list(re, im)` of `Ci(z)`.",
    },
    Entry {
        name: "sinh_integral_complex",
        params: "ss",
        eval: sinh_integral_complex,
        doc: "`sinh_integral_complex(re, im)`: `list(re, im)` of `Shi(z)`.",
    },
    Entry {
        name: "cosh_integral_complex",
        params: "ss",
        eval: cosh_integral_complex,
        doc: "`cosh_integral_complex(re, im)`: `list(re, im)` of `Chi(z)`.",
    },
    Entry {
        name: "bessel_y_complex",
        params: "sss",
        eval: bessel_y_complex,
        doc: "`bessel_y_complex(nu, re, im)`: `list(re, im)` of `Y_nu(z)`.",
    },
    Entry {
        name: "bessel_k_complex",
        params: "sss",
        eval: bessel_k_complex,
        doc: "`bessel_k_complex(nu, re, im)`: `list(re, im)` of `K_nu(z)`.",
    },
    Entry {
        name: "hankel_h1_complex",
        params: "sss",
        eval: hankel_h1_complex,
        doc: "`hankel_h1_complex(nu, re, im)`: `list(re, im)` of `H_nu^(1)(z)`.",
    },
    Entry {
        name: "hankel_h2_complex",
        params: "sss",
        eval: hankel_h2_complex,
        doc: "`hankel_h2_complex(nu, re, im)`: `list(re, im)` of `H_nu^(2)(z)`.",
    },
    Entry {
        name: "airy_ai_complex",
        params: "ss",
        eval: airy_ai_complex,
        doc: "`airy_ai_complex(re, im)`: `list(re, im)` of `Ai(z)`.",
    },
    Entry {
        name: "airy_ai_prime_complex",
        params: "ss",
        eval: airy_ai_prime_complex,
        doc: "`airy_ai_prime_complex(re, im)`: `list(re, im)` of `Ai'(z)`.",
    },
    Entry {
        name: "airy_bi_complex",
        params: "ss",
        eval: airy_bi_complex,
        doc: "`airy_bi_complex(re, im)`: `list(re, im)` of `Bi(z)`.",
    },
    Entry {
        name: "airy_bi_prime_complex",
        params: "ss",
        eval: airy_bi_prime_complex,
        doc: "`airy_bi_prime_complex(re, im)`: `list(re, im)` of `Bi'(z)`.",
    },
    Entry {
        name: "struve_h_complex",
        params: "sss",
        eval: struve_h_complex,
        doc: "`struve_h_complex(nu, re, im)`: `list(re, im)` of `H_nu(z)`.",
    },
    Entry {
        name: "struve_l_complex",
        params: "sss",
        eval: struve_l_complex,
        doc: "`struve_l_complex(nu, re, im)`: `list(re, im)` of `L_nu(z)`.",
    },
    Entry {
        name: "riemann_zeta_complex",
        params: "ss",
        eval: riemann_zeta_complex,
        doc: "`riemann_zeta_complex(re, im)`: `list(re, im)` of `zeta(s)`.",
    },
    Entry {
        name: "dirichlet_eta_complex",
        params: "ss",
        eval: dirichlet_eta_complex,
        doc: "`dirichlet_eta_complex(re, im)`: `list(re, im)` of `eta(s)`.",
    },
    Entry {
        name: "hurwitz_zeta_complex",
        params: "ssss",
        eval: hurwitz_zeta_complex,
        doc: "`hurwitz_zeta_complex(sr, si, ar, ai)`: `list(re, im)` of `zeta(s, a)`.",
    },
    Entry {
        name: "polylog_complex",
        params: "sss",
        eval: polylog_complex_op,
        doc: "`polylog_complex(s, re, im)`: `list(re, im)` of `Li_s(z)`.",
    },
    Entry {
        name: "lerch_phi_complex",
        params: "ssssss",
        eval: lerch_phi_complex,
        doc: "`lerch_phi_complex(zr, zi, sr, si, ar, ai)`: `list(re, im)` of `Phi(z, s, a)`.",
    },
    Entry {
        name: "elliptic_k_complex",
        params: "ss",
        eval: elliptic_k_complex,
        doc: "`elliptic_k_complex(re, im)`: `list(re, im)` of `K(m)`.",
    },
    Entry {
        name: "elliptic_e_complex",
        params: "ss",
        eval: elliptic_e_complex,
        doc: "`elliptic_e_complex(re, im)`: `list(re, im)` of `E(m)`.",
    },
    Entry {
        name: "jacobi_sn_complex",
        params: "ssss",
        eval: jacobi_sn_complex,
        doc: "`jacobi_sn_complex(ur, ui, mr, mi)`: `list(re, im)` of `sn(u | m)`.",
    },
    Entry {
        name: "jacobi_cn_complex",
        params: "ssss",
        eval: jacobi_cn_complex,
        doc: "`jacobi_cn_complex(ur, ui, mr, mi)`: `list(re, im)` of `cn(u | m)`.",
    },
    Entry {
        name: "jacobi_dn_complex",
        params: "ssss",
        eval: jacobi_dn_complex,
        doc: "`jacobi_dn_complex(ur, ui, mr, mi)`: `list(re, im)` of `dn(u | m)`.",
    },
    Entry {
        name: "jacobi_theta_complex",
        params: "sssss",
        eval: jacobi_theta_complex,
        doc: "`jacobi_theta_complex(k, zr, zi, qr, qi)`: `list(re, im)` of `theta_k(z, q)`.",
    },
    Entry {
        name: "weierstrass_p_complex",
        params: "ssssss",
        eval: weierstrass_p_complex,
        doc: "`weierstrass_p_complex(zr, zi, w1r, w1i, w2r, w2i)`: `list(re, im)` of `wp(z; w1, w2)`.",
    },
    Entry {
        name: "dedekind_eta_complex",
        params: "ss",
        eval: dedekind_eta_complex,
        doc: "`dedekind_eta_complex(re, im)`: `list(re, im)` of `eta(tau)`.",
    },
    Entry {
        name: "klein_j_complex",
        params: "ss",
        eval: klein_j_complex,
        doc: "`klein_j_complex(re, im)`: `list(re, im)` of `j(tau)`.",
    },
    Entry {
        name: "wright_omega_complex",
        params: "ss",
        eval: wright_omega_complex,
        doc: "`wright_omega_complex(re, im)`: `list(re, im)` of `omega(z)`.",
    },
    Entry {
        name: "whittaker_m_complex",
        params: "ssssss",
        eval: whittaker_m_complex,
        doc: "`whittaker_m_complex(kr, ki, mr, mi, zr, zi)`: `list(re, im)` of `M_{kappa, mu}(z)`.",
    },
    Entry {
        name: "whittaker_w_complex",
        params: "ssssss",
        eval: whittaker_w_complex,
        doc: "`whittaker_w_complex(kr, ki, mr, mi, zr, zi)`: `list(re, im)` of `W_{kappa, mu}(z)`.",
    },
    Entry {
        name: "parabolic_cylinder_d_complex",
        params: "ssss",
        eval: parabolic_cylinder_d_complex,
        doc: "`parabolic_cylinder_d_complex(nur, nui, zr, zi)`: `list(re, im)` of `D_nu(z)`.",
    },
    Entry {
        name: "hyp0f1_complex",
        params: "ssss",
        eval: hyp0f1_complex,
        doc: "`hyp0f1_complex(br, bi, zr, zi)`: `list(re, im)` of `0F1(; b; z)`.",
    },
    Entry {
        name: "hyp1f1_complex",
        params: "ssssss",
        eval: hyp1f1_complex,
        doc: "`hyp1f1_complex(ar, ai, br, bi, zr, zi)`: `list(re, im)` of `1F1(a; b; z)`.",
    },
    Entry {
        name: "hyp2f1_complex",
        params: "ssssssss",
        eval: hyp2f1_complex,
        doc: "`hyp2f1_complex(ar, ai, br, bi, cr, ci, zr, zi)`: `list(re, im)` of `2F1(a, b; c; z)`.",
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
