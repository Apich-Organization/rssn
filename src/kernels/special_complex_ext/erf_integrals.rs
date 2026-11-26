//! Complex error-function family, Fresnel auxiliary functions, exponential
//! integrals, trigonometric/hyperbolic integrals, and Spence's dilogarithm.

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::fresnel::erf_complex;
use crate::kernels::fresnel::faddeeva_w;
use crate::kernels::fresnel::fresnel_sc_complex;
use crate::kernels::special_ext::polylog_complex;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);
const EULER_GAMMA: f64 = 0.577_215_664_901_532_9;

/// Complex complementary error function `erfc(z) = 1 - erf(z)`.
#[must_use]
pub fn erfc_complex(z: Complex64) -> Complex64 {
    if z.re >= 0.0 {
        (-z * z).exp() * faddeeva_w(I * z)
    } else {
        Complex64::new(2.0, 0.0) - (-z * z).exp() * faddeeva_w(-I * z)
    }
}

/// Scaled complementary error function `erfcx(z) = exp(z^2) * erfc(z) = w(i*z)`.
#[must_use]
pub fn erfcx_complex(z: Complex64) -> Complex64 {
    if z.re >= 0.0 {
        faddeeva_w(I * z)
    } else {
        2.0 * (z * z).exp() - faddeeva_w(-I * z)
    }
}

/// Imaginary error function `erfi(z) = -i * erf(i * z)`.
#[must_use]
pub fn erfi_complex(z: Complex64) -> Complex64 {
    -I * erf_complex(I * z)
}

/// Faddeeva / Voigt plasma dispersion function `w(z) = exp(-z^2) * erfc(-i * z)`.
#[must_use]
pub fn voigt_complex(z: Complex64) -> Complex64 {
    faddeeva_w(z)
}

/// Complex Dawson integral `D+(z) = exp(-z^2) * int_0^z exp(t^2) dt`.
#[must_use]
pub fn dawson_complex(z: Complex64) -> Complex64 {
    if z.norm() < 0.25 {
        let z2 = z * z;
        let mut term = z;
        let mut sum = term;
        for k in 1..40 {
            term *= -2.0 * z2 / ((2 * k + 1) as f64);
            sum += term;
            if term.norm() <= 1e-16 * (1.0 + sum.norm()) {
                break;
            }
        }
        return sum;
    }
    if z.im.abs() > 25.0 {
        // In |Im(z)| > 25, exp(-z^2) overflows or underflows; use asymptotic expansion:
        // D+(z) ~ 1/(2z) + 1/(4z^3) + 3/(8z^5) + ... (plus i*sqrt(pi)/2 * exp(-z^2) when |arg| < pi/4)
        let inv_2z2 = 0.5 / (z * z);
        let mut term = 0.5 / z;
        let mut sum = term;
        for k in 1..15 {
            term *= ((2 * k - 1) as f64) * inv_2z2;
            sum += term;
            if term.norm() <= 1e-16 * sum.norm() {
                break;
            }
        }
        return sum;
    }
    let half_sqrt_pi = 0.5 * PI.sqrt();
    I * half_sqrt_pi * ((-z * z).exp() - faddeeva_w(z))
}

/// Auxiliary Fresnel functions `(f(z), g(z))` in `C` satisfying
/// `C(z) = 1/2 + f(z) sin(pi z^2 / 2) - g(z) cos(pi z^2 / 2)`.
#[must_use]
pub fn fresnel_fg_complex(z: Complex64) -> (Complex64, Complex64) {
    if z.norm() > 12.0 && z.arg().abs() < 0.22 * PI {
        // Asymptotic expansion for large Re(z) avoiding 0.5 - C(z) cancellation:
        let piz2 = PI * z * z;
        let inv_piz2_sq = ONE / (piz2 * piz2);
        let mut f_sum = ONE;
        let mut g_sum = ONE;
        let mut f_term = ONE;
        let mut g_term = ONE;
        for k in 1..=8 {
            let m = (4 * k) as f64;
            f_term *= -(m - 3.0) * (m - 1.0) * inv_piz2_sq;
            g_term *= -(m - 1.0) * (m + 1.0) * inv_piz2_sq;
            f_sum += f_term;
            g_sum += g_term;
        }
        let f = f_sum / (PI * z);
        let g = g_sum / (PI * piz2 * z);
        return (f, g);
    }
    let ((s_re, s_im), (c_re, c_im)) = fresnel_sc_complex(z.re, z.im);
    let s = Complex64::new(s_re, s_im);
    let c = Complex64::new(c_re, c_im);
    let half = Complex64::new(0.5, 0.0);
    let dc = c - half;
    let ds = s - half;
    let arg = (0.5 * PI) * z * z;
    let sin_a = arg.sin();
    let cos_a = arg.cos();
    let f = dc * sin_a - ds * cos_a;
    let g = -dc * cos_a - ds * sin_a;
    (f, g)
}

/// Complex exponential integral `E_1(z) = int_z^oo exp(-t)/t dt`, cut along `(-inf, 0]`.
#[must_use]
pub fn exp_integral_e1_complex(z: Complex64) -> Complex64 {
    if z == ZERO {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    let near_neg_cut = z.re < 0.0 && z.im.abs() < 4.0;
    if z.norm() <= 4.0 || (near_neg_cut && z.norm() <= 20.0) {
        // E_1(z) = -gamma - ln(z) - sum_{k>=1} (-1)^k z^k / (k * k!)
        let mut sum = ZERO;
        let mut term = ONE;
        for k in 1..250 {
            let kf = k as f64;
            term *= -z / kf;
            let add = term / kf;
            sum += add;
            if add.norm() <= 1e-16 * (1.0 + sum.norm()) {
                break;
            }
        }
        return -EULER_GAMMA - z.ln() - sum;
    }
    if near_neg_cut {
        // Optimal-truncation asymptotic series for |z| > 20 near (-inf, 0]:
        // E_1(z) ~ exp(-z)/z * sum_{k=0}^{N} (-1)^k k! / z^k +/- i*pi (for Im(z) >= 0 vs < 0)
        let n_opt = (z.norm().floor() as usize).clamp(10, 40);
        let mut sum = ONE;
        let mut term = ONE;
        for k in 1..=n_opt {
            term *= -(k as f64) / z;
            sum += term;
            if term.norm() <= 1e-16 * sum.norm() {
                break;
            }
        }
        let branch = if z.im >= 0.0 {
            -I * PI
        } else {
            I * PI
        };
        return (-z).exp() / z * sum + branch;
    }
    // Continued fraction for |z| > 4 off the negative real axis:
    // E_1(z) = exp(-z) * 1 / (z + 1 - 1^2 / (z + 3 - 2^2 / (z + 5 - ...)))
    let tiny = 1e-30_f64;
    let mut b = z + 1.0;
    let mut c = Complex64::new(1.0 / tiny, 0.0);
    let mut d = if b.norm() < tiny {
        Complex64::new(1.0 / tiny, 0.0)
    } else {
        ONE / b
    };
    let mut h = d;
    for k in 1..400 {
        let kf = k as f64;
        let an = -kf * kf;
        b += 2.0;
        d = an * d + b;
        if d.norm() < tiny {
            d = Complex64::new(tiny, 0.0);
        }
        c = b + an / c;
        if c.norm() < tiny {
            c = Complex64::new(tiny, 0.0);
        }
        d = ONE / d;
        let del = d * c;
        h *= del;
        if (del - ONE).norm() < 1e-15 {
            break;
        }
    }
    (-z).exp() * h
}

/// Generalized complex exponential integral `E_n(z)`.
#[must_use]
pub fn exp_integral_en_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    if n == 0 {
        return (-z).exp() / z;
    }
    if z == ZERO {
        return if n == 1 {
            Complex64::new(f64::INFINITY, 0.0)
        } else {
            Complex64::new(1.0 / f64::from(n - 1), 0.0)
        };
    }
    if n == 1 {
        return exp_integral_e1_complex(z);
    }
    // Lentz continued fraction for |z| > n off (-inf, 0]:
    if z.norm() > f64::from(n) && !(z.re < 0.0 && z.im.abs() < 4.0) {
        let tiny = 1e-30_f64;
        let nf = f64::from(n);
        let mut b = z + nf;
        let mut c = Complex64::new(1.0 / tiny, 0.0);
        let mut d = if b.norm() < tiny {
            Complex64::new(1.0 / tiny, 0.0)
        } else {
            ONE / b
        };
        let mut h = d;
        for k in 1..400 {
            let kf = k as f64;
            let an = -kf * (nf + kf - 1.0);
            b += 2.0;
            d = an * d + b;
            if d.norm() < tiny {
                d = Complex64::new(tiny, 0.0);
            }
            c = b + an / c;
            if c.norm() < tiny {
                c = Complex64::new(tiny, 0.0);
            }
            d = ONE / d;
            let del = d * c;
            h *= del;
            if (del - ONE).norm() < 1e-15 {
                break;
            }
        }
        return (-z).exp() * h;
    }
    let mut val = exp_integral_e1_complex(z);
    let exp_mz = (-z).exp();
    for k in 1..n {
        val = (exp_mz - z * val) / f64::from(k);
    }
    val
}

/// Complex exponential integral `Ei(z) = -E_1(-z) + i * pi * sign` (principal branch).
#[must_use]
pub fn exp_integral_ei_complex(z: Complex64) -> Complex64 {
    if z == ZERO {
        return Complex64::new(f64::NEG_INFINITY, 0.0);
    }
    let e1 = exp_integral_e1_complex(-z);
    let shift = if z.im > 0.0 || (z.im == 0.0 && z.re > 0.0) {
        I * PI
    } else {
        -I * PI
    };
    -e1 + shift
}

/// Complex logarithmic integral `li(z) = Ei(ln z)`.
#[must_use]
pub fn log_integral_complex(z: Complex64) -> Complex64 {
    if z == ZERO {
        return ZERO;
    }
    exp_integral_ei_complex(z.ln())
}

/// Complex sine and cosine integrals `(Si(z), Ci(z))`.
#[must_use]
pub fn sin_cos_integral_complex(z: Complex64) -> (Complex64, Complex64) {
    if z == ZERO {
        return (ZERO, Complex64::new(f64::NEG_INFINITY, 0.0));
    }
    if z.norm() <= 4.0 {
        let z2 = z * z;
        let mut si_term = z;
        let mut si = si_term;
        for k in 1..150 {
            let kf = (2 * k + 1) as f64;
            si_term *= -z2 / ((2 * k) as f64 * kf);
            let add = si_term / kf;
            si += add;
            if add.norm() <= 1e-16 * (1.0 + si.norm()) {
                break;
            }
        }
        let mut ci_sum = ZERO;
        let mut ci_term = ONE;
        for k in 1..150 {
            let kf = (2 * k) as f64;
            ci_term *= -z2 / ((2 * k - 1) as f64 * kf);
            let add = ci_term / kf;
            ci_sum += add;
            if add.norm() <= 1e-16 * (1.0 + ci_sum.norm()) {
                break;
            }
        }
        let ci = EULER_GAMMA + z.ln() + ci_sum;
        return (si, ci);
    }
    // Use E_1(i z) and E_1(-i z) for |z| > 4:
    // Ci(z) = -0.5 * (E_1(i z) + E_1(-i z)) (+ 2*pi*i if Re(z)<0...)
    // Si(z) = pi/2 + 0.5i * (E_1(i z) - E_1(-i z)) (adjusting for Re(z) < 0)
    let e_plus = exp_integral_e1_complex(I * z);
    let e_minus = exp_integral_e1_complex(-I * z);
    let si_base = 0.5 * I * (e_plus - e_minus);
    let si = if z.re >= 0.0 {
        Complex64::new(FRAC_PI_2, 0.0) + si_base
    } else {
        Complex64::new(-FRAC_PI_2, 0.0) + si_base
    };
    let mut ci = -0.5 * (e_plus + e_minus);
    if z.re < 0.0 {
        ci += if z.im >= 0.0 {
            I * PI
        } else {
            -I * PI
        };
    }
    (si, ci)
}

/// Complex hyperbolic sine and cosine integrals `(Shi(z), Chi(z))`.
#[must_use]
pub fn sinh_cosh_integral_complex(z: Complex64) -> (Complex64, Complex64) {
    let (si, ci) = sin_cos_integral_complex(I * z);
    let shi = -I * si;
    let chi = ci
        - if z.im >= 0.0 {
            I * FRAC_PI_2
        } else {
            -I * FRAC_PI_2
        };
    (shi, chi)
}

/// Spence's dilogarithm `spence(z) = Li_2(1 - z)`.
#[must_use]
pub fn spence_complex(z: Complex64) -> Complex64 {
    polylog_complex(Complex64::new(2.0, 0.0), ONE - z)
}
