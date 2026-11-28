//! Complex Dirichlet eta, Riemann zeta (full analytic continuation),
//! Hurwitz zeta, Lerch transcendent, Clausen functions, and Bernoulli/Euler
//! polynomials.

use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::quadrature::exp_sinh;
use crate::kernels::special_complex::gamma_complex;
use crate::kernels::special_complex::ln_gamma_complex;
use crate::kernels::special_ext::polylog_complex;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

/// Complex Dirichlet eta function `eta(s) = sum_{n>=1} (-1)^(n-1) / n^s`
/// using Borwein's accelerated alternating sum for small `|Im(s)|` and
/// `(1 - 2^(1-s)) * zeta(s)` elsewhere.
#[must_use]
pub fn dirichlet_eta_complex(s: Complex64) -> Complex64 {
    // Borwein weights d_k = M * sum_{i=0}^k (M + i - 1)! * 4^i / ((M - i)! * (2i)!)
    // With M = 20, every d_k fits EXACTLY in f64 (< 2^53) and truncation error is ~ 1 / 5.8^20 ~ 1e-16!
    const M: usize = 20;

    if (s - ONE).norm() < 1e-14 {
        return Complex64::new(std::f64::consts::LN_2, 0.0);
    }
    if s.im.abs() > 10.0 || s.re < -10.0 {
        return (ONE - Complex64::new(2.0, 0.0).powc(ONE - s)) * riemann_zeta_complex(s);
    }
    let mut d = [0.0_f64; M + 1];
    let mut acc = 0.0_f64;
    // Term i=0: (M-1)! / M! = 1/M, * M = 1.
    let mut term = 1.0_f64;
    acc += term;
    d[0] = acc;
    for i in 1..=M {
        let m_f = M as f64;
        let i_f = i as f64;
        term *= (m_f + i_f - 1.0) * (m_f - i_f + 1.0) * 4.0 / ((2.0 * i_f - 1.0) * (2.0 * i_f));
        acc += term;
        if let Some(slot) = d.get_mut(i) {
            *slot = acc;
        }
    }
    let d_m = d[M];
    let mut sum = ZERO;
    for k in 0..M {
        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
        let weight = d_m - d.get(k).copied().unwrap_or(0.0);
        let k1 = Complex64::new((k + 1) as f64, 0.0);
        sum += sign * weight / k1.powc(s);
    }
    sum / d_m
}

/// Complex Riemann zeta function `zeta(s)` over the entire complex plane.
///
///
/// Uses `eta(s) / (1 - 2^(1 - s))` for `Re(s) >= 0.5`, and the functional
/// equation `zeta(s) = 2^s * pi^(s - 1) * sin(pi * s / 2) * Gamma(1 - s) * zeta(1 - s)`
/// for `Re(s) < 0.5`.
#[must_use]
pub fn riemann_zeta_complex(s: Complex64) -> Complex64 {
    if (s - ONE).norm() < 1e-15 {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    if s.re >= 0.5 {
        if s.im.abs() > 10.0 {
            return hurwitz_zeta_complex(s, ONE);
        }
        let denom = ONE - Complex64::new(2.0, 0.0).powc(ONE - s);
        if denom.norm() < 1e-14 {
            return hurwitz_zeta_complex(s, ONE);
        }
        return dirichlet_eta_complex(s) / denom;
    }
    if s.im.abs() > 20.0 || s.re < -100.0 {
        let ln2 = std::f64::consts::LN_2;
        let ln_pi = PI.ln();
        let half_s = 0.5 * s;
        let ln_sin = if half_s.im.abs() > 20.0 {
            let y = half_s.im;
            if y > 0.0 {
                let t = Complex64::new(-2.0 * PI * y, 2.0 * PI * half_s.re);
                Complex64::new(PI * y - ln2, -PI * half_s.re + 0.5 * PI) + (ONE - t.exp()).ln()
            } else {
                let t = Complex64::new(2.0 * PI * y, -2.0 * PI * half_s.re);
                Complex64::new(-PI * y - ln2, PI * half_s.re - 0.5 * PI) + (ONE - t.exp()).ln()
            }
        } else {
            (PI * half_s).sin().ln()
        };
        let z_refl = riemann_zeta_complex(ONE - s);
        return (s * ln2 + (s - ONE) * ln_pi + ln_sin + ln_gamma_complex(ONE - s) + z_refl.ln())
            .exp();
    }
    let two_s = Complex64::new(2.0, 0.0).powc(s);
    let pi_s_m1 = Complex64::new(PI, 0.0).powc(s - ONE);
    let sin_term = (0.5 * PI * s).sin();
    let g_term = gamma_complex(ONE - s);
    let z_refl = riemann_zeta_complex(ONE - s);
    two_s * pi_s_m1 * sin_term * g_term * z_refl
}

/// Complex Hurwitz zeta function `zeta(s, a) = sum_{n>=0} (n + a)^(-s)`
/// via adaptive Euler-Maclaurin summation.
#[must_use]
pub fn hurwitz_zeta_complex(
    s: Complex64,
    mut a: Complex64,
) -> Complex64 {
    // B_{2k} / (2k)! for k = 1..=8
    const B_OVER_FACT: [f64; 8] = [
        1.0 / 12.0,
        -1.0 / 720.0,
        1.0 / 30_240.0,
        -1.0 / 1_209_600.0,
        1.0 / 47_900_160.0,
        -691.0 / 1_307_674_368_000.0,
        1.0 / 74_724_249_600.0,
        -3617.0 / 10_670_622_842_880_000.0,
    ];

    if (s - ONE).norm() < 1e-14 {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    if a.im == 0.0 && a.re <= 0.0 && (a.re - a.re.round()).abs() == 0.0 {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    let mut pre_sum = ZERO;
    while a.re <= 0.5 {
        pre_sum += a.powc(-s);
        a += ONE;
    }
    let n_steps = (16 + (0.7 * s.im.abs()).ceil() as usize).clamp(16, 512);
    let mut sum = pre_sum;
    for k in 0..n_steps {
        let ak = a + (k as f64);
        sum += ak.powc(-s);
    }
    let an = a + (n_steps as f64);
    let an_ms = an.powc(-s);
    let int_term = an.powc(ONE - s) / (s - ONE);
    let half_term = 0.5 * an_ms;

    let inv_an2 = ONE / (an * an);
    let mut poch = s;
    let mut pow_term = ONE / an;
    let mut em_sum = ZERO;
    for (idx, &coeff) in B_OVER_FACT.iter().enumerate() {
        if idx > 0 {
            let two_k = (2 * idx) as f64;
            poch *= (s + (two_k - 1.0)) * (s + two_k);
            pow_term *= inv_an2;
        }
        em_sum += coeff * poch * (an_ms * pow_term);
    }
    sum + int_term + half_term + em_sum
}

/// Complex Lerch transcendent `Phi(z, s, a) = sum_{n>=0} z^n / (n + a)^s`.
///
/// For positive real parts of `s` and `a`, the Euler integral continues the
/// series outside the unit disk. On the real cut `z > 1`, the value is the
/// limit from below, as for [`polylog_complex`].
#[must_use]
pub fn lerch_phi_complex(
    z: Complex64,
    s: Complex64,
    a: Complex64,
) -> Complex64 {
    if z == ONE {
        return hurwitz_zeta_complex(s, a);
    }
    if s == ZERO {
        return ONE / (ONE - z);
    }
    if z == ZERO {
        return a.powc(-s);
    }
    if z.im == 0.0 && s.im == 0.0 && a.im == 0.0 && z.re < 1.0 && a.re > 0.0 {
        let real = crate::kernels::special_ext::lerch(z.re, s.re, a.re);
        if !real.is_nan() {
            return Complex64::new(real, 0.0);
        }
    }
    if z.norm() < 0.9 {
        let mut sum = ZERO;
        let mut z_pow = ONE;
        for n in 0..400 {
            let term = z_pow * (a + f64::from(n)).powc(-s);
            sum += term;
            if term.norm() <= 1e-15 * (1.0 + sum.norm()) {
                break;
            }
            z_pow *= z;
        }
        return sum;
    }
    if a == ONE {
        let reduced = polylog_complex(s, z) / z;
        if reduced.is_finite() {
            return reduced;
        }
    }
    if s.re > 0.0 && a.re > 0.0 && z.is_finite() && s.is_finite() && a.is_finite() {
        let integral = lerch_euler_complex(z, s, a);
        if integral.is_finite() {
            return integral;
        }
    }
    // Euler series transformation for |z| >= 0.9 off [1, +inf):
    // sum_{n>=0} z^n u_n = (1/(1-z)) sum_{k>=0} (z/(z-1))^k (-1)^k Delta^k u_0
    // converges whenever |z/(z-1)| < 1 (i.e. Re(z) < 0.5).
    if z.re < 0.45 {
        const M: usize = 36;
        let mut diff = Vec::with_capacity(M);
        for n in 0..M {
            diff.push((a + (n as f64)).powc(-s));
        }
        let w = z / (z - ONE);
        let mut w_pow = ONE;
        let mut sum = ZERO;
        for _k in 0..M {
            let d0 = diff.first().copied().unwrap_or(ZERO);
            let add = w_pow * d0;
            sum += add;
            if add.norm() <= 1e-15 * (1.0 + sum.norm()) {
                break;
            }
            w_pow *= -w;
            for j in 0..(diff.len() - 1) {
                let next = diff.get(j + 1).copied().unwrap_or(ZERO);
                if let Some(cur) = diff.get_mut(j) {
                    *cur = next - *cur;
                }
            }
            diff.pop();
        }
        return sum / (ONE - z);
    }
    // The remaining orders and shifts use Wynn acceleration on the partial sums.
    let mut partials = Vec::with_capacity(48);
    let mut sum = ZERO;
    let mut z_pow = ONE;
    for n in 0..48 {
        sum += z_pow * (a + f64::from(n)).powc(-s);
        partials.push(sum);
        z_pow *= z;
    }
    let acc = crate::kernels::analytic_continuation::wynn_epsilon_complex(&partials);
    if acc.is_finite() { acc } else { sum }
}

// The Euler integral initially follows a unit-slope ray on the side opposite
// the pole log(z). The contour then approaches a horizontal line, keeping the
// weight bounded for complex shifts with positive real part. On the positive
// real cut the upper contour selects the limit from below.
// Euler representation: NIST DLMF 25.14.5.
#[allow(clippy::single_call_fn)]
fn lerch_euler_complex(
    z: Complex64,
    s: Complex64,
    a: Complex64,
) -> Complex64 {
    let direction = if z.re > 0.0 {
        if z.im > 0.0 { -1.0 } else { 1.0 }
    } else {
        0.0
    };
    // Integration by parts raises the endpoint power by one. It also avoids
    // the slow convergence of the original Euler integral for small Re(s).
    let by_parts = s.re < 1.0;
    let (order, log_gamma) = if by_parts {
        (s, ln_gamma_complex(s + ONE))
    } else {
        (s - ONE, ln_gamma_complex(s))
    };
    let integrand = |x: f64| {
        if x <= 0.0 {
            return ZERO;
        }
        let inv = (1.0 + x).recip();
        let y = direction * (x * inv);
        let t = Complex64::new(x, y);
        let log_weight = order * t.ln() - a * t - log_gamma;
        if log_weight.re < -745.0 {
            return ZERO;
        }
        let e = (-x).exp();
        let (sin, cos) = y.sin_cos();
        let q = z * Complex64::new(e * cos, -e * sin);
        let denominator = if x < 0.5 {
            // Preserve the small difference from 1 at the endpoint.
            let half_sin = (0.5 * y).sin();
            let em1 = Complex64::new((-x).exp_m1() - 2.0 * e * half_sin * half_sin, -e * sin);
            (ONE - z) - z * em1
        } else {
            ONE - q
        };
        let dt = Complex64::new(1.0, direction * inv * inv);
        let weight = log_weight.exp() * dt;
        let scale = denominator.re.abs().max(denominator.im.abs());
        // Complex division squares the denominator norm. Reduce it first
        // when q would overflow that intermediate.
        let (weight, q, denominator) = if scale > 1e150 {
            (
                weight.unscale(scale),
                q.unscale(scale),
                denominator.unscale(scale),
            )
        } else {
            (weight, q, denominator)
        };
        let value = weight / denominator;
        if by_parts {
            value * (a + q / denominator)
        } else {
            value
        }
    };
    Complex64::new(
        exp_sinh(|x| integrand(x).re, 0.0, 2e-14).value,
        exp_sinh(|x| integrand(x).im, 0.0, 2e-14).value,
    )
}

/// Complex Clausen function: `Im(Li_n(e^(i z)))` for even `n`, `Re(Li_n(e^(i z)))` for odd `n`.
#[must_use]
pub fn clausen_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    let eiz = (I * z).exp();
    let e_miz = (-I * z).exp();
    let s_c = Complex64::new(f64::from(n), 0.0);
    let lp = polylog_complex(s_c, eiz);
    let lm = polylog_complex(s_c, e_miz);
    if n.is_multiple_of(2) {
        -0.5 * I * (lp - lm)
    } else {
        0.5 * (lp + lm)
    }
}

/// Complex Bernoulli polynomial `B_n(z)`.
#[must_use]
pub fn bernoulli_poly_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    // Use recurrence / explicit binomial sum with Bernoulli numbers B_k.
    let n_usize = n as usize;
    let mut b = vec![0.0_f64; n_usize + 1];
    b[0] = 1.0;
    for m in 1..=n_usize {
        let mut s = 0.0_f64;
        let mut binom = 1.0_f64;
        for k in 0..m {
            s += binom * b.get(k).copied().unwrap_or(0.0) / ((m - k + 1) as f64);
            binom *= ((m - k) as f64) / ((k + 1) as f64);
        }
        if let Some(slot) = b.get_mut(m) {
            *slot = -s;
        }
    }
    let mut sum = ZERO;
    let mut binom = 1.0_f64;
    for k in 0..=n_usize {
        let bk = b.get(k).copied().unwrap_or(0.0);
        let power = i32::try_from(n_usize - k).unwrap_or(0);
        sum += binom * bk * z.powi(power);
        if k < n_usize {
            binom *= ((n_usize - k) as f64) / ((k + 1) as f64);
        }
    }
    sum
}

/// Complex Euler polynomial `E_n(z) = (2 / (n + 1)) * (B_{n+1}(z) - 2^(n+1) B_{n+1}(z / 2))`.
#[must_use]
pub fn euler_poly_complex(
    n: u32,
    z: Complex64,
) -> Complex64 {
    let np1 = n.saturating_add(1);
    let scale = 2.0_f64.powi(i32::try_from(np1).unwrap_or(1));
    let b1 = bernoulli_poly_complex(np1, z);
    let b2 = bernoulli_poly_complex(np1, 0.5 * z);
    (2.0 / f64::from(np1)) * (b1 - scale * b2)
}
