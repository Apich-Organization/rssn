//! Complex Dirichlet eta, Riemann zeta (full analytic continuation),
//! Hurwitz zeta, Lerch transcendent, Clausen functions, and Bernoulli/Euler
//! polynomials.

use std::f64::consts::PI;

use num_complex::Complex64;

use crate::kernels::special_complex::gamma_complex;
use crate::kernels::special_ext::polylog_complex;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);

/// Complex Dirichlet eta function `eta(s) = sum_{n>=1} (-1)^(n-1) / n^s`
/// using Borwein's accelerated alternating sum (converges rapidly for `Re(s) > -20`).
#[must_use]
pub fn dirichlet_eta_complex(s: Complex64) -> Complex64 {
    // Borwein weights d_k = M * sum_{i=0}^k (M + i - 1)! * 4^i / ((M - i)! * (2i)!)
    // With M = 20, every d_k fits EXACTLY in f64 (< 2^53) and truncation error is ~ 1 / 5.8^20 ~ 1e-16!
    const M: usize = 20;
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
    let two_s = Complex64::new(2.0, 0.0).powc(s);
    let pi_s_m1 = Complex64::new(PI, 0.0).powc(s - ONE);
    let sin_term = (0.5 * PI * s).sin();
    let g_term = gamma_complex(ONE - s);
    let z_refl = riemann_zeta_complex(ONE - s);
    two_s * pi_s_m1 * sin_term * g_term * z_refl
}

/// Complex Hurwitz zeta function `zeta(s, a) = sum_{n>=0} (n + a)^(-s)`
/// via Euler-Maclaurin summation.
#[must_use]
pub fn hurwitz_zeta_complex(
    s: Complex64,
    a: Complex64,
) -> Complex64 {
    if (s - ONE).norm() < 1e-14 {
        return Complex64::new(f64::INFINITY, 0.0);
    }
    const N: usize = 16;
    let mut sum = ZERO;
    for k in 0..N {
        let ak = a + (k as f64);
        sum += ak.powc(-s);
    }
    let an = a + (N as f64);
    let an_ms = an.powc(-s);
    let int_term = an.powc(ONE - s) / (s - ONE);
    let half_term = 0.5 * an_ms;

    // B_{2k} / (2k)! for k = 1..=6
    const B_OVER_FACT: [f64; 6] = [
        1.0 / 12.0,
        -1.0 / 720.0,
        1.0 / 30_240.0,
        -1.0 / 1_209_600.0,
        1.0 / 47_900_160.0,
        -691.0 / 1_307_674_368_000.0,
    ];
    let inv_an2 = ONE / (an * an);
    let mut poch = s; // s * (s+1) * ... * (s + 2k - 2)
    let mut pow_term = ONE / an; // an^(-s - (2k - 1)) = an_ms * an^(-(2k-1))
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
#[must_use]
pub fn lerch_phi_complex(
    z: Complex64,
    s: Complex64,
    a: Complex64,
) -> Complex64 {
    if z == ZERO {
        return a.powc(-s);
    }
    if z == ONE {
        return hurwitz_zeta_complex(s, a);
    }
    if z.norm() < 0.9 {
        let mut sum = ZERO;
        let mut z_pow = ONE;
        for n in 0..400 {
            let term = z_pow * (a + (n as f64)).powc(-s);
            sum += term;
            if term.norm() <= 1e-15 * (1.0 + sum.norm()) {
                break;
            }
            z_pow *= z;
        }
        return sum;
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
    // Fallback: split first few terms and use z-series up to convergence limit.
    let mut sum = ZERO;
    let mut z_pow = ONE;
    for n in 0..800 {
        let term = z_pow * (a + (n as f64)).powc(-s);
        sum += term;
        if term.norm() <= 1e-14 * (1.0 + sum.norm()) {
            break;
        }
        z_pow *= z;
    }
    sum
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
    if n % 2 == 0 {
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
