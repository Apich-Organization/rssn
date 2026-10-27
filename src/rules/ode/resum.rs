//! Borel–Padé resummation of formal power series, with lateral (median)
//! summation across Stokes directions.
//!
//! For `S(v) = Σ c_k v^k` the Borel transform is `B(s) = Σ_{k ≥ 1} c_k
//! s^{k-1} / (k-1)!` and `S(v) = c_0 + ∫_0^∞ e^{-s/v} B(s) ds` (along a ray
//! in the half plane where the exponential decays). Here `B` is continued
//! analytically by its Padé approximant `[M-1/M]` (the linear system is
//! solved in 640-bit fixed point, which keeps the poles accurate), the poles
//! of the approximant locate the Borel singularities, and the integral is
//! computed by composite Gauss–Legendre quadrature along a ray that avoids
//! them. When a pole lies near the ray (a Stokes direction, as for series
//! with non-alternating factorial growth) the contour is deformed to the two
//! sides and the median `(S₊ + S₋)/2` is returned, with half the jump as
//! part of the error estimate.

use num_bigint::BigInt;
use num_complex::Complex64;

use crate::rules::poly::galois::Fx;
use crate::rules::poly::galois::BITS;
use crate::rules::poly::galois_high::solve_system;

/// The result of a resummation.
#[derive(Clone, Copy, Debug)]
pub struct BorelSum {
    /// The sum.
    pub value: Complex64,
    /// An estimate of the error (Padé order change, and the lateral jump).
    pub error: f64,
    /// Whether a Borel singularity was near the integration ray (the value
    /// is then the median of the two lateral sums).
    pub stokes: bool,
}

fn to_fx(z: Complex64) -> Option<Fx> {
    Fx::from_f64(z.re, z.im)
}

fn from_fx(z: &Fx) -> Option<Complex64> {
    let scale = |v: &BigInt| -> Option<f64> {
        let q = num_rational::BigRational::new(v.clone(), crate::rules::poly::galois::one());
        num_traits::ToPrimitive::to_f64(&q)
    };
    let _ = BITS;
    Some(Complex64::new(scale(&z.re)?, scale(&z.im)?))
}

/// Padé approximant `[M-1/M]` of `Σ b_k s^k`: numerator and denominator
/// coefficients (ascending, denominator constant term 1).
fn pade(
    b: &[Complex64],
    m: usize,
) -> Option<(Vec<Complex64>, Vec<Complex64>)> {
    if m == 0 || b.len() < 2 * m {
        return None;
    }
    let l = m - 1;
    // Σ_{j=0}^{M} q_j b_{L+i-j} = 0, i = 1..M, q_0 = 1.
    let get = |index: i64| -> Complex64 { usize::try_from(index).ok().and_then(|i| b.get(i)).copied().unwrap_or_else(|| Complex64::new(0.0, 0.0)) };
    let mut matrix: Vec<Vec<Fx>> = Vec::new();
    let mut rhs: Vec<Vec<Fx>> = Vec::new();
    for i in 1..=m {
        let mut row = Vec::new();
        for j in 1..=m {
            row.push(to_fx(get(i64::try_from(l + i).ok()? - i64::try_from(j).ok()?))?);
        }
        matrix.push(row);
        rhs.push(vec![to_fx(-get(i64::try_from(l + i).ok()?))?]);
    }
    let solution = solve_system(matrix, rhs)?;
    let mut q = vec![Complex64::new(1.0, 0.0)];
    for row in &solution {
        q.push(from_fx(row.first()?)?);
    }
    let mut p = Vec::new();
    for i in 0..=l {
        let mut value = Complex64::new(0.0, 0.0);
        for (j, qj) in q.iter().enumerate() {
            if j <= i {
                value += qj * get(i64::try_from(i - j).ok()?);
            }
        }
        p.push(value);
    }
    Some((p, q))
}

fn horner(
    coefficients: &[Complex64],
    s: Complex64,
) -> Complex64 {
    coefficients.iter().rev().fold(Complex64::new(0.0, 0.0), |acc, c| acc * s + c)
}

/// Complex roots by the Durand–Kerner iteration.
pub(super) fn complex_roots(coefficients: &[Complex64]) -> Vec<Complex64> {
    let n = coefficients.len().saturating_sub(1);
    let Some(lead) = coefficients.last().copied().filter(|l| l.norm() > 0.0) else { return Vec::new() };
    if n == 0 {
        return Vec::new();
    }
    let monic: Vec<Complex64> = coefficients.iter().map(|c| c / lead).collect();
    let bound = 1.0 + monic.iter().take(n).map(|c| c.norm()).fold(0.0_f64, f64::max);
    let mut z: Vec<Complex64> = (0..n).map(|k| Complex64::from_polar(bound, 2.0 * std::f64::consts::PI * k as f64 / n as f64 + 0.4)).collect();
    for _ in 0..500 {
        let mut moved = 0.0_f64;
        for i in 0..n {
            let zi = z.get(i).copied().unwrap_or_default();
            let mut denominator = Complex64::new(1.0, 0.0);
            for (j, zj) in z.iter().enumerate() {
                if j != i {
                    denominator *= zi - zj;
                }
            }
            if denominator.norm() == 0.0 {
                continue;
            }
            let step = horner(&monic, zi) / denominator;
            if let Some(slot) = z.get_mut(i) {
                *slot = zi - step;
            }
            moved = moved.max(step.norm());
        }
        if moved < 1e-14 * bound {
            break;
        }
    }
    z
}

/// Gauss–Legendre nodes and weights on `[-1, 1]` (Newton on `P_n`).
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for i in 0..n {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut derivative = 1.0;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0, x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            derivative = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let step = p1 / derivative;
            x -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        out.push((x, 2.0 / ((1.0 - x * x) * derivative * derivative)));
    }
    out
}

/// `∫_0^∞ e^{-s/v} P(s)/Q(s) ds` along the ray `arg s = alpha`.
fn laplace_along(
    p: &[Complex64],
    q: &[Complex64],
    v: Complex64,
    alpha: f64,
) -> Complex64 {
    let direction = Complex64::from_polar(1.0, alpha);
    let decay = (direction / v).re;
    if decay <= 0.0 {
        return Complex64::new(f64::NAN, 0.0);
    }
    let r_max = 45.0 / decay;
    let rule = gauss_legendre(20);
    let panels = 120;
    let mut total = Complex64::new(0.0, 0.0);
    for panel in 0..panels {
        // Finer near the origin.
        let (a, b) = (r_max * (panel as f64 / panels as f64).powi(2), r_max * ((panel + 1) as f64 / panels as f64).powi(2));
        for (x, w) in &rule {
            let r = 0.5 * (b - a) * x + 0.5 * (b + a);
            let s = direction * r;
            let value = horner(p, s) / horner(q, s) * (-s / v).exp();
            total += value * (0.5 * (b - a) * w) * direction;
        }
    }
    total
}

/// Borel–Padé sum of `Σ c_k v^k` (see the module documentation).
#[must_use]
pub fn borel_pade(
    coefficients: &[Complex64],
    v: Complex64,
) -> Option<BorelSum> {
    // A series in v^g (gaps in the exponents) is summed in w = v^g.
    let stride = coefficients
        .iter()
        .enumerate()
        .skip(1)
        .filter(|(_, c)| c.norm() > 0.0)
        .map(|(k, _)| k)
        .fold(0_usize, |acc, k| {
            let (mut a, mut b) = (acc, k);
            while b != 0 {
                (a, b) = (b, a % b);
            }
            a
        });
    if stride > 1 {
        let reduced: Vec<Complex64> = coefficients.iter().step_by(stride).copied().collect();
        return borel_pade(&reduced, v.powi(i32::try_from(stride).ok()?));
    }
    let terms = coefficients.len();
    if terms < 5 {
        return None;
    }
    let c0 = coefficients.first().copied()?;
    // b_k = c_{k+1} / k!
    let mut b = Vec::new();
    let mut factorial = 1.0_f64;
    for k in 0..terms - 1 {
        if k > 0 {
            factorial *= k as f64;
        }
        b.push(coefficients.get(k + 1).copied()? / factorial);
    }
    let m = ((terms - 1) / 2).min(24);
    let sum_with = |m: usize| -> Option<(BorelSum, Vec<Complex64>)> {
        // An exactly rational Borel transform makes the system singular for large
        // orders: use the largest order that solves.
        let (p, q) = (1..=m).rev().find_map(|order| pade(&b, order))?;
        let poles = complex_roots(&q);
        let arg_v = v.arg();
        let near = poles.iter().filter(|z| z.norm() > 0.0).map(|z| ((z.arg() - arg_v + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI).abs()).fold(f64::INFINITY, f64::min);
        let (value, jump, stokes) = if near < 0.15 {
            let offset = 0.25;
            let (up, down) = (laplace_along(&p, &q, v, arg_v + offset), laplace_along(&p, &q, v, arg_v - offset));
            ((up + down) / 2.0, (up - down).norm() / 2.0, true)
        } else {
            (laplace_along(&p, &q, v, arg_v), 0.0, false)
        };
        Some((BorelSum { value: c0 + value, error: jump, stokes }, poles))
    };
    let (best, _) = sum_with(m)?;
    let error = match sum_with(m.saturating_sub(2).max(1)) {
        | Some((coarser, _)) => (best.value - coarser.value).norm().max(best.error),
        | None => best.error,
    };
    Some(BorelSum { error, ..best })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn euler(
        n: usize,
        sign: f64,
    ) -> Vec<Complex64> {
        let mut c = Vec::new();
        let mut factorial = 1.0_f64;
        for k in 0..n {
            if k > 0 {
                factorial *= k as f64;
            }
            c.push(Complex64::new(sign.powi(k as i32) * factorial, 0.0));
        }
        c
    }

    #[test]
    fn alternating_euler_series() {
        // Σ (-1)^k k! x^k = (1/x) e^{1/x} E1(1/x); at x = 0.1: 0.915633339...
        let sum = borel_pade(&euler(40, -1.0), Complex64::new(0.1, 0.0)).expect("sum");
        assert!((sum.value.re - 0.915_633_339_397_88).abs() < 1e-8, "{sum:?}");
        assert!(!sum.stokes);
        // The plain partial sum cannot do better than the smallest term, 10! 10^-10 ~ 4e-4.
        let truncated: f64 = (0..10).map(|k| (-1.0_f64).powi(k) * (1..=k).map(f64::from).product::<f64>() * 0.1_f64.powi(k)).sum();
        assert!((truncated - 0.915_633_339_397_88).abs() > 1e-5);
    }

    #[test]
    fn stokes_direction_gives_the_median() {
        // Σ k! x^k on the positive axis: the Borel transform 1/(1 - s) has its pole
        // on the integration ray; the median is the principal value
        // (1/x) e^{-1/x} Ei(1/x), at x = 0.05: 1.0559559...
        let sum = borel_pade(&euler(40, 1.0), Complex64::new(0.05, 0.0)).expect("sum");
        assert!(sum.stokes, "{sum:?}");
        assert!((sum.value.re - 1.055_955_905_593).abs() < 1e-6, "{sum:?}");
        // The imaginary part cancels in the median and the error estimate carries the jump.
        assert!(sum.value.im.abs() < 1e-6, "{sum:?}");
        assert!(sum.error > 1e-10 && sum.error < 1e-5, "{sum:?}");
    }
}
