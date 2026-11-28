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
//!
//! Independent Padé order intervals can run in parallel. Failed intervals
//! reuse the closest earlier successful approximant. Candidate sums are
//! compared in increasing requested order.

use std::sync::OnceLock;

use num_bigint::BigInt;
use num_complex::Complex64;

use crate::rules::poly::galois::Fx;
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
    Some(Complex64::new(scale(&z.re)?, scale(&z.im)?))
}

/// Padé approximant `[M-1/M]` of `Σ b_k s^k`: numerator and denominator
/// coefficients (ascending, denominator constant term 1).
fn pade(
    b: &[Complex64],
    fixed: &[Option<Fx>],
    m: usize,
) -> Option<(Vec<Complex64>, Vec<Complex64>)> {
    if m == 0 || b.len() < 2 * m {
        return None;
    }
    // Σ_{j=0}^{M} q_j b_{L+i-j} = 0, i = 1..M, q_0 = 1.
    let mut matrix = Vec::with_capacity(m);
    let mut rhs = Vec::with_capacity(m);
    for i in m..2 * m {
        matrix.push(
            fixed
                .get(i - m..i)?
                .iter()
                .rev()
                .cloned()
                .collect::<Option<Vec<_>>>()?,
        );
        let Fx { re, im } = fixed.get(i)?.as_ref()?;
        rhs.push(vec![Fx { re: -re, im: -im }]);
    }
    let solution = solve_system(matrix, rhs)?;
    let mut q = Vec::with_capacity(m + 1);
    q.push(Complex64::new(1.0, 0.0));
    for row in &solution {
        q.push(from_fx(row.first()?)?);
    }
    let mut p = Vec::with_capacity(m);
    for i in 0..m {
        let mut value = Complex64::new(0.0, 0.0);
        for (qj, bj) in q.iter().zip(b[..=i].iter().rev()) {
            value += qj * *bj;
        }
        p.push(value);
    }
    Some((p, q))
}

fn horner(
    coefficients: &[Complex64],
    s: Complex64,
) -> Complex64 {
    coefficients
        .iter()
        .rev()
        .fold(Complex64::new(0.0, 0.0), |acc, c| acc * s + c)
}

/// Complex roots by the Durand–Kerner iteration.
pub(super) fn complex_roots(coefficients: &[Complex64]) -> Vec<Complex64> {
    let n = coefficients.len().saturating_sub(1);
    let Some(lead) = coefficients.last().copied().filter(|l| l.norm() > 0.0) else {
        return Vec::new();
    };
    if n == 0 {
        return Vec::new();
    }
    let monic: Vec<Complex64> = coefficients.iter().map(|c| c / lead).collect();
    let bound = 1.0
        + monic
            .iter()
            .take(n)
            .map(|c| c.norm())
            .fold(0.0_f64, f64::max);
    let mut z: Vec<Complex64> = (0..n)
        .map(|k| {
            Complex64::from_polar(
                bound,
                2.0 * std::f64::consts::PI * k as f64 / n as f64 + 0.4,
            )
        })
        .collect();
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

/// Gauss–Legendre nodes and weights on `[-1, 1]`, backed by the shared reuse pool.
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    let rule = crate::kernels::reuse::global_pool().gauss_legendre_rule(n);
    rule.0.iter().copied().zip(rule.1.iter().copied()).collect()
}

/// `∫_0^∞ e^{-s/v} P(s)/Q(s) ds` along the ray `arg s = alpha`.
fn laplace_along(
    p: &[Complex64],
    q: &[Complex64],
    v: Complex64,
    alpha: f64,
) -> Complex64 {
    static RULE: OnceLock<Vec<(f64, f64)>> = OnceLock::new();
    let direction = Complex64::from_polar(1.0, alpha);
    let decay = (direction / v).re;
    if decay <= 0.0 {
        return Complex64::new(f64::NAN, 0.0);
    }
    let r_max = 45.0 / decay;
    let rule = RULE.get_or_init(|| gauss_legendre(20));
    let panels = 120_u32;
    let mut total = Complex64::new(0.0, 0.0);
    for panel in 0..panels {
        // Finer near the origin.
        let (a, b) = (
            r_max * (f64::from(panel) / f64::from(panels)).powi(2),
            r_max * (f64::from(panel + 1) / f64::from(panels)).powi(2),
        );
        for (x, w) in rule {
            let r = 0.5 * (b - a) * x + f64::midpoint(a, b);
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
    let top = ((terms - 1) / 2).min(24);
    // b_k = c_{k+1} / k!, through the highest Padé system's final entry.
    let mut b = Vec::with_capacity(2 * top);
    let mut factorial = 1.0_f64;
    for k in 0..2 * top {
        if k > 0 {
            factorial *= k as f64;
        }
        b.push(coefficients.get(k + 1).copied()? / factorial);
    }
    // Every Padé solve reads the same Borel coefficients. Keep unavailable
    // fixed-point entries so lower orders still succeed on clean prefixes.
    let fixed: Vec<_> = b.iter().copied().map(to_fx).collect();
    let orders: Vec<usize> = (2..=top)
        .step_by(2)
        .chain((!top.is_multiple_of(2)).then_some(top))
        .collect();
    // Search each degree at most once, within the interval for its requested
    // order. Each successful approximant needs just one contour evaluation.
    let sums = crate::kernels::par::map_indexed(orders.len(), 4, |i| {
        let lower = i.checked_sub(1).map_or(0, |j| orders[j]);
        (lower + 1..=orders[i])
            .rev()
            .find_map(|m| pade(&b, &fixed, m))
            .map(|(p, q)| sum_at_approximant(&p, &q, c0, v))
    });
    // Compare sums in requested order, including repeated approximants after
    // a singular interval. These repetitions still affect the error estimate.
    let mut previous: Option<BorelSum> = None;
    let mut best: Option<BorelSum> = None;
    for sum in sums {
        let Some(current) = sum.or(previous) else {
            continue;
        };
        if let Some(before) = previous {
            let change = (current.value - before.value).norm().max(current.error);
            if best.is_none_or(|b| change < b.error) {
                best = Some(BorelSum {
                    error: change,
                    ..current
                });
            }
        }
        previous = Some(current);
    }
    best.or(previous)
}

fn wrap_angle(angle: f64) -> f64 {
    (angle + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
}

/// The sum with the given Padé approximant. Away from the poles it integrates
/// along the original ray. Near a pole it rotates the ray without crossing
/// any poles. On a Stokes ray it takes the median of the two lateral sums,
/// using an offset adapted to the gap to the next pole.
fn sum_at_approximant(
    p: &[Complex64],
    q: &[Complex64],
    c0: Complex64,
    v: Complex64,
) -> BorelSum {
    let poles = complex_roots(q);
    let arg_v = v.arg();
    let relative: Vec<f64> = poles
        .iter()
        .filter(|z| z.norm() > 0.0)
        .map(|z| wrap_angle(z.arg() - arg_v))
        .collect();
    let nearest = relative
        .iter()
        .copied()
        .min_by(|a, b| a.abs().total_cmp(&b.abs()));
    let (value, jump, stokes) = match nearest {
        | Some(delta) if delta.abs() < 1e-3 => {
            // Poles on the ray: lateral sums, offset up to half the gap to the others.
            let gap = relative
                .iter()
                .map(|d| d.abs())
                .filter(|d| *d > 1e-3)
                .fold(f64::INFINITY, f64::min);
            let offset = (0.5 * gap).clamp(0.05, 0.25);
            let (up, down) = (
                laplace_along(p, q, v, arg_v + offset),
                laplace_along(p, q, v, arg_v - offset),
            );
            ((up + down) / 2.0, (up - down).norm() / 2.0, true)
        },
        | Some(delta) if delta.abs() < 0.15 => {
            // Rotate away from the near pole, staying clear of the others.
            let away = -delta.signum();
            let room = relative
                .iter()
                .filter(|d| (d.signum() - away).abs() < 0.5 && d.abs() > 1e-9)
                .map(|d| d.abs())
                .fold(0.4, f64::min);
            let rotation = away * (0.6 * room).clamp(0.02, 0.2);
            (laplace_along(p, q, v, arg_v + rotation), 0.0, false)
        },
        | _ => (laplace_along(p, q, v, arg_v), 0.0, false),
    };
    BorelSum {
        value: c0 + value,
        error: jump,
        stokes,
    }
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
        assert!(
            (sum.value.re - 0.915_633_339_397_88).abs() < 1e-8,
            "{sum:?}"
        );
        assert!(!sum.stokes);
        // The plain partial sum cannot do better than the smallest term, 10! 10^-10 ~ 4e-4.
        let truncated: f64 = (0..10)
            .map(|k| (-1.0_f64).powi(k) * (1..=k).map(f64::from).product::<f64>() * 0.1_f64.powi(k))
            .sum();
        assert!((truncated - 0.915_633_339_397_88).abs() > 1e-5);
    }

    #[test]
    fn sparse_exponents_reuse_the_reduced_borel_series() {
        let reduced = euler(16, -1.0);
        let mut sparse = vec![Complex64::default(); reduced.len() * 3];
        for (target, coefficient) in sparse.iter_mut().step_by(3).zip(&reduced) {
            *target = *coefficient;
        }
        let v = Complex64::new(0.45, 0.15);
        let expected = borel_pade(&reduced, v.powi(3)).expect("reduced sum");
        let sum = borel_pade(&sparse, v).expect("sparse sum");
        assert_eq!(sum.value, expected.value);
        assert_eq!(sum.error.to_bits(), expected.error.to_bits());
        assert_eq!(sum.stokes, expected.stokes);
    }

    #[test]
    fn capped_borel_orders_keep_the_same_result_across_worker_counts() {
        let serial = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let parallel = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let v = Complex64::new(0.1, 0.0);
        let mut coefficients = euler(49, -1.0);
        let expected = serial
            .install(|| borel_pade(&coefficients, v))
            .expect("capped sum");
        coefficients.extend([Complex64::new(f64::INFINITY, 0.0); 100]);
        let actual = parallel
            .install(|| borel_pade(&coefficients, v))
            .expect("parallel sum");
        assert_eq!(actual.value, expected.value);
        assert_eq!(actual.error.to_bits(), expected.error.to_bits());
        assert_eq!(actual.stokes, expected.stokes);
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

    #[test]
    fn near_stokes_pole_does_not_need_the_median() {
        // c_k = k! e^{-ik δ}, δ = 0.1: the Borel pole is 0.1 rad off the ray of v = 0.05,
        // so the Borel sum is the integral ∫_0^∞ e^{-t}/(1 - z t) dt with z = v e^{-iδ}
        // along the real axis (the pole at 1/z lies above it).
        let delta = 0.1_f64;
        let coefficients: Vec<Complex64> = euler(40, 1.0)
            .iter()
            .enumerate()
            .map(|(k, c)| c * Complex64::from_polar(1.0, -(k as f64) * delta))
            .collect();
        let sum = borel_pade(&coefficients, Complex64::new(0.05, 0.0)).expect("sum");
        assert!(!sum.stokes, "{sum:?}");
        // Reference by composite Gauss-Legendre of the Euler integral.
        let z = Complex64::from_polar(0.05, -delta);
        let rule = gauss_legendre(40);
        let mut reference = Complex64::new(0.0, 0.0);
        for panel in 0..200 {
            let (a, b) = (f64::from(panel) * 0.25, f64::from(panel + 1) * 0.25);
            for (x, w) in &rule {
                let t = 0.5 * (b - a) * x + f64::midpoint(a, b);
                reference += (-t).exp() / (Complex64::new(1.0, 0.0) - z * t) * (0.5 * (b - a) * w);
            }
        }
        assert!(
            (sum.value - reference).norm() < 1e-6,
            "{sum:?} vs {reference}"
        );
    }
}
