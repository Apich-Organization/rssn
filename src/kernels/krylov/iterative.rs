//! Scale-safe linear iterations. The public wrapper solves a normalized
//! residual correction and verifies convergence against the original system.
use super::IterResult;

// A signed mantissa and an unbounded (for f64 products) binary exponent.
// Dot products and coefficient ratios stay here until applied to a vector;
// a harmless underflowing r'z never masquerades as Krylov breakdown.
#[derive(Clone, Copy, Debug)]
struct Scaled {
    fraction: f64,
    exponent: i32,
}
impl Scaled {
    fn from(value: f64) -> Self {
        if value == 0.0 || !value.is_finite() {
            return Self {
                fraction: value,
                exponent: 0,
            };
        }
        let (value, shift) = if value.abs() < f64::MIN_POSITIVE {
            (value * 18014398509481984.0, -54)
        } else {
            (value, 0)
        };
        let bits = value.to_bits();
        Self {
            fraction: f64::from_bits(
                (bits & ((1u64 << 63) | ((1u64 << 52) - 1))) | (1023u64 << 52),
            ),
            exponent: ((bits >> 52) & 2047) as i32 - 1023 + shift,
        }
    }

    fn product(
        self,
        other: Self,
    ) -> Self {
        let mut result = Self::from(self.fraction * other.fraction);
        result.exponent += self.exponent + other.exponent;
        result
    }

    fn ratio(
        self,
        other: Self,
    ) -> Self {
        let mut result = Self::from(self.fraction / other.fraction);
        result.exponent += self.exponent - other.exponent;
        result
    }

    fn value(self) -> f64 {
        let mut value = self.fraction;
        let mut exponent = self.exponent;
        if value == 0.0 || !value.is_finite() {
            return value;
        }
        while exponent > 512 {
            value *= 2.0_f64.powi(512);
            exponent -= 512;
        }
        while exponent < -512 {
            value *= 2.0_f64.powi(-512);
            exponent += 512;
        }
        value * 2.0_f64.powi(exponent)
    }

    fn times(
        self,
        value: f64,
    ) -> f64 {
        self.product(Self::from(value)).value()
    }

    const fn negated(self) -> Self {
        Self {
            fraction: -self.fraction,
            exponent: self.exponent,
        }
    }

    fn zero(self) -> bool {
        self.fraction == 0.0
    }

    const fn finite(self) -> bool {
        self.fraction.is_finite()
    }
}
fn max_abs(values: &[f64]) -> f64 {
    values.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
}
fn compensated(values: impl Iterator<Item = f64>) -> f64 {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}
fn dot(
    a: &[f64],
    b: &[f64],
) -> Scaled {
    let (fraction, exponent) =
        super::accurate_dot::product_sum(a.iter().copied().zip(b.iter().copied()));
    Scaled { fraction, exponent }
}

fn norm(values: &[f64]) -> Scaled {
    let max = max_abs(values);
    if max == 0.0 {
        return Scaled::from(0.0);
    }
    Scaled::from(max).product(Scaled::from(
        compensated(values.iter().map(|v| (v / max).powi(2))).sqrt(),
    ))
}
// A nonzero residual must stay distinguishable from zero even when its
// relative magnitude lies below the smallest representable positive float.
fn relative_parts(
    rnorm: Scaled,
    bnorm: Scaled,
) -> f64 {
    if rnorm.zero() {
        0.0
    } else if bnorm.zero() {
        f64::INFINITY
    } else {
        rnorm.ratio(bnorm).value().max(f64::from_bits(1))
    }
}
fn relative(
    r: &[f64],
    bnorm: Scaled,
) -> f64 {
    relative_parts(norm(r), bnorm)
}
pub(super) fn original_residual_ratio(
    b: &[f64],
    ax: &[f64],
) -> Option<f64> {
    if b.len() != ax.len() || b.iter().chain(ax).any(|v| !v.is_finite()) {
        return None;
    }
    let (r, scale) = difference(b, ax);
    Some(relative_parts(scale.product(norm(&r)), norm(b)))
}
fn evaluate(
    op: &impl Fn(&[f64], &mut [f64]),
    x: &[f64],
    out: &mut [f64],
) -> bool {
    if x.iter().any(|v| !v.is_finite()) {
        return false;
    }
    out.fill(f64::NAN);
    op(x, out);
    out.iter().all(|v| v.is_finite())
}
// Subtract first when representable: normalizing by |b| before cancellation
// could erase an unsolved tiny component beside an already solved huge one.
// Only a genuinely overflowing difference needs common operand scaling.
fn difference(
    b: &[f64],
    ax: &[f64],
) -> (Vec<f64>, Scaled) {
    let mut r: Vec<f64> = b.iter().zip(ax).map(|(b, a)| b - a).collect();
    let factor = if r.iter().all(|value| value.is_finite()) {
        Scaled::from(1.0)
    } else {
        let max = max_abs(b).max(max_abs(ax));
        for ((value, b), a) in r.iter_mut().zip(b).zip(ax) {
            *value = if value.is_finite() {
                *value / max
            } else {
                b / max - a / max
            };
        }
        Scaled::from(max)
    };
    let rm = max_abs(&r);
    if rm == 0.0 {
        return (r, Scaled::from(0.0));
    }
    for value in &mut r {
        *value /= rm;
    }
    (r, factor.product(Scaled::from(rm)))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn products_and_ratios_survive_both_extremes() {
        for s in [1e-300, 1e300] {
            let d = dot(&[s, 2.0 * s], &[s, 2.0 * s]);
            assert!(d.finite() && !d.zero());
            assert!((d.ratio(d).value() - 1.0).abs() < 1e-15);
            let n = norm(&[s, 2.0 * s]).ratio(Scaled::from(s)).value();
            assert!((n - 5.0_f64.sqrt()).abs() < 1e-14);
        }
    }
    #[test]
    fn ratios_are_not_materialized_before_vector_updates() {
        let a = Scaled::from(1e300).product(Scaled::from(1e300));
        assert!(a.value().is_infinite());
        assert!((a.times(1e-300) / 1e300 - 1.0).abs() < 1e-14);
        let b = Scaled::from(1e-300).product(Scaled::from(1e-300));
        assert_eq!(b.value(), 0.0);
        assert!((b.times(1e300) / 1e-300 - 1.0).abs() < 1e-14);
    }
    #[test]
    fn scaled_updates_keep_underflowing_and_subnormal_coefficients() {
        let tiny = Scaled::from(1e-300).product(Scaled::from(1e-300));
        assert!((add_scaled(0.0, tiny, 1e300) / 1e-300 - 1.0).abs() < 1e-14);
        let subnormal = Scaled::from(1.25).product(Scaled::from(f64::from_bits(1)));
        assert_eq!(add_scaled(0.0, subnormal, 1e300), subnormal.times(1e300));
    }
    #[test]
    fn original_residual_ratios_keep_both_exponent_extremes() {
        let tiny = f64::from_bits(1);
        assert_eq!(
            original_residual_ratio(&[1e308, tiny], &[1e308, 0.0]),
            Some(tiny)
        );
        let b = [1e-308; 4];
        let ax = [-3.0, 1e-308, 1e-308, 1e-308];
        let actual = original_residual_ratio(&b, &ax).unwrap();
        assert!(actual.is_finite() && (actual / 1.5e308 - 1.0).abs() < 1e-14);
        assert_eq!(original_residual_ratio(&[1e308], &[-1e308]), Some(2.0));
        assert!(original_residual_ratio(&[1.0], &[f64::NAN]).is_none());
    }
    #[test]
    fn subnormal_parts_round_trip() {
        for bits in [1, 2, 12345, (1u64 << 52) - 1] {
            let n = f64::from_bits(bits);
            assert_eq!(Scaled::from(n).value(), n);
        }
    }
}

#[derive(Clone, Copy)]
pub(super) enum Method {
    Cg,
    Gmres(usize),
    Bicgstab,
}
struct Inner {
    x: Vec<f64>,
    iterations: usize,
    failed: bool,
}
const fn inner_result(
    x: Vec<f64>,
    iterations: usize,
    failed: bool,
) -> Inner {
    Inner {
        x,
        iterations,
        failed,
    }
}
#[derive(Clone, Copy)]
enum PreparedScale {
    Fast(f64),
    Wide(Scaled),
}
impl PreparedScale {
    fn new(scale: Scaled) -> Self {
        let value = scale.value();
        if value.is_normal() || scale.zero() {
            Self::Fast(value)
        } else {
            Self::Wide(scale)
        }
    }

    fn add(
        self,
        x: f64,
        y: f64,
    ) -> f64 {
        match self {
            | Self::Fast(value) => value.mul_add(y, x),
            | Self::Wide(scale) => x + scale.times(y),
        }
    }
}
#[cfg(test)]
fn add_scaled(
    x: f64,
    coefficient: Scaled,
    y: f64,
) -> f64 {
    PreparedScale::new(coefficient).add(x, y)
}

fn normalize(values: &mut [f64]) -> Scaled {
    let max = max_abs(values);
    if max != 0.0 {
        for value in values {
            *value /= max;
        }
    }
    Scaled::from(max)
}
fn residual(
    apply: &impl Fn(&[f64], &mut [f64]),
    x: &[f64],
    b: &[f64],
    out: &mut [f64],
) -> bool {
    if !evaluate(apply, x, out) {
        return false;
    }
    for (value, rhs) in out.iter_mut().zip(b) {
        *value = rhs - *value;
    }
    out.iter().all(|v| v.is_finite())
}

// Select a binary operator gain and a safe product-input magnitude. A tiny
// matrix and equally tiny RHS can have an ordinary finite solution even when
// solving A*y = normalized_rhs would require an unrepresentable y.
fn calibrate(
    apply: &impl Fn(&[f64], &mut [f64]),
    rhs: &[f64],
    input: &mut [f64],
    output: &mut [f64],
) -> Option<(Scaled, Option<f64>)> {
    output.fill(f64::NAN);
    apply(rhs, output);
    if output.iter().any(|v| v.is_nan()) {
        return None;
    }
    let observed = max_abs(output);
    let amplifier_exponent = if !observed.is_finite() || observed > 2.0_f64.powi(500) {
        -512
    } else if observed < 2.0_f64.powi(-500) {
        512
    } else {
        0
    };
    if amplifier_exponent == 0 {
        return Some((Scaled::from(1.0), None));
    }
    let amplifier = 2.0_f64.powi(amplifier_exponent);
    if amplifier_exponent != 0 {
        for (input, rhs) in input.iter_mut().zip(rhs) {
            *input = rhs * amplifier;
        }
        if !evaluate(apply, input, output) {
            return None;
        }
    }
    let mut gain = Scaled::from(max_abs(output)).ratio(Scaled::from(amplifier));
    if gain.zero() {
        gain = Scaled::from(1.0);
    }
    // Binary gains preserve significands when moving between systems.
    gain.fraction = 1.0;
    Some((gain, Some(amplifier)))
}
fn normalized_apply(
    apply: &impl Fn(&[f64], &mut [f64]),
    gain: Scaled,
    amplifier: f64,
    input: &mut [f64],
    x: &[f64],
    output: &mut [f64],
) {
    let mut xscale = Scaled::from(max_abs(x));
    if xscale.zero() {
        output.fill(0.0);
        return;
    }
    xscale.fraction = 1.0;
    for (input, value) in input.iter_mut().zip(x) {
        *input = Scaled::from(*value).ratio(xscale).times(amplifier);
    }
    if !evaluate(apply, input, output) {
        output.fill(f64::NAN);
        return;
    }
    let restore = xscale.ratio(Scaled::from(amplifier)).ratio(gain);
    for value in output {
        *value = restore.times(*value);
    }
}

// A normalized solve may discard components beyond the f64 dynamic range.
// Refine against the original residual using a fresh scale, sharing one strict
// iteration budget. This also makes zero tolerance mean an exactly zero
// original residual rather than a ratio that rounded down to zero.
pub(super) fn solve(
    apply: impl Fn(&[f64], &mut [f64]),
    precond: impl Fn(&[f64], &mut [f64]),
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
    method: Method,
) -> IterResult {
    let mut guess = x0.map(<[f64]>::to_vec);
    let mut total = 0;
    loop {
        let (mut result, failed) = match solve_correction(
            &apply,
            &precond,
            b,
            guess.as_deref(),
            tol,
            max_iter - total,
            method,
        ) {
            | Ok(result) => (result, false),
            | Err(result) => (result, true),
        };
        total += result.iterations;
        let stopped = failed
            || result.converged
            || result.iterations == 0
            || total >= max_iter
            || guess
                .as_deref()
                .is_some_and(|old| old == result.x.as_slice());
        if stopped {
            result.iterations = total;
            return result;
        }
        guess = Some(result.x);
    }
}
fn solve_correction(
    apply: impl Fn(&[f64], &mut [f64]),
    precond: impl Fn(&[f64], &mut [f64]),
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
    method: Method,
) -> Result<IterResult, IterResult> {
    let n = b.len();
    let valid_guess = x0.is_none_or(|x| x.len() == n && x.iter().all(|v| v.is_finite()));
    let base = if valid_guess {
        x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec)
    } else {
        vec![0.0; n]
    };
    let failure = |x, iterations| {
        Err(IterResult {
            x,
            iterations,
            residual: f64::INFINITY,
            converged: false,
        })
    };
    if !valid_guess || !tol.is_finite() || tol < 0.0 || b.iter().any(|v| !v.is_finite()) {
        return failure(base, 0);
    }
    let mut ax = vec![0.0; n];
    if base.iter().any(|v| *v != 0.0) && !evaluate(&apply, &base, &mut ax) {
        return failure(base, 0);
    }
    let (rhs, scale) = difference(b, &ax);
    let bn = norm(b);
    let initial = relative_parts(scale.product(norm(&rhs)), bn);
    if initial <= tol || max_iter == 0 {
        return Ok(IterResult {
            x: base,
            residual: initial,
            iterations: 0,
            converged: initial <= tol,
        });
    }
    let target = if bn.zero() {
        0.0
    } else {
        Scaled::from(tol)
            .product(bn)
            .ratio(scale.product(norm(&rhs)))
            .value()
    };
    let mut scratch = vec![0.0; n];
    let Some((gain, amplifier)) = calibrate(&apply, &rhs, &mut scratch, &mut ax) else {
        return failure(base, 0);
    };
    let scratch = std::cell::RefCell::new(scratch);
    let adjusted = |x: &[f64], out: &mut [f64]| {
        if let Some(amplifier) = amplifier {
            normalized_apply(&apply, gain, amplifier, &mut scratch.borrow_mut(), x, out);
        } else {
            apply(x, out);
        }
    };
    let inner = match method {
        | Method::Cg => cg(&adjusted, &precond, &rhs, target, max_iter),
        | Method::Gmres(restart) => gmres(&adjusted, &precond, &rhs, restart, target, max_iter),
        | Method::Bicgstab => bicgstab(&adjusted, &precond, &rhs, target, max_iter),
    };
    let correction_scale = PreparedScale::new(scale.ratio(gain));
    let x: Vec<_> = base
        .iter()
        .zip(&inner.x)
        .map(|(base, correction)| correction_scale.add(*base, *correction))
        .collect();
    if x.iter().any(|v| !v.is_finite()) {
        return failure(base, inner.iterations);
    }
    if inner.failed || !evaluate(&apply, &x, &mut ax) {
        return failure(x, inner.iterations);
    }
    let (actual, actual_scale) = difference(b, &ax);
    let actual_norm = actual_scale.product(norm(&actual));
    let result = relative_parts(actual_norm, bn);
    Ok(IterResult {
        x,
        residual: result,
        iterations: inner.iterations,
        converged: result.is_finite() && result <= tol,
    })
}
fn cg(
    apply: &impl Fn(&[f64], &mut [f64]),
    precond: &impl Fn(&[f64], &mut [f64]),
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> Inner {
    let n = b.len();
    let bn = norm(b);
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut z = vec![0.0; n];
    let mut ap = vec![0.0; n];
    if !evaluate(precond, &r, &mut z) {
        return inner_result(x, 0, true);
    }
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut it = 0;
    while it < max_iter {
        let ps = normalize(&mut p);
        rz = rz.ratio(ps);
        if ps.zero() || rz.zero() || !rz.finite() || !evaluate(apply, &p, &mut ap) {
            return inner_result(x, it, true);
        }
        let pap = dot(&p, &ap);
        if pap.zero() || !pap.finite() {
            return inner_result(x, it, false);
        }
        let alpha = rz.ratio(pap);
        let alpha_step = PreparedScale::new(alpha);
        let negative_alpha_step = PreparedScale::new(alpha.negated());
        let candidate: Vec<_> = x
            .iter()
            .zip(&p)
            .map(|(x, p)| alpha_step.add(*x, *p))
            .collect();
        if candidate.iter().any(|v| !v.is_finite()) {
            return inner_result(x, it, true);
        }
        x = candidate;
        it += 1;
        for i in 0..n {
            r[i] = negative_alpha_step.add(r[i], ap[i]);
        }
        if r.iter().any(|v| !v.is_finite()) {
            return inner_result(x, it, true);
        }
        let reset = if relative(&r, bn) <= tol {
            if !residual(apply, &x, b, &mut r) {
                return inner_result(x, it, true);
            }
            if relative(&r, bn) <= tol {
                break;
            }
            true
        } else {
            false
        };
        if !evaluate(precond, &r, &mut z) {
            return inner_result(x, it, true);
        }
        let next = dot(&r, &z);
        let beta = PreparedScale::new(next.ratio(rz));
        for i in 0..n {
            p[i] = if reset {
                z[i]
            } else {
                beta.add(z[i], p[i])
            };
        }
        rz = next;
    }
    inner_result(x, it, false)
}
// Back-substitution and reconstruction are kept separate from Arnoldi so
// a singular or unrepresentable projected system leaves the last finite iterate.
fn projected_update(
    x: &[f64],
    h: &[Vec<f64>],
    g: &[f64],
    z: &[Vec<f64>],
    used: usize,
) -> Option<Vec<f64>> {
    let mut y = vec![0.0; used];
    for i in (0..used).rev() {
        let mut value = g[i];
        for j in i + 1..used {
            value = (-h[i][j]).mul_add(y[j], value);
        }
        if h[i][i] == 0.0 {
            return None;
        }
        y[i] = value / h[i][i];
        if !y[i].is_finite() {
            return None;
        }
    }
    let mut candidate = x.to_vec();
    for j in 0..used {
        for q in 0..x.len() {
            candidate[q] = y[j].mul_add(z[j][q], candidate[q]);
        }
    }
    if candidate.iter().any(|v| !v.is_finite()) {
        return None;
    }
    Some(candidate)
}
fn gmres(
    apply: &impl Fn(&[f64], &mut [f64]),
    precond: &impl Fn(&[f64], &mut [f64]),
    b: &[f64],
    restart: usize,
    tol: f64,
    max_iter: usize,
) -> Inner {
    let n = b.len();
    let m = restart.clamp(1, n.max(1));
    let bn = norm(b);
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut total = 0;
    while total < max_iter && relative(&r, bn) > tol {
        let beta = norm(&r).value();
        if !beta.is_finite() || beta == 0.0 {
            return inner_result(x, total, true);
        }
        let mut v = vec![r.iter().map(|e| e / beta).collect::<Vec<_>>()];
        let mut z = Vec::with_capacity(m);
        let mut h = vec![vec![0.0; m]; m + 1];
        let (mut cs, mut sn) = (vec![0.0_f64; m], vec![0.0_f64; m]);
        let mut g = vec![0.0; m + 1];
        g[0] = beta;
        let mut used = 0;
        for j in 0..m {
            let mut zj = vec![0.0; n];
            let mut w = vec![0.0; n];
            if !evaluate(precond, &v[j], &mut zj) {
                return inner_result(x, total, true);
            }
            normalize(&mut zj);
            if !evaluate(apply, &zj, &mut w) {
                return inner_result(x, total, true);
            }
            z.push(zj);
            // Twice-modified Gram-Schmidt retains orthogonality near breakdown.
            for _ in 0..2 {
                for i in 0..=j {
                    let projection = dot(&w, &v[i]).value();
                    if !projection.is_finite() {
                        return inner_result(x, total, true);
                    }
                    h[i][j] += projection;
                    for q in 0..n {
                        w[q] = (-projection).mul_add(v[i][q], w[q]);
                    }
                }
            }
            h[j + 1][j] = norm(&w).value();
            if !h[j + 1][j].is_finite() {
                return inner_result(x, total, true);
            }
            let breakdown = h[j + 1][j] == 0.0;
            if !breakdown {
                v.push(w.iter().map(|e| e / h[j + 1][j]).collect());
            }
            for i in 0..j {
                let t = cs[i].mul_add(h[i][j], sn[i] * h[i + 1][j]);
                h[i + 1][j] = (-sn[i]).mul_add(h[i][j], cs[i] * h[i + 1][j]);
                h[i][j] = t;
            }
            let d = h[j][j].hypot(h[j + 1][j]);
            if !d.is_finite() {
                return inner_result(x, total, true);
            }
            if d == 0.0 {
                cs[j] = 1.0;
                sn[j] = 0.0;
            } else {
                cs[j] = h[j][j] / d;
                sn[j] = h[j + 1][j] / d;
            }
            h[j][j] = d;
            h[j + 1][j] = 0.0;
            g[j + 1] = -sn[j] * g[j];
            g[j] *= cs[j];
            total += 1;
            used = j + 1;
            if g[j + 1].abs() / bn.value() <= tol || breakdown || total >= max_iter {
                break;
            }
        }
        let Some(candidate) = projected_update(&x, &h, &g, &z, used) else {
            return inner_result(x, total, false);
        };
        x = candidate;
        if !residual(apply, &x, b, &mut r) {
            return inner_result(x, total, true);
        }
    }
    inner_result(x, total, false)
}
fn bicgstab(
    apply: &impl Fn(&[f64], &mut [f64]),
    precond: &impl Fn(&[f64], &mut [f64]),
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> Inner {
    let n = b.len();
    let bn = norm(b);
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut shadow = r.clone();
    let mut p = vec![0.0; n];
    let mut v = vec![0.0; n];
    let mut ph = vec![0.0; n];
    let mut sh = vec![0.0; n];
    let mut t = vec![0.0; n];
    let mut ph_scale = Scaled::from(1.0);
    let (mut rho, mut alpha, mut omega) = (Scaled::from(1.0), Scaled::from(1.0), Scaled::from(1.0));
    let mut reset = true;
    let mut it = 0;
    while it < max_iter {
        let next = dot(&shadow, &r);
        if next.zero() || !next.finite() {
            return inner_result(x, it, false);
        }
        let beta = next.ratio(rho).product(alpha.ratio(omega));
        rho = next;
        let beta_step = PreparedScale::new(beta);
        let previous_omega_step = PreparedScale::new(omega.negated().product(ph_scale));
        for i in 0..n {
            p[i] = if reset {
                r[i]
            } else {
                beta_step.add(r[i], previous_omega_step.add(p[i], v[i]))
            };
        }
        reset = false;
        if !evaluate(precond, &p, &mut ph) {
            return inner_result(x, it, true);
        }
        ph_scale = normalize(&mut ph);
        if !evaluate(apply, &ph, &mut v) {
            return inner_result(x, it, true);
        }
        let denom = dot(&shadow, &v).product(ph_scale);
        if denom.zero() || !denom.finite() {
            return inner_result(x, it, false);
        }
        alpha = rho.ratio(denom);
        let alpha_step = alpha.product(ph_scale);
        let positive_alpha = PreparedScale::new(alpha_step);
        let negative_alpha = PreparedScale::new(alpha_step.negated());
        let s: Vec<_> = r
            .iter()
            .zip(&v)
            .map(|(r, v)| negative_alpha.add(*r, *v))
            .collect();
        it += 1;
        if relative(&s, bn) <= tol {
            for i in 0..n {
                x[i] = positive_alpha.add(x[i], ph[i]);
            }
            r = s;
        } else {
            if !evaluate(precond, &s, &mut sh) {
                return inner_result(x, it, true);
            }
            let sh_scale = normalize(&mut sh);
            if !evaluate(apply, &sh, &mut t) {
                return inner_result(x, it, true);
            }
            let tt = dot(&t, &t);
            if tt.zero() || !tt.finite() {
                return inner_result(x, it, false);
            }
            omega = dot(&t, &s).ratio(tt).ratio(sh_scale);
            let omega_step = omega.product(sh_scale);
            let positive_omega = PreparedScale::new(omega_step);
            let negative_omega = PreparedScale::new(omega_step.negated());
            for i in 0..n {
                x[i] = positive_omega.add(positive_alpha.add(x[i], ph[i]), sh[i]);
                r[i] = negative_omega.add(s[i], t[i]);
            }
        }
        if x.iter().chain(&r).any(|v| !v.is_finite()) {
            return inner_result(x, it, true);
        }
        if relative(&r, bn) <= tol {
            if !residual(apply, &x, b, &mut r) {
                return inner_result(x, it, true);
            }
            if relative(&r, bn) <= tol {
                break;
            }
            shadow.clone_from(&r);
            reset = true;
        }
        if omega.zero() || !omega.finite() {
            break;
        }
    }
    inner_result(x, it, false)
}
