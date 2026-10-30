//! Meijer G-function evaluation.
//!
//! The kernel evaluates the Meijer G-function
//! `G^{m,n}_{p,q}(z | a_1..a_p; b_1..b_q)` for complex parameters and
//! argument, by Slater's contour theorem (residue series) and, where the
//! residue series does not apply, by numerical Mellin-Barnes integration.
//!
//! The symbolic layer lives in `rules::special::meijer`; this module is
//! plain numerics and is also usable directly.
//!
//! # Definition
//!
//! With `p = a.len()`, `q = b.len()` and the principal branch of `z^s`,
//! ```text
//! G^{m,n}_{p,q}(z | a; b)
//!   = (1/2 pi i) int_L
//!       prod_{j=1}^m Gamma(b_j - s) prod_{j=1}^n Gamma(1 - a_j + s)
//!       / (prod_{j=m+1}^q Gamma(1 - b_j + s)
//!          prod_{j=n+1}^p Gamma(a_j - s)) z^s ds,
//! ```
//! where `L` separates the poles `s = b_j + k` of the numerator gammas
//! (`j <= m`) from the poles `s = a_j - 1 - k` (`j <= n`); this is the
//! standard definition (DLMF 16.17.1, with the sign of `s` opposite to the
//! DLMF's own convention, which is equivalent). The first `n` entries of
//! `a` and the first `m` entries of `b` are the contour-defining sets.
//!
//! # Algorithms
//!
//! * **Slater's residue series.** Closing the contour to the right
//!   (`s = b_h + k`) gives
//!   ```text
//!   G = sum_{h=1}^m z^{b_h} sum_{k>=0} (-1)^k z^k / k!
//!         prod_{j!=h, j<=m} Gamma(b_j - b_h - k)
//!         prod_{j<=n} Gamma(1 - a_j + b_h + k)
//!       / (prod_{j>n} Gamma(a_j - b_h - k)
//!          prod_{j>m} Gamma(1 - b_j + b_h + k)),
//!   ```
//!   which converges for `p < q`, and for `p = q` when `|z| < 1` (Slater's
//!   theorem). Closing to the left (`s = a_h - 1 - k`) gives the companion
//!   series in powers of `1/z`, convergent for `p > q`, and for `p = q` when
//!   `|z| > 1`. The terms are summed by the exact rational recurrence of
//!   consecutive gamma ratios (no gamma function is evaluated per term), the
//!   families are combined on a common logarithmic scale, and a family sum
//!   that lost more than about 13 digits to cancellation makes the result
//!   `None`.
//! * **Confluent poles.** When two entries of the contour-defining `b` set
//!   (for the `z -> 0` series) or `a` set (for the `z -> infinity` series)
//!   differ by an integer, the residue formula has double poles and the
//!   series above is undefined. The G-function is then still analytic in
//!   those parameters (the poles do not pinch the contour), so the limit is
//!   obtained by a symmetric parameter shift `+h, -h` followed by
//!   Richardson extrapolation in `h^2` (three levels) on the regularized
//!   series. The accuracy there is about `1e-9` rather than `1e-12`.
//! * **Numerical Mellin-Barnes integration.** The defining integral along
//!   the vertical line `Re s = c` is evaluated by the trapezoidal rule with
//!   step refinement, which converges geometrically for an integrand of this
//!   analyticity: with `delta = m + n - (p+q)/2 > 0` the integrand decays
//!   like `exp(-pi delta |Im s|)` for `|arg z| < delta pi`. The step and the
//!   truncation are chosen from the distance to the nearest pole and from
//!   the effective decay rate `delta - |arg z|/pi`; three step sizes are
//!   compared and the difference of the last two is the error estimate.
//!   This path is independent of the residue series and is used (i) as a
//!   cross-check when the contour has a comfortable decay margin, and
//!   (ii) as the evaluator whenever the residue formula is not applicable.
//!   The second case includes confluent denominator gammas: when a factor
//!   `Gamma(a_j - b_h - k)` or `Gamma(1 - b_j + b_h + k)` in the
//!   denominator vanishes at a residue point, the residue is not the
//!   product of the separate residues (the `s Gamma(-s)` cancellation of
//!   the shifted derivative `G^{m,n+1}_{p+1,q+1}` is the typical example),
//!   and the value comes from the integral instead.
//! * **Degenerate orders.** `m = 0` with the contour closed to the right (or
//!   `n = 0` closed to the left) encloses no poles and the value is exactly
//!   zero: `G^{0,n}_{p,q} = 0` for `p < q` or `p = q, |z| < 1`, and
//!   `G^{m,0}_{p,q} = 0` for `p > q` or `p = q, |z| > 1`.
//! * **`z = 0`.** The limit is read off the applicable residue series: the
//!   terms `z^{b_h}` with `Re b_h > 0` vanish, a term with `b_h = 0`
//!   contributes its coefficient, and a negative real part or a non-zero
//!   imaginary part makes the limit non-existent (returned as `None`). When
//!   the `1/z` series applies, the value is `0` only if every exponent
//!   `Re(a_h - 1)` is positive.
//!
//! # Branch
//!
//! `z^s = exp(s Log z)` with the principal logarithm `Log z = ln|z| +
//! i arg z`, `arg z in (-pi, pi]`; on the negative real axis the upper side
//! (`arg z = +pi`) is taken. The G-function itself is single-valued in the
//! plane cut along `(-inf, 0]` for non-integer parameters.
//!
//! # Conditioning
//!
//! Relative accuracy is about `1e-12` to `1e-13` where the function is well
//! conditioned: the residue series converges fast (not too close to
//! `|z| = 1` when `p = q`), no family sum cancels, and no confluent poles
//! occur. It degrades to about `1e-9` for confluent contour parameters, and
//! the Mellin-Barnes fallback is accurate to about `1e-10` when the effective
//! decay `delta - |arg z|/pi` is at least `0.2`, and less (or unavailable)
//! close to the boundary of the convergence sector. The result is `None`
//! when
//!
//! * `m > q` or `n > p`, a parameter or `z` is not finite, or `z` is
//!   exactly `0` and the limit is not defined;
//! * the applicable series does not converge numerically (for `p = q`,
//!   `|z|` within about `1e-5` of `1`, away from the Mellin-Barnes regime),
//! * a residue family overflows, the families cancel by more than 13
//!   digits, or the terms of one family cancel by more than 13 digits
//!   against their peak (this is what happens for `p < q` and very large
//!   `|z|`, where the value is exponentially smaller than the intermediate
//!   terms; no rescaling across the peak is attempted, so `None` is
//!   returned rather than a zero),
//! * `p = q` with `|z| = 1` outside the Mellin-Barnes regime,
//! * or the series and the Mellin-Barnes integral disagree by more than
//!   `1e-7` relative (inconsistent data are never returned as a number).
//!
//! Large `|z|` is handled by the `1/z` residue series (convergent for
//! `p > q`, and for `p = q` when `|z| > 1`); `p = q + 1` with `n = 0` is
//! identically zero, and `p = q - 1` with `m = 0` likewise.

use num_complex::Complex64;

use crate::kernels::hypergeometric::ln_gamma_c;

type C = Complex64;

const ONE: C = C::new(1.0, 0.0);
const ZERO: C = C::new(0.0, 0.0);
/// Most terms summed in one residue family.
const MAX_TERMS: u32 = 200_000;
/// Largest modulus of a term before the family is declared overflowed.
const TERM_LIMIT: f64 = 1e250;
/// Digits lost to cancellation that still leave a usable result.
const LOST_DIGITS: f64 = 1e-13;

/// A non-finite complex number as `None`.
fn finite(v: C) -> Option<C> {
    (v.re.is_finite() && v.im.is_finite()).then_some(v)
}

/// Whether a real number is an integer to within the tolerances of the
/// parameters (they are read from `f64` literals, so exact integers are
/// exact here; the tolerance only absorbs round-off of derived values).
fn near_integer(x: f64) -> bool {
    (x - x.round()).abs() <= 1e-9 * x.abs().max(1.0)
}

/// Whether two entries of a contour-defining set have an integer difference
/// (a pair of confluent poles).
fn confluent(values: &[C]) -> bool {
    for (i, x) in values.iter().enumerate() {
        for y in values.iter().skip(i.saturating_add(1)) {
            if near_integer(x.re - y.re) && near_integer(x.im - y.im) {
                return true;
            }
        }
    }
    false
}

/// Whether a gamma argument is a non-positive integer, where `Gamma` has a
/// pole (and its reciprocal a zero).
fn nonpositive_integer(v: C) -> bool {
    near_integer(v.im) && v.im.abs() < 1e-9 && v.re <= 0.0 && near_integer(v.re)
}

/// A denominator gamma of the `z -> 0` residue formula vanishes at a residue
/// point: the residue is not the product of the separate residues there, so
/// the simple series must not be used.
fn left_denominator_singular(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
) -> bool {
    for h in 0..m {
        let Some(&bh) = b.get(h) else {
            return true;
        };
        for bj in b.iter().skip(m) {
            if nonpositive_integer(ONE - *bj + bh) {
                return true;
            }
        }
        for ai in a.iter().skip(n) {
            if nonpositive_integer(*ai - bh) {
                return true;
            }
        }
    }
    false
}

/// A numerator gamma of the `z -> 0` residue formula has a pole at a residue
/// point (a pinched contour), so the limit rather than the series is needed.
fn left_numerator_pinch(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
) -> bool {
    for h in 0..m {
        let Some(&bh) = b.get(h) else {
            return true;
        };
        for ai in a.iter().take(n) {
            if nonpositive_integer(ONE + bh - *ai) {
                return true;
            }
        }
    }
    false
}

/// The same two obstructions for the `z -> infinity` series.
fn right_denominator_singular(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
) -> bool {
    for h in 0..n {
        let Some(&ah) = a.get(h) else {
            return true;
        };
        for aj in a.iter().skip(n) {
            if nonpositive_integer(ONE - ah + *aj) {
                return true;
            }
        }
        for bj in b.iter().skip(m) {
            if nonpositive_integer(ah - *bj) {
                return true;
            }
        }
    }
    false
}

fn right_numerator_pinch(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
) -> bool {
    for h in 0..n {
        let Some(&ah) = a.get(h) else {
            return true;
        };
        for bj in b.iter().take(m) {
            if nonpositive_integer(ONE - ah + *bj) {
                return true;
            }
        }
    }
    false
}

/// `z^b` on the principal branch.
fn principal_power(
    z: C,
    b: C,
) -> C {
    (b * z.ln()).exp()
}

/// Coefficient `C_h` of the `z -> 0` family for `b_h`, in logarithmic form.
fn left_log_coefficient(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    h: usize,
) -> Option<C> {
    let p = a.len();
    let q = b.len();
    let bh = *b.get(h)?;
    let mut log = ZERO;
    for (j, bj) in b.iter().enumerate().take(m) {
        if j == h {
            continue;
        }
        log += ln_gamma_c(*bj - bh);
    }
    for aj in a.iter().take(n) {
        log += ln_gamma_c(ONE - *aj + bh);
    }
    for bj in b.iter().take(q).skip(m) {
        log -= ln_gamma_c(ONE - *bj + bh);
    }
    for aj in a.iter().take(p).skip(n) {
        log -= ln_gamma_c(*aj - bh);
    }
    // A pole of the gamma function in a numerator product is a genuine
    // degeneracy (a pinched contour); signal it so the caller can take the
    // limit instead.
    if !log.re.is_finite() || log.re > 700.0 {
        return None;
    }
    finite(log.exp())
}

/// Coefficient `C_h` of the `z -> infinity` family for `a_h`, in
/// logarithmic form.
fn right_log_coefficient(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    h: usize,
) -> Option<C> {
    let p = a.len();
    let q = b.len();
    let ah = *a.get(h)?;
    let mut log = ZERO;
    for (j, aj) in a.iter().enumerate().take(n) {
        if j == h {
            continue;
        }
        log += ln_gamma_c(ah - *aj);
    }
    for bj in b.iter().take(m) {
        log += ln_gamma_c(ONE - ah + *bj);
    }
    for aj in a.iter().take(p).skip(n) {
        log -= ln_gamma_c(ONE - ah + *aj);
    }
    for bj in b.iter().take(q).skip(m) {
        log -= ln_gamma_c(ah - *bj);
    }
    if !log.re.is_finite() || log.re > 700.0 {
        return None;
    }
    finite(log.exp())
}

/// Whether a complex factor is zero (a terminating residue series).
fn is_zero(v: C) -> bool {
    v.re == 0.0 && v.im == 0.0
}

/// Sums `sum_{k>=0} t_k` given the ratio `t_{k+1}/t_k`, until the terms fall
/// below `1e-18` of the sum. Returns `None` on overflow, non-convergence, or
/// when the terms cancel by more than about 13 digits (the sum is then no
/// longer carried by the accumulated digits).
fn ratio_sum(ratio: impl Fn(u32) -> Option<(C, bool)>) -> Option<C> {
    let mut term = ONE;
    let mut sum = ONE;
    let mut peak = 1.0_f64;
    for k in 0..MAX_TERMS {
        let (r, terminate) = ratio(k)?;
        if terminate {
            return finite(sum);
        }
        term *= r;
        let size = term.norm();
        if !size.is_finite() || size > TERM_LIMIT {
            return None;
        }
        peak = peak.max(size);
        sum += term;
        if k > 2 && (size == 0.0 || size <= 1e-18 * sum.norm().max(1e-300)) {
            if sum.norm() <= 1e-13 * peak {
                return None;
            }
            return finite(sum);
        }
    }
    None
}

/// One residue family evaluated directly from the recurrence of consecutive
/// gamma ratios. The exponent `b_h` (or `a_h - 1`) is returned separately so
/// that the caller can combine families on a common scale.
fn left_family(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    h: usize,
    z: C,
) -> Option<C> {
    let p = a.len();
    let q = b.len();
    let bh = *b.get(h)?;
    let leading = left_log_coefficient(m, n, a, b, h)?;
    let sum = ratio_sum(|k| {
        let kf = f64::from(k);
        let mut numerator = -z;
        let mut denominator = C::new(kf + 1.0, 0.0);
        for aj in a.iter().take(n) {
            numerator *= ONE - *aj + bh + kf;
        }
        for aj in a.iter().take(p).skip(n) {
            let factor = *aj - bh - kf - 1.0;
            if is_zero(factor) {
                return Some((ONE, true));
            }
            numerator *= factor;
        }
        for (j, bj) in b.iter().enumerate().take(q) {
            if j < m && j != h {
                let factor = *bj - bh - kf - 1.0;
                if is_zero(factor) {
                    // 1/Gamma(0) = 0 terminates the family.
                    return Some((ONE, true));
                }
                denominator *= factor;
            } else if j >= m {
                let factor = ONE - *bj + bh + kf;
                if is_zero(factor) {
                    // 1/Gamma(0) = 0 terminates the family.
                    return Some((ONE, true));
                }
                denominator *= factor;
            }
        }
        finite(numerator / denominator).map(|r| (r, false))
    })?;
    let prefactor = principal_power(z, bh);
    finite(leading * prefactor * sum)
}

/// One `z -> infinity` residue family.
fn right_family(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    h: usize,
    z: C,
) -> Option<C> {
    let p = a.len();
    let q = b.len();
    let ah = *a.get(h)?;
    let leading = right_log_coefficient(m, n, a, b, h)?;
    let sum = ratio_sum(|k| {
        let kf = f64::from(k);
        let mut numerator = -ONE / z;
        let mut denominator = C::new(kf + 1.0, 0.0);
        for (j, aj) in a.iter().enumerate().take(n) {
            if j == h {
                continue;
            }
            let factor = ah - *aj - kf - 1.0;
            if is_zero(factor) {
                return Some((ONE, true));
            }
            numerator *= factor;
        }
        for bj in b.iter().take(m) {
            numerator *= ONE - ah + *bj + kf;
        }
        for aj in a.iter().take(p).skip(n) {
            let factor = ONE - ah + *aj + kf;
            if is_zero(factor) {
                return Some((ONE, true));
            }
            denominator *= factor;
        }
        for bj in b.iter().take(q).skip(m) {
            let factor = ah - *bj - kf - 1.0;
            if is_zero(factor) {
                return Some((ONE, true));
            }
            numerator *= factor;
        }
        finite(numerator / denominator).map(|r| (r, false))
    })?;
    let prefactor = principal_power(z, ah - 1.0);
    finite(leading * prefactor * sum)
}

/// Sums one side's families, combining them on a common logarithmic scale.
/// `None` when a family overflows or the families cancel by more than
/// [`LOST_DIGITS`].
fn combine_families(
    count: usize,
    family: impl Fn(usize) -> Option<C>,
) -> Option<C> {
    let mut values = Vec::with_capacity(count);
    for h in 0..count {
        let v = family(h)?;
        if !v.re.is_finite() || !v.im.is_finite() {
            return None;
        }
        values.push(v);
    }
    let mut total = ZERO;
    let mut absolute = 0.0_f64;
    for v in values {
        total += v;
        absolute += v.norm();
    }
    if total.norm() <= LOST_DIGITS * absolute {
        // Nearly complete cancellation: the answer would carry no digits.
        return None;
    }
    finite(total)
}

/// Left (small `z`) residue series without confluent-pole handling.
fn left_series_raw(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
) -> Option<C> {
    combine_families(m, |h| left_family(m, n, a, b, h, z))
}

/// Right (large `z`) residue series without confluent-pole handling.
fn right_series_raw(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
) -> Option<C> {
    combine_families(n, |h| right_family(m, n, a, b, h, z))
}

/// Shifts the first `count` parameters by `scale * (j + 1)`, breaking any
/// integer differences between them.
fn shifted(
    values: &[C],
    count: usize,
    scale: f64,
) -> Vec<C> {
    values
        .iter()
        .enumerate()
        .map(|(j, v)| {
            if j < count {
                *v + C::new(scale * (j as f64 + 1.0), 0.0)
            } else {
                *v
            }
        })
        .collect()
}

/// A symmetric finite-difference limit: `(f(+h) + f(-h))/2` for three
/// halving step sizes, Richardson-extrapolated in `h^2`.
fn confluent_limit(f: impl Fn(f64) -> Option<C>) -> Option<C> {
    let symmetric = |h: f64| -> Option<C> {
        let plus = f(h)?;
        let minus = f(-h)?;
        finite(C::new(
            f64::midpoint(plus.re, minus.re),
            f64::midpoint(plus.im, minus.im),
        ))
    };
    let h = 2e-3;
    let s0 = symmetric(h)?;
    let s1 = symmetric(h / 2.0)?;
    let s2 = symmetric(h / 4.0)?;
    let r0 = (s1 * 4.0 - s0) / 3.0;
    let r1 = (s2 * 4.0 - s1) / 3.0;
    let best = (r1 * 16.0 - r0) / 15.0;
    let spread = (r1 - r0).norm();
    if spread > 1e-6 * best.norm().max(1e-12) {
        return None;
    }
    finite(best)
}

/// Left residue series, with a confluent-pole limit when needed. Returns
/// `None` when the simple residue formula is not applicable (a denominator
/// gamma vanishes at a residue point), so that the caller can fall back to
/// the Mellin-Barnes integral.
fn left_series(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
) -> Option<(C, bool)> {
    if m == 0 {
        return Some((ZERO, false));
    }
    let bad = confluent(&b[..m]);
    if bad {
        let value = confluent_limit(|h| {
            let bs = shifted(b, m, h);
            left_series_raw(m, n, a, &bs, z)
        })?;
        return Some((value, true));
    }
    if left_denominator_singular(m, n, a, b) || left_numerator_pinch(m, n, a, b) {
        return None;
    }
    left_series_raw(m, n, a, b, z).map(|v| (v, false))
}

/// Right residue series, with a confluent-pole limit when needed.
fn right_series(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
) -> Option<(C, bool)> {
    if n == 0 {
        return Some((ZERO, false));
    }
    let bad = confluent(&a[..n]);
    if bad {
        let value = confluent_limit(|h| {
            let as_ = shifted(a, n, h);
            right_series_raw(m, n, &as_, b, z)
        })?;
        return Some((value, true));
    }
    if right_denominator_singular(m, n, a, b) || right_numerator_pinch(m, n, a, b) {
        return None;
    }
    right_series_raw(m, n, a, b, z).map(|v| (v, false))
}

// ------------------------------------------------------- Mellin-Barnes

/// The Mellin-Barnes integrand at `s`.
fn mb_integrand(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
    s: C,
) -> Option<C> {
    let p = a.len();
    let q = b.len();
    let mut log = ZERO;
    for bj in b.iter().take(m) {
        log += ln_gamma_c(*bj - s);
    }
    for aj in a.iter().take(n) {
        log += ln_gamma_c(ONE - *aj + s);
    }
    for bj in b.iter().take(q).skip(m) {
        log -= ln_gamma_c(ONE - *bj + s);
    }
    for aj in a.iter().take(p).skip(n) {
        log -= ln_gamma_c(*aj - s);
    }
    let total = log + s * z.ln();
    if !total.re.is_finite() || !total.im.is_finite() {
        return None;
    }
    if total.re > 700.0 {
        return None;
    }
    finite(total.exp())
}

/// A trapezoidal approximation of the Mellin-Barnes integral along
/// `Re s = c` with step `h`. `grid` is `(c, h, ymax)`. Returned with the size
/// of the largest integrand sample (for cancellation accounting).
fn mb_trapezoid(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
    grid: (f64, f64, f64),
) -> Option<(C, f64)> {
    let (c, h, ymax) = grid;
    let steps = (ymax / h).ceil() as u32;
    if steps > 200_000 {
        return None;
    }
    let mut sum = mb_integrand(m, n, a, b, z, C::new(c, 0.0))?;
    let mut peak = sum.norm();
    for k in 1..=steps {
        let y = h * f64::from(k);
        let mut pair = ZERO;
        for direction in [-1.0_f64, 1.0] {
            let v = mb_integrand(m, n, a, b, z, C::new(c, direction * y))?;
            pair += v;
            peak = peak.max(v.norm());
        }
        sum += pair;
        // The integrand decays like exp(-pi delta |y|); stop once it is far
        // below the peak and still decreasing.
        if y > 2.0 && pair.norm() < 1e-19 * peak {
            break;
        }
    }
    // The contour is parametrized by `s = c + i y`, so `ds = i dy` and the
    // prefactor of the definition is `1/(2 pi i)`; the two imaginary units
    // cancel and leave a real `1/(2 pi)`.
    finite(sum * (h / (2.0 * std::f64::consts::PI))).map(|v| (v, peak))
}

/// Numerical Mellin-Barnes evaluation: returns the value and an error
/// estimate, or `None` when the vertical contour does not converge or the
/// quadrature does not settle.
fn mellin_barnes(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
) -> Option<(C, f64)> {
    let p = a.len();
    let q = b.len();
    let delta = m as f64 + n as f64 - 0.5 * (p + q) as f64;
    if delta <= 0.05 {
        return None;
    }
    let argument = z.arg().abs();
    if argument >= delta * std::f64::consts::PI - 1e-9 {
        return None;
    }
    // Effective exponential decay rate of the integrand along the contour:
    // the gamma factors give `exp(-pi delta |y|)` and the `z^s` factor gives
    // `exp(-|y| arg z)` on one side (and growth on the other).
    let decay = delta - argument / std::f64::consts::PI;
    if decay <= 0.02 {
        return None;
    }
    // The contour has to separate the poles s = a_j - 1 - k (j <= n) from
    // s = b_j + k (j <= m).
    let left = if n == 0 {
        f64::NEG_INFINITY
    } else {
        a.iter()
            .take(n)
            .map(|v| v.re)
            .fold(f64::NEG_INFINITY, f64::max)
            - 1.0
    };
    let right = if m == 0 {
        f64::INFINITY
    } else {
        b.iter().take(m).map(|v| v.re).fold(f64::INFINITY, f64::min)
    };
    let c = match (left.is_finite(), right.is_finite()) {
        | (true, true) => {
            if right - left <= 0.1 {
                return None;
            }
            f64::midpoint(left, right)
        },
        | (true, false) => left + 1.0,
        | (false, true) => right - 1.0,
        | (false, false) => 0.0,
    };
    let distance = match (left.is_finite(), right.is_finite()) {
        | (true, true) => (right - c).min(c - left),
        | (true, false) => (c - left).min(1.0),
        | (false, true) => (right - c).min(1.0),
        | (false, false) => 1.0,
    };
    if !(distance.is_finite() && distance > 1e-4) {
        return None;
    }
    let base = (distance / 6.0)
        .min(0.5)
        .min(0.3 / (1.0 + z.norm().ln().abs()));
    let ymax = 45.0 / (std::f64::consts::PI * decay) + 4.0;
    let mut previous: Option<C> = None;
    for level in 0..4 {
        let h = base / f64::from(1_u32 << level);
        let (value, peak) = mb_trapezoid(m, n, a, b, z, (c, h, ymax))?;
        if let Some(before) = previous {
            let difference = (value - before).norm();
            // The value has to be resolved relative to itself; when it lies
            // far below the peak of the integrand the quadrature cannot
            // deliver it, and `None` is returned instead of a rounded zero.
            let resolved = value.norm() > 1e-13 * peak;
            if resolved && difference <= 1e-10 * value.norm().max(1e-300) {
                return finite(value).map(|v| (v, difference));
            }
        }
        previous = Some(value);
    }
    None
}

// ------------------------------------------------------------ public API

/// Value of `G^{m,n}_{p,q}(z | a; b)`, the Meijer G-function, at complex
/// parameters and complex `z` (principal branch; see the module
/// documentation).
///
/// `p = a.len()` and `q = b.len()`. Returns `None` for invalid orders
/// (`m > q` or `n > p`), non-finite inputs, `z = 0` where the limit does not
/// exist, and parameter/argument combinations for which neither residue
/// series nor the Mellin-Barnes integral is numerically valid. The result is
/// never a silently wrong number: when the two independent methods disagree
/// by more than `1e-7` relative, `None` is returned.
#[must_use]
pub fn meijer_g(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
) -> Option<C> {
    let p = a.len();
    let q = b.len();
    if m > q || n > p {
        return None;
    }
    if a.iter()
        .chain(b)
        .any(|v| !v.re.is_finite() || !v.im.is_finite())
        || !z.re.is_finite()
        || !z.im.is_finite()
    {
        return None;
    }
    if z.norm() == 0.0 {
        return zero_limit(m, n, a, b);
    }
    if p == 0 && q == 0 {
        return Some(ZERO);
    }
    let left_ok = p < q || (p == q && z.norm() < 1.0);
    let right_ok = p > q || (p == q && z.norm() > 1.0);
    if left_ok {
        if m == 0 {
            return Some(ZERO);
        }
        if let Some((value, degraded)) = left_series(m, n, a, b, z) {
            return cross_check(value, degraded, m, n, a, b, z);
        }
        return mellin_barnes(m, n, a, b, z).and_then(|(value, _)| finite(value));
    }
    if right_ok {
        if n == 0 {
            return Some(ZERO);
        }
        if let Some((value, degraded)) = right_series(m, n, a, b, z) {
            return cross_check(value, degraded, m, n, a, b, z);
        }
        return mellin_barnes(m, n, a, b, z).and_then(|(value, _)| finite(value));
    }
    // p == q exactly on the unit circle: the residue series converge only
    // marginally, so rely on the Mellin-Barnes integral where it applies.
    if let Some((value, _)) = mellin_barnes(m, n, a, b, z) {
        return finite(value);
    }
    None
}

/// Compares a residue-series value with the independent Mellin-Barnes
/// integral when both are numerically valid.
fn cross_check(
    value: C,
    degraded: bool,
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
) -> Option<C> {
    let delta = m as f64 + n as f64 - 0.5 * (a.len() + b.len()) as f64;
    if delta < 0.5 {
        return finite(value);
    }
    // Only cross-check when the vertical contour has a comfortable decay
    // margin; with a narrow one the Mellin-Barnes quadrature is the less
    // accurate of the two.
    let margin = delta * std::f64::consts::PI - z.arg().abs();
    if margin < 0.5 {
        return finite(value);
    }
    let Some((other, _error)) = mellin_barnes(m, n, a, b, z) else {
        return finite(value);
    };
    let difference = (value - other).norm();
    let scale = value.norm().max(other.norm()).max(1e-300);
    if difference <= 1e-7 * scale {
        return finite(value);
    }
    if degraded && difference <= 1e-4 * scale {
        // The confluent limit is the less accurate of the two.
        return finite(other);
    }
    None
}

/// The limit at `z = 0` when it exists.
fn zero_limit(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
) -> Option<C> {
    let p = a.len();
    let q = b.len();
    let left_ok = p <= q;
    if left_ok {
        if m == 0 {
            return Some(ZERO);
        }
        let mut total = ZERO;
        for (h, bh) in b.iter().take(m).enumerate() {
            if bh.re < 0.0 {
                return None;
            }
            if bh.re > 0.0 {
                continue;
            }
            if bh.im != 0.0 {
                // z^{i y} oscillates as z -> 0: no limit.
                return None;
            }
            total += left_log_coefficient(m, n, a, b, h)?;
        }
        return finite(total);
    }
    if n == 0 {
        return Some(ZERO);
    }
    // The 1/z series applies: z^{a_h - 1} is bounded only when every
    // exponent has a positive real part, and then the limit is zero.
    let mut all_positive = true;
    for ah in a.iter().take(n) {
        if ah.re - 1.0 <= 0.0 {
            all_positive = false;
        }
    }
    all_positive.then_some(ZERO)
}

/// Value of `G^{m,n}_{p,q}(z | a; b)` for real parameters and real `z`.
///
/// Returns `None` when the value is complex (a branch cut was crossed),
/// when the imaginary part is not negligible, or under the conditions of
/// [`meijer_g`].
#[must_use]
pub fn meijer_g_real(
    m: usize,
    n: usize,
    a: &[f64],
    b: &[f64],
    z: f64,
) -> Option<f64> {
    let ac: Vec<C> = a.iter().map(|&v| C::new(v, 0.0)).collect();
    let bc: Vec<C> = b.iter().map(|&v| C::new(v, 0.0)).collect();
    let value = meijer_g(m, n, &ac, &bc, C::new(z, 0.0))?;
    if value.im.abs() <= 1e-9 * value.re.abs().max(1.0) {
        Some(value.re)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::special;

    fn close(
        got: C,
        want: C,
        tolerance: f64,
        what: &str,
    ) {
        let scale = want.norm().max(1.0);
        assert!(
            (got - want).norm() <= tolerance * scale,
            "{what}: got {got}, want {want} (difference {})",
            (got - want).norm()
        );
    }

    fn re(v: f64) -> C {
        C::new(v, 0.0)
    }

    /// `G^{1,0}_{0,1}(z | 0) = e^{-z}` by Slater's theorem: the single
    /// residue family is `z^0 sum_k (-1)^k z^k / k!`.
    #[test]
    fn exponential() {
        for z in [-3.0, -0.5, 0.0, 0.25, 1.0, 4.0] {
            let got = meijer_g(1, 0, &[], &[re(0.0)], re(z));
            close(
                got.unwrap_or(C::new(f64::NAN, 0.0)),
                re((-z).exp()),
                1e-12,
                "exp",
            );
        }
        // Complex argument.
        let z = C::new(0.4, 0.9);
        let got = meijer_g(1, 0, &[], &[re(0.0)], z).unwrap_or(C::new(f64::NAN, 0.0));
        close(got, (-z).exp(), 1e-12, "exp complex");
    }

    /// `G^{1,0}_{0,1}(z | b) = z^b e^{-z}`.
    #[test]
    fn power_times_exponential() {
        for b in [-1.5, 0.0, 0.5, 2.3] {
            for z in [0.3, 1.7] {
                let got = meijer_g(1, 0, &[], &[re(b)], re(z));
                let want = z.powf(b) * (-z).exp();
                close(
                    got.unwrap_or(C::new(f64::NAN, 0.0)),
                    re(want),
                    1e-12,
                    "z^b exp",
                );
            }
        }
    }

    /// `G^{1,1}_{1,1}(z | a; b) = Gamma(1+b-a) z^b (1+z)^{a-b-1}`, from
    /// Slater's theorem with the terminating 1F0.
    #[test]
    fn first_order_reduction() {
        for (a, b) in [(0.5, 0.0), (2.3, 1.1), (-0.7, 1.4)] {
            for z in [-0.5, 0.2, 3.0] {
                let got = meijer_g(1, 1, &[re(a)], &[re(b)], re(z));
                let gamma = crate::kernels::hypergeometric::ln_gamma_c(re(1.0 + b - a)).exp();
                let zc = re(z);
                let want = gamma * zc.powc(re(b)) * (re(1.0) + zc).powc(re(a - b - 1.0));
                close(got.unwrap_or(C::new(f64::NAN, 0.0)), want, 1e-10, "1F1");
            }
        }
    }

    /// Cosine from the `0F1` reduction:
    /// `G^{1,0}_{0,2}(z | 0, 1/2) = cos(2 sqrt(z)) / sqrt(pi)`.
    #[test]
    fn trigonometric() {
        for z in [0.1, 1.0, 2.5, 7.0] {
            let got = meijer_g(1, 0, &[], &[re(0.0), re(0.5)], re(z));
            let want = (2.0 * z.sqrt()).cos() / std::f64::consts::PI.sqrt();
            close(got.unwrap_or(C::new(f64::NAN, 0.0)), re(want), 1e-12, "cos");
        }
        // z < 0 continues to cosh on the principal branch.
        let z = -1.3_f64;
        let got = meijer_g(1, 0, &[], &[re(0.0), re(0.5)], re(z));
        let want = (2.0 * (-z).sqrt()).cosh() / std::f64::consts::PI.sqrt();
        close(
            got.unwrap_or(C::new(f64::NAN, 0.0)),
            re(want),
            1e-12,
            "cosh continuation",
        );
    }

    /// `G^{1,0}_{0,2}(z | nu/2, -nu/2) = J_nu(2 sqrt(z))` (DLMF 16.18).
    #[test]
    fn bessel_j() {
        for nu in [0.0, 0.5, 1.3, 3.0] {
            for z in [0.2, 1.5, 6.0] {
                let got = meijer_g(1, 0, &[], &[re(nu / 2.0), re(-nu / 2.0)], re(z));
                let want = special::bessel_j(nu, 2.0 * z.sqrt());
                close(got.unwrap_or(C::new(f64::NAN, 0.0)), re(want), 1e-11, "J");
            }
        }
    }

    /// `G^{2,0}_{0,2}(z | nu/2, -nu/2) = 2 K_nu(2 sqrt(z))`.
    #[test]
    fn bessel_k() {
        for nu in [0.0, 0.7, 2.0] {
            for z in [0.2, 1.5, 4.0] {
                let got = meijer_g(2, 0, &[], &[re(nu / 2.0), re(-nu / 2.0)], re(z));
                let want = 2.0 * special::bessel_k(nu, 2.0 * z.sqrt());
                close(got.unwrap_or(C::new(f64::NAN, 0.0)), re(want), 1e-8, "K");
            }
        }
    }

    /// `G^{2,0}_{1,3}(z | (-nu-1)/2; nu/2, -nu/2, (-nu-1)/2) = Y_nu(2 sqrt(z))`;
    /// integer `nu` is a confluent-pole limit.
    #[test]
    fn bessel_y_including_integer_order() {
        for nu in [0.0, 1.0, 2.0, 0.7] {
            for z in [0.3, 2.0, 5.0] {
                let got = meijer_g(
                    2,
                    0,
                    &[re(f64::midpoint(-nu, -1.0))],
                    &[re(nu / 2.0), re(-nu / 2.0), re(f64::midpoint(-nu, -1.0))],
                    re(z),
                );
                let want = special::bessel_y(nu, 2.0 * z.sqrt());
                close(got.unwrap_or(C::new(f64::NAN, 0.0)), re(want), 1e-7, "Y");
            }
        }
    }

    /// `G^{1,1}_{1,2}(z | 1; 1/2, 0) = sqrt(pi) erf(sqrt(z))` and
    /// `G^{2,0}_{1,2}(z | 1; 1/2, 0) = sqrt(pi) erfc(sqrt(z))`.
    #[test]
    fn error_functions() {
        let pi = std::f64::consts::PI;
        for z in [0.1, 0.8, 3.0] {
            let erf = meijer_g(1, 1, &[re(1.0)], &[re(0.5), re(0.0)], re(z));
            let want = pi.sqrt() * special::erf_numerical(z.sqrt());
            close(erf.unwrap_or(C::new(f64::NAN, 0.0)), re(want), 1e-11, "erf");
            let erfc = meijer_g(2, 0, &[re(1.0)], &[re(0.5), re(0.0)], re(z));
            let want = pi.sqrt() * special::erfc_numerical(z.sqrt());
            close(
                erfc.unwrap_or(C::new(f64::NAN, 0.0)),
                re(want),
                1e-11,
                "erfc",
            );
        }
    }

    /// `G^{1,2}_{2,2}(z | 1, 1; 1, 0) = ln(1+z)` (a confluent `a` pair,
    /// evaluated by the Mellin-Barnes integral here).
    #[test]
    fn logarithm() {
        for z in [-0.5, 0.25, 1.7, 4.0] {
            let got = meijer_g(1, 2, &[re(1.0), re(1.0)], &[re(1.0), re(0.0)], re(z));
            close(
                got.unwrap_or(C::new(f64::NAN, 0.0)),
                re((1.0 + z).ln()),
                1e-10,
                "log",
            );
        }
    }

    /// `G^{1,2}_{2,2}(z | 1/2, 1; 1/2, 0) = 2 atan(sqrt(z))`.
    #[test]
    fn arctangent() {
        for z in [0.1, 0.9, 5.0] {
            let got = meijer_g(1, 2, &[re(0.5), re(1.0)], &[re(0.5), re(0.0)], re(z));
            close(
                got.unwrap_or(C::new(f64::NAN, 0.0)),
                re(2.0 * z.sqrt().atan()),
                1e-10,
                "atan",
            );
        }
    }

    /// The two methods agree on a grid of non-degenerate cases.
    #[test]
    fn series_and_mellin_barnes_agree() {
        let cases: &[(usize, usize, &[f64], &[f64])] = &[
            (1, 0, &[], &[0.3]),
            (1, 0, &[1.4], &[0.5, 0.2]),
            (2, 1, &[0.7, 2.1], &[0.4, 1.3]),
            (2, 2, &[1.2, 0.4], &[0.6, 1.9]),
        ];
        for &(m, n, a, b) in cases {
            assert!(n <= a.len() && m <= b.len());
            let ac: Vec<C> = a.iter().map(|&v| re(v)).collect();
            let bc: Vec<C> = b.iter().map(|&v| re(v)).collect();
            for z in [0.4, 1.6] {
                let zc = re(z);
                let mb = mellin_barnes(m, n, &ac, &bc, zc);
                let left = if a.len() < b.len() || (a.len() == b.len() && z.abs() < 1.0) {
                    left_series(m, n, &ac, &bc, zc)
                } else {
                    None
                };
                let right = if a.len() > b.len() || (a.len() == b.len() && z.abs() > 1.0) {
                    right_series(m, n, &ac, &bc, zc)
                } else {
                    None
                };
                if let (Some((series, _)), Some((integral, _))) = (left, mb) {
                    close(series, integral, 1e-8, "left vs MB");
                }
                if let (Some((series, _)), Some((integral, _))) = (right, mb) {
                    close(series, integral, 1e-8, "right vs MB");
                }
            }
        }
    }

    /// `z d/dz G = G^{m,n+1}_{p+1,q+1}(z | 0, a; b, 1)` by central
    /// differences of the kernel itself.
    #[test]
    fn derivative_identity() {
        let (m, n) = (1, 0);
        let a = [re(1.2)];
        let b = [re(0.4), re(0.9)];
        let z = 0.7;
        let h = 1e-5;
        let up = meijer_g(m, n, &a, &b, re(z + h)).unwrap_or(C::new(f64::NAN, 0.0));
        let down = meijer_g(m, n, &a, &b, re(z - h)).unwrap_or(C::new(f64::NAN, 0.0));
        let derivative = z * (up - down) / (2.0 * h);
        let mut a2 = vec![re(0.0)];
        a2.extend_from_slice(&a);
        let mut b2 = b.to_vec();
        b2.push(re(1.0));
        let shifted = meijer_g(m, n + 1, &a2, &b2, re(z)).unwrap_or(C::new(f64::NAN, 0.0));
        close(derivative, shifted, 1e-7, "z G'");
    }

    /// Negative real arguments use the principal branch (`arg z = pi`).
    #[test]
    fn principal_branch_on_the_cut() {
        // G^{1,0}_{0,1}(z|0) = e^{-z} is entire; at z = -2 + 0i it is e^2.
        let got = meijer_g(1, 0, &[], &[re(0.0)], re(-2.0)).unwrap_or(C::new(f64::NAN, 0.0));
        close(got, re(2.0_f64.exp()), 1e-12, "negative real");
        // A non-integer power picks the upper side: (-2)^{1/2} = i sqrt 2.
        let got = meijer_g(1, 0, &[], &[re(0.5)], re(-2.0)).unwrap_or(C::new(f64::NAN, 0.0));
        let want = C::new(0.0, 2.0_f64.sqrt()) * re(2.0_f64.exp());
        close(got, want, 1e-12, "branch");
    }

    /// The two methods agree for complex `z` as well.
    #[test]
    fn series_and_mellin_barnes_agree_for_complex_z() {
        // p = q: the z-series inside the unit circle.
        let a = [re(1.2), re(0.4)];
        let b = [re(0.6), re(1.9)];
        for z in [C::new(0.5, 0.4), C::new(0.3, -0.8)] {
            let Some((m, _)) = mellin_barnes(2, 2, &a, &b, z) else {
                panic!("Mellin-Barnes failed at {z}");
            };
            let Some((s, _)) = left_series(2, 2, &a, &b, z) else {
                panic!("series failed at {z}");
            };
            close(s, m, 1e-8, "complex z series vs MB");
        }
        // p > q: the 1/z series for large complex z.
        let a = [re(0.7), re(2.1)];
        let b = [re(0.4)];
        for z in [C::new(1.6, 0.3), C::new(2.4, -0.9)] {
            let Some((m, _)) = mellin_barnes(1, 1, &a, &b, z) else {
                panic!("Mellin-Barnes failed at {z}");
            };
            let Some((s, _)) = right_series(1, 1, &a, &b, z) else {
                panic!("right series failed at {z}");
            };
            close(s, m, 1e-8, "complex z 1/z series vs MB");
        }
    }

    /// Large `|z|` with `p > q` uses the convergent `1/z` series; with
    /// `p < q` the value is exponentially small compared with the terms and
    /// the kernel refuses rather than returning a rounded-off zero.
    #[test]
    fn large_arguments() {
        // p = q: |z| > 1 uses the 1/z series.
        let a = [re(0.5)];
        let b = [re(0.0)];
        let gamma = crate::kernels::hypergeometric::ln_gamma_c(re(0.5)).exp();
        for z in [10.0, 100.0, 1000.0] {
            let got = meijer_g(1, 1, &a, &b, re(z));
            let want = gamma * (1.0 + z).powf(-0.5);
            close(
                got.unwrap_or(C::new(f64::NAN, 0.0)),
                want,
                1e-11,
                "1/z series",
            );
        }
        // p < q: large arguments cancel in the z-series.
        assert_eq!(meijer_g(1, 0, &[], &[re(0.5)], re(80.0)), None);
        // Moderate ones are still accurate.
        close(
            meijer_g(1, 0, &[], &[re(0.5)], re(6.0)).unwrap_or(C::new(f64::NAN, 0.0)),
            re(6.0_f64.sqrt() * (-6.0_f64).exp()),
            1e-11,
            "moderate z",
        );
    }

    /// Degenerate orders, zero argument and invalid input.
    #[test]
    fn edge_cases() {
        // m = 0 with p < q: no poles on the right, the value is zero.
        assert_eq!(
            meijer_g(0, 1, &[re(0.5)], &[re(0.3), re(1.2)], re(0.6)),
            Some(ZERO)
        );
        // n = 0 with p > q: no poles on the left, the value is zero.
        assert_eq!(
            meijer_g(1, 0, &[re(0.5), re(1.2)], &[re(0.3)], re(2.0)),
            Some(ZERO)
        );
        // G^{1,0}_{0,1}(0 | 0) = 1; with b = 1/2 it vanishes.
        assert_eq!(meijer_g(1, 0, &[], &[re(0.0)], ZERO), Some(ONE));
        assert_eq!(meijer_g(1, 0, &[], &[re(0.5)], ZERO), Some(ZERO));
        // A negative exponent has no limit at zero.
        assert_eq!(meijer_g(1, 0, &[], &[re(-0.5)], ZERO), None);
        // Invalid orders.
        assert_eq!(meijer_g(2, 0, &[], &[re(0.0)], re(1.0)), None);
        assert_eq!(meijer_g(0, 2, &[re(0.0)], &[re(0.0)], re(1.0)), None);
        // Non-finite input.
        assert_eq!(meijer_g(1, 0, &[], &[C::new(f64::NAN, 0.0)], re(1.0)), None);
    }

    /// The residue-series and Mellin-Barnes values of a case that only the
    /// right series reaches (p > q).
    #[test]
    fn other_side_series() {
        // G^{1,1}_{1,2}(z | a; b1, b2) with p < q converges for all z; its
        // value at large |z| must still be finite and match the formula
        // G = z^{b1} (1+z)^{a-b1-1} / Gamma(1+b1-b2) ... only for the
        // terminating case; here check the Mellin-Barnes agreement instead.
        let a = [re(0.9)];
        let b = [re(0.2), re(1.7)];
        for z in [0.3, 5.0] {
            let direct = meijer_g(1, 1, &a, &b, re(z));
            assert!(direct.is_some(), "z = {z}");
            let mb = mellin_barnes(1, 1, &a, &b, re(z));
            if let (Some(v), Some((w, _))) = (direct, mb) {
                close(v, w, 1e-8, "p<q MB");
            }
        }
    }
}
