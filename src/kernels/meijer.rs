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
//! * **The diagonal `p = q`.** The two residue series cover the two sides of
//!   the unit circle: the `z` series for `0 < |z| < 1` and the `1/z` series
//!   for `|z| > 1`. They are used directly, without the Mellin-Barnes
//!   condition `delta = m + n - (p+q)/2 > 0`; in particular the common
//!   `m + n <= q` cases, including residue points where a denominator gamma
//!   vanishes, are evaluated rather than falling through to a divergent
//!   integral. An empty side (`m = 0` inside, `n = 0` outside) is the zero
//!   function there, which is the value selected by the corresponding
//!   contour: `G^{1,0}_{1,1}` is non-zero for `|z| < 1` and zero for
//!   `|z| > 1`, while `G^{0,1}_{1,1}` is zero inside and non-zero outside.
//!   Exactly `|z| = 1` is the common boundary of the two series; the value is
//!   taken from the Mellin-Barnes integral when `delta > 0`, and otherwise
//!   the request stays unevaluated because the radial limits from the two
//!   sides need not agree.
//! * **Confluent poles.** When two entries of the contour-defining `b` set
//!   (for the `z -> 0` series) or `a` set (for the `z -> infinity` series)
//!   differ by an integer, the poles of the integrand merge and the
//!   separate-residue recurrence is undefined. The residue at a pole of
//!   order `r` is then the coefficient of `t^{r-1}` in the Taylor expansion
//!   of `t^r f(p+t)`. The gamma factors are differentiated exactly through
//!   the logarithmic derivatives of `Gamma` (the Lanczos formula
//!   differentiated term by term, i.e. the digamma/polygamma hierarchy), and
//!   `t Gamma(-k + t)` is evaluated from the reflection formula instead of as
//!   `0 * infinity`. This is an exact confluent evaluation rather than a
//!   finite-difference limit, and it has the same accuracy as the ordinary
//!   residue path. The same expansion covers residue points at which a
//!   denominator gamma vanishes: the reciprocal gamma is entire, so its zero
//!   is a genuine zero of that residue term.
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
//!   (ii) when the residue formula is genuinely inapplicable: a pinched
//!   contour (a left and a right pole coincide, so the G-function itself is
//!   singular), or `p = q` exactly on `|z| = 1`.
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
//! `|z| = 1` when `p = q`) and no family sum cancels. Confluent contour
//! parameters (repeated or integer-spaced `a`/`b`) do not degrade this: they
//! are evaluated through the exact residue expansion described above, whose
//! only error comes from the gamma function itself. The Mellin-Barnes
//! fallback is accurate to about `1e-10` when the effective decay
//! `delta - |arg z|/pi` is at least `0.2`, and less (or unavailable) close to
//! the boundary of the convergence sector. The result is `None` when
//!
//! * `m > q` or `n > p`, a parameter or `z` is not finite, or `z` is
//!   exactly `0` and the limit is not defined;
//! * the contour is pinched (some `a_j - b_k` is a positive integer with
//!   `j <= n`, `k <= m`), so the G-function has no finite value at those
//!   parameters;
//! * a denominator gamma vanishes identically along a whole residue family
//!   (`a_j - b_h` a non-positive integer for some `j > n`, or symmetrically
//!   `a_h - b_j` for some `j > m`). The residue series is then identically
//!   zero, while the G-function is the limit of the non-degenerate parameter
//!   values; no value is returned unless the Mellin-Barnes integral
//!   converges;
//! * the applicable series does not converge numerically (for `p = q`,
//!   `|z|` within about `1e-5` of `1`, away from the Mellin-Barnes regime);
//! * a residue family overflows, the families cancel by more than 13
//!   digits, or the terms of one family cancel by more than 13 digits
//!   against their peak (this is what happens for `p < q` and very large
//!   `|z|`, where the value is exponentially smaller than the intermediate
//!   terms; no rescaling across the peak is attempted, so `None` is
//!   returned rather than a zero);
//! * `p = q` with `|z| = 1` outside the Mellin-Barnes regime: both residue
//!   series are then on their circle of convergence and the radial limits
//!   from the two sides need not exist, or need not agree;
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

// -------------------------------------------------- exact confluent residues
//
// When contour parameters coincide (or are an integer apart), several gamma
// factors have poles at the same point and the residue is that of the whole
// integrand, not the product of the separate residues.  The residue at a pole
// `p` of order `r` is the coefficient of `t^(r-1)` in the Taylor expansion of
// `t^r f(p+t)`; the gamma factors are expanded exactly with the logarithmic
// derivatives of `Gamma` (the Lanczos formula differentiated term by term,
// i.e. the digamma hierarchy).  This replaces any finite-difference or
// Richardson extrapolation and is accurate to the precision of the underlying
// gamma evaluation.

/// Lanczos coefficients; kept identical to `hypergeometric::ln_gamma_c` so
/// that the confluent Taylor series and the ordinary gamma values agree.
const LANCZOS: [f64; 9] = [
    0.999_999_999_999_809_9,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_1,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

/// Largest pole order handled by the exact residue expansion.
const MAX_POLE_ORDER: usize = 64;

/// A truncated power series in `t`; all series passed to [`ser_mul`] have the
/// same length.
type Ser = Vec<C>;

/// Truncated product of two series.
fn ser_mul(
    x: &Ser,
    y: &Ser,
) -> Ser {
    let n = x.len().min(y.len());
    let mut out = vec![ZERO; n];
    for i in 0..n {
        if is_zero(x[i]) {
            continue;
        }
        for j in 0..(n - i) {
            out[i + j] += x[i] * y[j];
        }
    }
    out
}

fn ser_neg(x: &Ser) -> Ser {
    x.iter().map(|v| -*v).collect()
}

fn ser_sub(
    x: &Ser,
    y: &Ser,
) -> Ser {
    x.iter().zip(y).map(|(a, b)| *a - *b).collect()
}

fn ser_scale(
    x: &Ser,
    factor: C,
) -> Ser {
    x.iter().map(|v| *v * factor).collect()
}

/// `exp(x)` as a series (the recurrence `n w_n = sum_{k=1}^n k x_k w_{n-k}`).
fn ser_exp(x: &Ser) -> Option<Ser> {
    let n = x.len();
    let mut w = vec![ZERO; n];
    w[0] = x[0].exp();
    finite(w[0])?;
    for k in 1..n {
        let mut acc = ZERO;
        for j in 1..=k {
            acc += x[j] * w[k - j] * (j as f64);
        }
        w[k] = acc * (1.0 / k as f64);
    }
    w.iter()
        .all(|v| v.re.is_finite() && v.im.is_finite())
        .then_some(w)
}

/// `1/x` as a series; `x[0]` must be non-zero.
fn ser_recip(x: &Ser) -> Option<Ser> {
    let n = x.len();
    if x[0].norm() == 0.0 {
        return None;
    }
    let inv0 = ONE / x[0];
    let mut out = vec![ZERO; n];
    out[0] = inv0;
    for k in 1..n {
        let mut acc = ZERO;
        for j in 1..=k {
            acc += x[j] * out[k - j];
        }
        out[k] = -acc * inv0;
    }
    out.iter()
        .all(|v| v.re.is_finite() && v.im.is_finite())
        .then_some(out)
}

/// `ln(1 + x)` as a series, for `x[0] = 0`.
fn ser_log1p(x: &Ser) -> Ser {
    let n = x.len();
    let mut out = vec![ZERO; n];
    for m in 0..n.saturating_sub(1) {
        let mut acc = ZERO;
        for k in 1..=m {
            acc += x[k] * out[m - k + 1] * ((m - k + 1) as f64);
        }
        out[m + 1] = x[m + 1] - acc * (1.0 / (m + 1) as f64);
    }
    out
}

/// `ln(c + sign * t)` as a series.
fn ser_ln_linear(
    c: C,
    sign: f64,
    n: usize,
) -> Option<Ser> {
    if c.norm() == 0.0 {
        return None;
    }
    let mut u = vec![ZERO; n];
    if n > 1 {
        u[1] = C::new(sign, 0.0) / c;
    }
    let mut out = ser_log1p(&u);
    out[0] += c.ln();
    Some(out)
}

/// Taylor coefficients of `ln Gamma(beta + sign * t)` from the Lanczos
/// formula, for `Re beta >= 2`.  Differentiating the closed form gives the
/// digamma/polygamma hierarchy exactly; no finite differences are involved.
fn ln_gamma_lanczos_series(
    beta: C,
    sign: f64,
    n: usize,
) -> Option<Ser> {
    if beta.re < 2.0 {
        return None;
    }
    // x(t) = a_0 + sum_{i=1}^{8} a_i / (c_i + sign t), c_i = beta - 1 + i.
    let mut x = vec![ZERO; n];
    x[0] = C::new(LANCZOS[0], 0.0);
    for i in 1..LANCZOS.len() {
        let c = beta - 1.0 + i as f64;
        let inv = ONE / c;
        let mut term = C::new(LANCZOS[i], 0.0) * inv;
        for coefficient in &mut x {
            *coefficient += term;
            term = term * inv * C::new(-sign, 0.0);
        }
    }
    let x0 = x[0];
    if x0.norm() < 1e-250 {
        return None;
    }
    let mut u = vec![ZERO; n];
    for (l, coefficient) in u.iter_mut().enumerate().skip(1) {
        *coefficient = x[l] / x0;
    }
    let mut out = ser_log1p(&u);
    out[0] += x0.ln();

    // ln(t0 + sign t), t0 = beta + 6.5.
    let t0 = beta + 6.5;
    let mut ln_t = vec![ZERO; n];
    ln_t[0] = t0.ln();
    let ratio = C::new(sign, 0.0) / t0;
    let mut power = ONE;
    for l in 1..n {
        power *= ratio;
        let magnitude = if l % 2 == 1 {
            1.0 / l as f64
        } else {
            -1.0 / l as f64
        };
        ln_t[l] = power * magnitude;
    }

    // ln Gamma = (beta - 1/2 + sign t) ln(t0 + sign t) - (t0 + sign t)
    //            + ln x + (1/2) ln(2 pi).
    for (l, value) in out.iter_mut().enumerate() {
        let mut term = (beta - 0.5) * ln_t[l];
        if l > 0 {
            term += ln_t[l - 1] * sign;
        }
        *value += term;
    }
    out[0] -= t0;
    if n > 1 {
        out[1] -= sign;
    }
    out[0] += 0.5 * (2.0 * std::f64::consts::PI).ln();
    Some(out)
}

/// Taylor coefficients of `ln Gamma(alpha + sign * t)`.
///
/// For `Re alpha < 3` the recurrence `Gamma(z+1) = z Gamma(z)` shifts the
/// argument to a region where the Lanczos series is accurate; the shift is
/// analytic, so the Taylor coefficients are exact up to the Lanczos error.
fn ln_gamma_series(
    alpha: C,
    sign: f64,
    n: usize,
) -> Option<Ser> {
    let shift = (3.0 - alpha.re).ceil().max(0.0);
    if !shift.is_finite() || shift > 5000.0 {
        return None;
    }
    let count = shift as usize;
    let beta = alpha + C::new(shift, 0.0);
    let mut out = ln_gamma_lanczos_series(beta, sign, n)?;
    for k in 0..count {
        let c = alpha + C::new(k as f64, 0.0);
        let term = ser_ln_linear(c, sign, n)?;
        out = ser_sub(&out, &term);
    }
    Some(out)
}

/// Taylor coefficients of `Gamma(alpha + sign * t)`; `None` at a pole.
fn gamma_series(
    alpha: C,
    sign: f64,
    n: usize,
) -> Option<Ser> {
    if nonpositive_integer(alpha) {
        return None;
    }
    ser_exp(&ln_gamma_series(alpha, sign, n)?)
}

/// The non-positive integer `-l` equal to `v`, if any.
fn nonpositive_integer_index(v: C) -> Option<usize> {
    if v.im.abs() > 1e-9 || !near_integer(v.re) || v.re > 1e-9 {
        return None;
    }
    let l = -v.re.round();
    if !(0.0..=1.0e6).contains(&l) {
        return None;
    }
    usize::try_from(l as i64).ok()
}

/// Taylor coefficients of `1 / Gamma(alpha + sign * t)`.  The reciprocal
/// gamma is entire, so its zeros at non-positive integers are handled by the
/// recurrence rather than by the (singular) logarithmic derivative.
fn recip_gamma_series(
    alpha: C,
    sign: f64,
    n: usize,
) -> Option<Ser> {
    if let Some(l) = nonpositive_integer_index(alpha) {
        if l > 512 {
            return None;
        }
        let mut product = vec![ZERO; n];
        product[0] = ONE;
        for i in 0..=l {
            let c = i as f64 - l as f64;
            let mut linear = vec![ZERO; n];
            linear[0] = C::new(c, 0.0);
            if n > 1 {
                linear[1] = C::new(sign, 0.0);
            }
            product = ser_mul(&product, &linear);
        }
        let log_gamma = ln_gamma_series(ONE, sign, n)?;
        let inverse = ser_exp(&ser_neg(&log_gamma))?;
        return Some(ser_mul(&product, &inverse));
    }
    let log_gamma = ln_gamma_series(alpha, sign, n)?;
    ser_exp(&ser_neg(&log_gamma))
}

/// Taylor coefficients of `t Gamma(-k + sign * t)` (the residue factor of a
/// `k`-th order pole).  From the recurrence and the reflection formula
/// ```text
/// t Gamma(-k + sign t) = sign (-1)^k (pi t / sin(pi t))
///                        / Gamma(k + 1 - sign t),
/// ```
/// which is analytic at `t = 0` and avoids any `0 * infinity`.
fn pole_factor_series(
    k: usize,
    sign: f64,
    n: usize,
) -> Option<Ser> {
    // pi t / sin(pi t) = 1 / sinc(pi t).
    let mut sinc = vec![ZERO; n];
    let pi_squared = std::f64::consts::PI * std::f64::consts::PI;
    let mut term = 1.0_f64;
    let mut degree = 0usize;
    while degree < n {
        sinc[degree] = C::new(term, 0.0);
        // The next even coefficient of sin(pi t)/(pi t).
        let l = degree / 2;
        term *= -pi_squared / (((2 * l + 2) * (2 * l + 3)) as f64);
        degree += 2;
    }
    let quotient = ser_recip(&sinc)?;
    let gamma = recip_gamma_series(C::new((k + 1) as f64, 0.0), -sign, n)?;
    let parity = if k.is_multiple_of(2) {
        1.0
    } else {
        -1.0
    };
    Some(ser_scale(
        &ser_mul(&quotient, &gamma),
        C::new(sign * parity, 0.0),
    ))
}

/// `v` as a non-negative integer, with the tolerance of the parameters.
fn nonnegative_integer(v: C) -> Option<usize> {
    if v.im.abs() > 1e-9 * v.re.abs().max(1.0) || !near_integer(v.re) {
        return None;
    }
    let k = v.re.round();
    if !(-0.5..=1.0e6).contains(&k) {
        return None;
    }
    usize::try_from(k as i64).ok()
}

/// Residue of the Mellin-Barnes integrand at `s = p`, or `None` if `p` is not
/// a pole or if the contour is pinched there (a left and a right pole
/// coincide, so the G-function itself is singular).
fn residue_at(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    p: C,
    z: C,
) -> Option<C> {
    let mut b_poles: Vec<usize> = Vec::new();
    for &bj in b.iter().take(m) {
        if let Some(k) = nonnegative_integer(p - bj) {
            b_poles.push(k);
        }
    }
    let mut a_poles: Vec<usize> = Vec::new();
    for &aj in a.iter().take(n) {
        if let Some(k) = nonnegative_integer(aj - ONE - p) {
            a_poles.push(k);
        }
    }
    if b_poles.is_empty() && a_poles.is_empty() {
        return None;
    }
    if !b_poles.is_empty() && !a_poles.is_empty() {
        // Pinched contour: the value is not analytic at these parameters.
        return None;
    }
    let order = b_poles.len() + a_poles.len();
    if order > MAX_POLE_ORDER {
        return None;
    }
    let ncoef = order;
    // z^{p+t} = z^p sum_l (ln z)^l t^l / l!.
    let ln_z = z.ln();
    let mut series = vec![ZERO; ncoef];
    let mut term = (p * ln_z).exp();
    let mut factorial = 1.0_f64;
    for (l, coefficient) in series.iter_mut().enumerate() {
        *coefficient = term / factorial;
        term *= ln_z;
        factorial *= (l + 1) as f64;
    }
    for &k in &b_poles {
        series = ser_mul(&series, &pole_factor_series(k, -1.0, ncoef)?);
    }
    for &k in &a_poles {
        series = ser_mul(&series, &pole_factor_series(k, 1.0, ncoef)?);
    }
    for (j, &bj) in b.iter().enumerate() {
        if j < m {
            if nonnegative_integer(p - bj).is_some() {
                continue;
            }
            series = ser_mul(&series, &gamma_series(bj - p, -1.0, ncoef)?);
        } else {
            let factor = recip_gamma_series(ONE - bj + p, 1.0, ncoef)?;
            series = ser_mul(&series, &factor);
        }
    }
    for (j, &aj) in a.iter().enumerate() {
        if j < n {
            if nonnegative_integer(aj - ONE - p).is_some() {
                continue;
            }
            series = ser_mul(&series, &gamma_series(ONE - aj + p, 1.0, ncoef)?);
        } else {
            let factor = recip_gamma_series(aj - p, -1.0, ncoef)?;
            series = ser_mul(&series, &factor);
        }
    }
    finite(series[ncoef - 1])
}

/// Groups parameters that differ by integers, returning the lowest member of
/// each class and the largest integer offset from it.
fn integer_classes(values: &[C]) -> Vec<(C, usize)> {
    let mut used = vec![false; values.len()];
    let mut out = Vec::new();
    for i in 0..values.len() {
        if used[i] {
            continue;
        }
        used[i] = true;
        let mut group = vec![i];
        for j in (i + 1)..values.len() {
            if used[j] {
                continue;
            }
            let difference = values[j] - values[i];
            if near_integer(difference.re) && near_integer(difference.im) {
                used[j] = true;
                group.push(j);
            }
        }
        let mut base = values[group[0]];
        for &g in &group {
            if values[g].re < base.re {
                base = values[g];
            }
        }
        let mut highest = 0usize;
        for &g in &group {
            let offset = (values[g].re - base.re).round();
            if offset > highest as f64 {
                highest = offset as usize;
            }
        }
        out.push((base, highest));
    }
    out
}

/// Sums the residues at one integer class of pole locations.
fn class_residues(
    mn: (usize, usize),
    a: &[C],
    b: &[C],
    z: C,
    base: C,
    highest: usize,
    descending: bool,
) -> Option<C> {
    let (m, n) = mn;
    let start = if descending {
        base + (highest as f64) - 1.0
    } else {
        base
    };
    let step = if descending { -1.0 } else { 1.0 };
    let mut total = ZERO;
    let mut absolute = 0.0_f64;
    let mut peak = 0.0_f64;
    for k in 0..MAX_TERMS {
        let p = start + C::new(step * f64::from(k), 0.0);
        let term = residue_at(m, n, a, b, p, z)?;
        total += term;
        let size = term.norm();
        absolute += size;
        peak = peak.max(size);
        if k > 2 && peak > 0.0 && size <= 1e-18 * peak {
            if total.norm() <= LOST_DIGITS * absolute {
                return None;
            }
            return finite(total);
        }
    }
    None
}

/// Exact residue series.  `descending` selects the `1/z` series (poles of the
/// `a` set, closed to the left); otherwise the `z` series (poles of the `b`
/// set, closed to the right).
///
/// The contour orientation contributes a sign: closing to the right is
/// clockwise and gives `G = -sum residues`, closing to the left is
/// counter-clockwise and gives `G = +sum residues`.
fn exact_residue_series(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
    z: C,
    descending: bool,
) -> Option<C> {
    let values = if descending {
        &a[..n]
    } else {
        &b[..m]
    };
    if values.is_empty() {
        return Some(ZERO);
    }
    let classes = integer_classes(values);
    let sum = combine_families(classes.len(), |index| {
        let (base, highest) = classes[index];
        class_residues((m, n), a, b, z, base, highest, descending)
    })?;
    Some(if descending { sum } else { -sum })
}

/// A denominator gamma of the `z -> 0` formula vanishes at *every* residue
/// point of a family (`a_j - b_h` a non-positive integer with `j > n`).  The
/// whole family's residues are then zero, while the G-function itself is the
/// limit of the non-degenerate parameter values and need not vanish; the
/// residue series is not a representation of that limit.  Callers must fall
/// back to the Mellin-Barnes integral (or report the value as unavailable)
/// rather than return the zero sum.
fn left_family_identically_zero(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
) -> bool {
    for bh in b.iter().take(m) {
        for aj in a.iter().skip(n) {
            if nonpositive_integer(*aj - *bh) {
                return true;
            }
        }
    }
    false
}

/// The same obstruction for the `z -> infinity` series.
fn right_family_identically_zero(
    m: usize,
    n: usize,
    a: &[C],
    b: &[C],
) -> bool {
    for ah in a.iter().take(n) {
        for bj in b.iter().skip(m) {
            if nonpositive_integer(*ah - *bj) {
                return true;
            }
        }
    }
    false
}

/// Left residue series, with an exact confluent evaluation when needed.
/// Returns `None` when the residue formula is not applicable (a numerator
/// gamma is pinched by the contour, or a whole family of residues vanishes
/// identically), so that the caller can fall back to the Mellin-Barnes
/// integral.
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
    if left_family_identically_zero(m, n, a, b) {
        return None;
    }
    // A merging of contour poles (integer-spaced `b`), or a denominator
    // gamma that vanishes at an isolated residue point, makes the
    // separate-residue recurrence invalid.  Evaluate the actual residues of
    // the integrand instead; this is exact and also covers `p = q`, where
    // the Mellin-Barnes fallback is unavailable.
    if confluent(&b[..m]) || left_denominator_singular(m, n, a, b) {
        return exact_residue_series(m, n, a, b, z, false).map(|v| (v, true));
    }
    if left_numerator_pinch(m, n, a, b) {
        return None;
    }
    left_series_raw(m, n, a, b, z)
        .map(|v| (v, false))
        .or_else(|| exact_residue_series(m, n, a, b, z, false).map(|v| (v, true)))
}

/// Right residue series, with an exact confluent evaluation when needed.
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
    if right_family_identically_zero(m, n, a, b) {
        return None;
    }
    if confluent(&a[..n]) || right_denominator_singular(m, n, a, b) {
        return exact_residue_series(m, n, a, b, z, true).map(|v| (v, true));
    }
    if right_numerator_pinch(m, n, a, b) {
        return None;
    }
    right_series_raw(m, n, a, b, z)
        .map(|v| (v, false))
        .or_else(|| exact_residue_series(m, n, a, b, z, true).map(|v| (v, true)))
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

    /// Hard-coded references from mpmath 1.3.0 at 60 digits.  These cover
    /// repeated `b`, integer-spaced `b`, triple poles, repeated `a`, the
    /// integer-order Bessel `Y` case, and the `p = q` denominator-gamma
    /// singularity.  The exact confluent evaluation must reproduce them to a
    /// relative error of `1e-11`; the former Richardson shim was only good to
    /// `1e-9`--`1e-7`.
    #[test]
    fn confluent_poles_match_mpmath() {
        let cases: &[(usize, usize, &[f64], &[f64], f64, f64)] = &[
            // Repeated b: G^{2,0}_{0,2}(z | 0, 0) = 2 K_0(2 sqrt z).
            (2, 0, &[], &[0.0, 0.0], 0.25, 0.842_048_876_481_416_666_7),
            (2, 0, &[], &[0.0, 0.0], 2.0, 0.084_783_547_996_802_991_543),
            // Repeated b at b = 1/2.
            (2, 0, &[], &[0.5, 0.5], 0.3, 0.403_050_315_801_433_244_41),
            (2, 0, &[], &[0.5, 0.5], 1.5, 0.162_155_520_296_842_533_77),
            // Integer-spaced but distinct b.
            (2, 0, &[], &[0.0, 1.0], 0.3, 0.562_569_856_924_219_936_87),
            (2, 0, &[], &[0.0, 1.0], 2.0, 0.139_667_474_015_293_142_86),
            // A third-order pole.
            (
                3,
                0,
                &[],
                &[0.0, 0.0, 0.0],
                0.3,
                0.635_055_911_879_765_542_83,
            ),
            (
                3,
                0,
                &[],
                &[0.0, 0.0, 0.0],
                2.0,
                0.060_771_033_962_088_616_491,
            ),
            // Integer-order Bessel Y is a confluent b pair.
            (
                2,
                0,
                &[-0.5],
                &[0.0, 0.0, -0.5],
                0.3,
                0.158_975_073_307_206_066_11,
            ),
            (
                2,
                0,
                &[-0.5],
                &[0.0, 0.0, -0.5],
                2.0,
                0.428_287_398_117_322_672_15,
            ),
            (
                2,
                0,
                &[-1.0],
                &[0.5, -0.5, -1.0],
                0.3,
                -0.701_755_262_444_025_362_85,
            ),
            (
                2,
                0,
                &[-1.0],
                &[0.5, -0.5, -1.0],
                2.0,
                0.273_118_770_266_153_836_43,
            ),
            (
                2,
                0,
                &[-1.5],
                &[1.0, -1.0, -1.5],
                0.3,
                -1.440_199_030_269_041_488_6,
            ),
            (
                2,
                0,
                &[-1.5],
                &[1.0, -1.0, -1.5],
                2.0,
                -0.235_163_263_592_794_489_16,
            ),
            // Repeated a: G^{1,2}_{2,2}(z | 1, 1; 1, 0) = ln(1 + z).
            (
                1,
                2,
                &[1.0, 1.0],
                &[1.0, 0.0],
                0.5,
                0.405_465_108_108_164_381_98,
            ),
            (
                1,
                2,
                &[1.0, 1.0],
                &[1.0, 0.0],
                2.0,
                1.098_612_288_668_109_691_4,
            ),
            (
                1,
                2,
                &[1.0, 1.0],
                &[1.0, 0.0],
                4.0,
                1.609_437_912_434_100_374_6,
            ),
            // Repeated a in the left contour set for p > q.
            (
                1,
                2,
                &[0.5, 0.5, 1.7],
                &[0.3],
                0.3,
                0.959_091_063_085_557_961_69,
            ),
            (
                1,
                2,
                &[0.5, 0.5, 1.7],
                &[0.3],
                2.0,
                1.119_095_137_359_783_299_3,
            ),
            // Repeated a for p = q.
            (
                2,
                2,
                &[1.0, 1.0],
                &[0.3, 2.2],
                0.4,
                6.368_721_541_069_509_005_0,
            ),
            (
                2,
                2,
                &[1.0, 1.0],
                &[0.3, 2.2],
                2.0,
                9.763_775_708_768_645_719_8,
            ),
            // A denominator gamma vanishes at a residue point: b_2 - b_1 is
            // an integer, and the simple residue product would be 0 * inf.
            (
                1,
                1,
                &[0.5],
                &[0.0, 1.0],
                0.25,
                -0.184_043_558_896_726_299_82,
            ),
            (
                1,
                1,
                &[0.5],
                &[0.0, 1.0],
                2.0,
                -0.457_025_793_745_633_966_82,
            ),
            (2, 1, &[0.5], &[0.0, 1.0], 0.25, 1.411_754_207_277_880_155_8),
            (2, 1, &[0.5], &[0.0, 1.0], 2.0, 0.871_498_554_376_879_410_74),
        ];
        for &(m, n, a, b, z, want) in cases {
            let ac: Vec<C> = a.iter().map(|&v| re(v)).collect();
            let bc: Vec<C> = b.iter().map(|&v| re(v)).collect();
            let got = meijer_g(m, n, &ac, &bc, re(z)).unwrap_or(C::new(f64::NAN, 0.0));
            let scale = want.abs().max(1e-300);
            assert!(
                (got - re(want)).norm() <= 1e-11 * scale,
                "G^{{{m},{n}}} a={a:?} b={b:?} z={z}: got {got}, want {want}"
            );
        }
    }

    /// `p = q` on both sides of the unit circle, including the four decades
    /// required of the boundary fix.  `G^{1,1}_{1,1}` has `delta = 1`, while
    /// `G^{1,1}_{2,2}` with `b = (0, 1)` has `delta = 0` and a vanishing
    /// denominator gamma; the latter used to fall through to a Mellin-Barnes
    /// integral that cannot converge.
    #[test]
    fn diagonal_p_q_matches_mpmath() {
        let cases: &[(usize, usize, &[f64], &[f64], f64, f64)] = &[
            (1, 1, &[0.5], &[0.2], 0.25, 0.841_480_775_090_993_030_31),
            (1, 1, &[0.5], &[0.2], 2.0, 0.691_057_725_642_243_145_75),
            (1, 1, &[0.5], &[0.2], 1e-4, 0.205_713_506_285_645_853_95),
            (1, 1, &[0.5], &[0.2], 1e4, 0.012_979_644_764_970_066_240),
            (
                1,
                1,
                &[0.5],
                &[0.0, 1.0],
                0.25,
                -0.184_043_558_896_726_299_82,
            ),
            (
                1,
                1,
                &[0.5],
                &[0.0, 1.0],
                2.0,
                -0.457_025_793_745_633_966_82,
            ),
            (
                1,
                1,
                &[0.5],
                &[0.0, 1.0],
                1e-4,
                -0.000_088_616_046_120_272_746_723,
            ),
            (
                1,
                1,
                &[0.5],
                &[0.0, 1.0],
                1e4,
                -0.005_000_375_070_333_015_891_5,
            ),
            (
                2,
                2,
                &[1.0, 1.0],
                &[0.3, 2.2],
                0.4,
                6.368_721_541_069_509_005_0,
            ),
            (
                2,
                2,
                &[1.0, 1.0],
                &[0.3, 2.2],
                2.0,
                9.763_775_708_768_645_719_8,
            ),
            (
                2,
                2,
                &[1.0, 1.0],
                &[0.3, 2.2],
                1e-4,
                0.543_079_051_556_264_931_37,
            ),
            (
                2,
                2,
                &[1.0, 1.0],
                &[0.3, 2.2],
                1e4,
                36.305_913_650_456_759_948,
            ),
            // The empty-side cases follow the contour: G^{1,0}_{1,1} is
            // non-zero inside the circle and zero outside (mpmath 1.3.0 and
            // SymPy agree on this branch).
            (1, 0, &[0.5], &[0.2], 0.25, 0.309_846_143_397_960_091_42),
            (1, 0, &[0.5], &[0.2], 2.0, 0.0),
            (1, 0, &[0.5], &[0.2], 1e-4, 0.052_982_369_817_943_530_816),
            (1, 0, &[0.5], &[0.2], 1e4, 0.0),
        ];
        for &(m, n, a, b, z, want) in cases {
            let ac: Vec<C> = a.iter().map(|&v| re(v)).collect();
            let bc: Vec<C> = b.iter().map(|&v| re(v)).collect();
            let got = meijer_g(m, n, &ac, &bc, re(z)).unwrap_or(C::new(f64::NAN, 0.0));
            let scale = want.abs().max(1e-300);
            assert!(
                (got - re(want)).norm() <= 1e-11 * scale,
                "G^{{{m},{n}}} a={a:?} b={b:?} z={z}: got {got}, want {want}"
            );
        }
    }

    /// `G^{1,1}_{1,1} = Gamma(1+b-a) z^b (1+z)^{a-b-1}` on both sides of the
    /// branch point at `z = -1`, on both sides of `|z| = 1`, and at large
    /// `|z|`.
    #[test]
    fn first_order_on_both_sides_of_the_branch() {
        let (a, b) = (0.5, 0.2);
        let gamma = crate::kernels::hypergeometric::ln_gamma_c(re(1.0 + b - a)).exp();
        let formula = |z: f64| gamma * re(z).powc(re(b)) * (re(1.0) + re(z)).powc(re(a - b - 1.0));
        // mpmath 1.3.0 reference at z = -2 (upper side of the cut).
        let cut_reference = C::new(
            1.889_976_881_964_643_485_7e-54,
            -1.491_074_025_307_378_697_9,
        );
        for z in [2.0, 1e4, -0.9] {
            let got = meijer_g(1, 1, &[re(a)], &[re(b)], re(z)).unwrap_or(C::new(f64::NAN, 0.0));
            close(got, formula(z), 1e-12, "branch sides");
        }
        let got = meijer_g(1, 1, &[re(a)], &[re(b)], re(-2.0)).unwrap_or(C::new(f64::NAN, 0.0));
        close(got, cut_reference, 1e-12, "branch cut");
        close(formula(-2.0), cut_reference, 1e-12, "formula on the cut");
    }

    /// `G^{1,2}_{2,2}(z | 1, 1; 1, 0) = ln(1+z)` (a confluent `a` pair).
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
