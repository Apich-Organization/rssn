//! Laurent series via FFT, Padé analytic continuation, Wynn-epsilon sequence
//! acceleration, Borel-Padé resummation, and D-finite (holonomic) ODE analytic
//! continuation with monodromy matrices.

use num_complex::Complex64;

use crate::kernels::complex::circle_spectrum;
use crate::kernels::complex_roots::aberth_ehrlich;
use crate::kernels::complex_scalar::Complex;
use crate::kernels::contour_deform::ContourSegment;
use crate::kernels::contour_deform::integrate_segment_checked;
use crate::kernels::reuse::global_pool;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);

/// Extracts Laurent coefficients `c_k` for `k in -neg_order..=pos_order` around
/// `center` on the circle `|z - center| = radius` using a single FFT.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "sample count bounds all indices and floating normalization does not panic"
)]
pub fn laurent_coeffs_fft<F>(
    f: F,
    center: Complex64,
    radius: f64,
    neg_order: usize,
    pos_order: usize,
    n_samples: usize,
) -> Vec<(i32, Complex64)>
where
    F: Fn(Complex64) -> Complex64,
{
    let Some(min_n) = neg_order
        .checked_add(pos_order)
        .and_then(|order| order.checked_add(2))
    else {
        return Vec::new();
    };
    let n = n_samples.max(min_n).max(8);
    if !center.is_finite()
        || !radius.is_finite()
        || radius <= 0.0
        || i32::try_from(n).is_err()
        || crate::cancellation::is_cancelled()
    {
        return Vec::new();
    }
    let sample = |point| {
        if crate::cancellation::is_cancelled() {
            Complex64::new(f64::NAN, f64::NAN)
        } else {
            f(point)
        }
    };
    let buffer = circle_spectrum(sample, center, radius, n);
    if buffer.iter().any(|value| !value.is_finite()) {
        return Vec::new();
    }

    let nf = n as f64;
    let mut out = Vec::with_capacity(neg_order + pos_order + 1);
    for k in -(neg_order as i32)..=(pos_order as i32) {
        let idx = k.rem_euclid(n as i32) as usize;
        let ck = buffer.get(idx).copied().unwrap_or_default();
        let coeff = (ck / nf) / radius.powi(k);
        out.push((k, coeff));
    }
    out
}

/// Solves a dense complex linear system `A x = b` via Gaussian elimination with partial pivoting.
fn solve_complex_linear(
    mut a: Vec<Vec<Complex64>>,
    mut b: Vec<Complex64>,
) -> Option<Vec<Complex64>> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|r| r.len() != n) {
        return None;
    }
    for col in 0..n {
        let mut piv = col;
        let mut best = 0.0_f64;
        for (r, row) in a.iter().enumerate().skip(col) {
            let v = row.get(col)?.norm();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-30 || !best.is_finite() {
            return None;
        }
        if piv != col {
            a.swap(col, piv);
            b.swap(col, piv);
        }
        let diag = *a.get(col)?.get(col)?;
        for r in (col + 1)..n {
            let factor = *a.get(r)?.get(col)? / diag;
            if factor == ZERO {
                continue;
            }
            for j in col..n {
                let piv_val = *a.get(col)?.get(j)?;
                if let Some(cell) = a.get_mut(r).and_then(|row| row.get_mut(j)) {
                    *cell -= factor * piv_val;
                }
            }
            let b_col = *b.get(col)?;
            if let Some(slot) = b.get_mut(r) {
                *slot -= factor * b_col;
            }
        }
    }
    let mut x = vec![ZERO; n];
    for i in (0..n).rev() {
        let mut s = *b.get(i)?;
        for j in (i + 1)..n {
            s -= *a.get(i)?.get(j)? * *x.get(j)?;
        }
        let diag = *a.get(i)?.get(i)?;
        if let Some(slot) = x.get_mut(i) {
            *slot = s / diag;
        }
    }
    Some(x)
}

/// Constructs the `[m/n]` complex Padé approximant `P(z) / Q(z)` from Taylor
/// coefficients `c_0, c_1, ..., c_{m+n}` (with `Q(0) = 1`).
#[must_use]
pub fn pade_approximant_complex(
    taylor_coeffs: &[Complex64],
    m: usize,
    n: usize,
) -> Option<(Vec<Complex64>, Vec<Complex64>)> {
    if taylor_coeffs.len() <= m.checked_add(n)? {
        return None;
    }
    let mut q = vec![ONE];
    if n > 0 {
        let mut mat = vec![vec![ZERO; n]; n];
        let mut rhs = vec![ZERO; n];
        for i in 0..n {
            let row_idx = m + 1 + i;
            if let Some(slot) = rhs.get_mut(i) {
                *slot = -*taylor_coeffs.get(row_idx)?;
            }
            for j in 0..n {
                let c_idx = row_idx.checked_sub(j + 1);
                if let Some(idx) = c_idx
                    && let Some(cell) = mat.get_mut(i).and_then(|r| r.get_mut(j))
                {
                    *cell = *taylor_coeffs.get(idx)?;
                }
            }
        }
        let Some(b_sol) = solve_complex_linear(mat, rhs) else {
            return pade_approximant_complex(taylor_coeffs, m, n - 1);
        };
        q.extend(b_sol);
    }
    let mut p = vec![ZERO; m + 1];
    for (k, p_slot) in p.iter_mut().enumerate() {
        let mut s = ZERO;
        for j in 0..=k.min(n) {
            s += *q.get(j)? * *taylor_coeffs.get(k - j)?;
        }
        *p_slot = s;
    }
    Some((p, q))
}

/// Evaluates rational function `P(z) / Q(z)` via Horner's method.
#[must_use]
pub fn eval_pade_complex(
    p: &[Complex64],
    q: &[Complex64],
    z: Complex64,
) -> Complex64 {
    let horner = |coeffs: &[Complex64]| -> Complex64 {
        let mut acc = ZERO;
        for &c in coeffs.iter().rev() {
            acc = acc * z + c;
        }
        acc
    };
    let numerator = horner(p);
    let denominator = horner(q);
    let ratio =
        Complex::new(numerator.re, numerator.im).div(Complex::new(denominator.re, denominator.im));
    Complex64::new(ratio.re, ratio.im)
}

/// Wynn's epsilon algorithm for accelerating convergence of complex sequences.
///
/// Returns the last finite even-column estimate, or the last input if no
/// estimate is available. An empty sequence returns zero.
#[must_use]
pub fn wynn_epsilon_complex(partial_sums: &[Complex64]) -> Complex64 {
    let Some(&last) = partial_sums.last() else {
        return ZERO;
    };
    if partial_sums.len() < 2 {
        return last;
    }
    let mut best = last;
    let mut work = crate::kernels::reuse::fft_scratch_for(partial_sums.len().saturating_mul(2));
    work.resize(partial_sums.len(), ZERO);
    work.extend_from_slice(partial_sums);
    let (previous, current) = work.split_at_mut(partial_sums.len());
    super::convergence::wynn_epsilon_columns(
        previous,
        current,
        |delta| {
            let small = delta.re.abs() < 1e-30 && delta.im.abs() < 1e-30 && delta.norm() < 1e-30;
            (!small).then(|| ONE / delta)
        },
        |value| {
            if value.is_finite() {
                best = value;
            }
        },
    );
    best
}

// Divide complex values without squaring an extreme denominator first.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "nonfinite complex quotients are checked by the caller"
)]
fn complex_ratio(
    numerator: Complex64,
    denominator: Complex64,
) -> Complex64 {
    let scale = denominator.re.abs().max(denominator.im.abs());
    (numerator / scale) / (denominator / scale)
}

/// Borel–Padé resummation using a Laplace integral of the rational Borel transform.
///
/// A positive-axis pole selects the upper lateral integral. Returns `None`
/// when the Padé model or the integration cannot be resolved in binary64.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "coefficient orders and quadrature are bounded; complex failures propagate as nonfinite status"
)]
pub fn borel_pade_resum(
    coeffs: &[Complex64],
    z: Complex64,
    pade_m: usize,
    pade_n: usize,
) -> Option<Complex64> {
    let used = pade_m.checked_add(pade_n)?.checked_add(1)?;
    if used > 171 || !z.is_finite() || crate::cancellation::is_cancelled() {
        return None;
    }
    let coeffs = coeffs.get(..used)?;
    if coeffs.iter().any(|coefficient| !coefficient.is_finite()) {
        return None;
    }
    if z == ZERO {
        return coeffs.first().copied();
    }
    let mut borel = global_pool().scratch_complex64(used);
    let mut fact = 1.0;
    for (k, &coefficient) in coeffs.iter().enumerate() {
        if k > 0 {
            fact *= k as f64;
        }
        borel.push(coefficient / fact);
    }
    let (p, mut q) = pade_approximant_complex(&borel, pade_m, pade_n)?;
    while q.len() > 1 && q.last() == Some(&ZERO) {
        q.pop();
    }
    let Some(p_degree) = p.iter().rposition(|&coefficient| coefficient != ZERO) else {
        return Some(ZERO);
    };
    let p = p.get(..=p_degree)?;
    let degree_delta = p_degree as f64 + 1.0 - q.len() as f64;
    let direction = borel_direction(&q, z)?;
    let evaluate = |u: Complex64| {
        let remaining = 1.0 - u.re;
        let t = direction * (u.re / remaining);
        let argument = z * t;
        let (value, exponent) = if argument.norm() > 1.0 {
            // Reverse the polynomials at large arguments and combine their net
            // power with the Laplace exponent before either can overflow.
            let inverse = complex_ratio(ONE, argument);
            let reverse = |polynomial: &[Complex64]| {
                polynomial
                    .iter()
                    .fold(ZERO, |sum, coefficient| sum * inverse + coefficient)
            };
            (
                complex_ratio(reverse(p), reverse(&q)),
                degree_delta * argument.ln() - t,
            )
        } else {
            (eval_pade_complex(p, &q, argument), -t)
        };
        direction * exponent.exp() * value / (remaining * remaining)
    };
    let segment = ContourSegment::Line {
        start: ZERO,
        end: ONE,
    };
    let result = integrate_segment_checked(evaluate, segment, 1e-11, 4096).ok()?;
    (result.converged && result.value.is_finite()).then_some(result.value)
}

#[allow(
    clippy::arithmetic_side_effects,
    clippy::single_call_fn,
    reason = "select a pole-free lateral direction before the finite integration"
)]
fn borel_direction(
    denominator: &[Complex64],
    z: Complex64,
) -> Option<Complex64> {
    if denominator.len() < 2 {
        return Some(ONE);
    }
    let poles = aberth_ehrlich(denominator, 256, 1e-12).ok()?;
    let mut axial = false;
    let mut angle = 0.2_f64;
    for root in poles {
        let pole = complex_ratio(root, z);
        if !pole.is_finite() {
            return None;
        }
        if pole.re > 0.0 {
            let pole_angle = pole.arg();
            if pole_angle.abs() <= 1e-10 {
                axial = true;
            } else if pole_angle > 0.0 {
                angle = angle.min(0.5 * pole_angle);
            }
        }
    }
    Some(if axial {
        Complex64::from_polar(1.0, angle)
    } else {
        ONE
    })
}

/// Continues a Taylor series along a chain of overlapping disks.
///
/// The series `sum c_k (z - a_0)^k` is re-centered at each successive point in `centers`,
/// returning the coefficients at the last center.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "complex coordinate differences use the supplied floating-point scale"
)]
pub fn weierstrass_chain_continue(
    initial_coeffs: &[Complex64],
    centers: &[Complex64],
) -> Option<Vec<Complex64>> {
    if initial_coeffs.is_empty() || centers.is_empty() {
        return None;
    }
    let mut coeffs = initial_coeffs.to_vec();
    if centers.len() == 1 {
        return Some(coeffs);
    }
    let n = coeffs.len();
    let mut next = vec![ZERO; n];
    for window in centers.windows(2) {
        let &[previous, center] = window else {
            return None;
        };
        translate_into(&coeffs, center - previous, &mut next);
        std::mem::swap(&mut coeffs, &mut next);
    }
    Some(coeffs)
}

/// Linear homogeneous ODE `sum_{j=0}^r p_j(z) y^(j)(z) = 0` with polynomial
/// coefficients `coeffs[j]` in ascending powers of `z`.
#[derive(Clone, Debug)]
pub struct HolonomicOde {
    /// Polynomial coefficients `[p_0(z), p_1(z), ..., p_r(z)]`.
    pub coeffs: Vec<Vec<Complex64>>,
}

impl HolonomicOde {
    /// Order `r` of the differential equation.
    #[must_use]
    pub const fn order(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }
}

/// A step evaluates all continued solution columns against one shifted polynomial
/// table. Falling factorial factors do not depend on the path or basis column.
struct TaylorPlan {
    order: usize,
    degree: usize,
    cutoff: usize,
    falling: Vec<Vec<f64>>,
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "a bounded Taylor table and finite floating-point recurrence handle every accepted step"
)]
impl TaylorPlan {
    #[allow(
        clippy::single_call_fn,
        reason = "bounds and precomputes the Taylor recurrence before any path is stepped"
    )]
    fn new(
        ode: &HolonomicOde,
        order: usize,
    ) -> Option<Self> {
        let degree = ode.order();
        let max_power = ode
            .coeffs
            .iter()
            .filter_map(|polynomial| {
                polynomial
                    .iter()
                    .rposition(|&coefficient| coefficient != ZERO)
            })
            .max()
            .unwrap_or(0);
        // A polynomial factor can delay the next nonzero Taylor term by
        // its degree plus derivative order. Cover that recurrence width
        // both in the expansion and in the tail comparison.
        let recurrence_width = max_power.checked_add(degree)?.checked_add(1)?;
        let order = order.max(recurrence_width.checked_mul(2)?.checked_add(6)?);
        let cutoff = order.saturating_sub(recurrence_width.max(4));
        if degree == 0
            || degree > 64
            || order > 512
            || ode.coeffs.iter().flatten().any(|value| !value.is_finite())
            || ode.coeffs.last()?.iter().all(|&value| value == ZERO)
        {
            return None;
        }
        let falling = (0..=order)
            .map(|n| {
                let mut product = 1.0;
                (0..=degree)
                    .map(|derivative| {
                        if derivative > 0 {
                            product *= n.saturating_sub(derivative - 1) as f64;
                        }
                        product
                    })
                    .collect()
            })
            .collect();
        Some(Self {
            order,
            degree,
            cutoff,
            falling,
        })
    }

    fn expand(
        &self,
        shifted: &[Vec<Complex64>],
        state: &[Complex64],
        series: &mut [Complex64],
    ) -> Option<()> {
        let lead = *shifted.last()?.first()?;
        if lead == ZERO || !lead.is_finite() {
            return None;
        }
        let lead_scale = lead.re.abs().max(lead.im.abs());
        let normalized_lead = lead / lead_scale;
        let inverse = normalized_lead.conj() / normalized_lead.norm_sqr();
        for (j, (slot, &initial)) in series.iter_mut().zip(state).enumerate() {
            *slot = initial / self.falling.get(j)?.get(j)?;
        }
        for n in 0..=(self.order - self.degree) {
            let mut rhs = ZERO;
            for (j, polynomial) in shifted.iter().enumerate() {
                let start = usize::from(j == self.degree);
                for (m, &coefficient) in polynomial.iter().enumerate().take(n + 1).skip(start) {
                    let k = n - m + j;
                    let falling = *self.falling.get(k)?.get(j)?;
                    rhs += coefficient * falling * *series.get(k)?;
                }
            }
            let k = n + self.degree;
            let factor = *self.falling.get(k)?.get(self.degree)?;
            let value = -(rhs / lead_scale) * inverse / factor;
            if !value.is_finite() {
                return None;
            }
            *series.get_mut(k)? = value;
        }
        Some(())
    }

    fn evaluate(
        &self,
        series: &[Complex64],
        step: Complex64,
        powers: &[f64],
        state: &mut [Complex64],
    ) -> Option<bool> {
        let mut accurate = true;
        let max_component = |z: Complex64| z.re.abs().max(z.im.abs());
        let cutoff = self.cutoff + 1;
        let (low, high) = series.split_at_checked(cutoff)?;
        let reference = series.iter().take(state.len()).enumerate().try_fold(
            f64::MIN_POSITIVE,
            |largest, (j, value)| {
                Some(largest.max(max_component(*value) * self.falling.get(j)?.get(j)?))
            },
        )?;
        // Accumulate the tail from its highest degree, then continue the
        // same Horner chain through the lower polynomial.
        for (j, slot) in state.iter_mut().enumerate() {
            let mut full = ZERO;
            for (k, &coefficient) in high.iter().enumerate().rev() {
                full = full * step + coefficient * *self.falling.get(cutoff + k)?.get(j)?;
            }
            let tail = max_component(full) * powers.get(cutoff - j)?;
            for (k, &coefficient) in low.iter().enumerate().skip(j).rev() {
                full = full * step + coefficient * *self.falling.get(k)?.get(j)?;
            }
            if !full.is_finite() || !tail.is_finite() {
                return Some(false);
            }
            let tolerance = 1.4e-13 * max_component(full).max(reference);
            accurate &= tail <= tolerance;
            *slot = full;
        }
        Some(accurate)
    }
}

#[allow(
    clippy::single_call_fn,
    reason = "compute only the final leading polynomial when validating the endpoint"
)]
fn translate_polynomial(
    polynomial: &[Complex64],
    center: Complex64,
) -> Vec<Complex64> {
    let mut shifted = vec![ZERO; polynomial.len()];
    translate_into(polynomial, center, &mut shifted);
    shifted
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "the fixed coefficient slice bounds the reverse Horner shift"
)]
fn translate_into(
    polynomial: &[Complex64],
    center: Complex64,
    shifted: &mut [Complex64],
) {
    shifted.fill(ZERO);
    for (degree, &coefficient) in polynomial.iter().rev().enumerate() {
        for j in (1..=degree).rev() {
            if let (Some(&previous), Some(&value)) = (shifted.get(j - 1), shifted.get(j))
                && let Some(slot) = shifted.get_mut(j)
            {
                *slot = value * center + previous;
            }
        }
        if let Some(first) = shifted.first_mut() {
            *first = *first * center + coefficient;
        }
    }
}

fn singularity_bound(leading: &[Complex64]) -> Option<f64> {
    let (&constant, tail) = leading.split_first()?;
    let base = constant.norm();
    if base == 0.0 || !base.is_finite() {
        return None;
    }
    let max_tail = tail
        .iter()
        .map(|coefficient| coefficient.norm())
        .fold(0.0_f64, f64::max);
    let bound = if max_tail == 0.0 {
        f64::INFINITY
    } else if tail.iter().skip(1).all(|&coefficient| coefficient == ZERO) {
        base / max_tail
    } else {
        // Reverse-polynomial Cauchy bound avoids computing any roots.
        1.0 / (1.0 + max_tail / base)
    };
    (bound > 0.0).then_some(bound)
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "finite Taylor steps are bounded by a shared path budget and checked endpoints"
)]
fn continue_columns(
    ode: &HolonomicOde,
    path: &[Complex64],
    states: &mut [Vec<Complex64>],
    taylor_order: usize,
) -> Option<()> {
    let plan = TaylorPlan::new(ode, taylor_order)?;
    if path.is_empty()
        || path.iter().any(|point| !point.is_finite())
        || states
            .iter()
            .any(|state| state.len() != plan.degree || state.iter().any(|value| !value.is_finite()))
    {
        return None;
    }
    let mut series = vec![vec![ZERO; plan.order + 1]; states.len()];
    let mut next_states = vec![vec![ZERO; plan.degree]; states.len()];
    let mut shifted: Vec<_> = ode
        .coeffs
        .iter()
        .map(|polynomial| vec![ZERO; polynomial.len()])
        .collect();
    let mut powers = vec![1.0; plan.cutoff + 2];
    for window in path.windows(2) {
        let &[mut center, target] = window else {
            return None;
        };
        for _ in 0..10_000 {
            if crate::cancellation::is_cancelled() {
                return None;
            }
            if center == target {
                break;
            }
            let remaining = target - center;
            let distance = remaining.norm();
            if !distance.is_finite() || distance <= 0.0 {
                return None;
            }
            for (polynomial, output) in ode.coeffs.iter().zip(&mut shifted) {
                translate_into(polynomial, center, output);
            }
            let limit = 0.45 * singularity_bound(shifted.last()?)?;
            let mut fraction = (limit / distance).min(1.0);
            for (state, out) in states.iter().zip(&mut series) {
                plan.expand(&shifted, state, out)?;
            }
            let mut accepted = None;
            for _ in 0..60 {
                let step = remaining * fraction;
                let next_center = if fraction < 1.0 {
                    center + step
                } else {
                    target
                };
                if next_center == center || !next_center.is_finite() {
                    return None;
                }
                let shift = next_center - center;
                let step_norm = shift.norm();
                let mut power = 1.0;
                for value in &mut powers {
                    *value = power;
                    power *= step_norm;
                }
                let mut accurate = true;
                for (coefficients, out) in series.iter().zip(&mut next_states) {
                    accurate &= plan.evaluate(coefficients, shift, &powers, out)?;
                }
                if accurate {
                    accepted = Some(next_center);
                    break;
                }
                fraction *= 0.5;
            }
            center = accepted?;
            for (slot, evaluated) in states.iter_mut().zip(&mut next_states) {
                std::mem::swap(slot, evaluated);
            }
        }
        if center != target {
            return None;
        }
    }
    // Single-point paths and the final endpoint also need an ordinary center.
    singularity_bound(&translate_polynomial(ode.coeffs.last()?, *path.last()?))?;
    Some(())
}

/// Continues the derivative state `[y, y', ..., y^(r-1)]` along a polygonal path.
/// Returns `None` for an invalid ODE, a singular path, cancellation or failed
/// Taylor refinement.
#[must_use]
pub fn continue_holonomic_along_path(
    ode: &HolonomicOde,
    path: &[Complex64],
    initial_state: &[Complex64],
    taylor_order: usize,
) -> Option<Vec<Complex64>> {
    let mut state = [initial_state.to_vec()];
    continue_columns(ode, path, &mut state, taylor_order)?;
    state.into_iter().next()
}

/// Computes the monodromy matrix by continuing the identity basis around the
/// supplied loop, closing it at its first point.
#[must_use]
pub fn holonomic_monodromy_matrix(
    ode: &HolonomicOde,
    loop_pts: &[Complex64],
    taylor_order: usize,
) -> Option<Vec<Vec<Complex64>>> {
    let r = ode.order();
    if r == 0 || r > 64 {
        return None;
    }
    let &first = loop_pts.first()?;
    let mut closed = loop_pts.to_vec();
    if closed.last() != Some(&first) {
        closed.push(first);
    }
    let mut identity = vec![vec![ZERO; r]; r];
    for (i, column) in identity.iter_mut().enumerate() {
        *column.get_mut(i)? = ONE;
    }
    continue_columns(ode, &closed, &mut identity, taylor_order)?;
    let mut matrix = vec![vec![ZERO; r]; r];
    for (row, out) in matrix.iter_mut().enumerate() {
        for (slot, column) in out.iter_mut().zip(&identity) {
            *slot = *column.get(row)?;
        }
    }
    Some(matrix)
}
