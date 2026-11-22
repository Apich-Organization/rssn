//! # Numerical Polynomial Operations
//!
//! This module provides a `Polynomial` struct and associated functions for numerical
//! operations on polynomials with `f64` coefficients. It supports evaluation,
//! batch evaluation, interpolation, composition, exponentiation, numerical GCD,
//! differentiation, arithmetic (addition, subtraction, fast Karatsuba/FFT multiplication,
//! division), and finding real roots.

use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Sub;

use num_complex::Complex64;
use rustfft::FftDirection;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::dd::Dd;
use crate::kernels::real_roots;
use crate::kernels::reuse;

/// Minimum shorter-operand length (`min(deg A, deg B) + 1`) for Karatsuba multiplication.
///
/// Measured crossover on this machine (16 cores, `--release`, coefficient
/// magnitudes around 25, interleaved A/B, minimum of nine trials): at a
/// shorter length of 36 Karatsuba takes 1363 ns against the schoolbook loop's
/// 1300 ns, and from 40 upwards it is consistently 7-14% ahead (40: 1519 vs
/// 1751 ns; 44: 1689 vs 1969 ns; 64: 4208 vs 4681 ns). The switch sits at 40.
/// The integer coefficient paths are exact, so a different summation order
/// cannot change their result; the `f64` path keeps the same association as
/// the schoolbook loop because Karatsuba only re-associates *additions of
/// partial products*, which is why the differential test compares against the
/// schoolbook oracle with a tolerance instead of bit equality.
const KARATSUBA_MIN_LEN: usize = 40;
/// Minimum shorter-operand length (`min(deg A, deg B) + 1`) for FFT multiplication.
///
/// Measured crossover on this machine, same harness: the packed FFT and
/// Karatsuba are within noise of each other at 44 and 45 (1.06x and 1.00x of
/// Karatsuba's time) and the FFT then wins every size above (48: 0.73x, 52:
/// 0.63x, 64: 0.45x, 96: 0.50x, 128: 0.31x). 48 is the first size where the
/// FFT is clearly ahead, so the switch sits there.
const FFT_MIN_LEN: usize = 48;
/// Degree threshold from which `compose` switches from Horner to Brent-Kung baby-step giant-step.
const BRENT_KUNG_MIN_DEGREE: usize = 16;

/// Represents a polynomial with f64 coefficients for numerical operations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Polynomial {
    /// Coefficients of the polynomial, from highest degree to lowest.
    pub coeffs: Vec<f64>,
}

/// Why checked polynomial division could not produce finite coefficients.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolynomialDivisionError {
    /// The divisor is empty or every coefficient is zero.
    ZeroDivisor,
    /// An input coefficient is nonfinite.
    NonFiniteInput,
    /// A quotient or remainder coefficient exceeded the available range.
    NonFiniteResult,
}

impl std::fmt::Display for PolynomialDivisionError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        f.write_str(match self {
            | Self::ZeroDivisor => "polynomial divisor must be nonzero",
            | Self::NonFiniteInput => "polynomial division requires finite coefficients",
            | Self::NonFiniteResult => "polynomial division produced a nonfinite coefficient",
        })
    }
}

impl std::error::Error for PolynomialDivisionError {}

#[inline]
fn mul_schoolbook_add(
    a: &[f64],
    b: &[f64],
    out: &mut [f64],
) {
    for (i, &c1) in a.iter().enumerate() {
        let dst = &mut out[i..i + b.len()];
        let mut b_chunks = b.chunks_exact(4);
        let mut dst_chunks = dst.chunks_exact_mut(4);
        for (dc, bc) in dst_chunks.by_ref().zip(b_chunks.by_ref()) {
            if let ([d0, d1, d2, d3], &[b0, b1, b2, b3]) = (dc, bc) {
                *d0 += c1 * b0;
                *d1 += c1 * b1;
                *d2 += c1 * b2;
                *d3 += c1 * b3;
            }
        }
        for (d, &c2) in dst_chunks
            .into_remainder()
            .iter_mut()
            .zip(b_chunks.remainder())
        {
            *d += c1 * c2;
        }
    }
}

fn karatsuba_scratch_len(max_len: usize) -> usize {
    max_len.saturating_mul(8).saturating_add(128)
}

fn mul_karatsuba_rec(
    a: &[f64],
    b: &[f64],
    out: &mut [f64],
    scratch: &mut [f64],
) {
    let na = a.len();
    let nb = b.len();
    if na == 0 || nb == 0 {
        return;
    }
    if na.min(nb) < KARATSUBA_MIN_LEN {
        if na <= nb {
            mul_schoolbook_add(a, b, out);
        } else {
            mul_schoolbook_add(b, a, out);
        }
        return;
    }

    // Ensure `a` is the shorter or equal slice.
    if na > nb {
        mul_karatsuba_rec(b, a, out, scratch);
        return;
    }

    let m = nb / 2;
    if na <= m {
        // Unbalanced split: divide `b` into `b0` and `b1`.
        let (b0, b1) = b.split_at(m);
        mul_karatsuba_rec(a, b0, &mut out[..na + m - 1], scratch);
        mul_karatsuba_rec(a, b1, &mut out[m..na + nb - 1], scratch);
        return;
    }

    let (a0, a1) = a.split_at(m);
    let (b0, b1) = b.split_at(m);
    let a1_len = a1.len();
    let b1_len = b1.len();
    let sum_a_len = m.max(a1_len);
    let sum_b_len = m.max(b1_len);
    let z0_len = 2 * m - 1;
    let z2_len = a1_len + b1_len - 1;
    let z1_len = sum_a_len + sum_b_len - 1;
    let needed = sum_a_len + sum_b_len + z0_len + z2_len + z1_len;

    if scratch.len() < needed {
        mul_schoolbook_add(a, b, out);
        return;
    }

    let (sum_a, rest) = scratch.split_at_mut(sum_a_len);
    let (sum_b, rest) = rest.split_at_mut(sum_b_len);
    let (z0, rest) = rest.split_at_mut(z0_len);
    let (z2, rest) = rest.split_at_mut(z2_len);
    let (z1, child_scratch) = rest.split_at_mut(z1_len);

    z0.fill(0.0);
    z2.fill(0.0);
    z1.fill(0.0);

    mul_karatsuba_rec(a0, b0, z0, child_scratch);
    mul_karatsuba_rec(a1, b1, z2, child_scratch);

    sum_a[..m].copy_from_slice(a0);
    for slot in &mut sum_a[m..] {
        *slot = 0.0;
    }
    for (slot, &val) in sum_a.iter_mut().zip(a1) {
        *slot += val;
    }

    sum_b[..m].copy_from_slice(b0);
    for slot in &mut sum_b[m..] {
        *slot = 0.0;
    }
    for (slot, &val) in sum_b.iter_mut().zip(b1) {
        *slot += val;
    }

    mul_karatsuba_rec(sum_a, sum_b, z1, child_scratch);

    for (slot, &val) in z1.iter_mut().zip(z0.iter()) {
        *slot -= val;
    }
    for (slot, &val) in z1.iter_mut().zip(z2.iter()) {
        *slot -= val;
    }

    for (slot, &val) in out[..z0_len].iter_mut().zip(z0.iter()) {
        *slot += val;
    }
    for (slot, &val) in out[m..m + z1_len].iter_mut().zip(z1.iter()) {
        *slot += val;
    }
    for (slot, &val) in out[2 * m..2 * m + z2_len].iter_mut().zip(z2.iter()) {
        *slot += val;
    }
}

/// Evaluates `coeffs` at every point of `xs`, writing the results into `out`.
///
/// `out` must have the same length as `xs`. Each result is bit-identical to
/// [`Polynomial::eval`]: the accumulator for a point follows the same Horner
/// chain, and only four points are advanced at once so that four independent
/// multiply-add chains overlap.
fn eval_into(
    coeffs: &[f64],
    xs: &[f64],
    out: &mut [f64],
) {
    let mut x_chunks = xs.chunks_exact(4);
    let mut out_chunks = out.chunks_exact_mut(4);
    for (oc, xc) in out_chunks.by_ref().zip(x_chunks.by_ref()) {
        if let ([o0, o1, o2, o3], &[x0, x1, x2, x3]) = (oc, xc) {
            let mut acc0 = 0.0_f64;
            let mut acc1 = 0.0_f64;
            let mut acc2 = 0.0_f64;
            let mut acc3 = 0.0_f64;
            for &c in coeffs {
                acc0 = acc0.mul_add(x0, c);
                acc1 = acc1.mul_add(x1, c);
                acc2 = acc2.mul_add(x2, c);
                acc3 = acc3.mul_add(x3, c);
            }
            *o0 = acc0;
            *o1 = acc1;
            *o2 = acc2;
            *o3 = acc3;
        }
    }
    for (o, &x) in out_chunks
        .into_remainder()
        .iter_mut()
        .zip(x_chunks.remainder())
    {
        *o = coeffs.iter().fold(0.0, |acc, &c| acc.mul_add(x, c));
    }
}

fn trim_leading_zeros(mut coeffs: Vec<f64>) -> Vec<f64> {
    let first_nonzero = coeffs.iter().position(|&c| c != 0.0);
    match first_nonzero {
        | Some(0) => coeffs,
        | Some(idx) => {
            coeffs.drain(..idx);
            coeffs
        },
        | None => {
            if coeffs.is_empty() {
                coeffs
            } else {
                vec![0.0]
            }
        },
    }
}

fn trim_leading_tol(
    coeffs: &[f64],
    tol: f64,
) -> Vec<f64> {
    let max_abs = coeffs.iter().fold(0.0_f64, |m, &c| m.max(c.abs()));
    if max_abs == 0.0 || max_abs <= tol {
        return vec![0.0];
    }
    let cutoff = tol * max_abs.max(1.0);
    let first = coeffs.iter().position(|&c| c.abs() > cutoff);
    match first {
        | Some(idx) => coeffs[idx..].to_vec(),
        | None => vec![0.0],
    }
}

impl Polynomial {
    /// Creates a new `Polynomial` from coefficients (highest degree first).
    #[must_use]
    pub const fn new(coeffs: Vec<f64>) -> Self {
        Self { coeffs }
    }

    /// Evaluates the polynomial at a given point `x` using Horner's method.
    ///
    /// Horner's method is an efficient algorithm for evaluating polynomials.
    /// For a polynomial `P(x) = c_n*x^n + c_{n-1}*x^{n-1} + ... + c_1*x + c_0`,
    /// it computes `P(x) = (...((c_n*x + c_{n-1})*x + c_{n-2})*x + ... + c_1)*x + c_0`.
    ///
    /// # Arguments
    /// * `x` - The point at which to evaluate the polynomial.
    ///
    /// # Returns
    /// The value of the polynomial at `x` as an `f64`.
    #[must_use]
    pub fn eval(
        &self,
        x: f64,
    ) -> f64 {
        self.coeffs.iter().fold(0.0, |acc, &c| acc.mul_add(x, c))
    }

    /// Evaluates the polynomial at a slice of points `xs` using 4-way unrolled
    /// Horner evaluation and deterministic multi-core parallelism when the
    /// workload exceeds the execution threshold.
    ///
    /// The points are split into disjoint contiguous chunks; every chunk runs
    /// the same unrolled kernel, so each output element is produced by exactly
    /// one call and the result is bit-identical to [`Self::eval`] for every
    /// thread count.
    #[must_use]
    pub fn eval_batch(
        &self,
        xs: &[f64],
    ) -> Vec<f64> {
        if xs.is_empty() {
            return Vec::new();
        }
        if self.coeffs.is_empty() {
            return vec![0.0; xs.len()];
        }
        let mut out = vec![0.0; xs.len()];
        // `should_split_work` sees the polynomial degree, which a plain chunk
        // count cannot; when it says no, one chunk keeps the call sequential.
        let chunk = if crate::execution::should_split_work(xs.len(), self.coeffs.len()) {
            crate::execution::choose_chunk_size(xs.len(), 8)
        } else {
            xs.len()
        };
        crate::kernels::par::for_chunks_mut(&mut out, chunk, 2, |index, block| {
            let start = index.saturating_mul(chunk);
            eval_into(&self.coeffs, &xs[start..start + block.len()], block);
        });
        out
    }

    /// Multiplies two polynomials using the 4-way unrolled $O(n m)$ schoolbook algorithm.
    #[must_use]
    pub fn mul_schoolbook(
        &self,
        rhs: &Self,
    ) -> Self {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return Self { coeffs: Vec::new() };
        }
        let mut new_coeffs = vec![0.0; self.coeffs.len() + rhs.coeffs.len() - 1];
        if self.coeffs.len() <= rhs.coeffs.len() {
            mul_schoolbook_add(&self.coeffs, &rhs.coeffs, &mut new_coeffs);
        } else {
            mul_schoolbook_add(&rhs.coeffs, &self.coeffs, &mut new_coeffs);
        }
        Self { coeffs: new_coeffs }
    }

    /// Multiplies two polynomials using recursive Karatsuba multiplication
    /// backed by a reusable scratch buffer from [`crate::kernels::reuse::f64_scratch`].
    #[must_use]
    pub fn mul_karatsuba(
        &self,
        rhs: &Self,
    ) -> Self {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return Self { coeffs: Vec::new() };
        }
        let na = self.coeffs.len();
        let nb = rhs.coeffs.len();
        if na.min(nb) < KARATSUBA_MIN_LEN {
            return self.mul_schoolbook(rhs);
        }
        let mut out = vec![0.0; na + nb - 1];
        let scratch_len = karatsuba_scratch_len(na.max(nb));
        let mut scratch = reuse::f64_scratch(scratch_len);
        scratch.resize(scratch_len, 0.0);
        mul_karatsuba_rec(&self.coeffs, &rhs.coeffs, &mut out, &mut scratch);
        Self { coeffs: out }
    }

    /// Multiplies two polynomials using a packed 2-for-1 real-input complex FFT
    /// backed by [`crate::kernels::reuse::fft_plan`] and [`crate::kernels::reuse::fft_scratch_for`].
    ///
    /// Falls back to [`Self::mul_karatsuba`] if any input coefficient is non-finite.
    ///
    /// # Accuracy
    ///
    /// The FFT accumulates a relative error of a few `f64` epsilons times
    /// `log n`, where the schoolbook loop stays at one or two epsilons and
    /// Karatsuba at a handful. That is a deliberate trade for the 1.4–3x
    /// speedup above the threshold; a caller that needs the tighter
    /// association should call [`Self::mul_karatsuba`] directly, and
    /// [`Self::compose`] sidesteps the question by accumulating in
    /// double-double. `pow` inherits whatever tier its operands select, which
    /// is why `p.pow(3)` and `p * p * p` agree exactly: both go through this
    /// same dispatch.
    #[must_use]
    pub fn mul_fft(
        &self,
        rhs: &Self,
    ) -> Self {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return Self { coeffs: Vec::new() };
        }
        if self
            .coeffs
            .iter()
            .chain(&rhs.coeffs)
            .any(|c| !c.is_finite())
        {
            return self.mul_karatsuba(rhs);
        }
        let out_len = self.coeffs.len() + rhs.coeffs.len() - 1;
        let fft_len = out_len.next_power_of_two().max(2);
        let fwd = reuse::fft_plan(fft_len, FftDirection::Forward);
        let inv = reuse::fft_plan(fft_len, FftDirection::Inverse);
        let inplace_scratch_len = fwd
            .get_inplace_scratch_len()
            .max(inv.get_inplace_scratch_len());

        let mut data = reuse::fft_scratch_for(fft_len);
        let mut freq = reuse::fft_scratch_for(fft_len);
        let mut fft_work = reuse::fft_scratch_for(inplace_scratch_len);

        data.resize(fft_len, Complex64::new(0.0, 0.0));
        for (i, &a) in self.coeffs.iter().enumerate() {
            data[i].re = a;
        }
        for (i, &b) in rhs.coeffs.iter().enumerate() {
            data[i].im = b;
        }

        fft_work.resize(fwd.get_inplace_scratch_len(), Complex64::new(0.0, 0.0));
        fwd.process_with_scratch(&mut data, &mut fft_work);

        // Extract A_k * B_k = -0.25 * i * (Z_k^2 - conj(Z_{-k mod N})^2).
        freq.resize(fft_len, Complex64::new(0.0, 0.0));
        let minus_quarter_i = Complex64::new(0.0, -0.25);
        for k in 0..fft_len {
            let mirror = if k == 0 { 0 } else { fft_len - k };
            let zk = data[k];
            let z_mirror_conj = data[mirror].conj();
            freq[k] = (zk * zk - z_mirror_conj * z_mirror_conj) * minus_quarter_i;
        }

        fft_work.resize(inv.get_inplace_scratch_len(), Complex64::new(0.0, 0.0));
        inv.process_with_scratch(&mut freq, &mut fft_work);

        let scale = 1.0 / (fft_len as f64);
        let coeffs: Vec<f64> = freq[..out_len].iter().map(|c| c.re * scale).collect();
        Self { coeffs }
    }

    /// Builds the monic polynomial $\prod_{i} (x - r_i)$ from its roots using a
    /// balanced divide-and-conquer subproduct tree.
    #[must_use]
    pub fn from_roots(roots: &[f64]) -> Self {
        fn build_tree(slice: &[f64]) -> Polynomial {
            match slice {
                | [] => Polynomial::new(vec![1.0]),
                | &[r] => Polynomial::new(vec![1.0, -r]),
                | &[r0, r1] => Polynomial::new(vec![1.0, -(r0 + r1), r0 * r1]),
                | _ => {
                    let mid = slice.len() / 2;
                    let left = build_tree(&slice[..mid]);
                    let right = build_tree(&slice[mid..]);
                    left * right
                },
            }
        }
        build_tree(roots)
    }

    /// Computes the unique interpolating polynomial of degree `< xs.len()` passing
    /// through `(xs[i], ys[i])` using Newton divided differences converted to the
    /// descending power basis.
    ///
    /// # Errors
    /// Returns an error if `xs` and `ys` have different lengths, are empty,
    /// contain non-finite values, or contain duplicate `x` nodes.
    pub fn interpolate_lagrange(
        xs: &[f64],
        ys: &[f64],
    ) -> Result<Self, String> {
        if xs.is_empty() {
            return Err("interpolation requires at least one point".to_string());
        }
        if xs.len() != ys.len() {
            return Err("xs and ys must have the same length".to_string());
        }
        if xs.iter().chain(ys.iter()).any(|v| !v.is_finite()) {
            return Err("interpolation points must be finite".to_string());
        }
        let n = xs.len();
        let mut dd = ys.to_vec();
        for j in 1..n {
            for i in (j..n).rev() {
                let dx = xs[i] - xs[i - j];
                if dx == 0.0 {
                    return Err("interpolation nodes must be distinct".to_string());
                }
                dd[i] = (dd[i] - dd[i - 1]) / dx;
                if !dd[i].is_finite() {
                    return Err("interpolation divided difference overflowed".to_string());
                }
            }
        }

        // Convert Newton form: P(x) = (...((dd[n-1]*(x - xs[n-2]) + dd[n-2])*(x - xs[n-3]) + ...)+ dd[0])
        // Maintain descending-degree coefficients.
        let mut coeffs = vec![dd[n - 1]];
        for i in (0..n - 1).rev() {
            let xi = xs[i];
            let mut next = vec![0.0; coeffs.len() + 1];
            for (k, &c) in coeffs.iter().enumerate() {
                next[k] += c;
                next[k + 1] -= c * xi;
            }
            if let Some(constant) = next.last_mut() {
                *constant += dd[i];
            }
            coeffs = next;
        }
        Ok(Self {
            coeffs: trim_leading_zeros(coeffs),
        })
    }

    /// Evaluates the Lagrange interpolating polynomial through `(xs[i], ys[i])`
    /// directly at `x` using the second (true) barycentric interpolation formula.
    ///
    /// # Errors
    /// Returns an error if `xs` and `ys` have different lengths, are empty,
    /// contain non-finite values, or contain duplicate `x` nodes.
    pub fn interpolate_barycentric_eval(
        xs: &[f64],
        ys: &[f64],
        x: f64,
    ) -> Result<f64, String> {
        if xs.is_empty() {
            return Err("barycentric evaluation requires at least one point".to_string());
        }
        if xs.len() != ys.len() {
            return Err("xs and ys must have the same length".to_string());
        }
        if !x.is_finite() || xs.iter().chain(ys.iter()).any(|v| !v.is_finite()) {
            return Err("barycentric evaluation requires finite inputs".to_string());
        }
        let n = xs.len();
        if n == 1 {
            return Ok(ys[0]);
        }

        // Normalize node spacing to avoid overflow/underflow in weight products.
        let min_x = xs.iter().copied().fold(f64::INFINITY, f64::min);
        let max_x = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let span = (max_x - min_x).max(1.0);
        let inv_scale = 4.0 / span;

        let mut weights = vec![1.0_f64; n];
        for i in 0..n {
            let xi = xs[i];
            let mut prod = 1.0_f64;
            for (j, &xj) in xs.iter().enumerate() {
                if i == j {
                    continue;
                }
                let diff = (xi - xj) * inv_scale;
                if diff == 0.0 {
                    return Err("barycentric interpolation nodes must be distinct".to_string());
                }
                prod *= diff;
            }
            if prod == 0.0 || !prod.is_finite() {
                return Err("barycentric weight overflowed or underflowed".to_string());
            }
            weights[i] = 1.0 / prod;
        }

        let mut numer = 0.0_f64;
        let mut denom = 0.0_f64;
        for i in 0..n {
            let diff = x - xs[i];
            if diff.abs() <= f64::EPSILON * xs[i].abs().max(1.0) {
                return Ok(ys[i]);
            }
            let term = weights[i] / diff;
            numer = term.mul_add(ys[i], numer);
            denom += term;
        }
        Ok(numer / denom)
    }

    /// Computes the polynomial composition $P(\text{inner}(x))$ where `self` is $P$.
    ///
    /// The arithmetic runs in double-double and the result is rounded to `f64`
    /// once, at the end. Composition is ill-conditioned whenever the
    /// intermediate terms $c_k Q(x)^k$ are much larger than the result: a `f64`
    /// Horner and a `f64` Brent-Kung composition then disagree in the last
    /// digits, not because either is wrong but because they round at different
    /// points. Carrying the same two algorithms over [`Dd`] drops the
    /// accumulated error to about `1e-32` relative, so both round to the same
    /// `f64` and the answer is the correctly-rounded one rather than whichever
    /// strategy happened to run.
    #[must_use]
    pub fn compose(
        &self,
        inner: &Self,
    ) -> Self {
        if self.coeffs.is_empty() {
            return Self { coeffs: Vec::new() };
        }
        if self.coeffs.len() == 1 {
            return self.clone();
        }
        if inner.coeffs.is_empty() {
            return Self {
                coeffs: vec![*self.coeffs.last().unwrap_or(&0.0)],
            };
        }

        // Coefficients are stored highest-degree first; the composition reads
        // better constant-term first, so both operands are reversed on the way
        // in and the result is reversed on the way out.
        let p_asc: Vec<f64> = self.coeffs.iter().rev().copied().collect();
        let q: Vec<Dd> = inner
            .coeffs
            .iter()
            .rev()
            .map(|&c| Dd::from_f64(c))
            .collect();
        let acc = if p_asc.len() - 1 < BRENT_KUNG_MIN_DEGREE {
            compose_horner_dd(&p_asc, &q)
        } else {
            compose_brent_kung_dd(&p_asc, &q)
        };
        let coeffs: Vec<f64> = trim_leading_zeros(acc.iter().rev().map(|c| c.to_f64()).collect());
        Self { coeffs }
    }

    /// Raises the polynomial to the integer power `exp` by binary exponentiation.
    #[must_use]
    pub fn pow(
        &self,
        mut exp: usize,
    ) -> Self {
        if exp == 0 {
            return Self::new(vec![1.0]);
        }
        if self.coeffs.is_empty() {
            return Self::new(Vec::new());
        }
        if exp == 1 {
            return self.clone();
        }
        let mut result = Self::new(vec![1.0]);
        let mut base = self.clone();
        while exp > 0 {
            if exp & 1 == 1 {
                result = &result * &base;
            }
            exp >>= 1;
            if exp > 0 {
                base = &base * &base;
            }
        }
        result
    }

    /// Computes the numerical monic greatest common divisor of `self` and `other`
    /// using the Euclidean algorithm with tolerance `tol` for leading-term and
    /// remainder cleanup.
    #[must_use]
    pub fn gcd(
        &self,
        other: &Self,
        tol: f64,
    ) -> Self {
        let tol = tol.abs().max(f64::EPSILON * 16.0);
        let mut a = trim_leading_tol(&self.coeffs, tol);
        let mut b = trim_leading_tol(&other.coeffs, tol);

        let is_effectively_zero = |v: &[f64], t: f64| -> bool { v.iter().all(|&c| c.abs() <= t) };

        if is_effectively_zero(&a, tol) && is_effectively_zero(&b, tol) {
            return Self::new(vec![0.0]);
        }
        if is_effectively_zero(&a, tol) {
            a = b;
            b = vec![0.0];
        }

        let scale_init = a
            .iter()
            .chain(b.iter())
            .fold(0.0_f64, |m, &c| m.max(c.abs()))
            .max(1.0);

        while !is_effectively_zero(&b, tol * scale_init) {
            // Normalize `b` by its leading coefficient for numerical stability.
            let lead_b = b[0];
            if lead_b.abs() <= tol * scale_init {
                break;
            }
            for c in &mut b {
                *c /= lead_b;
            }
            if a.len() < b.len() {
                std::mem::swap(&mut a, &mut b);
                continue;
            }
            let norm_a = a.iter().fold(0.0_f64, |m, &c| m.max(c.abs())).max(1.0);
            // Only the remainder is needed, and `b` stays live as the next
            // dividend, so the elimination runs in place on `a` without
            // cloning the divisor or materializing the quotient.
            let (_, rem) = Self::new(std::mem::take(&mut a)).divide_by::<false>(&b);
            let mut r = rem.coeffs;
            for c in &mut r {
                if c.abs() <= tol * norm_a {
                    *c = 0.0;
                }
            }
            let r_trimmed = trim_leading_tol(&r, tol);
            a = b;
            b = if is_effectively_zero(&r_trimmed, tol * norm_a) {
                vec![0.0]
            } else {
                r_trimmed
            };
        }

        let a = trim_leading_tol(&a, tol);
        let lead = a.first().copied().unwrap_or(0.0);
        if lead.abs() <= tol {
            return Self::new(vec![0.0]);
        }
        let monic_coeffs: Vec<f64> = a
            .into_iter()
            .map(|c| {
                let v = c / lead;
                if v.abs() <= tol { 0.0 } else { v }
            })
            .collect();
        Self::new(monic_coeffs)
    }

    /// Finds the real roots of the polynomial.
    ///
    /// This method combines Sturm's theorem for root isolation with Newton's method
    /// for refining the roots. Sturm's theorem provides disjoint intervals, each
    /// containing exactly one real root, which are then used as starting points
    /// for Newton's method to converge to the root.
    ///
    /// # Returns
    /// A `Result` containing a `Vec<f64>` of the real roots found, or an error string
    /// if root isolation or refinement fails.
    ///
    /// # Errors
    /// Returns an error if root isolation or refinement fails.
    pub fn find_roots(&self) -> Result<Vec<f64>, String> {
        let derivative = self.derivative();

        let isolating_intervals = real_roots::isolate_real_roots(self, 1e-9)?;

        let mut roots = Vec::new();

        for (a, b) in isolating_intervals {
            let mut guess = f64::midpoint(a, b);

            for _ in 0..30 {
                let f_val = self.eval(guess);

                let f_prime_val = derivative.eval(guess);

                if f_prime_val.abs() < 1e-12 {
                    break;
                }

                let next_guess = guess - f_val / f_prime_val;

                if (next_guess - guess).abs() < 1e-12 {
                    guess = next_guess;

                    break;
                }

                guess = next_guess;
            }

            roots.push(guess);
        }

        Ok(roots)
    }

    /// Returns the derivative of the polynomial.
    ///
    /// The derivative is computed by applying the power rule to each term.
    /// For a term `c*x^n`, its derivative is `(c*n)*x^(n-1)`.
    ///
    /// # Returns
    /// A new `Polynomial` representing the derivative.
    #[must_use]
    pub fn derivative(&self) -> Self {
        if self.coeffs.len() <= 1 {
            return Self { coeffs: vec![0.0] };
        }

        let mut new_coeffs = Vec::with_capacity(self.coeffs.len() - 1);

        let n = (self.coeffs.len() - 1) as f64;

        for (i, &c) in self.coeffs.iter().enumerate().take(self.coeffs.len() - 1) {
            new_coeffs.push(c * (n - i as f64));
        }

        Self { coeffs: new_coeffs }
    }

    /// Performs polynomial long division.
    ///
    /// The divisor must have a nonzero leading coefficient. Use
    /// [`Self::try_long_division`] to validate the inputs and results.
    ///
    /// # Returns
    /// A tuple `(quotient, remainder)` as `Polynomial`s.
    ///
    /// # Panics
    /// Panics if the divisor has no coefficients.
    #[must_use]
    pub fn long_division(
        self,
        divisor: &Self,
    ) -> (Self, Self) {
        self.divide_by::<false>(&divisor.coeffs)
    }

    /// Divides finite coefficients by a nonzero descending-coefficient polynomial.
    ///
    /// Leading zeros in the divisor are skipped. A fused update can retain a
    /// finite coefficient when a product overflows before it is subtracted.
    ///
    /// # Errors
    /// Returns an error for zero divisors, nonfinite inputs or nonfinite results.
    pub fn try_long_division(
        self,
        divisor: &Self,
    ) -> Result<(Self, Self), PolynomialDivisionError> {
        if self
            .coeffs
            .iter()
            .chain(&divisor.coeffs)
            .any(|x| !x.is_finite())
        {
            return Err(PolynomialDivisionError::NonFiniteInput);
        }
        let start = divisor
            .coeffs
            .iter()
            .position(|&x| x != 0.0)
            .ok_or(PolynomialDivisionError::ZeroDivisor)?;
        let (q, r) = self.divide_by::<true>(&divisor.coeffs[start..]);
        if q.coeffs.iter().chain(&r.coeffs).any(|x| !x.is_finite()) {
            return Err(PolynomialDivisionError::NonFiniteResult);
        }
        Ok((q, r))
    }

    fn divide_by<const CHECKED: bool>(
        mut self,
        divisor: &[f64],
    ) -> (Self, Self) {
        let (&lead, tail) = divisor.split_first().expect("polynomial divisor is empty");
        let eliminated = self
            .coeffs
            .len()
            .checked_sub(divisor.len())
            .map_or(0, |d| d + 1);
        let mut quotient = Vec::with_capacity(eliminated.max(1));
        let safe_factor = if CHECKED {
            (f64::MAX * 0.5) / tail.iter().fold(0.0_f64, |max, term| max.max(term.abs()))
        } else {
            f64::INFINITY
        };
        // Keep the original allocation during elimination. Only the trailing
        // remainder needs to move when the calculation is complete.
        for offset in 0..eliminated {
            let (leading, remaining) = self.coeffs[offset..]
                .split_first_mut()
                .expect("each quotient term eliminates one coefficient");
            let factor = *leading / lead;
            quotient.push(factor);
            if CHECKED && factor.abs() > safe_factor {
                for (value, &term) in remaining.iter_mut().zip(tail) {
                    let product = term * factor;
                    *value = if product.is_finite() {
                        *value - product
                    } else {
                        (-factor).mul_add(term, *value)
                    };
                }
            } else {
                for (value, &term) in remaining.iter_mut().zip(tail) {
                    *value -= term * factor;
                }
            }
        }
        if quotient.is_empty() {
            quotient.push(0.0);
        }

        let remainder = &self.coeffs[eliminated..];
        let scale = remainder.iter().fold(0.0_f64, |max, c| max.max(c.abs()));
        let trim = remainder
            .iter()
            .take(remainder.len().saturating_sub(1))
            .take_while(|c| (!CHECKED || c.is_finite()) && c.abs() <= f64::EPSILON * 16.0 * scale)
            .count();

        if eliminated + trim != 0 {
            self.coeffs.drain(..eliminated + trim);
        }

        (Self { coeffs: quotient }, self)
    }

    /// Returns the degree of the polynomial.
    #[must_use]
    pub const fn degree(&self) -> usize {
        if self.coeffs.is_empty() {
            return 0;
        }

        self.coeffs.len() - 1
    }

    /// Checks if the polynomial is zero (within epsilon).
    #[must_use]
    pub fn is_zero(
        &self,
        epsilon: f64,
    ) -> bool {
        self.coeffs.iter().all(|&c| c.abs() < epsilon)
    }

    /// Returns the indefinite integral of the polynomial (constant of integration = 0).
    #[must_use]
    pub fn integral(&self) -> Self {
        if self.coeffs.is_empty() {
            return Self { coeffs: vec![0.0] };
        }

        let mut new_coeffs = Vec::with_capacity(self.coeffs.len() + 1);

        let d = self.degree() as f64;

        for (i, &c) in self.coeffs.iter().enumerate() {
            new_coeffs.push(c / (d - i as f64 + 1.0));
        }

        new_coeffs.push(0.0);

        Self { coeffs: new_coeffs }
    }
}

impl Add for Polynomial {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {
        // Coefficients are stored highest degree first, so the shorter operand
        // is aligned against the tail. Keeping the operand order of the sum
        // (`self + rhs`) preserves the exact floating-point result, including
        // signed zeros and NaN payloads.
        if self.coeffs.len() >= rhs.coeffs.len() {
            let mut long = self.coeffs;
            let offset = long.len() - rhs.coeffs.len();
            for (slot, &value) in long[offset..].iter_mut().zip(&rhs.coeffs) {
                *slot += value;
            }
            Self { coeffs: long }
        } else {
            let mut long = rhs.coeffs;
            let offset = long.len() - self.coeffs.len();
            for (slot, &value) in long[offset..].iter_mut().zip(&self.coeffs) {
                *slot = value + *slot;
            }
            Self { coeffs: long }
        }
    }
}

impl Sub for Polynomial {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        // Same tail alignment as `Add`, preserving `self - rhs` coefficient by
        // coefficient.
        if self.coeffs.len() >= rhs.coeffs.len() {
            let mut long = self.coeffs;
            let offset = long.len() - rhs.coeffs.len();
            for (slot, &value) in long[offset..].iter_mut().zip(&rhs.coeffs) {
                *slot -= value;
            }
            Self { coeffs: long }
        } else {
            let mut long = rhs.coeffs;
            let offset = long.len() - self.coeffs.len();
            for (slot, &value) in long[offset..].iter_mut().zip(&self.coeffs) {
                *slot = value - *slot;
            }
            Self { coeffs: long }
        }
    }
}

impl Mul for Polynomial {
    type Output = Self;

    fn mul(
        self,
        rhs: Self,
    ) -> Self {
        (&self).mul(&rhs)
    }
}

impl Mul<&Polynomial> for &Polynomial {
    type Output = Polynomial;

    fn mul(
        self,
        rhs: &Polynomial,
    ) -> Polynomial {
        if self.coeffs.is_empty() || rhs.coeffs.is_empty() {
            return Polynomial { coeffs: Vec::new() };
        }
        let min_len = self.coeffs.len().min(rhs.coeffs.len());
        if min_len < KARATSUBA_MIN_LEN {
            self.mul_schoolbook(rhs)
        } else if min_len < FFT_MIN_LEN {
            self.mul_karatsuba(rhs)
        } else {
            self.mul_fft(rhs)
        }
    }
}

/// `a * b` over double-double coefficients, schoolbook order.
///
/// Schoolbook is the right choice here: the composition needs the product to be
/// accurate, not fast, and a Karatsuba split would re-associate the sums for no
/// measurable gain at these degrees.
#[expect(
    clippy::single_call_fn,
    reason = "one multiplication kernel shared by both composition strategies"
)]
fn dd_poly_mul(
    a: &[Dd],
    b: &[Dd],
) -> Vec<Dd> {
    if a.is_empty() || b.is_empty() {
        return vec![Dd::ZERO];
    }
    let mut out = vec![Dd::ZERO; a.len() + b.len() - 1];
    for (i, &ai) in a.iter().enumerate() {
        if ai.hi == 0.0 && ai.lo == 0.0 {
            continue;
        }
        for (j, &bj) in b.iter().enumerate() {
            if let Some(slot) = out.get_mut(i + j) {
                *slot += ai * bj;
            }
        }
    }
    dd_trim(out)
}

/// `acc += factor * addend`, extending `acc` to the addend's length.
fn dd_poly_add_scaled(
    acc: &mut Vec<Dd>,
    addend: &[Dd],
    factor: Dd,
) {
    if acc.len() < addend.len() {
        acc.resize(addend.len(), Dd::ZERO);
    }
    for (slot, &term) in acc.iter_mut().zip(addend) {
        *slot += factor * term;
    }
}

/// Drops trailing exact zeros, keeping at least one coefficient.
fn dd_trim(mut coeffs: Vec<Dd>) -> Vec<Dd> {
    while coeffs.len() > 1 && coeffs.last().is_some_and(|c| c.hi == 0.0 && c.lo == 0.0) {
        coeffs.pop();
    }
    if coeffs.is_empty() {
        coeffs.push(Dd::ZERO);
    }
    coeffs
}

/// Horner composition `P(Q)`, both given constant-term first.
#[expect(
    clippy::single_call_fn,
    reason = "one of the two composition strategies"
)]
fn compose_horner_dd(
    p_asc: &[f64],
    q: &[Dd],
) -> Vec<Dd> {
    let Some((&lead, rest)) = p_asc.split_last() else {
        return vec![Dd::ZERO];
    };
    let mut acc = vec![Dd::from_f64(lead)];
    for &c in rest.iter().rev() {
        acc = dd_poly_mul(&acc, q);
        dd_poly_add_scaled(&mut acc, &[Dd::from_f64(c)], Dd::ONE);
    }
    acc
}

/// Brent-Kung baby-step giant-step composition `P(Q)`, constant-term first.
///
/// Writes `P(y) = sum_j P_j(y) * (y^m)^j` with `m = ceil(sqrt(deg P + 1))`, so
/// the number of polynomial multiplications falls from `deg P` to about
/// `2 * sqrt(deg P)`.
#[expect(
    clippy::single_call_fn,
    reason = "one of the two composition strategies"
)]
fn compose_brent_kung_dd(
    p_asc: &[f64],
    q: &[Dd],
) -> Vec<Dd> {
    let n = p_asc.len();
    let m = (n as f64).sqrt().ceil() as usize;
    let mut powers: Vec<Vec<Dd>> = Vec::with_capacity(m + 1);
    powers.push(vec![Dd::ONE]);
    for i in 1..=m {
        let previous = powers.get(i - 1).map_or(&[][..], Vec::as_slice);
        powers.push(dd_poly_mul(previous, q));
    }
    let giant_step = powers.pop().unwrap_or_else(|| vec![Dd::ZERO]);

    let num_chunks = n.div_ceil(m);
    let mut blocks: Vec<Vec<Dd>> = Vec::with_capacity(num_chunks);
    for chunk in 0..num_chunks {
        let low = chunk * m;
        let high = ((chunk + 1) * m).min(n);
        let mut block = vec![Dd::ZERO];
        for d in low..high {
            let Some(&c) = p_asc.get(d) else {
                continue;
            };
            if c == 0.0 {
                continue;
            }
            let Some(power) = powers.get(d - low) else {
                continue;
            };
            dd_poly_add_scaled(&mut block, power, Dd::from_f64(c));
        }
        blocks.push(dd_trim(block));
    }

    let mut acc = blocks.pop().unwrap_or_else(|| vec![Dd::ZERO]);
    while let Some(block) = blocks.pop() {
        let mut product = dd_poly_mul(&acc, &giant_step);
        dd_poly_add_scaled(&mut product, &block, Dd::ONE);
        acc = dd_trim(product);
    }
    acc
}

impl Div for Polynomial {
    type Output = Self;

    fn div(
        self,
        rhs: Self,
    ) -> Self {
        self.long_division(&rhs).0
    }
}

impl Mul<f64> for Polynomial {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {
        let new_coeffs = self.coeffs.iter().map(|&c| c * rhs).collect();

        Self { coeffs: new_coeffs }
    }
}

impl Div<f64> for Polynomial {
    type Output = Self;

    fn div(
        self,
        rhs: f64,
    ) -> Self {
        let new_coeffs = self.coeffs.iter().map(|&c| c / rhs).collect();

        Self { coeffs: new_coeffs }
    }
}

impl Polynomial {
    /// Divides the polynomial by a scalar.
    ///
    /// # Arguments
    /// * `rhs` - The scalar to divide by.
    ///
    /// # Returns
    /// A `Result` containing the new `Polynomial`, or an error if the scalar is zero.
    ///
    /// # Errors
    /// Returns an error if the divisor scalar is zero.
    pub fn div_scalar(
        self,
        rhs: f64,
    ) -> Result<Self, String> {
        if rhs == 0.0 {
            return Err("Division by \
                        zero scalar"
                .to_string());
        }

        let new_coeffs = self.coeffs.iter().map(|&c| c / rhs).collect();

        Ok(Self { coeffs: new_coeffs })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_poly(
        len: usize,
        seed: u64,
    ) -> Polynomial {
        let mut s = seed;
        let coeffs = (0..len)
            .map(|_| {
                s = s
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let v = ((s >> 33) as i64 % 201) - 100;
                (v as f64) * 0.25
            })
            .collect();
        Polynomial::new(coeffs)
    }

    #[test]
    fn schoolbook_karatsuba_and_fft_agree_across_tier_boundaries() {
        for (na, nb) in [(17, 29), (33, 33), (64, 97), (129, 131), (256, 300)] {
            let a = sample_poly(na, (na as u64) * 17 + 1);
            let b = sample_poly(nb, (nb as u64) * 31 + 7);
            let sb = a.mul_schoolbook(&b);
            let ka = a.mul_karatsuba(&b);
            let fft = a.mul_fft(&b);
            let auto = a.clone() * b.clone();

            assert_eq!(sb.coeffs.len(), na + nb - 1);
            assert_eq!(ka.coeffs.len(), sb.coeffs.len());
            assert_eq!(fft.coeffs.len(), sb.coeffs.len());
            assert_eq!(auto.coeffs.len(), sb.coeffs.len());

            for i in 0..sb.coeffs.len() {
                let ref_v = sb.coeffs[i];
                let scale = ref_v.abs().max(1.0);
                assert!(
                    (ka.coeffs[i] - ref_v).abs() <= 1e-11 * scale,
                    "Karatsuba mismatch at {i} for ({na}, {nb}): got {}, expected {ref_v}",
                    ka.coeffs[i]
                );
                assert!(
                    (fft.coeffs[i] - ref_v).abs() <= 1e-10 * scale,
                    "FFT mismatch at {i} for ({na}, {nb}): got {}, expected {ref_v}",
                    fft.coeffs[i]
                );
                assert!(
                    (auto.coeffs[i] - ref_v).abs() <= 1e-10 * scale,
                    "Auto mul mismatch at {i} for ({na}, {nb}): got {}, expected {ref_v}",
                    auto.coeffs[i]
                );
            }
        }
    }

    #[test]
    fn mul_fft_falls_back_cleanly_on_non_finite_inputs() {
        let mut a = sample_poly(140, 42);
        let b = sample_poly(140, 99);
        a.coeffs[10] = f64::NAN;
        let res = a.clone() * b.clone();
        let expected = a.mul_karatsuba(&b);
        assert_eq!(res.coeffs.len(), expected.coeffs.len());
        assert!(res.coeffs[10].is_nan());
    }

    #[test]
    fn eval_batch_is_bit_identical_to_single_eval_across_thread_pools() {
        let p = sample_poly(120, 777);
        let xs: Vec<f64> = (0..513).map(|i| -2.0 + 4.0 * (i as f64) / 512.0).collect();
        let expected: Vec<f64> = xs.iter().map(|&x| p.eval(x)).collect();

        let serial = crate::execution::with_scoped_pool(Some(1), || p.eval_batch(&xs));
        let parallel = crate::execution::with_scoped_pool(Some(4), || p.eval_batch(&xs));

        assert_eq!(serial.len(), expected.len());
        for i in 0..xs.len() {
            assert_eq!(serial[i].to_bits(), expected[i].to_bits());
            assert_eq!(parallel[i].to_bits(), expected[i].to_bits());
        }
    }

    #[test]
    fn from_roots_and_lagrange_and_barycentric_interpolation() {
        let roots = [-3.0, -1.5, 0.5, 2.0, 4.25];
        let p = Polynomial::from_roots(&roots);
        assert_eq!(p.degree(), roots.len());
        for &r in &roots {
            assert!(p.eval(r).abs() < 1e-12, "p({r}) = {}", p.eval(r));
        }

        let xs = [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
        let ys: Vec<f64> = xs.iter().map(|&x| p.eval(x)).collect();
        let interp = Polynomial::interpolate_lagrange(&xs, &ys).unwrap();
        assert_eq!(interp.coeffs.len(), p.coeffs.len());
        for (c1, c2) in interp.coeffs.iter().zip(&p.coeffs) {
            assert!((c1 - c2).abs() < 1e-10, "{c1} vs {c2}");
        }

        for test_x in [-2.0, -0.75, 0.5, 1.25, 3.0, 3.5] {
            let bary = Polynomial::interpolate_barycentric_eval(&xs, &ys, test_x).unwrap();
            let exact = p.eval(test_x);
            assert!(
                (bary - exact).abs() < 1e-10,
                "bary({test_x}) = {bary}, exact = {exact}"
            );
        }
    }

    #[test]
    fn auto_mul_matches_schoolbook_across_retuned_thresholds() {
        // Sizes straddling KARATSUBA_MIN_LEN (40) and FFT_MIN_LEN (48), plus
        // the power-of-two boundaries where the FFT padding jumps.
        let cases = [
            (39, 39),
            (40, 40),
            (44, 44),
            (47, 47),
            (48, 48),
            (49, 49),
            (64, 64),
            (65, 65),
            (128, 128),
            (129, 129),
            (200, 41),
            (48, 400),
        ];
        for (na, nb) in cases {
            let a = sample_poly(na, (na as u64) * 13 + 3);
            let b = sample_poly(nb, (nb as u64) * 29 + 5);
            let oracle = a.mul_schoolbook(&b);
            let auto = a.clone() * b.clone();
            assert_eq!(auto.coeffs.len(), na + nb - 1, "length for ({na}, {nb})");
            for i in 0..oracle.coeffs.len() {
                let reference = oracle.coeffs[i];
                let scale = reference.abs().max(1.0);
                assert!(
                    (auto.coeffs[i] - reference).abs() <= 1e-10 * scale,
                    "auto mul mismatch at {i} for ({na}, {nb}): got {}, expected {reference}",
                    auto.coeffs[i]
                );
            }
        }
    }

    #[test]
    fn eval_batch_is_bit_identical_across_thread_pools_at_parallel_sizes() {
        // 4096 points at degree 200 clears `should_split_work`, so the chunked
        // parallel path is exercised rather than the single-chunk fallback.
        let p = sample_poly(201, 424_242);
        let xs: Vec<f64> = (0..4096)
            .map(|i| -3.0 + 6.0 * (i as f64) / 4095.0)
            .collect();
        let expected: Vec<f64> = xs.iter().map(|&x| p.eval(x)).collect();

        let one = crate::execution::with_scoped_pool(Some(1), || p.eval_batch(&xs));
        let four = crate::execution::with_scoped_pool(Some(4), || p.eval_batch(&xs));

        assert_eq!(one.len(), expected.len());
        assert_eq!(four.len(), expected.len());
        for i in 0..xs.len() {
            assert_eq!(one[i].to_bits(), expected[i].to_bits(), "serial at {i}");
            assert_eq!(four[i].to_bits(), expected[i].to_bits(), "parallel at {i}");
        }

        // `mul` and `compose` take no part in the parallel split, so their
        // results must be bit-identical across thread counts; the product is
        // compared with the schoolbook oracle to the FFT tolerance because at
        // this length the auto path selects the FFT.
        let a = sample_poly(70, 11);
        let b = sample_poly(70, 12);
        let oracle = a.mul_schoolbook(&b);
        let prod1 = crate::execution::with_scoped_pool(Some(1), || a.clone() * b.clone());
        let prod4 = crate::execution::with_scoped_pool(Some(4), || a.clone() * b.clone());
        assert_eq!(prod1.coeffs.len(), oracle.coeffs.len());
        for i in 0..oracle.coeffs.len() {
            assert_eq!(
                prod1.coeffs[i].to_bits(),
                prod4.coeffs[i].to_bits(),
                "mul differs across thread counts at {i}"
            );
            let scale = oracle.coeffs[i].abs().max(1.0);
            assert!(
                (prod1.coeffs[i] - oracle.coeffs[i]).abs() <= 1e-10 * scale,
                "mul mismatch at {i}: got {}, expected {}",
                prod1.coeffs[i],
                oracle.coeffs[i]
            );
        }
    }

    #[test]
    fn compose_matches_pointwise_oracle_across_both_paths() {
        // `inner` is contractive on the sampled interval (|inner(x)| <= 0.653
        // for |x| <= 1.9), which keeps the composed coefficients bounded. With
        // a non-contractive inner polynomial the power-basis coefficients grow
        // like |inner|^deg (1e18 at deg 90) and the pointwise comparison
        // becomes meaningless for any implementation.
        let inner = Polynomial::new(vec![0.25, 0.0, -0.25]);
        // 15 is the last Horner degree, 16 the first Brent-Kung degree.
        for deg in [15, 16, 37, 90] {
            let outer = sample_poly(deg + 1, (deg as u64) * 5 + 11);
            let composed = outer.compose(&inner);
            // The leading coefficient of the composition is
            // `outer[0] * inner[0]^deg`, so the degree is exact whenever the
            // sample's leading coefficient is nonzero.
            if outer.coeffs[0] != 0.0 {
                assert_eq!(composed.degree(), deg * inner.degree());
            }
            for step in 0..40 {
                let x = -1.9 + 3.8 * (step as f64) / 39.0;
                let got = composed.eval(x);
                let expected = outer.eval(inner.eval(x));
                let scale = expected.abs().max(1.0);
                assert!(
                    (got - expected).abs() <= 1e-9 * scale,
                    "compose mismatch at deg={deg}, x={x}: got {got}, expected {expected}"
                );
            }
        }
    }

    /// The two composition strategies must round to the same `f64`.
    ///
    /// In `f64` they do not: Horner and Brent-Kung accumulate in different
    /// orders, so on an ill-conditioned composition (a non-contractive inner
    /// polynomial, where the intermediate terms dwarf the result) they differ
    /// in the last digits. In double-double the accumulated error is about
    /// `1e-32`, far below the rounding of the final result, so the two agree
    /// bit for bit — which is the property this test pins.
    #[test]
    fn compose_strategies_agree_bit_for_bit() {
        let inner = Polynomial::new(vec![3.0, -2.0, 1.0]);
        let q: Vec<Dd> = inner
            .coeffs
            .iter()
            .rev()
            .map(|&c| Dd::from_f64(c))
            .collect();
        for deg in [4_usize, 15, 16, 17, 40, 90] {
            let outer = sample_poly(deg + 1, (deg as u64) * 7 + 3);
            let p_asc: Vec<f64> = outer.coeffs.iter().rev().copied().collect();
            let horner = compose_horner_dd(&p_asc, &q);
            let brent_kung = compose_brent_kung_dd(&p_asc, &q);
            let round = |v: &Vec<Dd>| -> Vec<f64> {
                trim_leading_zeros(v.iter().rev().map(|c| c.to_f64()).collect())
            };
            assert_eq!(round(&horner), round(&brent_kung), "degree {deg}");
        }
    }

    /// The FFT tier must stay inside its stated error bound.
    ///
    /// This is a regression guard, not a proof: if a future change makes the
    /// packed transform lose more precision (a different packing, a missed
    /// normalisation), the relative deviation from the schoolbook oracle grows
    /// and this fails.
    #[test]
    fn fft_tier_stays_within_its_error_bound() {
        for len in [FFT_MIN_LEN, 64, 97, 128, 200] {
            let a = sample_poly(len, len as u64 * 3 + 1);
            let b = sample_poly(len, len as u64 * 5 + 2);
            let fast = a.mul_fft(&b);
            let oracle = a.mul_schoolbook(&b);
            assert_eq!(fast.coeffs.len(), oracle.coeffs.len(), "len {len}");
            let scale = oracle
                .coeffs
                .iter()
                .fold(0.0_f64, |m, &c| m.max(c.abs()))
                .max(1.0);
            let worst = fast
                .coeffs
                .iter()
                .zip(&oracle.coeffs)
                .fold(0.0_f64, |m, (&x, &y)| m.max((x - y).abs()));
            assert!(
                worst <= 1e-11 * scale,
                "len {len}: worst absolute deviation {worst} against scale {scale}"
            );
        }
    }

    #[test]
    fn compose_pow_and_numerical_gcd() {
        // Test compose across both Horner (< 16) and Brent-Kung (>= 16) degrees.
        let inner = Polynomial::new(vec![2.0, -1.0, 0.5]);
        for deg in [5, 20] {
            let outer = sample_poly(deg + 1, deg as u64 + 123);
            let comp = outer.compose(&inner);
            for x in [-1.5, -0.5, 0.25, 1.0, 2.0] {
                let got = comp.eval(x);
                let expected = outer.eval(inner.eval(x));
                let scale = expected.abs().max(1.0);
                assert!(
                    (got - expected).abs() <= 1e-9 * scale,
                    "compose mismatch at deg={deg}, x={x}: got {got}, expected {expected}"
                );
            }
        }

        // Test pow
        let base = Polynomial::new(vec![1.0, -2.0, 3.0]);
        let p5 = base.pow(5);
        for x in [-1.0, 0.5, 2.0] {
            assert!((p5.eval(x) - base.eval(x).powi(5)).abs() < 1e-10);
        }

        // Test numerical GCD: gcd((x - 1)(x + 2)(x - 3), (x - 1)(x - 3)(x + 4)) = (x - 1)(x - 3) = x^2 - 4x + 3
        let g_true = Polynomial::from_roots(&[1.0, 3.0]);
        let p1 = g_true.clone() * Polynomial::from_roots(&[-2.0, 0.5]);
        let p2 = g_true.clone() * Polynomial::from_roots(&[-4.0, 2.5]);
        let g = p1.gcd(&p2, 1e-9);
        assert_eq!(g.coeffs.len(), g_true.coeffs.len());
        for (c1, c2) in g.coeffs.iter().zip(&g_true.coeffs) {
            assert!((c1 - c2).abs() < 1e-7, "gcd coeff {c1} vs {c2}");
        }
    }
}
