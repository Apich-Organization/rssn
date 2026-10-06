//! # Interpolation and approximation
//!
//! * cubic splines (natural, clamped, not-a-knot) with derivatives and
//!   integrals;
//! * monotone shape-preserving PCHIP (Fritsch-Carlson / Fritsch-Butland);
//! * B-splines: Cox-de Boor basis, evaluation by de Boor's algorithm,
//!   interpolation with clamped knots;
//! * Chebyshev approximation: coefficients, Clenshaw evaluation,
//!   differentiation, integration, adaptive "chebfun-style" construction
//!   and roots (colleague matrix);
//! * AAA rational approximation (real data, barycentric form) and Pade
//!   approximants from Taylor coefficients.
#![allow(
    clippy::manual_midpoint,
    clippy::missing_const_for_fn,
    clippy::struct_field_names,
    clippy::or_fun_call,
    clippy::manual_swap,
    clippy::if_not_else,
    clippy::unnecessary_sort_by,
    clippy::while_float,
    clippy::too_long_first_doc_paragraph,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::needless_pass_by_value,
    clippy::manual_map,
    clippy::unnecessary_map_or,
    clippy::suboptimal_flops,
    clippy::similar_names,
    clippy::unreadable_literal,
    clippy::excessive_precision,
    clippy::needless_range_loop,
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::cognitive_complexity,
    clippy::option_if_let_else,
    clippy::many_single_char_names,
    clippy::type_complexity,
    clippy::too_many_arguments,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::map_unwrap_or,
    clippy::cast_lossless
)]

use std::f64::consts::PI;

use crate::kernels::dense::{self, Mat};

/// Errors of the approximation routines.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApproxError {
    /// Inputs have inconsistent lengths, are too short or not increasing.
    Invalid,
    /// An underlying linear solve failed.
    Singular,
}

fn check_grid(x: &[f64], y: &[f64], min: usize) -> Result<(), ApproxError> {
    if x.len() != y.len() || x.len() < min || x.windows(2).any(|w| w[1] <= w[0]) {
        return Err(ApproxError::Invalid);
    }
    Ok(())
}

/// Boundary condition of a cubic spline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SplineBc {
    /// Zero second derivative at both ends.
    Natural,
    /// Prescribed first derivatives `(left, right)`.
    Clamped(f64, f64),
    /// Third derivative continuous across the second and last-but-one knots.
    NotAKnot,
}

/// Piecewise cubic spline `y = a + b t + c t^2 + d t^3` with `t = x - x_i`.
#[derive(Debug, Clone)]
pub struct CubicSpline {
    x: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

impl CubicSpline {
    /// Builds an interpolating spline through `(x, y)`.
    ///
    /// # Errors
    /// [`ApproxError::Invalid`] for fewer than 3 points (4 for not-a-knot)
    /// or non-increasing abscissae.
    pub fn new(x: &[f64], y: &[f64], bc: SplineBc) -> Result<Self, ApproxError> {
        let n = x.len();
        check_grid(x, y, if bc == SplineBc::NotAKnot { 4 } else { 3 })?;
        let h: Vec<f64> = x.windows(2).map(|w| w[1] - w[0]).collect();
        // Solve for second-derivative-like coefficients c_i (c = M/2).
        let mut m = Mat::zeros(n, n);
        let mut r = vec![0.0; n];
        for i in 1..n - 1 {
            m.set(i, i - 1, h[i - 1]);
            m.set(i, i, 2.0 * (h[i - 1] + h[i]));
            m.set(i, i + 1, h[i]);
            r[i] = 3.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
        }
        match bc {
            SplineBc::Natural => {
                m.set(0, 0, 1.0);
                m.set(n - 1, n - 1, 1.0);
            }
            SplineBc::Clamped(l, rr) => {
                m.set(0, 0, 2.0 * h[0]);
                m.set(0, 1, h[0]);
                r[0] = 3.0 * ((y[1] - y[0]) / h[0] - l);
                m.set(n - 1, n - 2, h[n - 2]);
                m.set(n - 1, n - 1, 2.0 * h[n - 2]);
                r[n - 1] = 3.0 * (rr - (y[n - 1] - y[n - 2]) / h[n - 2]);
            }
            SplineBc::NotAKnot => {
                m.set(0, 0, h[1]);
                m.set(0, 1, -(h[0] + h[1]));
                m.set(0, 2, h[0]);
                m.set(n - 1, n - 3, h[n - 2]);
                m.set(n - 1, n - 2, -(h[n - 3] + h[n - 2]));
                m.set(n - 1, n - 1, h[n - 3]);
            }
        }
        let c = dense::solve(&m, &r).map_err(|_| ApproxError::Singular)?;
        let mut b = vec![0.0; n];
        let mut d = vec![0.0; n];
        for i in 0..n - 1 {
            b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0;
            d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
        }
        // last knot: derivative / third-derivative continuation
        let hl = h[n - 2];
        b[n - 1] = b[n - 2] + 2.0 * c[n - 2] * hl + 3.0 * d[n - 2] * hl * hl;
        Ok(Self { x: x.to_vec(), a: y.to_vec(), b, c, d })
    }

    fn seg(&self, x: f64) -> usize {
        self.x.partition_point(|&v| v <= x).saturating_sub(1).min(self.x.len() - 2)
    }

    /// Spline value at `x` (extrapolates with the end cubics).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        let i = self.seg(x);
        let t = x - self.x[i];
        self.a[i] + t * (self.b[i] + t * (self.c[i] + t * self.d[i]))
    }

    /// Derivative of order 1, 2 or 3 at `x` (order 0 is the value).
    #[must_use]
    pub fn derivative(&self, x: f64, order: u8) -> f64 {
        let i = self.seg(x);
        let t = x - self.x[i];
        match order {
            0 => self.eval(x),
            1 => self.b[i] + t * (2.0 * self.c[i] + 3.0 * self.d[i] * t),
            2 => 2.0 * self.c[i] + 6.0 * self.d[i] * t,
            3 => 6.0 * self.d[i],
            _ => 0.0,
        }
    }

    /// Exact integral of the spline over `[lo, hi]`.
    #[must_use]
    pub fn integrate(&self, lo: f64, hi: f64) -> f64 {
        if lo > hi {
            return -self.integrate(hi, lo);
        }
        let anti = |i: usize, t: f64| {
            t * (self.a[i] + t * (self.b[i] / 2.0 + t * (self.c[i] / 3.0 + t * self.d[i] / 4.0)))
        };
        let (il, ih) = (self.seg(lo), self.seg(hi));
        if il == ih {
            return anti(il, hi - self.x[il]) - anti(il, lo - self.x[il]);
        }
        let mut s = anti(il, self.x[il + 1] - self.x[il]) - anti(il, lo - self.x[il]);
        for i in il + 1..ih {
            s += anti(i, self.x[i + 1] - self.x[i]);
        }
        s + anti(ih, hi - self.x[ih])
    }
}

/// Monotone piecewise cubic Hermite interpolant (PCHIP).
///
/// Preserves the monotonicity of the data: no overshoot between knots.
#[derive(Debug, Clone)]
pub struct Pchip {
    x: Vec<f64>,
    y: Vec<f64>,
    d: Vec<f64>,
}

impl Pchip {
    /// Builds the interpolant.
    ///
    /// # Errors
    /// [`ApproxError::Invalid`] for fewer than 2 points or non-increasing `x`.
    pub fn new(x: &[f64], y: &[f64]) -> Result<Self, ApproxError> {
        check_grid(x, y, 2)?;
        let n = x.len();
        let h: Vec<f64> = x.windows(2).map(|w| w[1] - w[0]).collect();
        let del: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / h[i]).collect();
        let mut d = vec![0.0; n];
        if n == 2 {
            d.fill(del[0]);
        } else {
            for i in 1..n - 1 {
                if del[i - 1] * del[i] > 0.0 {
                    let w1 = 2.0 * h[i] + h[i - 1];
                    let w2 = h[i] + 2.0 * h[i - 1];
                    d[i] = (w1 + w2) / (w1 / del[i - 1] + w2 / del[i]);
                }
            }
            let end = |h0: f64, h1: f64, d0: f64, d1: f64| {
                let mut v = ((2.0 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
                if v.signum() != d0.signum() {
                    v = 0.0;
                } else if d0.signum() != d1.signum() && v.abs() > 3.0 * d0.abs() {
                    v = 3.0 * d0;
                }
                v
            };
            d[0] = end(h[0], h[1], del[0], del[1]);
            d[n - 1] = end(h[n - 2], h[n - 3], del[n - 2], del[n - 3]);
        }
        Ok(Self { x: x.to_vec(), y: y.to_vec(), d })
    }

    /// Value at `x`.
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        let i = self.x.partition_point(|&v| v <= x).saturating_sub(1).min(self.x.len() - 2);
        let h = self.x[i + 1] - self.x[i];
        let s = (x - self.x[i]) / h;
        let (h00, h10, h01, h11) = (
            (1.0 + 2.0 * s) * (1.0 - s) * (1.0 - s),
            s * (1.0 - s) * (1.0 - s),
            s * s * (3.0 - 2.0 * s),
            s * s * (s - 1.0),
        );
        h00 * self.y[i] + h10 * h * self.d[i] + h01 * self.y[i + 1] + h11 * h * self.d[i + 1]
    }

    /// Derivative at `x`.
    #[must_use]
    pub fn derivative(&self, x: f64) -> f64 {
        let i = self.x.partition_point(|&v| v <= x).saturating_sub(1).min(self.x.len() - 2);
        let h = self.x[i + 1] - self.x[i];
        let s = (x - self.x[i]) / h;
        let (d00, d10, d01, d11) = (
            6.0 * s * (s - 1.0),
            (1.0 - s) * (1.0 - 3.0 * s),
            6.0 * s * (1.0 - s),
            s * (3.0 * s - 2.0),
        );
        (d00 * self.y[i] + d01 * self.y[i + 1]) / h + d10 * self.d[i] + d11 * self.d[i + 1]
    }
}

/// Value of the `i`-th B-spline basis function of the given `degree` on
/// `knots` at `x` (Cox-de Boor recursion; right-continuous, closed at the
/// last knot).
#[must_use]
pub fn bspline_basis(knots: &[f64], degree: usize, i: usize, x: f64) -> f64 {
    if degree == 0 {
        let last = knots.len() - 1;
        let inside = knots[i] <= x && x < knots[i + 1];
        let at_end = x == knots[last] && knots[i] < knots[i + 1] && knots[i + 1] == knots[last];
        return if inside || at_end { 1.0 } else { 0.0 };
    }
    let mut v = 0.0;
    let d1 = knots[i + degree] - knots[i];
    if d1 > 0.0 {
        v += (x - knots[i]) / d1 * bspline_basis(knots, degree - 1, i, x);
    }
    let d2 = knots[i + degree + 1] - knots[i + 1];
    if d2 > 0.0 {
        v += (knots[i + degree + 1] - x) / d2 * bspline_basis(knots, degree - 1, i + 1, x);
    }
    v
}

/// Evaluates a B-spline curve `sum c_i B_{i,p}(x)` by de Boor's algorithm.
///
/// Requires `knots.len() == coeffs.len() + degree + 1`; `x` is clamped to
/// the valid range.
#[must_use]
pub fn bspline_eval(knots: &[f64], coeffs: &[f64], degree: usize, x: f64) -> f64 {
    let n = coeffs.len();
    let lo = knots[degree];
    let hi = knots[n];
    let x = x.clamp(lo, hi);
    let mut k = degree;
    while k < n - 1 && knots[k + 1] <= x {
        k += 1;
    }
    let mut d: Vec<f64> = (0..=degree).map(|j| coeffs[j + k - degree]).collect();
    for r in 1..=degree {
        for j in (r..=degree).rev() {
            let i = j + k - degree;
            let den = knots[i + degree + 1 - r] - knots[i];
            let alpha = if den == 0.0 { 0.0 } else { (x - knots[i]) / den };
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
        }
    }
    d[degree]
}

/// Interpolates `(x, y)` with a B-spline of the given degree on a clamped
/// knot vector (knot averaging). Returns `(knots, coefficients)`.
///
/// # Errors
/// [`ApproxError::Invalid`] if `degree >= x.len()` or `x` not increasing,
/// [`ApproxError::Singular`] if the collocation matrix is singular.
pub fn bspline_interpolate(
    x: &[f64],
    y: &[f64],
    degree: usize,
) -> Result<(Vec<f64>, Vec<f64>), ApproxError> {
    check_grid(x, y, degree + 1)?;
    let n = x.len();
    let mut knots = vec![x[0]; degree + 1];
    for j in 1..n - degree {
        let s: f64 = x[j..j + degree].iter().sum();
        knots.push(s / degree as f64);
    }
    knots.extend(std::iter::repeat_n(x[n - 1], degree + 1));
    let mut a = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            a.set(i, j, bspline_basis(&knots, degree, j, x[i]));
        }
    }
    let c = dense::solve(&a, y).map_err(|_| ApproxError::Singular)?;
    Ok((knots, c))
}

/// Chebyshev points of the second kind (extrema)
/// `cos(pi k / n)` mapped to `[a, b]`, `k = 0..=n`.
#[must_use]
pub fn chebyshev_points(n: usize, a: f64, b: f64) -> Vec<f64> {
    (0..=n).map(|k| 0.5 * (a + b) + 0.5 * (b - a) * (PI * k as f64 / n as f64).cos()).collect()
}

/// A function on `[a, b]` represented by Chebyshev coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct ChebSeries {
    /// Coefficients `c_0 .. c_n` of `sum c_k T_k(t)`.
    pub coeffs: Vec<f64>,
    /// Left end of the interval.
    pub a: f64,
    /// Right end of the interval.
    pub b: f64,
}

/// Clenshaw evaluation of `sum c_k T_k(t)` at `t in [-1, 1]`.
#[must_use]
pub fn clenshaw(c: &[f64], t: f64) -> f64 {
    let (mut b1, mut b2) = (0.0, 0.0);
    for &ck in c.iter().skip(1).rev() {
        let b0 = ck + 2.0 * t * b1 - b2;
        b2 = b1;
        b1 = b0;
    }
    c.first().copied().unwrap_or(0.0) + t * b1 - b2
}

impl ChebSeries {
    /// Interpolates `f` at the `n + 1` Chebyshev extrema.
    pub fn from_fn<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, n: usize) -> Self {
        let pts = chebyshev_points(n, a, b);
        let v: Vec<f64> = pts.iter().map(|&x| f(x)).collect();
        let mut c = vec![0.0; n + 1];
        for j in 0..=n {
            let mut s = 0.0;
            for k in 0..=n {
                let w = if k == 0 || k == n { 0.5 } else { 1.0 };
                s += w * v[k] * (PI * (j * k) as f64 / n as f64).cos();
            }
            c[j] = 2.0 * s / n as f64;
        }
        c[0] *= 0.5;
        c[n] *= 0.5;
        Self { coeffs: c, a, b }
    }

    /// Adaptive construction: doubles the degree (up to `max_degree`)
    /// until the trailing coefficients fall below `tol * max|c|`, then
    /// chops the negligible tail.
    pub fn adaptive<F: Fn(f64) -> f64>(f: F, a: f64, b: f64, tol: f64, max_degree: usize) -> Self {
        let mut n = 8;
        loop {
            let mut s = Self::from_fn(&f, a, b, n);
            let big = s.coeffs.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
            let cut = tol * big.max(f64::MIN_POSITIVE);
            let tail_ok = s.coeffs.iter().rev().take(3.min(s.coeffs.len())).all(|v| v.abs() <= cut);
            if tail_ok || n >= max_degree {
                while s.coeffs.len() > 1 && s.coeffs.last().is_some_and(|v| v.abs() <= cut) {
                    s.coeffs.pop();
                }
                return s;
            }
            n *= 2;
        }
    }

    /// Degree of the series.
    #[must_use]
    pub fn degree(&self) -> usize {
        self.coeffs.len().saturating_sub(1)
    }

    /// Value at `x in [a, b]`.
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        clenshaw(&self.coeffs, (2.0 * x - self.a - self.b) / (self.b - self.a))
    }

    /// Derivative series.
    #[must_use]
    pub fn derivative(&self) -> Self {
        let n = self.coeffs.len();
        let mut d = vec![0.0; n.max(2) + 1];
        for k in (1..n).rev() {
            d[k - 1] = d[k + 1] + 2.0 * k as f64 * self.coeffs[k];
        }
        d.truncate(n.saturating_sub(1).max(1));
        d[0] *= 0.5;
        let s = 2.0 / (self.b - self.a);
        Self { coeffs: d.iter().map(|v| v * s).collect(), a: self.a, b: self.b }
    }

    /// Antiderivative series vanishing at `a`.
    #[must_use]
    pub fn integral(&self) -> Self {
        let n = self.coeffs.len();
        let scale = 0.5 * (self.b - self.a);
        let mut c = vec![0.0; n + 1];
        let get = |k: usize| self.coeffs.get(k).copied().unwrap_or(0.0);
        for k in 1..=n {
            let hi = get(k + 1);
            let lo = if k == 1 { 2.0 * get(0) } else { get(k - 1) };
            c[k] = (lo - hi) / (2.0 * k as f64) * scale;
        }
        let at_minus_one: f64 =
            c.iter().enumerate().skip(1).map(|(k, v)| if k % 2 == 0 { *v } else { -*v }).sum();
        c[0] = -at_minus_one;
        Self { coeffs: c, a: self.a, b: self.b }
    }

    /// Definite integral over `[a, b]`.
    #[must_use]
    pub fn definite_integral(&self) -> f64 {
        let s = 0.5 * (self.b - self.a);
        let mut sum = 0.0;
        for (k, &c) in self.coeffs.iter().enumerate().step_by(2) {
            sum += c * 2.0 / (1.0 - (k * k) as f64);
        }
        s * sum
    }

    /// Real roots in `[a, b]` from the eigenvalues of the colleague matrix.
    ///
    /// # Errors
    /// [`ApproxError::Singular`] if the eigenvalue iteration fails.
    pub fn roots(&self) -> Result<Vec<f64>, ApproxError> {
        let mut c = self.coeffs.clone();
        while c.len() > 1 && c.last().is_some_and(|v| *v == 0.0) {
            c.pop();
        }
        let n = c.len() - 1;
        if n == 0 {
            return Ok(Vec::new());
        }
        let mut m = Mat::zeros(n, n);
        for i in 0..n - 1 {
            m.set(i, i + 1, 0.5);
            m.set(i + 1, i, 0.5);
        }
        m.set(0, 1, if n > 1 { 1.0 } else { 0.0 });
        if n == 1 {
            m.set(0, 0, -c[0] / c[1]);
        } else {
            for j in 0..n {
                let v = m.at(n - 1, j) - 0.5 * c[j] / c[n];
                m.set(n - 1, j, v);
            }
        }
        let ev = dense::eigenvalues(&m).map_err(|_| ApproxError::Singular)?;
        let mut out: Vec<f64> = ev
            .iter()
            .filter(|&&(re, im)| im.abs() < 1e-8 && re.abs() <= 1.0 + 1e-10)
            .map(|&(re, _)| 0.5 * (self.a + self.b) + 0.5 * (self.b - self.a) * re)
            .collect();
        out.sort_by(f64::total_cmp);
        Ok(out)
    }
}

/// Barycentric rational approximant produced by [`aaa`].
#[derive(Debug, Clone)]
pub struct Aaa {
    /// Support points.
    pub z: Vec<f64>,
    /// Function values at the support points.
    pub f: Vec<f64>,
    /// Barycentric weights.
    pub w: Vec<f64>,
    /// Final maximal sample error.
    pub error: f64,
}

impl Aaa {
    /// Evaluates the rational function at `x`.
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        let (mut num, mut den) = (0.0, 0.0);
        for k in 0..self.z.len() {
            let d = x - self.z[k];
            if d == 0.0 {
                return self.f[k];
            }
            let t = self.w[k] / d;
            num += t * self.f[k];
            den += t;
        }
        num / den
    }

    /// Degree (`numerator = denominator = m - 1`).
    #[must_use]
    pub fn degree(&self) -> usize {
        self.z.len().saturating_sub(1)
    }
}

/// AAA algorithm (Nakatsukasa, Sete, Trefethen) on real samples
/// `(z_i, f_i)`: greedy support-point selection and a least-squares
/// weight vector from the SVD of the Loewner matrix.
///
/// # Errors
/// [`ApproxError::Invalid`] for fewer than 3 samples.
pub fn aaa(z: &[f64], f: &[f64], tol: f64, max_terms: usize) -> Result<Aaa, ApproxError> {
    let m = z.len();
    if m < 3 || f.len() != m {
        return Err(ApproxError::Invalid);
    }
    let fmax = f.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
    let mut rest: Vec<usize> = (0..m).collect();
    let mean = f.iter().sum::<f64>() / m as f64;
    let mut r = vec![mean; m];
    let (mut zj, mut fj): (Vec<f64>, Vec<f64>) = (Vec::new(), Vec::new());
    let mut w = Vec::new();
    let mut err = f64::INFINITY;
    for _ in 0..max_terms.min(m - 1) {
        let (pos, &j) = rest
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| (f[**a] - r[**a]).abs().total_cmp(&(f[**b] - r[**b]).abs()))
            .unwrap_or((0, &0));
        zj.push(z[j]);
        fj.push(f[j]);
        rest.remove(pos);
        let k = zj.len();
        if rest.len() < k {
            break;
        }
        let mut a = Mat::zeros(rest.len(), k);
        for (row, &i) in rest.iter().enumerate() {
            for c in 0..k {
                a.set(row, c, (f[i] - fj[c]) / (z[i] - zj[c]));
            }
        }
        let Ok(s) = dense::svd(&a) else { break };
        let last = k - 1;
        w = (0..k).map(|c| s.v.at(c, last)).collect();
        for i in 0..m {
            let (mut num, mut den) = (0.0, 0.0);
            let mut hit = None;
            for c in 0..k {
                let d = z[i] - zj[c];
                if d == 0.0 {
                    hit = Some(c);
                    break;
                }
                num += w[c] * fj[c] / d;
                den += w[c] / d;
            }
            r[i] = hit.map_or(num / den, |c| fj[c]);
        }
        err = rest.iter().map(|&i| (f[i] - r[i]).abs()).fold(0.0, f64::max);
        if err <= tol * fmax {
            break;
        }
    }
    Ok(Aaa { z: zj, f: fj, w, error: err })
}

/// Pade approximant `[m/n]` from Taylor coefficients `c_0, c_1, ...`
/// (needs `m + n + 1` coefficients). Returns `(numerator, denominator)`
/// with ascending coefficients and `denominator[0] = 1`.
///
/// # Errors
/// [`ApproxError::Invalid`] if too few coefficients, [`ApproxError::Singular`]
/// if the approximant does not exist.
pub fn pade(c: &[f64], m: usize, n: usize) -> Result<(Vec<f64>, Vec<f64>), ApproxError> {
    if c.len() < m + n + 1 {
        return Err(ApproxError::Invalid);
    }
    let cc = |i: isize| if i < 0 { 0.0 } else { c[i as usize] };
    let mut q = vec![1.0];
    if n > 0 {
        let mut a = Mat::zeros(n, n);
        let mut b = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                a.set(i, j, cc((m + i) as isize - j as isize));
            }
            b[i] = -cc((m + i + 1) as isize);
        }
        q.extend(dense::solve(&a, &b).map_err(|_| ApproxError::Singular)?);
    }
    let p: Vec<f64> = (0..=m)
        .map(|k| (0..=k.min(n)).map(|j| q[j] * cc((k - j) as isize)).sum())
        .collect();
    Ok((p, q))
}

/// Evaluates a polynomial with ascending coefficients by Horner's scheme.
#[must_use]
pub fn horner(c: &[f64], x: f64) -> f64 {
    c.iter().rev().fold(0.0, |acc, &v| acc * x + v)
}
