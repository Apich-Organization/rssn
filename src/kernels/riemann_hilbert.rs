//! Numerical solvers for Riemann–Hilbert problems.
//!
//! Includes scalar Cauchy and Plemelj–Sokhotski integrals, scalar arc factorization, and Olver’s
//! spectral collocation solver for `2 × 2` matrix problems on `(-1, 1)`.

use crate::kernels::complex_scalar::Complex as ScalarComplex;
use num_complex::Complex;
use std::f64::consts::PI;

#[path = "riemann_hilbert_scalar.rs"]
pub(crate) mod scalar;

pub(crate) use scalar::ScalarRhError;
pub(crate) use scalar::cauchy_transform_arc_checked;
pub(crate) use scalar::plemelj_boundary_values_checked;
pub(crate) use scalar::scalar_factorization_checked;

fn complex_quotient(
    numerator: Complex<f64>,
    denominator: Complex<f64>,
) -> Complex<f64> {
    let value = ScalarComplex::new(numerator.re, numerator.im)
        .div(ScalarComplex::new(denominator.re, denominator.im));
    Complex::new(value.re, value.im)
}

/// 2x2 complex matrix for Riemann–Hilbert jump problems.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix2 {
    /// Row-major entries `[m00, m01, m10, m11]`.
    pub m: [Complex<f64>; 4],
}

impl Matrix2 {
    /// Identity 2x2 matrix.
    #[must_use]
    pub const fn identity() -> Self {
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        Self {
            m: [one, zero, zero, one],
        }
    }

    /// Zero 2x2 matrix.
    #[must_use]
    pub const fn zero() -> Self {
        let z = Complex::new(0.0, 0.0);
        Self { m: [z, z, z, z] }
    }

    /// Matrix from 4 complex entries `[m00, m01, m10, m11]`.
    #[must_use]
    pub const fn new(
        m00: Complex<f64>,
        m01: Complex<f64>,
        m10: Complex<f64>,
        m11: Complex<f64>,
    ) -> Self {
        Self {
            m: [m00, m01, m10, m11],
        }
    }

    /// Matrix addition.
    #[must_use]
    #[expect(
        clippy::should_implement_trait,
        reason = "Explicit matrix methods are part of the borrowed-operand API"
    )]
    pub fn add(
        self,
        other: Self,
    ) -> Self {
        Self {
            m: [
                self.m[0] + other.m[0],
                self.m[1] + other.m[1],
                self.m[2] + other.m[2],
                self.m[3] + other.m[3],
            ],
        }
    }

    /// Matrix subtraction.
    #[must_use]
    #[expect(
        clippy::should_implement_trait,
        reason = "Explicit matrix methods are part of the borrowed-operand API"
    )]
    pub fn sub(
        self,
        other: Self,
    ) -> Self {
        Self {
            m: [
                self.m[0] - other.m[0],
                self.m[1] - other.m[1],
                self.m[2] - other.m[2],
                self.m[3] - other.m[3],
            ],
        }
    }

    /// Matrix multiplication.
    #[must_use]
    #[expect(
        clippy::should_implement_trait,
        reason = "Explicit matrix methods are part of the borrowed-operand API"
    )]
    pub fn mul(
        self,
        other: Self,
    ) -> Self {
        Self {
            m: [
                self.m[0] * other.m[0] + self.m[1] * other.m[2],
                self.m[0] * other.m[1] + self.m[1] * other.m[3],
                self.m[2] * other.m[0] + self.m[3] * other.m[2],
                self.m[2] * other.m[1] + self.m[3] * other.m[3],
            ],
        }
    }

    /// Scalar multiplication.
    #[must_use]
    pub fn scale(
        self,
        s: Complex<f64>,
    ) -> Self {
        Self {
            m: [self.m[0] * s, self.m[1] * s, self.m[2] * s, self.m[3] * s],
        }
    }

    /// Matrix inverse with row equilibration. Returns `None` for a numerically
    /// singular normalized matrix or when its inverse is outside binary64 range.
    #[must_use]
    pub fn inv(self) -> Option<Self> {
        let [a, b, c, d] = self.m;
        let [upper, lower] = [[a, b], [c, d]].map(|row| {
            row.into_iter()
                .flat_map(|z| [z.re.abs(), z.im.abs()])
                .fold(0.0_f64, f64::max)
        });
        if upper == 0.0 || lower == 0.0 || self.m.iter().any(|z| !z.is_finite()) {
            return None;
        }
        let [a, b, c, d] = [a / upper, b / upper, c / lower, d / lower];
        let det = a * d - b * c;
        if det.norm() < 1e-15 {
            return None;
        }
        let m = [
            complex_quotient(d, det) / upper,
            complex_quotient(-b, det) / lower,
            complex_quotient(-c, det) / upper,
            complex_quotient(a, det) / lower,
        ];
        m.iter().all(|z| z.is_finite()).then_some(Self { m })
    }

    /// Frobenius norm of the 2x2 matrix.
    #[must_use]
    pub fn norm(self) -> f64 {
        let squared = self.m.iter().map(Complex::norm_sqr).sum::<f64>();
        if squared.is_normal() {
            squared.sqrt()
        } else {
            self.m.iter().fold(0.0_f64, |size, z| size.hypot(z.norm()))
        }
    }
}

/// Evaluates the scalar Cauchy transform on the normalized arc, at a point off `[-1, 1]`.
/// Returns complex NaN when the estimate fails or does not meet `tol`.
pub fn cauchy_transform_arc<F>(
    phi: &F,
    z: Complex<f64>,
    tol: f64,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    cauchy_transform_arc_checked(phi, z, tol, scalar::DEFAULT_PANELS)
        .ok()
        .filter(|result| result.converged)
        .map_or(Complex::new(f64::NAN, f64::NAN), |result| result.value)
}

/// Computes the Plemelj boundary limits on the open arc. Their jump is `phi(x)`.
///
/// # Errors
/// Rejects invalid inputs, nonfinite results, or failure to meet the absolute tolerance.
pub fn plemelj_boundary_values<F>(
    phi: &F,
    x: f64,
    tol: f64,
) -> Result<(Complex<f64>, Complex<f64>), String>
where
    F: Fn(f64) -> Complex<f64>,
{
    let result = plemelj_boundary_values_checked(phi, x, tol, scalar::DEFAULT_PANELS)
        .map_err(|error| error.to_string())?;
    if result.converged {
        Ok((result.upper, result.lower))
    } else {
        Err("Plemelj quadrature did not meet its absolute tolerance".into())
    }
}

/// Factors a nonzero scalar jump on `(-1, 1)` with normalization at infinity.
///
/// The phase starts at the principal value at `-1`. The phase mesh must resolve the jump.
/// Returns complex NaN when the estimate fails or does not meet `tol`.
pub fn scalar_rh_factorize_arc<F>(
    g: &F,
    z: Complex<f64>,
    tol: f64,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    scalar_factorization_checked(g, z, tol, 128, scalar::DEFAULT_PANELS)
        .ok()
        .filter(|result| result.estimate.converged)
        .map_or(Complex::new(f64::NAN, f64::NAN), |result| {
            result.estimate.value
        })
}

/// Barycentric differentiation matrix for arbitrary distinct nodes `x`.
fn barycentric_diff_matrix(x: &[f64]) -> Result<Vec<Vec<f64>>, String> {
    let n = x.len();
    let bw = super::barycentric::weights(x)?;
    let mut d = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        let xi = *x.get(i).unwrap_or(&0.0);
        let bi = *bw.get(i).unwrap_or(&1.0);
        let mut row_sum = 0.0_f64;
        for j in 0..n {
            if i != j {
                let xj = *x.get(j).unwrap_or(&0.0);
                let bj = *bw.get(j).unwrap_or(&1.0);
                let val = (bj / bi) / (xi - xj);
                row_sum += val;
                if let Some(cell) = d.get_mut(i).and_then(|r| r.get_mut(j)) {
                    *cell = val;
                }
            }
        }
        if let Some(cell) = d.get_mut(i).and_then(|r| r.get_mut(i)) {
            *cell = -row_sum;
        }
    }
    Ok(d)
}

/// Solution to a 2x2 matrix Riemann–Hilbert problem $Y^+(t) = Y^-(t) G(t)$ on $(-1, 1)$
/// with $Y(\infty) = I$, represented by its jump density $U(t) = Y^+(t) - Y^-(t)$
/// at Gauss–Legendre collocation nodes.
#[derive(Debug, Clone)]
pub struct MatrixRhSolution {
    /// Collocation nodes on $(-1, 1)$.
    pub nodes: Vec<f64>,
    /// Quadrature weights on $(-1, 1)$.
    pub weights: Vec<f64>,
    /// Nodal jump density $`U(x_j)` = `Y^+(x_j)` - Y^-(x_j)$.
    pub density: Vec<Matrix2>,
}

impl MatrixRhSolution {
    fn eval_density_at(
        &self,
        x: f64,
    ) -> Matrix2 {
        if let Some((_, &value)) = self
            .nodes
            .iter()
            .zip(&self.density)
            .find(|(node, _)| (x - **node).abs() < 1e-14)
        {
            return value;
        }
        let Ok(bw) = super::barycentric::weights(&self.nodes) else {
            return Matrix2::zero().scale(Complex::new(f64::NAN, f64::NAN));
        };
        let mut num = Matrix2::zero();
        let mut den = 0.0_f64;
        for ((&xj, &bj), &uj) in self.nodes.iter().zip(bw.iter()).zip(self.density.iter()) {
            let diff = x - xj;
            let w = bj / diff;
            den += w;
            num = num.add(uj.scale(Complex::new(w, 0.0)));
        }
        num.scale(Complex::new(1.0 / den, 0.0))
    }

    /// Evaluates $Y(z) = I + \frac{1}{2\pi i} \int_{-1}^1 \frac{U(t)}{t - z} \, dt$ off the arc,
    /// using singularity subtraction near `[-1, 1]` for high accuracy up to the boundary.
    #[must_use]
    pub fn eval(
        &self,
        z: Complex<f64>,
    ) -> Matrix2 {
        let two_pi_i = Complex::new(0.0, 2.0 * PI);
        if z.re.abs() < 1.0 && z.im.abs() < 0.25 && z.im.abs() > 0.0 {
            let u_x = self.eval_density_at(z.re);
            let mut acc = Matrix2::zero();
            for ((&t, &w), &u) in self
                .nodes
                .iter()
                .zip(self.weights.iter())
                .zip(self.density.iter())
            {
                let difference = u.sub(u_x);
                if difference != Matrix2::zero() {
                    let factor = complex_quotient(Complex::new(w, 0.0), Complex::new(t, 0.0) - z);
                    acc = acc.add(difference.scale(factor));
                }
            }
            // int_{-1}^1 dt/(t - z) = ln(1 - z) - ln(-1 - z)
            let log_jump = (Complex::new(1.0, 0.0) - z).ln() - (Complex::new(-1.0, 0.0) - z).ln();
            let total = acc.add(u_x.scale(log_jump));
            return Matrix2::identity().add(total.scale(Complex::new(1.0, 0.0) / two_pi_i));
        }
        let mut acc = Matrix2::zero();
        for ((&t, &w), &u) in self
            .nodes
            .iter()
            .zip(self.weights.iter())
            .zip(self.density.iter())
        {
            let factor = complex_quotient(Complex::new(w, 0.0), Complex::new(t, 0.0) - z);
            acc = acc.add(u.scale(factor));
        }
        Matrix2::identity().add(acc.scale(Complex::new(1.0, 0.0) / two_pi_i))
    }
}

/// Solves a two-by-two matrix Riemann–Hilbert problem on `(-1, 1)`.
///
/// Solves the 2x2 matrix Riemann–Hilbert problem $Y^+(t) = Y^-(t) G(t)$ on $(-1, 1)$
/// with $Y(\infty) = I$ (assuming $G(\pm 1) = I$) via spectral collocation on `n_nodes`
/// Gauss–Legendre nodes using the singular integral equation $U - C_-\[U\](G - I) = G - I$.
///
/// # Errors
///
/// Returns an error for fewer than four collocation nodes or when the collocation system cannot be assembled or solved.
pub fn solve_matrix_rh_chebyshev<F>(
    jump_g: F,
    n_nodes: usize,
) -> Result<MatrixRhSolution, String>
where
    F: Fn(f64) -> Matrix2,
{
    if n_nodes < 4 {
        return Err("n_nodes must be >= 4".to_string());
    }
    let plan = super::reuse::global_pool().gauss_legendre_rule(n_nodes);
    let (nodes, weights) = plan.as_ref();
    let diff = barycentric_diff_matrix(nodes)?;
    let two_pi_i = Complex::new(0.0, 2.0 * PI);

    // Build scalar C_- matrix of size n_nodes x n_nodes:
    // (C_- u)_i = (1/(2*pi*i)) * (sum_{m!=i} (u_m - u_i)/(x_m - x_i) w_m + (D u)_i w_i + u_i ln((1-x_i)/(1+x_i))) - 0.5 u_i
    let mut c_minus = vec![vec![Complex::new(0.0, 0.0); n_nodes]; n_nodes];
    for i in 0..n_nodes {
        let xi = *nodes.get(i).ok_or_else(|| "Index".to_string())?;
        let wi = *weights.get(i).ok_or_else(|| "Index".to_string())?;
        let mut diag_sum = Complex::new(0.0, 0.0);
        for m in 0..n_nodes {
            let xm = *nodes.get(m).ok_or_else(|| "Index".to_string())?;
            let wm = *weights.get(m).ok_or_else(|| "Index".to_string())?;
            let d_im = *diff
                .get(i)
                .and_then(|r| r.get(m))
                .ok_or_else(|| "Index".to_string())?;
            let off = if i == m {
                0.0
            } else {
                wm / (xm - xi)
            };
            let entry = Complex::new(off + wi * d_im, 0.0) / two_pi_i;
            if let Some(cell) = c_minus.get_mut(i).and_then(|r| r.get_mut(m)) {
                *cell = entry;
            }
            if i != m {
                diag_sum += Complex::new(wm / (xm - xi), 0.0);
            }
        }
        let log_term = Complex::new(((1.0 - xi) / (1.0 + xi)).ln(), 0.0);
        if let Some(cell) = c_minus.get_mut(i).and_then(|r| r.get_mut(i)) {
            *cell += (log_term - diag_sum) / two_pi_i - Complex::new(0.5, 0.0);
        }
    }

    // Solve 2n x 2n complex linear system per column of (U_00, U_10) and (U_01, U_11):
    // Notice (U - C_-[U] * K)_i = K_i where K_i = G(x_i) - I:
    // Row 0 of U (u_00, u_01) and Row 1 of U (u_10, u_11) decouple!
    // For each row r in {0, 1} of U:
    // [u_{r,0}(x_i), u_{r,1}(x_i)] - sum_m C_{i,m} [u_{r,0}(x_m), u_{r,1}(x_m)] * K_i = [K_{r,0}(x_i), K_{r,1}(x_i)].
    let mut k_mats = Vec::with_capacity(n_nodes);
    for &x in nodes {
        if crate::cancellation::is_cancelled() {
            return Err("Riemann–Hilbert solve cancelled".into());
        }
        let jump = jump_g(x);
        if jump.m.iter().any(|z| !z.is_finite()) {
            return Err("Riemann–Hilbert jump is nonfinite".into());
        }
        k_mats.push(jump.sub(Matrix2::identity()));
    }

    let dim = 2 * n_nodes;
    let mut sys = vec![vec![Complex::new(0.0, 0.0); dim]; dim];
    for i in 0..n_nodes {
        let ki = *k_mats.get(i).ok_or_else(|| "Index".to_string())?;
        let r0 = 2 * i;
        let r1 = 2 * i + 1;
        for m in 0..n_nodes {
            let c_im = *c_minus
                .get(i)
                .and_then(|r| r.get(m))
                .ok_or_else(|| "Index".to_string())?;
            let c0 = 2 * m;
            let c1 = 2 * m + 1;
            // Subtract c_im * K_i^T from block (i, m)
            if let Some(cell) = sys.get_mut(r0).and_then(|r| r.get_mut(c0)) {
                *cell -= c_im * ki.m[0];
            }
            if let Some(cell) = sys.get_mut(r0).and_then(|r| r.get_mut(c1)) {
                *cell -= c_im * ki.m[2];
            }
            if let Some(cell) = sys.get_mut(r1).and_then(|r| r.get_mut(c0)) {
                *cell -= c_im * ki.m[1];
            }
            if let Some(cell) = sys.get_mut(r1).and_then(|r| r.get_mut(c1)) {
                *cell -= c_im * ki.m[3];
            }
        }
        if let Some(cell) = sys.get_mut(r0).and_then(|r| r.get_mut(r0)) {
            *cell += Complex::new(1.0, 0.0);
        }
        if let Some(cell) = sys.get_mut(r1).and_then(|r| r.get_mut(r1)) {
            *cell += Complex::new(1.0, 0.0);
        }
    }

    let rhs = k_mats
        .iter()
        .flat_map(|k| [[k.m[0], k.m[2]], [k.m[1], k.m[3]]])
        .collect();
    let solved = solve_linear_c(sys, rhs)?;
    let (pairs, _) = solved.as_chunks::<2>();
    let density = pairs
        .iter()
        .map(|[first, second]| Matrix2::new(first[0], second[0], first[1], second[1]))
        .collect();
    Ok(MatrixRhSolution {
        nodes: nodes.clone(),
        weights: weights.clone(),
        density,
    })
}

// Eliminate the common coefficient system once, updating both right-hand sides.
fn solve_linear_c(
    mut matrix: Vec<Vec<Complex<f64>>>,
    mut rhs: Vec<[Complex<f64>; 2]>,
) -> Result<Vec<[Complex<f64>; 2]>, String> {
    for col in 0..matrix.len() {
        if crate::cancellation::is_cancelled() {
            return Err("Riemann–Hilbert solve cancelled".into());
        }
        let (pivot, magnitude) = matrix
            .iter()
            .enumerate()
            .skip(col)
            .map(|(index, row)| (index, row[col].norm()))
            .max_by(|(a, va), (b, vb)| va.total_cmp(vb).then_with(|| b.cmp(a)))
            .unwrap_or((col, 0.0));
        if !magnitude.is_finite() || magnitude < 1e-14 {
            return Err("Riemann–Hilbert collocation pivot is singular or nonfinite".into());
        }
        matrix.swap(col, pivot);
        rhs.swap(col, pivot);
        let (top, remaining) = matrix.split_at_mut(col + 1);
        let pivot_row = &top[col];
        let (top_rhs, remaining_rhs) = rhs.split_at_mut(col + 1);
        let [first, second] = top_rhs[col];
        let pivot_value = pivot_row[col];
        for (row, current) in remaining.iter_mut().zip(remaining_rhs) {
            let factor = complex_quotient(row[col], pivot_value);
            for (entry, &value) in row[col + 1..].iter_mut().zip(&pivot_row[col + 1..]) {
                *entry -= factor * value;
            }
            current[0] -= factor * first;
            current[1] -= factor * second;
        }
    }
    for (index, row) in matrix.iter().enumerate().rev() {
        let (current, following) = rhs.split_at_mut(index + 1);
        let [mut first, mut second] = current[index];
        for (&entry, &[x, y]) in row[index + 1..].iter().zip(following.iter()) {
            first -= entry * x;
            second -= entry * y;
        }
        let solved = [
            complex_quotient(first, row[index]),
            complex_quotient(second, row[index]),
        ];
        if solved.iter().any(|z| !z.is_finite()) {
            return Err("Riemann–Hilbert density is outside the finite binary64 range".into());
        }
        current[index] = solved;
    }
    Ok(rhs)
}
