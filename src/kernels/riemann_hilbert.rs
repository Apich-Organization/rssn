//! Riemann–Hilbert problem (RHP) solvers:
//! scalar Cauchy/Plemelj–Sokhotski integrals, scalar arc factorization,
//! and Olver's Chebyshev/Gauss–Legendre spectral collocation solver for
//! 2x2 matrix Riemann–Hilbert problems on (-1, 1).

use crate::kernels::integrate::gauss_kronrod_any;
use num_complex::Complex;
use std::f64::consts::PI;

/// 2x2 complex matrix for Riemann–Hilbert jump problems.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Matrix2 {
    /// Row-major entries `[m00, m01, m10, m11]`.
    pub m: [Complex<f64>; 4],
}

impl Matrix2 {
    /// Identity 2x2 matrix.
    pub const fn identity() -> Self {
        let one = Complex::new(1.0, 0.0);
        let zero = Complex::new(0.0, 0.0);
        Self {
            m: [one, zero, zero, one],
        }
    }

    /// Zero 2x2 matrix.
    pub const fn zero() -> Self {
        let z = Complex::new(0.0, 0.0);
        Self { m: [z, z, z, z] }
    }

    /// Matrix from 4 complex entries `[m00, m01, m10, m11]`.
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
    pub fn scale(
        self,
        s: Complex<f64>,
    ) -> Self {
        Self {
            m: [self.m[0] * s, self.m[1] * s, self.m[2] * s, self.m[3] * s],
        }
    }

    /// Matrix inverse.
    pub fn inv(self) -> Option<Self> {
        let det = self.m[0] * self.m[3] - self.m[1] * self.m[2];
        if det.norm() < 1e-15 {
            return None;
        }
        Some(Self {
            m: [
                self.m[3] / det,
                -self.m[1] / det,
                -self.m[2] / det,
                self.m[0] / det,
            ],
        })
    }

    /// Frobenius norm of the 2x2 matrix.
    pub fn norm(self) -> f64 {
        (self.m[0].norm_sqr() + self.m[1].norm_sqr() + self.m[2].norm_sqr() + self.m[3].norm_sqr())
            .sqrt()
    }
}

/// Evaluates the scalar Cauchy transform $C[\phi](z) = \frac{1}{2\pi i} \int_{-1}^1 \frac{\phi(t)}{t - z} \, dt$
/// for $z \notin [-1, 1]$.
pub fn cauchy_transform_arc<F>(
    phi: &F,
    z: Complex<f64>,
    tol: f64,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    let re = gauss_kronrod_any(
        &|t| (phi(t) / (Complex::new(t, 0.0) - z)).re,
        -1.0,
        1.0,
        tol,
        300,
    )
    .value;
    let im = gauss_kronrod_any(
        &|t| (phi(t) / (Complex::new(t, 0.0) - z)).im,
        -1.0,
        1.0,
        tol,
        300,
    )
    .value;
    Complex::new(re, im) / Complex::new(0.0, 2.0 * PI)
}

/// Computes the Plemelj–Sokhotski boundary limits $(C_+[\phi](x), C_-[\phi](x))$
/// on $x \in (-1, 1)$ via singularity subtraction: $C_+ - C_- = \phi(x)$.
pub fn plemelj_boundary_values<F>(
    phi: &F,
    x: f64,
    tol: f64,
) -> Result<(Complex<f64>, Complex<f64>), String>
where
    F: Fn(f64) -> Complex<f64>,
{
    if x.abs() >= 1.0 {
        return Err("x must lie strictly inside (-1, 1)".to_string());
    }
    let px = phi(x);
    let diff_re = gauss_kronrod_any(
        &|t: f64| {
            let dt = t - x;
            if dt.abs() < 1e-12 {
                0.0
            } else {
                ((phi(t) - px) / dt).re
            }
        },
        -1.0,
        1.0,
        tol,
        300,
    )
    .value;
    let diff_im = gauss_kronrod_any(
        &|t: f64| {
            let dt = t - x;
            if dt.abs() < 1e-12 {
                0.0
            } else {
                ((phi(t) - px) / dt).im
            }
        },
        -1.0,
        1.0,
        tol,
        300,
    )
    .value;
    let pv = Complex::new(diff_re, diff_im) + px * ((1.0 - x) / (1.0 + x)).ln();
    let cauchy_pv = pv / Complex::new(0.0, 2.0 * PI);
    Ok((cauchy_pv + px * 0.5, cauchy_pv - px * 0.5))
}

/// Solves the scalar Riemann–Hilbert problem $X^+(t) = G(t) X^-(t)$ on $(-1, 1)$
/// with $X(\infty) = 1$ via continuous phase unwrapping of $\ln G(t)$.
pub fn scalar_rh_factorize_arc<F>(
    g: &F,
    z: Complex<f64>,
    tol: f64,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    // Sample [-1, 1] at 129 points to unwrap arg(G(t)) continuously
    let n_grid = 128_usize;
    let mut unwrapped_args = Vec::with_capacity(n_grid + 1);
    let mut prev_arg = g(-1.0).arg();
    unwrapped_args.push(prev_arg);
    for k in 1..=n_grid {
        let t = -1.0 + 2.0 * (k as f64) / (n_grid as f64);
        let mut a = g(t).arg();
        while a - prev_arg > PI {
            a -= 2.0 * PI;
        }
        while a - prev_arg < -PI {
            a += 2.0 * PI;
        }
        prev_arg = a;
        unwrapped_args.push(a);
    }
    cauchy_transform_arc(
        &|t| {
            let gt = g(t);
            let u = ((t + 1.0) * 0.5 * (n_grid as f64)).clamp(0.0, n_grid as f64);
            let idx = (u.floor() as usize).min(n_grid - 1);
            let frac = u - (idx as f64);
            let ref_arg = unwrapped_args[idx] * (1.0 - frac) + unwrapped_args[idx + 1] * frac;
            let mut a = gt.arg();
            while a - ref_arg > PI {
                a -= 2.0 * PI;
            }
            while a - ref_arg < -PI {
                a += 2.0 * PI;
            }
            Complex::new(gt.norm().ln(), a)
        },
        z,
        tol,
    )
    .exp()
}

/// Computes $n$-point Gauss–Legendre nodes and weights on $(-1, 1)$.
fn gauss_legendre_nodes_weights(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut nodes = vec![0.0_f64; n];
    let mut weights = vec![0.0_f64; n];
    let m = (n + 1) / 2;
    for i in 0..m {
        let mut z = (PI * ((i as f64) + 0.75) / ((n as f64) + 0.5)).cos();
        let mut dp = 0.0_f64;
        for _ in 0..40 {
            let mut p1 = 1.0_f64;
            let mut p2 = 0.0_f64;
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                let jf = j as f64;
                p1 = ((2.0 * jf - 1.0) * z * p2 - (jf - 1.0) * p3) / jf;
            }
            dp = (n as f64) * (z * p1 - p2) / (z * z - 1.0);
            let z_next = z - p1 / dp;
            if (z_next - z).abs() < 1e-15 {
                z = z_next;
                break;
            }
            z = z_next;
        }
        let w = 2.0 / ((1.0 - z * z) * dp * dp);
        if let Some(slot) = nodes.get_mut(i) {
            *slot = -z;
        }
        if let Some(slot) = nodes.get_mut(n - 1 - i) {
            *slot = z;
        }
        if let Some(slot) = weights.get_mut(i) {
            *slot = w;
        }
        if let Some(slot) = weights.get_mut(n - 1 - i) {
            *slot = w;
        }
    }
    (nodes, weights)
}

/// Barycentric weights for distinct nodes `x`.
fn barycentric_weights(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut bw = vec![1.0_f64; n];
    for j in 0..n {
        let xj = *x.get(j).unwrap_or(&0.0);
        let mut prod = 1.0_f64;
        for (k, &xk) in x.iter().enumerate() {
            if k != j {
                prod *= xj - xk;
            }
        }
        if let Some(slot) = bw.get_mut(j) {
            *slot = 1.0 / prod;
        }
    }
    bw
}

/// Barycentric differentiation matrix for arbitrary distinct nodes `x`.
fn barycentric_diff_matrix(x: &[f64]) -> Vec<Vec<f64>> {
    let n = x.len();
    let bw = barycentric_weights(x);
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
    d
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
    /// Nodal jump density $U(x_j) = Y^+(x_j) - Y^-(x_j)$.
    pub density: Vec<Matrix2>,
}

impl MatrixRhSolution {
    fn eval_density_at(
        &self,
        x: f64,
    ) -> Matrix2 {
        let bw = barycentric_weights(&self.nodes);
        let mut num = Matrix2::zero();
        let mut den = 0.0_f64;
        for ((&xj, &bj), &uj) in self.nodes.iter().zip(bw.iter()).zip(self.density.iter()) {
            let diff = x - xj;
            if diff.abs() < 1e-14 {
                return uj;
            }
            let w = bj / diff;
            den += w;
            num = num.add(uj.scale(Complex::new(w, 0.0)));
        }
        num.scale(Complex::new(1.0 / den, 0.0))
    }

    /// Evaluates $Y(z) = I + \frac{1}{2\pi i} \int_{-1}^1 \frac{U(t)}{t - z} \, dt$ off the arc,
    /// using singularity subtraction near `[-1, 1]` for high accuracy up to the boundary.
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
                let factor = Complex::new(w, 0.0) / (Complex::new(t, 0.0) - z);
                acc = acc.add(u.sub(u_x).scale(factor));
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
            let factor = Complex::new(w, 0.0) / (Complex::new(t, 0.0) - z);
            acc = acc.add(u.scale(factor));
        }
        Matrix2::identity().add(acc.scale(Complex::new(1.0, 0.0) / two_pi_i))
    }
}

/// Solves the 2x2 matrix Riemann–Hilbert problem $Y^+(t) = Y^-(t) G(t)$ on $(-1, 1)$
/// with $Y(\infty) = I$ (assuming $G(\pm 1) = I$) via spectral collocation on `n_nodes`
/// Gauss–Legendre nodes using the singular integral equation $U - C_-[U](G - I) = G - I$.
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
    let (nodes, weights) = gauss_legendre_nodes_weights(n_nodes);
    let diff = barycentric_diff_matrix(&nodes);
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
            let off = if i != m {
                wm / (xm - xi)
            } else {
                0.0
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
    let k_mats: Vec<Matrix2> = nodes
        .iter()
        .map(|&x| jump_g(x).sub(Matrix2::identity()))
        .collect();

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

    let mut rhs_row0 = vec![Complex::new(0.0, 0.0); dim];
    let mut rhs_row1 = vec![Complex::new(0.0, 0.0); dim];
    for (i, ki) in k_mats.iter().enumerate() {
        if let Some(s) = rhs_row0.get_mut(2 * i) {
            *s = ki.m[0];
        }
        if let Some(s) = rhs_row0.get_mut(2 * i + 1) {
            *s = ki.m[1];
        }
        if let Some(s) = rhs_row1.get_mut(2 * i) {
            *s = ki.m[2];
        }
        if let Some(s) = rhs_row1.get_mut(2 * i + 1) {
            *s = ki.m[3];
        }
    }

    let sol0 = solve_linear_c(&sys, &rhs_row0)?;
    let sol1 = solve_linear_c(&sys, &rhs_row1)?;
    let mut density = Vec::with_capacity(n_nodes);
    for i in 0..n_nodes {
        let u00 = *sol0.get(2 * i).ok_or_else(|| "Index".to_string())?;
        let u01 = *sol0.get(2 * i + 1).ok_or_else(|| "Index".to_string())?;
        let u10 = *sol1.get(2 * i).ok_or_else(|| "Index".to_string())?;
        let u11 = *sol1.get(2 * i + 1).ok_or_else(|| "Index".to_string())?;
        density.push(Matrix2::new(u00, u01, u10, u11));
    }

    Ok(MatrixRhSolution {
        nodes,
        weights,
        density,
    })
}

fn solve_linear_c(
    a: &[Vec<Complex<f64>>],
    b: &[Complex<f64>],
) -> Result<Vec<Complex<f64>>, String> {
    let n = b.len();
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        let mut max_v = mat
            .get(col)
            .and_then(|r| r.get(col))
            .map(|c| c.norm())
            .unwrap_or(0.0);
        for row in (col + 1)..n {
            let v = mat
                .get(row)
                .and_then(|r| r.get(col))
                .map(|c| c.norm())
                .unwrap_or(0.0);
            if v > max_v {
                max_v = v;
                piv = row;
            }
        }
        if max_v < 1e-14 {
            return Err("Singular matrix in RH collocation solver".to_string());
        }
        if piv != col {
            mat.swap(col, piv);
            rhs.swap(col, piv);
        }
        let diag = mat
            .get(col)
            .and_then(|r| r.get(col))
            .copied()
            .unwrap_or(Complex::new(1.0, 0.0));
        for row in (col + 1)..n {
            let factor = mat
                .get(row)
                .and_then(|r| r.get(col))
                .copied()
                .unwrap_or_default()
                / diag;
            for k in col..n {
                let p_val = mat
                    .get(col)
                    .and_then(|r| r.get(k))
                    .copied()
                    .unwrap_or_default();
                if let Some(cell) = mat.get_mut(row).and_then(|r| r.get_mut(k)) {
                    *cell -= factor * p_val;
                }
            }
            let rhs_col = *rhs.get(col).unwrap_or(&Complex::default());
            if let Some(r_slot) = rhs.get_mut(row) {
                *r_slot -= factor * rhs_col;
            }
        }
    }
    let mut x = vec![Complex::new(0.0, 0.0); n];
    for i in (0..n).rev() {
        let mut sum = *rhs.get(i).unwrap_or(&Complex::default());
        for j in (i + 1)..n {
            let a_ij = mat
                .get(i)
                .and_then(|r| r.get(j))
                .copied()
                .unwrap_or_default();
            let x_j = *x.get(j).unwrap_or(&Complex::default());
            sum -= a_ij * x_j;
        }
        let diag = mat
            .get(i)
            .and_then(|r| r.get(i))
            .copied()
            .unwrap_or(Complex::new(1.0, 0.0));
        if let Some(slot) = x.get_mut(i) {
            *slot = sum / diag;
        }
    }
    Ok(x)
}
