//! Semidefinite programming by the SDPAD alternating-direction method.
//!
//! The problem is `minimise c . x` subject to the linear matrix inequality
//! `F(x) = F0 + sum_i x_i F_i ⪰ 0` and optional linear constraints
//! `a_r . x (rel_r) b_r` (equalities, `<=`, and `>=`). Equality rows are
//! enforced directly in the quadratic `x`-update, while linear inequality rows
//! are written in non-negative slack form `G x + h = z_lin >= 0` and handled
//! alongside the PSD cone `F(x) = Z ⪰ 0` by the over-relaxed ADMM of Wen,
//! Goldfarb and Yin (SDPAD):
//!
//! ```text
//! x     <- argmin_x  c . x + <Y, F(x) - Z> + (rho/2) ||F(x) - Z||_F^2
//!                    + <y_lin, G x + h - z_lin> + (rho/2) ||G x + h - z_lin||_2^2
//!                    s.t. A_eq x = b_eq
//! Z     <- proj_{Z ⪰ 0} ( gamma F(x) + (1 - gamma) Z + Y / rho )
//! z_lin <- proj_{z >= 0} ( gamma (G x + h) + (1 - gamma) z_lin + y_lin / rho )
//! Y     <- Y + rho ( gamma F(x) + (1 - gamma) Z - Z )
//! y_lin <- y_lin + rho ( gamma (G x + h) + (1 - gamma) z_lin - z_lin )
//! ```
//!
//! with the over-relaxation factor `gamma = 1.6`, `rho = 1.0`, and an adaptive
//! `rho` adjustment every fifty iterations. The projection onto the positive
//! semidefinite cone uses a cyclic Jacobi symmetric eigensolver implemented in
//! this module, so no LAPACK binding is required.

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};

use super::lp::Relation;

/// A semidefinite program.
#[derive(Debug, Clone, PartialEq)]
pub struct SdpProblem {
    /// Objective coefficients, `minimise c . x` over `m` variables.
    pub c: Vec<f64>,
    /// Constant term `F0` of the linear matrix inequality.
    pub f0: Vec<Vec<f64>>,
    /// Matrices `F_i` of the linear matrix inequality (one per variable).
    pub fs: Vec<Vec<Vec<f64>>>,
    /// Optional linear constraint rows `a_r . x (rel_r) b_r`.
    pub a: Vec<Vec<f64>>,
    /// Sense of each linear constraint row (`Eq`, `Le`, or `Ge`).
    pub rel: Vec<Relation>,
    /// Right-hand side of each linear constraint row.
    pub b: Vec<f64>,
}

/// Result of [`solve_sdp`].
#[derive(Debug, Clone, PartialEq)]
pub struct SdpSolution {
    /// Primal point.
    pub x: Vec<f64>,
    /// Primal objective value `c . x`.
    pub objective: f64,
    /// Number of ADMM iterations.
    pub iterations: usize,
    /// Combined primal residual `sqrt(||F(x) - Z||_F^2 + ||G x + h - z_lin||_2^2)`.
    pub primal_residual: f64,
    /// Combined dual residual `rho sqrt(||Z - Z_prev||_F^2 + ||z_lin - z_lin_prev||_2^2)`.
    pub dual_residual: f64,
    /// Smallest eigenvalue of `F(x)` at the returned point.
    pub min_eigenvalue: f64,
}

/// Solves a semidefinite program by SDPAD.
///
/// # Errors
/// Returns an error string when the dimensions disagree, when an entry is not
/// finite, when `max_iter` is zero, or when the `x`-update system is singular.
pub fn solve_sdp(
    p: &SdpProblem,
    max_iter: usize,
) -> Result<SdpSolution, String> {
    let m = validate(p)?;
    if max_iter == 0 {
        return Err("iteration budget is zero".to_string());
    }
    let n = p.f0.len();

    // Partition linear rows into equalities `A_eq x = b_eq` and non-negative
    // linear cone constraints `G x + h_lin = z_lin >= 0`.
    let mut a_eq: Vec<Vec<f64>> = Vec::new();
    let mut b_eq: Vec<f64> = Vec::new();
    let mut g_lin: Vec<Vec<f64>> = Vec::new();
    let mut h_lin: Vec<f64> = Vec::new();
    for r in 0..p.a.len() {
        match p.rel[r] {
            | Relation::Eq => {
                a_eq.push(p.a[r].clone());
                b_eq.push(p.b[r]);
            },
            | Relation::Le => {
                // a_r . x <= b_r  <=>  (-a_r) . x + b_r >= 0
                g_lin.push(p.a[r].iter().map(|v| -v).collect());
                h_lin.push(p.b[r]);
            },
            | Relation::Ge => {
                // a_r . x >= b_r  <=>  a_r . x - b_r >= 0
                g_lin.push(p.a[r].clone());
                h_lin.push(-p.b[r]);
            },
        }
    }
    let me = a_eq.len();
    let mi = g_lin.len();

    // H_ij = <F_i, F_j> + sum_k g_ki g_kj.
    let mut h = vec![vec![0.0; m]; m];
    for i in 0..m {
        for j in i..m {
            let mut v = inner(&p.fs[i], &p.fs[j]);
            for k in 0..mi {
                v += g_lin[k][i] * g_lin[k][j];
            }
            h[i][j] = v;
            h[j][i] = v;
        }
    }

    let mut x = vec![0.0; m];
    let mut z = vec![vec![0.0; n]; n];
    let mut y = vec![vec![0.0; n]; n];
    let mut z_lin = vec![0.0; mi];
    let mut y_lin = vec![0.0; mi];

    let mut rho = 1.0f64;
    let gamma = 1.6f64;
    let mut iterations = 0usize;
    let mut primal_residual = f64::INFINITY;
    let mut dual_residual = f64::INFINITY;

    for it in 0..max_iter {
        // x-update: minimise the augmented Lagrangian over the affine set A_eq x = b_eq.
        let r = sub(&p.f0, &z);
        let mut d = vec![0.0; m];
        for i in 0..m {
            let mut acc = p.c[i] + inner(&y, &p.fs[i]) + rho * inner(&r, &p.fs[i]);
            for k in 0..mi {
                acc += g_lin[k][i] * (y_lin[k] + rho * (h_lin[k] - z_lin[k]));
            }
            d[i] = acc;
        }
        x = solve_x_update(&h, &a_eq, &d, &b_eq, rho, m, me)?;

        // Z-update with over-relaxation for the PSD cone.
        let fx = linear_map(&p.f0, &p.fs, &x);
        let mut fhat = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                fhat[i][j] = gamma * fx[i][j] + (1.0 - gamma) * z[i][j];
            }
        }
        let mut input = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in 0..n {
                input[i][j] = fhat[i][j] + y[i][j] / rho;
            }
        }
        let z_new = project_psd(&input);

        // z_lin-update with over-relaxation for the linear non-negative cone.
        let mut sl = vec![0.0; mi];
        let mut shat = vec![0.0; mi];
        let mut z_lin_new = vec![0.0; mi];
        for k in 0..mi {
            sl[k] = dot(&g_lin[k], &x) + h_lin[k];
            shat[k] = gamma * sl[k] + (1.0 - gamma) * z_lin[k];
            z_lin_new[k] = (shat[k] + y_lin[k] / rho).max(0.0);
        }

        // Y-update (PSD and linear duals).
        for i in 0..n {
            for j in 0..n {
                y[i][j] += rho * (fhat[i][j] - z_new[i][j]);
            }
        }
        for k in 0..mi {
            y_lin[k] += rho * (shat[k] - z_lin_new[k]);
        }

        // Combined primal and dual residuals.
        let mat_p = sub(&fx, &z_new);
        let mat_d = sub(&z_new, &z);
        let mut lin_p_sq = 0.0;
        let mut lin_d_sq = 0.0;
        for k in 0..mi {
            let dp = sl[k] - z_lin_new[k];
            let dd = z_lin_new[k] - z_lin[k];
            lin_p_sq += dp * dp;
            lin_d_sq += dd * dd;
        }
        primal_residual = (inner(&mat_p, &mat_p) + lin_p_sq).max(0.0).sqrt();
        dual_residual = rho * (inner(&mat_d, &mat_d) + lin_d_sq).max(0.0).sqrt();
        z = z_new;
        z_lin = z_lin_new;
        iterations = it + 1;

        // Adaptive rho every fifty iterations.
        if iterations.is_multiple_of(50) {
            if primal_residual > 10.0 * dual_residual {
                rho *= 2.0;
                for row in &mut y {
                    for v in row.iter_mut() {
                        *v *= 0.5;
                    }
                }
                for v in &mut y_lin {
                    *v *= 0.5;
                }
            } else if dual_residual > 10.0 * primal_residual {
                rho *= 0.5;
                for row in &mut y {
                    for v in row.iter_mut() {
                        *v *= 2.0;
                    }
                }
                for v in &mut y_lin {
                    *v *= 2.0;
                }
            }
        }

        let scale_p = 1.0 + (inner(&z, &z) + dot(&z_lin, &z_lin)).max(0.0).sqrt();
        let scale_d = 1.0 + (inner(&y, &y) + dot(&y_lin, &y_lin)).max(0.0).sqrt();
        if primal_residual < 1e-8 * scale_p && dual_residual < 1e-8 * scale_d {
            break;
        }
    }

    let fx = linear_map(&p.f0, &p.fs, &x);
    let min_eigenvalue = smallest_eigenvalue(&fx);

    Ok(SdpSolution {
        objective: dot(&p.c, &x),
        x,
        iterations,
        primal_residual,
        dual_residual,
        min_eigenvalue,
    })
}

/// Solves the quadratic x-update with the equality rows `A x = b`.
fn solve_x_update(
    h: &[Vec<f64>],
    a: &[Vec<f64>],
    d: &[f64],
    b: &[f64],
    rho: f64,
    m: usize,
    me: usize,
) -> Result<Vec<f64>, String> {
    let size = m + me;
    let mut k = vec![vec![0.0; size]; size];
    for i in 0..m {
        for j in 0..m {
            k[i][j] = rho * h[i][j];
        }
        k[i][i] += 1e-10;
    }
    for r in 0..me {
        for i in 0..m {
            k[i][m + r] = a[r][i];
            k[m + r][i] = a[r][i];
        }
    }
    let mut rhs = vec![0.0; size];
    for i in 0..m {
        rhs[i] = -d[i];
    }
    rhs[m..(me + m)].copy_from_slice(&b[..me]);
    let mat = Mat::from_rows(&k).map_err(|e| e.to_string())?;
    match dense::solve(&mat, &rhs) {
        | Ok(sol) => Ok(sol[0..m].to_vec()),
        | Err(_) => {
            for i in 0..m {
                k[i][i] += 1e-6;
            }
            let mat2 = Mat::from_rows(&k).map_err(|e| e.to_string())?;
            dense::solve(&mat2, &rhs)
                .map(|sol| sol[0..m].to_vec())
                .map_err(|_| "singular sdp x-update".to_string())
        },
    }
}

/// `F0 + sum_i x_i F_i`.
fn linear_map(
    f0: &[Vec<f64>],
    fs: &[Vec<Vec<f64>>],
    x: &[f64],
) -> Vec<Vec<f64>> {
    let n = f0.len();
    let mut out = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut v = f0[i][j];
            for (k, xk) in x.iter().enumerate() {
                if *xk != 0.0 {
                    v += xk * fs[k][i][j];
                }
            }
            out[i][j] = v;
        }
    }
    out
}

/// Projection of a symmetric matrix onto the positive semidefinite cone.
fn project_psd(s: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = s.len();
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            a[i][j] = f64::midpoint(s[i][j], s[j][i]);
        }
    }
    let (vals, vecs) = jacobi_eigen(&mut a);
    let mut z = vec![vec![0.0; n]; n];
    for k in 0..n {
        let lam = vals[k].max(0.0);
        if lam == 0.0 {
            continue;
        }
        for i in 0..n {
            for j in 0..n {
                z[i][j] += lam * vecs[i][k] * vecs[j][k];
            }
        }
    }
    z
}

/// Smallest eigenvalue of a symmetric matrix.
fn smallest_eigenvalue(s: &[Vec<f64>]) -> f64 {
    let n = s.len();
    if n == 0 {
        return 0.0;
    }
    let mut a = s.to_vec();
    let (vals, _) = jacobi_eigen(&mut a);
    vals.iter().copied().fold(f64::INFINITY, f64::min)
}

/// Cyclic Jacobi eigensolver for a symmetric matrix.
///
/// Returns the eigenvalues and the matrix whose columns are the corresponding
/// eigenvectors (eigenvalues are not sorted).
fn jacobi_eigen(a: &mut [Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..100 {
        let mut off = 0.0;
        for p in 0..n {
            for q in p + 1..n {
                off += a[p][q] * a[p][q];
            }
        }
        if off.sqrt() <= 1e-15 {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p][q];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + theta.mul_add(theta, 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / t.mul_add(t, 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (kp, kq) = (a[k][p], a[k][q]);
                    a[k][p] = c * kp - s * kq;
                    a[k][q] = s * kp + c * kq;
                }
                for k in 0..n {
                    let (pk, qk) = (a[p][k], a[q][k]);
                    a[p][k] = c * pk - s * qk;
                    a[q][k] = s * pk + c * qk;
                }
                for k in 0..n {
                    let (kp, kq) = (v[k][p], v[k][q]);
                    v[k][p] = c * kp - s * kq;
                    v[k][q] = s * kp + c * kq;
                }
            }
        }
    }
    let vals = (0..n).map(|i| a[i][i]).collect();
    (vals, v)
}

/// Frobenius inner product `<A, B>`.
fn inner(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> f64 {
    let mut s = 0.0;
    for i in 0..a.len() {
        for j in 0..a[i].len() {
            s += a[i][j] * b[i][j];
        }
    }
    s
}

fn sub(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut out = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            out[i][j] = a[i][j] - b[i][j];
        }
    }
    out
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Validates shapes and finiteness, returning the number of variables.
fn validate(p: &SdpProblem) -> Result<usize, String> {
    let m = p.c.len();
    let n = p.f0.len();
    if p.f0.iter().any(|row| row.len() != n) {
        return Err("f0 must be a square matrix".to_string());
    }
    if p.fs.len() != m {
        return Err("one matrix per objective variable is required".to_string());
    }
    for mat in &p.fs {
        if mat.len() != n || mat.iter().any(|row| row.len() != n) {
            return Err("every f_i must match the size of f0".to_string());
        }
        if mat.iter().any(|row| row.iter().any(|v| !v.is_finite())) {
            return Err("f_i must be finite".to_string());
        }
    }
    if p.f0.iter().any(|row| row.iter().any(|v| !v.is_finite())) {
        return Err("f0 must be finite".to_string());
    }
    if p.c.iter().any(|v| !v.is_finite()) {
        return Err("objective coefficients must be finite".to_string());
    }
    if p.rel.len() != p.a.len() || p.b.len() != p.a.len() {
        return Err("linear rows, senses and right-hand sides must agree".to_string());
    }
    for row in &p.a {
        if row.len() != m {
            return Err("linear rows must have one entry per variable".to_string());
        }
        if row.iter().any(|v| !v.is_finite()) {
            return Err("linear coefficients must be finite".to_string());
        }
    }
    if p.b.iter().any(|v| !v.is_finite()) {
        return Err("right-hand sides must be finite".to_string());
    }
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scalar_psd_constraint() {
        // min x  s.t.  [[x, 0], [0, 1]] ⪰ 0  ->  x = 0.
        let p = SdpProblem {
            c: vec![1.0],
            f0: vec![vec![0.0, 0.0], vec![0.0, 1.0]],
            fs: vec![vec![vec![1.0, 0.0], vec![0.0, 0.0]]],
            a: Vec::new(),
            rel: Vec::new(),
            b: Vec::new(),
        };
        let s = solve_sdp(&p, 20000).unwrap();
        assert!(s.x[0].abs() < 1e-4, "x = {}", s.x[0]);
        assert!(s.objective.abs() < 1e-4);
        assert!(s.min_eigenvalue > -1e-6);
    }

    #[test]
    fn two_by_two_lmi() {
        // min x1 + x2  s.t.  [[x1, 1], [1, x2]] ⪰ 0  ->  (1, 1), objective 2.
        let p = SdpProblem {
            c: vec![1.0, 1.0],
            f0: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            fs: vec![
                vec![vec![1.0, 0.0], vec![0.0, 0.0]],
                vec![vec![0.0, 0.0], vec![0.0, 1.0]],
            ],
            a: Vec::new(),
            rel: Vec::new(),
            b: Vec::new(),
        };
        let s = solve_sdp(&p, 20000).unwrap();
        assert!((s.x[0] - 1.0).abs() < 1e-3, "x0 = {}", s.x[0]);
        assert!((s.x[1] - 1.0).abs() < 1e-3, "x1 = {}", s.x[1]);
        assert!((s.objective - 2.0).abs() < 1e-3);
    }

    #[test]
    fn with_equality_row() {
        // min x1 + 2 x2  s.t.  diag(x1, x2) ⪰ 0,  x1 + x2 = 3
        // -> x = (3, 0), objective 3.
        let p = SdpProblem {
            c: vec![1.0, 2.0],
            f0: vec![vec![0.0, 0.0], vec![0.0, 0.0]],
            fs: vec![
                vec![vec![1.0, 0.0], vec![0.0, 0.0]],
                vec![vec![0.0, 0.0], vec![0.0, 1.0]],
            ],
            a: vec![vec![1.0, 1.0]],
            rel: vec![Relation::Eq],
            b: vec![3.0],
        };
        let s = solve_sdp(&p, 20000).unwrap();
        assert!((s.x[0] - 3.0).abs() < 1e-3, "x0 = {}", s.x[0]);
        assert!(s.x[1].abs() < 1e-3, "x1 = {}", s.x[1]);
        assert!((s.objective - 3.0).abs() < 1e-3);
    }

    #[test]
    fn with_le_inequality_row() {
        // min x1 + x2  s.t.  [[x1, 1], [1, x2]] ⪰ 0,  x1 <= 0.5
        // Since x1 * x2 >= 1 and x1 <= 0.5, opt is (0.5, 2.0), objective 2.5.
        let p = SdpProblem {
            c: vec![1.0, 1.0],
            f0: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            fs: vec![
                vec![vec![1.0, 0.0], vec![0.0, 0.0]],
                vec![vec![0.0, 0.0], vec![0.0, 1.0]],
            ],
            a: vec![vec![1.0, 0.0]],
            rel: vec![Relation::Le],
            b: vec![0.5],
        };
        let s = solve_sdp(&p, 20000).unwrap();
        assert!((s.x[0] - 0.5).abs() < 1e-3, "x0 = {}", s.x[0]);
        assert!((s.x[1] - 2.0).abs() < 1e-3, "x1 = {}", s.x[1]);
        assert!((s.objective - 2.5).abs() < 1e-3, "obj = {}", s.objective);
        assert!(s.min_eigenvalue > -1e-5, "min_eig = {}", s.min_eigenvalue);
        assert!(s.x[0] <= 0.5 + 1e-4, "x0 = {}", s.x[0]);
    }

    #[test]
    fn with_ge_inequality_row() {
        // min x1 + x2  s.t.  [[x1, 1], [1, x2]] ⪰ 0,  x1 >= 2.0
        // Since x1 * x2 >= 1 and x1 >= 2.0, opt is (2.0, 0.5), objective 2.5.
        let p = SdpProblem {
            c: vec![1.0, 1.0],
            f0: vec![vec![0.0, 1.0], vec![1.0, 0.0]],
            fs: vec![
                vec![vec![1.0, 0.0], vec![0.0, 0.0]],
                vec![vec![0.0, 0.0], vec![0.0, 1.0]],
            ],
            a: vec![vec![1.0, 0.0]],
            rel: vec![Relation::Ge],
            b: vec![2.0],
        };
        let s = solve_sdp(&p, 20000).unwrap();
        assert!((s.x[0] - 2.0).abs() < 1e-3, "x0 = {}", s.x[0]);
        assert!((s.x[1] - 0.5).abs() < 1e-3, "x1 = {}", s.x[1]);
        assert!((s.objective - 2.5).abs() < 1e-3, "obj = {}", s.objective);
        assert!(s.min_eigenvalue > -1e-5, "min_eig = {}", s.min_eigenvalue);
        assert!(s.x[0] >= 2.0 - 1e-4, "x0 = {}", s.x[0]);
    }

    #[test]
    fn with_mixed_eq_le_ge_rows() {
        // min x1 + 2 x2 + 3 x3
        // s.t. [[x1, 1, 0], [1, x2, 0], [0, 0, x3]] ⪰ 0,
        //      x1 + x2 + x3 = 5,  x1 <= 2,  x3 >= 1.
        // Opt: (2, 2, 1), objective 9.
        let p = SdpProblem {
            c: vec![1.0, 2.0, 3.0],
            f0: vec![
                vec![0.0, 1.0, 0.0],
                vec![1.0, 0.0, 0.0],
                vec![0.0, 0.0, 0.0],
            ],
            fs: vec![
                vec![
                    vec![1.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0],
                ],
                vec![
                    vec![0.0, 0.0, 0.0],
                    vec![0.0, 1.0, 0.0],
                    vec![0.0, 0.0, 0.0],
                ],
                vec![
                    vec![0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 0.0],
                    vec![0.0, 0.0, 1.0],
                ],
            ],
            a: vec![
                vec![1.0, 1.0, 1.0],
                vec![1.0, 0.0, 0.0],
                vec![0.0, 0.0, 1.0],
            ],
            rel: vec![Relation::Eq, Relation::Le, Relation::Ge],
            b: vec![5.0, 2.0, 1.0],
        };
        let s = solve_sdp(&p, 20000).unwrap();
        assert!((s.x[0] - 2.0).abs() < 1e-3, "x0 = {}", s.x[0]);
        assert!((s.x[1] - 2.0).abs() < 1e-3, "x1 = {}", s.x[1]);
        assert!((s.x[2] - 1.0).abs() < 1e-3, "x2 = {}", s.x[2]);
        assert!((s.objective - 9.0).abs() < 1e-3, "obj = {}", s.objective);
        assert!(s.min_eigenvalue > -1e-5, "min_eig = {}", s.min_eigenvalue);
        assert!((s.x[0] + s.x[1] + s.x[2] - 5.0).abs() < 1e-4);
        assert!(s.x[0] <= 2.0 + 1e-4);
        assert!(s.x[2] >= 1.0 - 1e-4);
    }
}
