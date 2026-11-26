//! Sequential quadratic programming with a sparse QP subproblem.
//!
//! The dense reference solver [`crate::kernels::optim::sqp`] forms and
//! factorises a dense `n x n` matrix for every QP subproblem, costing `O(n^3)`
//! per iteration. This module is the sparse alternative for problems whose
//! Lagrangian Hessian and constraint Jacobian are sparse: each subproblem
//!
//! ```text
//! minimise  0.5 p'Hp + g'p
//! subject to  J_eq p = -c_eq,   J_ineq p <= -c_ineq
//! ```
//!
//! is solved by an active-set method whose equality-constrained core is a
//! **projected, diagonally preconditioned conjugate-gradient iteration on the
//! working set**. Every CG step uses only sparse matrix–vector products
//! (`O(nnz)` work), so the subproblem costs `O(nnz * iterations)` and never
//! builds a dense `n x n` factorisation. The only dense factorisation is the
//! small `m x m` system `A_W A_W'`, where `m` is the number of active rows
//! (equalities plus active inequalities). [`crate::kernels::optim::sqp`] remains
//! the reference implementation for small dense problems.

use crate::kernels::dense::Lu;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};

/// A symmetric sparse matrix stored by upper-triangular triplets.
///
/// Duplicate entries are summed by the matrix–vector product.
#[derive(Debug, Clone, PartialEq)]
pub struct SparseSym {
    /// Number of rows (equal to the number of columns).
    pub n: usize,
    /// Row index of every entry.
    pub rows: Vec<usize>,
    /// Column index of every entry.
    pub cols: Vec<usize>,
    /// Value of every entry.
    pub vals: Vec<f64>,
}

/// A general sparse matrix stored by triplets.
#[derive(Debug, Clone, PartialEq)]
pub struct SparseMat {
    /// Number of rows.
    pub rows_n: usize,
    /// Number of columns.
    pub cols_n: usize,
    /// Row index of every entry.
    pub rows: Vec<usize>,
    /// Column index of every entry.
    pub cols: Vec<usize>,
    /// Value of every entry.
    pub vals: Vec<f64>,
}

/// The nonlinear program seen by [`solve_sqp`].
///
/// Constraints are returned with the equality rows first and the inequality
/// rows (`g(x) <= 0`) after them.
pub trait NlpOracle {
    /// Objective value `f(x)`.
    fn objective(
        &self,
        x: &[f64],
    ) -> f64;
    /// Gradient `grad f(x)`.
    fn gradient(
        &self,
        x: &[f64],
    ) -> Vec<f64>;
    /// Hessian of the Lagrangian `grad^2 f + sum lambda_i grad^2 c_i`.
    fn hessian(
        &self,
        x: &[f64],
        mult_eq: &[f64],
        mult_ineq: &[f64],
    ) -> SparseSym;
    /// Constraint values, equalities first, then inequalities.
    fn constraints(
        &self,
        x: &[f64],
    ) -> Vec<f64>;
    /// Jacobian of the constraints, rows ordered as in [`NlpOracle::constraints`].
    fn jacobian(
        &self,
        x: &[f64],
    ) -> SparseMat;
    /// Number of equality constraints.
    fn n_eq(&self) -> usize;
    /// Number of inequality constraints.
    fn n_ineq(&self) -> usize;
}

/// Options of [`solve_sqp`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SqpOptions {
    /// Maximum number of SQP iterations.
    pub max_iter: usize,
    /// KKT tolerance for convergence.
    pub tol: f64,
    /// Maximum step norm per iteration (trust radius).
    pub trust_radius: f64,
}

/// Result of [`solve_sqp`].
#[derive(Debug, Clone, PartialEq)]
pub struct SqpResult {
    /// Final iterate.
    pub x: Vec<f64>,
    /// Objective value at `x`.
    pub objective: f64,
    /// Equality multipliers.
    pub mult_eq: Vec<f64>,
    /// Inequality multipliers.
    pub mult_ineq: Vec<f64>,
    /// Number of SQP iterations.
    pub iterations: usize,
    /// Final KKT residual (stationarity plus constraint violation).
    pub kkt_residual: f64,
}

/// Solves a sparse nonlinear program by sequential quadratic programming.
///
/// # Errors
/// Returns an error string when the dimensions reported by the oracle are
/// inconsistent or a working-set system is singular.
pub fn solve_sqp(
    oracle: &dyn NlpOracle,
    x0: &[f64],
    opts: &SqpOptions,
) -> Result<SqpResult, String> {
    let n = x0.len();
    let n_eq = oracle.n_eq();
    let n_ineq = oracle.n_ineq();
    let n_con = n_eq + n_ineq;
    let tol = if opts.tol > 0.0 {
        opts.tol
    } else {
        1e-8
    };
    let trust = if opts.trust_radius > 0.0 {
        opts.trust_radius
    } else {
        1e3
    };

    let mut x = x0.to_vec();
    let mut mult_eq = vec![0.0; n_eq];
    let mut mult_ineq = vec![0.0; n_ineq];
    let mut iterations = 0usize;
    let mut kkt_residual = kkt(oracle, &x, &mult_eq, &mult_ineq)?;

    for _ in 0..opts.max_iter {
        if kkt_residual < tol {
            break;
        }
        let g = oracle.gradient(&x);
        if g.len() != n {
            return Err("gradient length does not match the iterate".to_string());
        }
        let cvals = oracle.constraints(&x);
        if cvals.len() != n_con {
            return Err("constraint length does not match the declared counts".to_string());
        }
        let jac = oracle.jacobian(&x);
        let rows = build_rows(&jac, n_con, n)?;
        let h = oracle.hessian(&x, &mult_eq, &mult_ineq);
        if h.n != n {
            return Err("hessian size does not match the iterate".to_string());
        }

        let (mut p, lam_eq, lam_ineq) =
            solve_qp_subproblem(&h, &g, &rows, &cvals, n_eq, n_ineq, n)?;

        // Trust region.
        let pn = norm(&p);
        if pn > trust {
            let s = trust / pn;
            for v in &mut p {
                *v *= s;
            }
        }

        // Backtracking line search on the L1 merit function.
        let pen = 1.0 + max_multiplier(&lam_eq, &lam_ineq);
        let merit0 = merit(oracle, &x, pen);
        let mut alpha = 1.0;
        let mut trial = x.clone();
        let mut accepted = false;
        for _ in 0..30 {
            trial = x.iter().zip(&p).map(|(xi, pi)| xi + alpha * pi).collect();
            let mt = merit(oracle, &trial, pen);
            if mt <= merit0 - 1e-4 * alpha * merit0.abs().max(1e-12) {
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            break;
        }

        x = trial;
        mult_eq = lam_eq;
        mult_ineq = lam_ineq;
        iterations += 1;
        kkt_residual = kkt(oracle, &x, &mult_eq, &mult_ineq)?;
    }

    Ok(SqpResult {
        objective: oracle.objective(&x),
        x,
        mult_eq,
        mult_ineq,
        iterations,
        kkt_residual,
    })
}

/// Solves the sparse QP subproblem by an active-set projected-CG method.
fn solve_qp_subproblem(
    h: &SparseSym,
    g: &[f64],
    rows: &[Vec<(usize, f64)>],
    cvals: &[f64],
    n_eq: usize,
    n_ineq: usize,
    n: usize,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), String> {
    let n_con = n_eq + n_ineq;
    let mut active: Vec<usize> = (0..n_eq).collect();
    let tol_m = 1e-9;
    let mut p = vec![0.0; n];
    let mut lambda: Vec<f64> = Vec::new();

    for _ in 0..(n_con + n + 10) {
        let (pp, lam) = solve_eq_qp(h, g, rows, &active, cvals, n)?;
        p = pp;
        lambda = lam;

        // Drop the most negative inequality multiplier.
        let mut worst: Option<(usize, f64)> = None;
        for (pos, &ai) in active.iter().enumerate() {
            if ai >= n_eq {
                let val = lambda.get(pos).copied().unwrap_or(0.0);
                if val < -tol_m && worst.is_none_or(|(_, v)| val < v) {
                    worst = Some((pos, val));
                }
            }
        }
        if let Some((pos, _)) = worst {
            active.remove(pos);
            continue;
        }

        // Add the most violated inactive inequality.
        let mut add: Option<(usize, f64)> = None;
        for i in n_eq..n_con {
            if active.contains(&i) {
                continue;
            }
            let mut val = cvals[i];
            for &(c, v) in &rows[i] {
                val += v * p[c];
            }
            if val > 1e-9 && add.is_none_or(|(_, vv)| val > vv) {
                add = Some((i, val));
            }
        }
        if let Some((i, _)) = add {
            active.push(i);
            continue;
        }
        break;
    }

    let mut mult_eq = vec![0.0; n_eq];
    let mut mult_ineq = vec![0.0; n_ineq];
    for (pos, &ai) in active.iter().enumerate() {
        let val = lambda.get(pos).copied().unwrap_or(0.0);
        if ai < n_eq {
            mult_eq[ai] = val;
        } else {
            mult_ineq[ai - n_eq] = val;
        }
    }
    Ok((p, mult_eq, mult_ineq))
}

/// Equality-constrained QP `min 0.5 p'Hp + g'p s.t. A_W p = -c_W` by
/// projected preconditioned CG.
#[allow(unused_assignments)]
fn solve_eq_qp(
    h: &SparseSym,
    g: &[f64],
    rows: &[Vec<(usize, f64)>],
    active: &[usize],
    cvals: &[f64],
    n: usize,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let k = active.len();

    // Diagonal of H for preconditioning, plus a tiny shift for definiteness.
    let mut hdiag = vec![0.0; n];
    let mut dmax = 0.0f64;
    for t in 0..h.rows.len() {
        if h.rows[t] == h.cols[t] && h.rows[t] < n {
            hdiag[h.rows[t]] += h.vals[t];
        }
    }
    for d in &hdiag {
        dmax = dmax.max(d.abs());
    }
    let delta = 1e-10 * (1.0 + dmax);

    let hvec = |d: &[f64]| -> Vec<f64> {
        let mut y = sym_matvec(h, d);
        for i in 0..n {
            y[i] += delta * d[i];
        }
        y
    };
    let diag_floor = 1e-6 * (1.0 + dmax);
    let prec = |r: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|i| {
                // Jacobi scaling where the diagonal is meaningful, identity
                // where it vanishes (so a singular direction is never
                // amplified).
                if hdiag[i] > diag_floor {
                    r[i] / hdiag[i]
                } else {
                    r[i]
                }
            })
            .collect()
    };

    let lu_opt: Option<Lu>;
    let p0: Vec<f64>;
    if k == 0 {
        lu_opt = None;
        p0 = vec![0.0; n];
    } else {
        let mut aa = vec![vec![0.0; k]; k];
        for i in 0..k {
            for j in 0..k {
                let mut s = 0.0;
                for &(c1, v1) in &rows[active[i]] {
                    for &(c2, v2) in &rows[active[j]] {
                        if c1 == c2 {
                            s += v1 * v2;
                        }
                    }
                }
                aa[i][j] = s;
            }
        }
        let lu = factor_small(&aa)?;
        let cw: Vec<f64> = active.iter().map(|&i| cvals[i]).collect();
        let w = lu.solve(&cw);
        let mut pp = vec![0.0; n];
        for (idx, &ai) in active.iter().enumerate() {
            let wk = w.get(idx).copied().unwrap_or(0.0);
            for &(c, v) in &rows[ai] {
                pp[c] -= v * wk;
            }
        }
        p0 = pp;
        lu_opt = Some(lu);
    }

    let proj = |u: &[f64]| -> Vec<f64> {
        match &lu_opt {
            | Some(lu) => project(u, active, rows, lu),
            | None => u.to_vec(),
        }
    };

    let mut p = p0;
    let hp = hvec(&p);
    let grad: Vec<f64> = (0..n).map(|i| hp[i] + g[i]).collect();
    let mut r = proj(&grad);
    let mut z = proj(&prec(&r));
    let mut d: Vec<f64> = z.iter().map(|v| -v).collect();
    let mut rz = dot(&r, &z);
    let r0 = norm(&r).max(1.0);
    let max_cg = 5 * n + 100;

    for _ in 0..max_cg {
        if rz.abs() < 1e-30 {
            break;
        }
        let hd = hvec(&d);
        let mut dhd = dot(&d, &hd);
        if dhd <= 1e-16 {
            dhd = 1e-8 * dot(&d, &d);
            if dhd <= 1e-30 {
                break;
            }
        }
        let alpha = rz / dhd;
        for i in 0..n {
            p[i] += alpha * d[i];
        }
        let hp = hvec(&p);
        let grad2: Vec<f64> = (0..n).map(|i| hp[i] + g[i]).collect();
        let r_new = proj(&grad2);
        if norm(&r_new) <= 1e-10 * r0 {
            r = r_new;
            break;
        }
        let z_new = proj(&prec(&r_new));
        let rz_new = dot(&r_new, &z_new);
        let beta = rz_new / rz;
        for i in 0..n {
            d[i] = -z_new[i] + beta * d[i];
        }
        r = r_new;
        z = z_new;
        rz = rz_new;
    }

    // Multipliers from the final stationarity residual.
    let lambda = if k == 0 {
        Vec::new()
    } else {
        let hp = hvec(&p);
        let gg: Vec<f64> = (0..n).map(|i| hp[i] + g[i]).collect();
        let mut au = vec![0.0; k];
        for (idx, &ai) in active.iter().enumerate() {
            let mut s = 0.0;
            for &(c, v) in &rows[ai] {
                s += v * gg[c];
            }
            au[idx] = s;
        }
        match &lu_opt {
            | Some(lu) => lu.solve(&au).iter().map(|v| -v).collect(),
            | None => Vec::new(),
        }
    };

    Ok((p, lambda))
}

/// Orthogonal projection onto `null(A_W)`.
fn project(
    u: &[f64],
    active: &[usize],
    rows: &[Vec<(usize, f64)>],
    lu: &Lu,
) -> Vec<f64> {
    let k = active.len();
    let mut au = vec![0.0; k];
    for (idx, &ai) in active.iter().enumerate() {
        let mut s = 0.0;
        for &(c, v) in &rows[ai] {
            s += v * u[c];
        }
        au[idx] = s;
    }
    let w = lu.solve(&au);
    let mut out = u.to_vec();
    for (idx, &ai) in active.iter().enumerate() {
        let wk = w.get(idx).copied().unwrap_or(0.0);
        if wk != 0.0 {
            for &(c, v) in &rows[ai] {
                out[c] -= v * wk;
            }
        }
    }
    out
}

/// Factors a small dense matrix, retrying with a diagonal shift.
fn factor_small(a: &[Vec<f64>]) -> Result<Lu, String> {
    let mat = Mat::from_rows(a).map_err(|e| e.to_string())?;
    match dense::lu_factor(&mat) {
        | Ok(lu) => Ok(lu),
        | Err(_) => {
            let mut ar = a.to_vec();
            let mut scale = 1.0f64;
            for (i, row) in ar.iter().enumerate() {
                scale = scale.max(row[i].abs());
            }
            for (i, row) in ar.iter_mut().enumerate() {
                row[i] += 1e-10 * scale;
            }
            let mat2 = Mat::from_rows(&ar).map_err(|e| e.to_string())?;
            dense::lu_factor(&mat2).map_err(|_| "singular qp working set".to_string())
        },
    }
}

/// Sparse symmetric matrix–vector product, summing duplicate entries.
fn sym_matvec(
    h: &SparseSym,
    x: &[f64],
) -> Vec<f64> {
    let mut y = vec![0.0; h.n];
    for t in 0..h.rows.len() {
        let i = h.rows[t];
        let j = h.cols[t];
        let v = h.vals[t];
        if i >= h.n || j >= h.n {
            continue;
        }
        if i == j {
            y[i] += v * x[i];
        } else {
            y[i] += v * x[j];
            y[j] += v * x[i];
        }
    }
    y
}

/// Splits a sparse Jacobian into per-row entry lists.
fn build_rows(
    jac: &SparseMat,
    n_con: usize,
    n: usize,
) -> Result<Vec<Vec<(usize, f64)>>, String> {
    if jac.rows_n != n_con || jac.cols_n != n {
        return Err("jacobian has the wrong shape".to_string());
    }
    let mut rows = vec![Vec::new(); n_con];
    for t in 0..jac.rows.len() {
        let r = jac.rows[t];
        let c = jac.cols[t];
        if r >= n_con || c >= n {
            return Err("jacobian index out of range".to_string());
        }
        rows[r].push((c, jac.vals[t]));
    }
    Ok(rows)
}

/// L1 merit function `f + pen * (sum |eq| + sum max(0, ineq))`.
fn merit(
    oracle: &dyn NlpOracle,
    x: &[f64],
    pen: f64,
) -> f64 {
    let f = oracle.objective(x);
    let c = oracle.constraints(x);
    let neq = oracle.n_eq();
    let mut v = 0.0;
    for (i, ci) in c.iter().enumerate() {
        if i < neq {
            v += ci.abs();
        } else {
            v += ci.max(0.0);
        }
    }
    f + pen * v
}

fn max_multiplier(
    eq: &[f64],
    ineq: &[f64],
) -> f64 {
    let mut m = 0.0f64;
    for v in eq.iter().chain(ineq.iter()) {
        m = m.max(v.abs());
    }
    m
}

/// KKT residual: stationarity inf-norm plus constraint violation.
fn kkt(
    oracle: &dyn NlpOracle,
    x: &[f64],
    mult_eq: &[f64],
    mult_ineq: &[f64],
) -> Result<f64, String> {
    let g = oracle.gradient(x);
    let c = oracle.constraints(x);
    let jac = oracle.jacobian(x);
    let n_eq = oracle.n_eq();
    let n_ineq = oracle.n_ineq();
    let n_con = n_eq + n_ineq;
    if g.len() != x.len() || c.len() != n_con || jac.rows_n != n_con || jac.cols_n != x.len() {
        return Err("oracle returned inconsistent dimensions".to_string());
    }
    let mut mu = vec![0.0; n_con];
    for i in 0..n_eq {
        mu[i] = mult_eq.get(i).copied().unwrap_or(0.0);
    }
    for i in 0..n_ineq {
        mu[n_eq + i] = mult_ineq.get(i).copied().unwrap_or(0.0);
    }
    let mut jt = vec![0.0; x.len()];
    for t in 0..jac.rows.len() {
        let r = jac.rows[t];
        let cc = jac.cols[t];
        if r < n_con && cc < x.len() {
            jt[cc] += jac.vals[t] * mu[r];
        }
    }
    let mut stat = 0.0f64;
    for i in 0..x.len() {
        stat = stat.max((g[i] + jt[i]).abs());
    }
    let mut viol = 0.0f64;
    for i in 0..n_con {
        if i < n_eq {
            viol = viol.max(c[i].abs());
        } else {
            viol = viol.max(c[i].max(0.0));
        }
    }
    Ok(stat + viol)
}

fn norm(v: &[f64]) -> f64 {
    crate::kernels::simd::norm2_sq_f64(v).sqrt()
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    crate::kernels::simd::dot_f64(a, b)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hock–Schittkowski 1 style unconstrained problem: `(x0-1)^2 + x1^2`.
    struct Hs1;

    impl NlpOracle for Hs1 {
        fn objective(
            &self,
            x: &[f64],
        ) -> f64 {
            (x[0] - 1.0).powi(2) + x[1] * x[1]
        }

        fn gradient(
            &self,
            x: &[f64],
        ) -> Vec<f64> {
            vec![2.0 * (x[0] - 1.0), 2.0 * x[1]]
        }

        fn hessian(
            &self,
            _x: &[f64],
            _me: &[f64],
            _mi: &[f64],
        ) -> SparseSym {
            SparseSym {
                n: 2,
                rows: vec![0, 1],
                cols: vec![0, 1],
                vals: vec![2.0, 2.0],
            }
        }

        fn constraints(
            &self,
            _x: &[f64],
        ) -> Vec<f64> {
            Vec::new()
        }

        fn jacobian(
            &self,
            _x: &[f64],
        ) -> SparseMat {
            SparseMat {
                rows_n: 0,
                cols_n: 2,
                rows: Vec::new(),
                cols: Vec::new(),
                vals: Vec::new(),
            }
        }

        fn n_eq(&self) -> usize {
            0
        }

        fn n_ineq(&self) -> usize {
            0
        }
    }

    /// A quadratic with optimum `(a, a)` and value `-0.5141`, `a = sqrt(1.5)`,
    /// under the active inequality `x0 + x1 <= 2a`.
    struct Hs2 {
        a: f64,
        b: f64,
    }

    impl NlpOracle for Hs2 {
        fn objective(
            &self,
            x: &[f64],
        ) -> f64 {
            0.5 * (x[0] * x[0] + x[1] * x[1]) - self.a * (x[0] + x[1]) + self.b
        }

        fn gradient(
            &self,
            x: &[f64],
        ) -> Vec<f64> {
            vec![x[0] - self.a, x[1] - self.a]
        }

        fn hessian(
            &self,
            _x: &[f64],
            _me: &[f64],
            _mi: &[f64],
        ) -> SparseSym {
            SparseSym {
                n: 2,
                rows: vec![0, 1],
                cols: vec![0, 1],
                vals: vec![1.0, 1.0],
            }
        }

        fn constraints(
            &self,
            x: &[f64],
        ) -> Vec<f64> {
            vec![x[0] + x[1] - 2.0 * self.a]
        }

        fn jacobian(
            &self,
            _x: &[f64],
        ) -> SparseMat {
            SparseMat {
                rows_n: 1,
                cols_n: 2,
                rows: vec![0, 0],
                cols: vec![0, 1],
                vals: vec![1.0, 1.0],
            }
        }

        fn n_eq(&self) -> usize {
            0
        }

        fn n_ineq(&self) -> usize {
            1
        }
    }

    /// Hock–Schittkowski 6: `min (1-x0)^2 s.t. 10(x1 - x0^2) = 0`.
    struct Hs6;

    impl NlpOracle for Hs6 {
        fn objective(
            &self,
            x: &[f64],
        ) -> f64 {
            (1.0 - x[0]).powi(2)
        }

        fn gradient(
            &self,
            x: &[f64],
        ) -> Vec<f64> {
            vec![-2.0 * (1.0 - x[0]), 0.0]
        }

        fn hessian(
            &self,
            _x: &[f64],
            me: &[f64],
            _mi: &[f64],
        ) -> SparseSym {
            let lam = me.first().copied().unwrap_or(0.0);
            SparseSym {
                n: 2,
                rows: vec![0, 1],
                cols: vec![0, 1],
                vals: vec![2.0 - 20.0 * lam, 0.0],
            }
        }

        fn constraints(
            &self,
            x: &[f64],
        ) -> Vec<f64> {
            vec![10.0 * (x[1] - x[0] * x[0])]
        }

        fn jacobian(
            &self,
            x: &[f64],
        ) -> SparseMat {
            SparseMat {
                rows_n: 1,
                cols_n: 2,
                rows: vec![0, 0],
                cols: vec![0, 1],
                vals: vec![-20.0 * x[0], 10.0],
            }
        }

        fn n_eq(&self) -> usize {
            1
        }

        fn n_ineq(&self) -> usize {
            0
        }
    }

    /// Tridiagonal quadratic `0.5 x'Hx - b'x` with optimum `x = 1`.
    struct Tridiag {
        n: usize,
    }

    impl NlpOracle for Tridiag {
        fn objective(
            &self,
            x: &[f64],
        ) -> f64 {
            let mut f = 0.0;
            for i in 0..self.n {
                f += x[i] * x[i];
                if i + 1 < self.n {
                    f -= x[i] * x[i + 1];
                }
            }
            // -b.x with b = H * ones.
            f - bval(0, self.n) * x[0]
                - if self.n > 1 {
                    bval(self.n - 1, self.n) * x[self.n - 1]
                } else {
                    0.0
                }
        }

        fn gradient(
            &self,
            x: &[f64],
        ) -> Vec<f64> {
            let mut g = vec![0.0; self.n];
            for i in 0..self.n {
                let mut v = 2.0 * x[i];
                if i > 0 {
                    v -= x[i - 1];
                }
                if i + 1 < self.n {
                    v -= x[i + 1];
                }
                v -= bval(i, self.n);
                g[i] = v;
            }
            g
        }

        fn hessian(
            &self,
            _x: &[f64],
            _me: &[f64],
            _mi: &[f64],
        ) -> SparseSym {
            let mut rows = Vec::new();
            let mut cols = Vec::new();
            let mut vals = Vec::new();
            for i in 0..self.n {
                rows.push(i);
                cols.push(i);
                vals.push(2.0);
                if i + 1 < self.n {
                    rows.push(i);
                    cols.push(i + 1);
                    vals.push(-1.0);
                }
            }
            SparseSym {
                n: self.n,
                rows,
                cols,
                vals,
            }
        }

        fn constraints(
            &self,
            _x: &[f64],
        ) -> Vec<f64> {
            Vec::new()
        }

        fn jacobian(
            &self,
            _x: &[f64],
        ) -> SparseMat {
            SparseMat {
                rows_n: 0,
                cols_n: self.n,
                rows: Vec::new(),
                cols: Vec::new(),
                vals: Vec::new(),
            }
        }

        fn n_eq(&self) -> usize {
            0
        }

        fn n_ineq(&self) -> usize {
            0
        }
    }

    /// Right-hand side `b = H * ones` of the tridiagonal problem.
    fn bval(
        i: usize,
        n: usize,
    ) -> f64 {
        let mut v = 2.0;
        if i > 0 {
            v -= 1.0;
        }
        if i + 1 < n {
            v -= 1.0;
        }
        v
    }

    #[test]
    fn hs1_unconstrained() {
        let opts = SqpOptions {
            max_iter: 50,
            tol: 1e-10,
            trust_radius: 10.0,
        };
        let r = solve_sqp(&Hs1, &[0.0, 0.0], &opts).unwrap();
        assert!((r.x[0] - 1.0).abs() < 1e-8, "x0 = {}", r.x[0]);
        assert!(r.x[1].abs() < 1e-8, "x1 = {}", r.x[1]);
        assert!(r.objective.abs() < 1e-10);
        assert!(r.kkt_residual < 1e-8);
    }

    #[test]
    fn hs2_active_inequality() {
        let a = 1.5_f64.sqrt();
        let b = 1.5 - 0.5141;
        let oracle = Hs2 { a, b };
        let opts = SqpOptions {
            max_iter: 50,
            tol: 1e-10,
            trust_radius: 10.0,
        };
        let r = solve_sqp(&oracle, &[0.0, 0.0], &opts).unwrap();
        assert!((r.x[0] - a).abs() < 1e-8, "x0 = {}", r.x[0]);
        assert!((r.x[1] - a).abs() < 1e-8, "x1 = {}", r.x[1]);
        assert!((r.objective + 0.5141).abs() < 1e-6, "f = {}", r.objective);
    }

    #[test]
    fn hs6_equality_constrained() {
        let opts = SqpOptions {
            max_iter: 50,
            tol: 1e-10,
            trust_radius: 10.0,
        };
        let r = solve_sqp(&Hs6, &[0.5, 0.5], &opts).unwrap();
        assert!((r.x[0] - 1.0).abs() < 1e-6, "x0 = {}", r.x[0]);
        assert!((r.x[1] - 1.0).abs() < 1e-6, "x1 = {}", r.x[1]);
        assert!(r.objective.abs() < 1e-8);
    }

    #[test]
    fn sparse_tridiagonal_200() {
        let n = 200;
        let oracle = Tridiag { n };
        let x0 = vec![0.0; n];
        let opts = SqpOptions {
            max_iter: 20,
            tol: 1e-8,
            trust_radius: 1e3,
        };
        let r = solve_sqp(&oracle, &x0, &opts).unwrap();
        let max_err = r.x.iter().fold(0.0f64, |acc, v| acc.max((v - 1.0).abs()));
        assert!(max_err < 1e-5, "max error = {}", max_err);
        assert!(r.kkt_residual < 1e-5);
    }
}
