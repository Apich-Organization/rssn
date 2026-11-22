//! Convex quadratic programming by the Goldfarb–Idnani dual active-set method.
//!
//! The problem is `minimise 0.5 x'Qx + c'x` subject to general linear rows and
//! finite lower bounds with optional finite upper bounds. Every constraint is
//! rewritten in the normal form `a_i . x >= b_i` (equalities contribute both
//! `a_i . x >= b_i` and `-a_i . x >= -b_i`; upper bounds become `-x_j >= -u_j`).
//!
//! Goldfarb–Idnani starts from the unconstrained minimiser `x = -Q^{-1} c`
//! with an empty working set, then repeatedly adds a violated constraint. For a
//! candidate constraint `p` and working set `W` the search direction `dx` and
//! the multiplier changes `du` solve the block system
//!
//! ```text
//! [ Q   -A_W' ] [ dx ]   [ a_p ]
//! [ A_W   0   ] [ du ] = [  0  ]
//! ```
//!
//! and the multiplier of `p` is raised from zero until either `p` becomes
//! active (it joins `W`) or an active multiplier would turn negative (that
//! constraint leaves `W`). The iteration is finite for a positive definite
//! `Q`; a positive semidefinite `Q` is regularised on the diagonal by
//! `1e-12 * max(1, ||Q||_inf)` only when the Cholesky factorisation fails.

use std::collections::BTreeSet;

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};

use super::lp::Relation;

/// A convex quadratic program `minimise 0.5 x'Qx + c'x`.
#[derive(Debug, Clone, PartialEq)]
pub struct QpProblem {
    /// Symmetric positive semidefinite `n x n` matrix `Q`.
    pub q: Vec<Vec<f64>>,
    /// Linear objective coefficients `c`.
    pub c: Vec<f64>,
    /// Constraint rows; each row has one entry per variable.
    pub a: Vec<Vec<f64>>,
    /// Sense of each constraint row.
    pub rel: Vec<Relation>,
    /// Right-hand side of each constraint row.
    pub b: Vec<f64>,
    /// Finite lower bound of each variable.
    pub lower: Vec<f64>,
    /// Optional finite upper bound of each variable (`None` means `+inf`).
    pub upper: Vec<Option<f64>>,
}

/// Result of [`solve_qp`].
#[derive(Debug, Clone, PartialEq)]
pub struct QpSolution {
    /// Optimal point.
    pub x: Vec<f64>,
    /// Optimal objective value `0.5 x'Qx + c'x`.
    pub objective: f64,
    /// One multiplier per row of `a` (equalities carry a signed multiplier).
    pub multipliers: Vec<f64>,
    /// Number of active-set steps taken.
    pub iterations: usize,
    /// Row indices of `a` that are active at the solution, in ascending order.
    pub active: Vec<usize>,
}

/// A constraint in normal form `a . x >= b`.
struct Cons {
    a: Vec<f64>,
    b: f64,
    orig: Option<usize>,
}

/// Solves a convex quadratic program.
///
/// # Errors
/// Returns an error string when the dimensions disagree, when an entry is not
/// finite, when `q` has a negative eigenvalue (a non-convex problem), when the
/// problem is infeasible, or when the active-set iteration exceeds `max_iter`.
pub fn solve_qp(
    p: &QpProblem,
    max_iter: usize,
) -> Result<QpSolution, String> {
    let n = validate(p)?;
    if max_iter == 0 {
        return Err("iteration budget is zero".to_string());
    }

    // Symmetric part of Q.
    let mut q_orig = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            q_orig[i][j] = 0.5 * (p.q[i][j] + p.q[j][i]);
        }
    }

    // Factorise Q, regularising a semidefinite matrix only when needed.
    let qmat = Mat::from_rows(&q_orig).map_err(|e| e.to_string())?;
    let mut q_alg = q_orig.clone();
    let chol = match dense::cholesky(&qmat) {
        | Ok(l) => l,
        | Err(_) => {
            let (vals, _) = dense::eigen_symmetric(&qmat).map_err(|e| e.to_string())?;
            let scale = inf_norm(&q_orig).max(1.0);
            let min_eig = vals.first().copied().unwrap_or(0.0);
            if min_eig < -1e-8 * scale {
                return Err("quadratic matrix is not positive semidefinite".to_string());
            }
            let reg = 1e-12 * scale;
            for (i, row) in q_alg.iter_mut().enumerate() {
                row[i] += reg;
            }
            let qreg = Mat::from_rows(&q_alg).map_err(|e| e.to_string())?;
            dense::cholesky(&qreg)
                .map_err(|_| "quadratic matrix is not positive definite".to_string())?
        },
    };

    // Inverse of the (regularised) Q via the Cholesky factor.
    let ginv = inverse_from_cholesky(&chol, n);

    // Build the normal-form constraint list.
    let mut cons: Vec<Cons> = Vec::new();
    let mut row_map: Vec<Vec<(usize, f64)>> = vec![Vec::new(); p.a.len()];
    for i in 0..p.a.len() {
        let row = p.a[i].clone();
        let bi = p.b[i];
        match p.rel[i] {
            | Relation::Le => {
                let neg: Vec<f64> = row.iter().map(|v| -v).collect();
                let idx = cons.len();
                cons.push(Cons {
                    a: neg,
                    b: -bi,
                    orig: Some(i),
                });
                row_map[i].push((idx, 1.0));
            },
            | Relation::Ge => {
                let idx = cons.len();
                cons.push(Cons {
                    a: row,
                    b: bi,
                    orig: Some(i),
                });
                row_map[i].push((idx, 1.0));
            },
            | Relation::Eq => {
                let idx = cons.len();
                cons.push(Cons {
                    a: row.clone(),
                    b: bi,
                    orig: Some(i),
                });
                row_map[i].push((idx, 1.0));
                let neg: Vec<f64> = row.iter().map(|v| -v).collect();
                let idx2 = cons.len();
                cons.push(Cons {
                    a: neg,
                    b: -bi,
                    orig: Some(i),
                });
                row_map[i].push((idx2, -1.0));
            },
        }
    }
    for j in 0..n {
        let mut row = vec![0.0; n];
        row[j] = 1.0;
        cons.push(Cons {
            a: row,
            b: p.lower[j],
            orig: None,
        });
    }
    for j in 0..n {
        if let Some(u) = p.upper[j] {
            let mut row = vec![0.0; n];
            row[j] = -1.0;
            cons.push(Cons {
                a: row,
                b: -u,
                orig: None,
            });
        }
    }

    // Start at the unconstrained minimiser x = -Q^{-1} c.
    let cvec: Vec<f64> = p.c.iter().map(|v| -v).collect();
    let mut x = mat_vec(&ginv, &cvec);

    let mut in_working = vec![false; cons.len()];
    let mut w: Vec<usize> = Vec::new();
    let mut u: Vec<f64> = Vec::new();
    let mut iterations = 0usize;
    let tol = 1e-9;
    let inner_cap = 100 * (cons.len() + n + 1);
    let mut converged = false;

    'outer: for _ in 0..max_iter {
        // Most violated constraint outside the working set.
        let mut worst: Option<(usize, f64)> = None;
        for i in 0..cons.len() {
            if in_working[i] {
                continue;
            }
            let viol = cons[i].b - dot(&cons[i].a, &x);
            if viol > tol && worst.is_none_or(|(_, v)| viol > v) {
                worst = Some((i, viol));
            }
        }
        let Some((p_idx, _)) = worst else {
            converged = true;
            break 'outer;
        };

        let mut u_p = 0.0;
        let mut inner = 0usize;
        loop {
            inner += 1;
            iterations += 1;
            if inner > inner_cap {
                return Err("qp solver failed to converge".to_string());
            }
            let a_p = cons[p_idx].a.clone();
            let gap_p = cons[p_idx].b - dot(&a_p, &x);
            if gap_p <= tol {
                if u_p > tol {
                    w.push(p_idx);
                    u.push(u_p);
                    in_working[p_idx] = true;
                }
                break;
            }
            let (dx, du) = direction(&ginv, &cons, &w, &a_p, n)?;
            let ap_dx = dot(&a_p, &dx);
            if ap_dx <= 1e-14 {
                // The violated constraint cannot be improved along this
                // direction, so it conflicts with the working set.
                return Err("infeasible".to_string());
            }
            let alpha_feas = gap_p / ap_dx;

            let mut alpha_block = f64::INFINITY;
            let mut block_pos: Option<usize> = None;
            for (j, &du_j) in du.iter().enumerate() {
                if du_j < -1e-14 {
                    let t = (-u.get(j).copied().unwrap_or(0.0) / du_j).max(0.0);
                    if t < alpha_block {
                        alpha_block = t;
                        block_pos = Some(j);
                    }
                }
            }

            let alpha = alpha_feas.min(alpha_block);
            for (k, xk) in x.iter_mut().enumerate() {
                *xk += alpha * dx[k];
            }
            for j in 0..u.len() {
                u[j] += alpha * du[j];
            }
            u_p += alpha;

            if alpha_block < alpha_feas {
                if let Some(j) = block_pos {
                    let removed = w.remove(j);
                    u.remove(j);
                    in_working[removed] = false;
                }
                u_p = u_p.max(0.0);
                continue;
            }
            u_p = u_p.max(0.0);
            w.push(p_idx);
            u.push(u_p);
            in_working[p_idx] = true;
            break;
        }
    }

    if !converged {
        return Err("qp solver failed to converge".to_string());
    }

    // Multipliers of the original rows and the active-row set.
    let mut umult = vec![0.0; cons.len()];
    for (j, &ci) in w.iter().enumerate() {
        umult[ci] = u.get(j).copied().unwrap_or(0.0);
    }
    let mut multipliers = vec![0.0; p.a.len()];
    for i in 0..p.a.len() {
        let mut m = 0.0;
        for &(ci, sign) in &row_map[i] {
            m += sign * umult[ci];
        }
        multipliers[i] = m;
    }
    let mut active_set: BTreeSet<usize> = BTreeSet::new();
    for &ci in &w {
        if let Some(i) = cons[ci].orig {
            active_set.insert(i);
        }
    }

    let objective = 0.5 * quad_form(&q_orig, &x) + dot(&p.c, &x);
    Ok(QpSolution {
        x,
        objective,
        multipliers,
        iterations,
        active: active_set.into_iter().collect(),
    })
}

/// Solves the block system for `(dx, du)` given the working set `w`.
fn direction(
    ginv: &Mat,
    cons: &[Cons],
    w: &[usize],
    a_p: &[f64],
    _n: usize,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let g_a = mat_vec(ginv, a_p);
    let k = w.len();
    if k == 0 {
        return Ok((g_a, Vec::new()));
    }
    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(k);
    for &ci in w {
        cols.push(mat_vec(ginv, &cons[ci].a));
    }
    let mut m = vec![vec![0.0; k]; k];
    let mut rhs = vec![0.0; k];
    for i in 0..k {
        for j in 0..k {
            m[i][j] = dot(&cons[w[i]].a, &cols[j]);
        }
        rhs[i] = -dot(&cons[w[i]].a, &g_a);
    }
    let du = solve_small(&m, &rhs)?;
    let mut dx = g_a;
    for j in 0..k {
        let dj = du[j];
        if dj != 0.0 {
            for (t, dxt) in dx.iter_mut().enumerate() {
                *dxt += dj * cols[j][t];
            }
        }
    }
    Ok((dx, du))
}

/// Solves a small dense linear system, retrying with a tiny diagonal shift.
fn solve_small(
    m: &[Vec<f64>],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let mat = Mat::from_rows(m).map_err(|e| e.to_string())?;
    match dense::solve(&mat, rhs) {
        | Ok(x) => Ok(x),
        | Err(_) => {
            let mut mr = m.to_vec();
            let mut scale = 1.0f64;
            for (i, row) in mr.iter().enumerate() {
                scale = scale.max(row[i].abs());
            }
            for (i, row) in mr.iter_mut().enumerate() {
                row[i] += 1e-10 * scale;
            }
            let mat2 = Mat::from_rows(&mr).map_err(|e| e.to_string())?;
            dense::solve(&mat2, rhs).map_err(|_| "singular active set".to_string())
        },
    }
}

fn inverse_from_cholesky(
    l: &Mat,
    n: usize,
) -> Mat {
    let mut inv = Mat::zeros(n, n);
    for j in 0..n {
        let mut e = vec![0.0; n];
        e[j] = 1.0;
        let col = chol_solve(l, &e);
        for i in 0..n {
            inv.data[i * n + j] = col[i];
        }
    }
    inv
}

/// Solves `L L' x = b` for a lower-triangular Cholesky factor `L`.
fn chol_solve(
    l: &Mat,
    b: &[f64],
) -> Vec<f64> {
    let n = b.len();
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l.at(i, k) * y[k];
        }
        let d = l.at(i, i);
        y[i] = if d == 0.0 { 0.0 } else { s / d };
    }
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in i + 1..n {
            s -= l.at(k, i) * y[k];
        }
        let d = l.at(i, i);
        y[i] = if d == 0.0 { 0.0 } else { s / d };
    }
    y
}

fn mat_vec(
    m: &Mat,
    x: &[f64],
) -> Vec<f64> {
    (0..m.rows)
        .map(|i| {
            let mut s = 0.0;
            for j in 0..m.cols {
                s += m.at(i, j) * x.get(j).copied().unwrap_or(0.0);
            }
            s
        })
        .collect()
}

fn quad_form(
    q: &[Vec<f64>],
    x: &[f64],
) -> f64 {
    let mut s = 0.0;
    for i in 0..q.len() {
        let mut row = 0.0;
        for j in 0..q[i].len() {
            row += q[i][j] * x.get(j).copied().unwrap_or(0.0);
        }
        s += x.get(i).copied().unwrap_or(0.0) * row;
    }
    s
}

fn inf_norm(q: &[Vec<f64>]) -> f64 {
    let mut best = 0.0f64;
    for row in q {
        let mut s = 0.0;
        for v in row {
            s += v.abs();
        }
        best = best.max(s);
    }
    best
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Validates shapes and finiteness, returning the number of variables.
fn validate(p: &QpProblem) -> Result<usize, String> {
    let n = p.c.len();
    if p.q.len() != n || p.q.iter().any(|row| row.len() != n) {
        return Err("q must be a square matrix matching the objective".to_string());
    }
    if p.lower.len() != n || p.upper.len() != n {
        return Err("bound vectors must match the number of variables".to_string());
    }
    if p.rel.len() != p.a.len() || p.b.len() != p.a.len() {
        return Err("constraint rows, senses and right-hand sides must agree".to_string());
    }
    for row in &p.q {
        if row.iter().any(|v| !v.is_finite()) {
            return Err("q must be finite".to_string());
        }
    }
    if p.c.iter().any(|v| !v.is_finite()) {
        return Err("objective coefficients must be finite".to_string());
    }
    if p.lower.iter().any(|v| !v.is_finite()) {
        return Err("lower bounds must be finite".to_string());
    }
    for u in &p.upper {
        if let Some(v) = u {
            if !v.is_finite() {
                return Err("upper bounds must be finite when present".to_string());
            }
        }
    }
    for row in &p.a {
        if row.len() != n {
            return Err("constraint rows must have one entry per variable".to_string());
        }
        if row.iter().any(|v| !v.is_finite()) {
            return Err("constraint coefficients must be finite".to_string());
        }
    }
    if p.b.iter().any(|v| !v.is_finite()) {
        return Err("right-hand sides must be finite".to_string());
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(
        a: f64,
        b: f64,
    ) -> bool {
        (a - b).abs() <= 1e-8
    }

    #[test]
    fn halfspace_qp() {
        // min x^2 + y^2  s.t.  x + y >= 1  ->  (0.5, 0.5), objective 0.5.
        let p = QpProblem {
            q: vec![vec![2.0, 0.0], vec![0.0, 2.0]],
            c: vec![0.0, 0.0],
            a: vec![vec![1.0, 1.0]],
            rel: vec![Relation::Ge],
            b: vec![1.0],
            lower: vec![-1000.0, -1000.0],
            upper: vec![None, None],
        };
        let s = solve_qp(&p, 1000).unwrap();
        assert!(close(s.x[0], 0.5), "x0 = {}", s.x[0]);
        assert!(close(s.x[1], 0.5), "x1 = {}", s.x[1]);
        assert!(close(s.objective, 0.5));
        assert_eq!(s.active, vec![0]);
    }

    #[test]
    fn bounded_box_qp() {
        // min (x-2)^2 + (y-3)^2  s.t.  0 <= x,y <= 1  ->  (1, 1).
        // The reported objective is 0.5 x'Qx + c'x = x^2 + y^2 - 4x - 6y,
        // which equals -8 at the solution (the +13 constant is dropped).
        let p = QpProblem {
            q: vec![vec![2.0, 0.0], vec![0.0, 2.0]],
            c: vec![-4.0, -6.0],
            a: Vec::new(),
            rel: Vec::new(),
            b: Vec::new(),
            lower: vec![0.0, 0.0],
            upper: vec![Some(1.0), Some(1.0)],
        };
        let s = solve_qp(&p, 1000).unwrap();
        assert!(close(s.x[0], 1.0), "x0 = {}", s.x[0]);
        assert!(close(s.x[1], 1.0), "x1 = {}", s.x[1]);
        assert!(close(s.objective, -8.0), "objective = {}", s.objective);
    }

    #[test]
    fn equality_row_qp() {
        // min x^2 + y^2  s.t.  x + y = 2  ->  (1, 1), objective 2.
        let p = QpProblem {
            q: vec![vec![2.0, 0.0], vec![0.0, 2.0]],
            c: vec![0.0, 0.0],
            a: vec![vec![1.0, 1.0]],
            rel: vec![Relation::Eq],
            b: vec![2.0],
            lower: vec![-1000.0, -1000.0],
            upper: vec![None, None],
        };
        let s = solve_qp(&p, 1000).unwrap();
        assert!(close(s.x[0], 1.0), "x0 = {}", s.x[0]);
        assert!(close(s.x[1], 1.0), "x1 = {}", s.x[1]);
        assert!(close(s.objective, 2.0));
        assert_eq!(s.multipliers.len(), 1);
    }

    #[test]
    fn nonconvex_is_rejected() {
        // Q with a negative eigenvalue.
        let p = QpProblem {
            q: vec![vec![-2.0, 0.0], vec![0.0, 2.0]],
            c: vec![0.0, 0.0],
            a: Vec::new(),
            rel: Vec::new(),
            b: Vec::new(),
            lower: vec![-1000.0, -1000.0],
            upper: vec![None, None],
        };
        assert!(solve_qp(&p, 100).is_err());
    }

    #[test]
    fn dimension_mismatch_is_error() {
        let p = QpProblem {
            q: vec![vec![2.0]],
            c: vec![0.0, 0.0],
            a: Vec::new(),
            rel: Vec::new(),
            b: Vec::new(),
            lower: vec![-1.0, -1.0],
            upper: vec![None, None],
        };
        assert!(solve_qp(&p, 100).is_err());
    }
}
