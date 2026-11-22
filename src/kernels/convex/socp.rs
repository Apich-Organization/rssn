//! Second-order cone programming by a primal–dual interior-point method.
//!
//! A cone constraint is given in the form
//! `||A_i x + b_i||_2 <= c_i . x + d_i`, whose slack `s_i = (t_i, z_i)` with
//! `t_i = c_i . x + d_i` and `z_i = A_i x + b_i` belongs to the second-order
//! cone `Q = { s : s_0 >= ||s_{1:}||_2 }`. Linear rows and variable bounds are
//! handled as one-dimensional cones, so the whole problem is
//!
//! ```text
//! minimise  c . x   subject to   G x + h = s,   s in K,   E x = f,
//! ```
//!
//! with `K` a product of second-order cones. The solver is an infeasible-start
//! primal–dual method: it keeps `s` and the dual `y` strictly inside `K` and
//! drives the residuals `Gx + h - s`, `Ex - f` and `c - G'y + E'nu` to zero
//! while reducing the duality gap `s . y`. Complementarity is linearised
//! through the symmetric scaling given by the Hessian of the cone barrier
//! `phi(s) = -log det(s)`,
//!
//! ```text
//! dy + mu * Hess phi(s) ds = -y - mu * grad phi(s),
//! ```
//!
//! which keeps the primal and dual iterates balanced (it is the barrier-based
//! member of the Nesterov–Todd family). The resulting symmetric positive
//! definite system is solved densely with [`crate::kernels::dense::lu_factor`].
//! The centring parameter is `mu = sigma * gap / deg` with `sigma = 0.3` and a
//! `0.5` fraction-to-boundary step rule; this keeps the slack well away from
//! the cone apex, where the barrier Hessian would otherwise blow up.

use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};

use super::lp::Relation;

/// A second-order cone constraint `||A x + b||_2 <= c . x + d`.
#[derive(Debug, Clone, PartialEq)]
pub struct SocCone {
    /// Rows of the affine map inside the norm.
    pub a: Vec<Vec<f64>>,
    /// Constant term inside the norm.
    pub b: Vec<f64>,
    /// Coefficients of the right-hand side `c . x + d`.
    pub c: Vec<f64>,
    /// Constant of the right-hand side `c . x + d`.
    pub d: f64,
}

/// A second-order cone program.
#[derive(Debug, Clone, PartialEq)]
pub struct SocpProblem {
    /// Objective coefficients, `minimise c . x`.
    pub c: Vec<f64>,
    /// Linear constraint rows.
    pub a: Vec<Vec<f64>>,
    /// Sense of each linear constraint row.
    pub rel: Vec<Relation>,
    /// Right-hand side of each linear constraint row.
    pub b: Vec<f64>,
    /// Second-order cone constraints.
    pub cones: Vec<SocCone>,
    /// Finite lower bound of each variable.
    pub lower: Vec<f64>,
    /// Optional finite upper bound of each variable (`None` means `+inf`).
    pub upper: Vec<Option<f64>>,
}

/// Result of [`solve_socp`].
#[derive(Debug, Clone, PartialEq)]
pub struct SocpSolution {
    /// Primal point.
    pub x: Vec<f64>,
    /// Primal objective value `c . x`.
    pub objective: f64,
    /// Number of interior-point iterations.
    pub iterations: usize,
    /// Norm of the primal residual `(Gx + h - s, Ex - f)`.
    pub primal_residual: f64,
    /// Norm of the dual residual `c - G'y + E'nu`.
    pub dual_residual: f64,
    /// Duality gap `s . y`.
    pub gap: f64,
}

/// One block of the product cone.
struct Block {
    start: usize,
    dim: usize,
    soc: bool,
}

/// Solves a second-order cone program.
///
/// # Errors
/// Returns an error string when the dimensions disagree, when an entry is not
/// finite, or when the KKT system is singular.
pub fn solve_socp(
    p: &SocpProblem,
    max_iter: usize,
) -> Result<SocpSolution, String> {
    let n = validate(p)?;
    if max_iter == 0 {
        return Err("iteration budget is zero".to_string());
    }

    // Assemble the conic rows G x + h = s and the equalities E x = f.
    let mut g: Vec<Vec<f64>> = Vec::new();
    let mut h: Vec<f64> = Vec::new();
    let mut blocks: Vec<Block> = Vec::new();
    let mut eq_rows: Vec<Vec<f64>> = Vec::new();
    let mut f: Vec<f64> = Vec::new();

    for i in 0..p.a.len() {
        match p.rel[i] {
            | Relation::Le => {
                let row: Vec<f64> = p.a[i].iter().map(|v| -v).collect();
                blocks.push(Block {
                    start: g.len(),
                    dim: 1,
                    soc: false,
                });
                g.push(row);
                h.push(p.b[i]);
            },
            | Relation::Ge => {
                blocks.push(Block {
                    start: g.len(),
                    dim: 1,
                    soc: false,
                });
                g.push(p.a[i].clone());
                h.push(-p.b[i]);
            },
            | Relation::Eq => {
                eq_rows.push(p.a[i].clone());
                f.push(p.b[i]);
            },
        }
    }
    for j in 0..n {
        let mut row = vec![0.0; n];
        row[j] = 1.0;
        blocks.push(Block {
            start: g.len(),
            dim: 1,
            soc: false,
        });
        g.push(row);
        h.push(-p.lower[j]);
    }
    for j in 0..n {
        if let Some(u) = p.upper[j] {
            let mut row = vec![0.0; n];
            row[j] = -1.0;
            blocks.push(Block {
                start: g.len(),
                dim: 1,
                soc: false,
            });
            g.push(row);
            h.push(u);
        }
    }
    for cone in &p.cones {
        let dim = 1 + cone.a.len();
        let soc = dim > 1;
        let start = g.len();
        g.push(cone.c.clone());
        h.push(cone.d);
        for (k, row) in cone.a.iter().enumerate() {
            g.push(row.clone());
            h.push(cone.b[k]);
        }
        blocks.push(Block { start, dim, soc });
    }

    let m_tot = g.len();
    let me = eq_rows.len();
    let n_lin = blocks.iter().filter(|b| !b.soc).count();
    let n_soc = blocks.iter().filter(|b| b.soc).count();
    let bdeg = (n_lin + 2 * n_soc).max(1) as f64;

    // Strictly interior starting point: the identity of each cone.
    let mut x = vec![0.0; n];
    let mut s = vec![0.0; m_tot];
    let mut y = vec![0.0; m_tot];
    for b in &blocks {
        s[b.start] = 1.0;
        y[b.start] = 1.0;
    }
    let mut nu = vec![0.0; me];

    let mut iterations = 0usize;
    let mut primal_residual: f64;
    let mut dual_residual: f64;
    let mut gap: f64;
    let sigma = 0.3f64;

    for _ in 0..max_iter {
        let gx = g_matvec(&g, &x, m_tot, n);
        let mut rp = vec![0.0; m_tot];
        for r in 0..m_tot {
            rp[r] = gx[r] + h[r] - s[r];
        }
        let mut re = vec![0.0; me];
        for k in 0..me {
            let mut acc = -f[k];
            for j in 0..n {
                acc += eq_rows[k][j] * x[j];
            }
            re[k] = acc;
        }
        let gty = gt_matvec(&g, &y, m_tot, n);
        let mut rd = vec![0.0; n];
        for j in 0..n {
            let mut acc = p.c[j] - gty[j];
            for k in 0..me {
                acc += eq_rows[k][j] * nu[k];
            }
            rd[j] = acc;
        }
        gap = s.iter().zip(&y).map(|(a, b)| a * b).sum();
        primal_residual = norm(&rp) + norm(&re);
        dual_residual = norm(&rd);

        if gap < 1e-9 && primal_residual < 1e-7 && dual_residual < 1e-7 {
            break;
        }

        // Centring parameter of the target central-path point.
        let mu = (sigma * gap / bdeg).max(0.0);

        // Scaled Hessian of the complementarity: M = mu G' Hess phi(s) G.
        let mut mm = vec![vec![0.0; n]; n];
        for j in 0..n {
            let mut v = vec![0.0; m_tot];
            for r in 0..m_tot {
                v[r] = g[r][j];
            }
            let hv = apply_hess(&blocks, &s, &v);
            for k in 0..n {
                let mut acc = 0.0;
                for r in 0..m_tot {
                    acc += g[r][k] * hv[r];
                }
                mm[k][j] = mu * acc;
            }
        }

        let (dx, ds, dy, dnu) = direction(
            &blocks, &s, &y, &rp, &rd, &re, &mm, &g, &eq_rows, n, m_tot, me, mu,
        )?;

        let alpha = step_length(&blocks, &s, &y, &ds, &dy);
        for j in 0..n {
            x[j] += alpha * dx[j];
        }
        for r in 0..m_tot {
            s[r] += alpha * ds[r];
            y[r] += alpha * dy[r];
        }
        for k in 0..me {
            nu[k] += alpha * dnu[k];
        }
        iterations += 1;
    }

    // Final residuals at the returned iterate.
    let gx = g_matvec(&g, &x, m_tot, n);
    let mut rp = vec![0.0; m_tot];
    for r in 0..m_tot {
        rp[r] = gx[r] + h[r] - s[r];
    }
    let mut re = vec![0.0; me];
    for k in 0..me {
        let mut acc = -f[k];
        for j in 0..n {
            acc += eq_rows[k][j] * x[j];
        }
        re[k] = acc;
    }
    let gty = gt_matvec(&g, &y, m_tot, n);
    let mut rd = vec![0.0; n];
    for j in 0..n {
        let mut acc = p.c[j] - gty[j];
        for k in 0..me {
            acc += eq_rows[k][j] * nu[k];
        }
        rd[j] = acc;
    }
    gap = s.iter().zip(&y).map(|(a, b)| a * b).sum();
    primal_residual = norm(&rp) + norm(&re);
    dual_residual = norm(&rd);

    Ok(SocpSolution {
        objective: dot(&p.c, &x),
        x,
        iterations,
        primal_residual,
        dual_residual,
        gap,
    })
}

/// Computes one Newton direction targeting the central path at `mu`.
#[allow(clippy::too_many_arguments)]
fn direction(
    blocks: &[Block],
    s: &[f64],
    y: &[f64],
    rp: &[f64],
    rd: &[f64],
    re: &[f64],
    mm: &[Vec<f64>],
    g: &[Vec<f64>],
    eq_rows: &[Vec<f64>],
    n: usize,
    m_tot: usize,
    me: usize,
    mu: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), String> {
    let gphi = barrier_grad(blocks, s);
    let mut ypg = vec![0.0; m_tot];
    for r in 0..m_tot {
        ypg[r] = y[r] + mu * gphi[r];
    }
    let gt_ypg = gt_matvec(g, &ypg, m_tot, n);
    let hrp = apply_hess(blocks, s, rp);
    let gt_hrp = gt_matvec(g, &hrp, m_tot, n);

    let mut rhs = vec![0.0; n + me];
    for j in 0..n {
        rhs[j] = -rd[j] - gt_ypg[j] - mu * gt_hrp[j];
    }
    for k in 0..me {
        rhs[n + k] = -re[k];
    }

    let sol = solve_kkt(mm, eq_rows, &rhs, n, me)?;
    let dx = sol[0..n].to_vec();
    let dnu = sol[n..].to_vec();

    let gdx = g_matvec(g, &dx, m_tot, n);
    let mut ds = vec![0.0; m_tot];
    for r in 0..m_tot {
        ds[r] = gdx[r] + rp[r];
    }

    let hds = apply_hess(blocks, s, &ds);
    let mut dy = vec![0.0; m_tot];
    for r in 0..m_tot {
        dy[r] = -y[r] - mu * gphi[r] - mu * hds[r];
    }

    Ok((dx, ds, dy, dnu))
}

/// Gradient of the cone barrier `phi` at `s`, block by block.
fn barrier_grad(
    blocks: &[Block],
    s: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; s.len()];
    for b in blocks {
        let sb = &s[b.start..b.start + b.dim];
        if b.dim == 1 {
            let d = sb[0].max(f64::MIN_POSITIVE);
            out[b.start] = -1.0 / d;
        } else {
            let s0 = sb[0];
            let mut s2 = 0.0;
            for k in 1..b.dim {
                s2 += sb[k] * sb[k];
            }
            let det = (s0 * s0 - s2).max(1e-300);
            let c = 2.0 / det;
            out[b.start] = -c * s0;
            for k in 1..b.dim {
                out[b.start + k] = c * sb[k];
            }
        }
    }
    out
}

/// Applies the cone barrier Hessian `Hess phi(s)` to `v`, block by block.
fn apply_hess(
    blocks: &[Block],
    s: &[f64],
    v: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; v.len()];
    for b in blocks {
        let sb = &s[b.start..b.start + b.dim];
        let vb = &v[b.start..b.start + b.dim];
        if b.dim == 1 {
            let d = sb[0].max(f64::MIN_POSITIVE);
            out[b.start] = vb[0] / (d * d);
        } else {
            let s0 = sb[0];
            let mut s2 = 0.0;
            for k in 1..b.dim {
                s2 += sb[k] * sb[k];
            }
            let det = (s0 * s0 - s2).max(1e-300);
            let d2 = det * det;
            let mut sv = 0.0;
            for k in 1..b.dim {
                sv += sb[k] * vb[k];
            }
            let v0 = vb[0];
            out[b.start] = (2.0 * (s0 * s0 + s2) / d2) * v0 - (4.0 * s0 / d2) * sv;
            let coef = 2.0 / det;
            let sc = 4.0 * sv / d2;
            for k in 1..b.dim {
                out[b.start + k] = -(4.0 * s0 * sb[k] / d2) * v0 + coef * vb[k] + sc * sb[k];
            }
        }
    }
    out
}

/// Largest step keeping `s + alpha ds` inside every block of the cone.
fn max_step(
    s: &[f64],
    ds: &[f64],
) -> f64 {
    if s.len() == 1 {
        if ds[0] < -1e-300 {
            let r = -s[0] / ds[0];
            return if r > 0.0 { r } else { f64::INFINITY };
        }
        return f64::INFINITY;
    }
    let dim = s.len();
    let mut a0 = s[0] * s[0];
    let mut a1 = 2.0 * s[0] * ds[0];
    let mut a2 = ds[0] * ds[0];
    for k in 1..dim {
        a0 -= s[k] * s[k];
        a1 -= 2.0 * s[k] * ds[k];
        a2 -= ds[k] * ds[k];
    }
    if a0 <= 0.0 {
        return 0.0;
    }
    if a2.abs() < 1e-300 {
        if a1 < 0.0 {
            let r = -a0 / a1;
            return if r > 0.0 { r } else { f64::INFINITY };
        }
        return f64::INFINITY;
    }
    let disc = a1 * a1 - 4.0 * a2 * a0;
    if disc < 0.0 {
        return f64::INFINITY;
    }
    let sq = disc.sqrt();
    let r1 = (-a1 - sq) / (2.0 * a2);
    let r2 = (-a1 + sq) / (2.0 * a2);
    let mut best = f64::INFINITY;
    for r in [r1, r2] {
        if r > 1e-15 && r < best {
            best = r;
        }
    }
    best
}

/// Common step length for the primal and dual blocks.
fn step_length(
    blocks: &[Block],
    s: &[f64],
    y: &[f64],
    ds: &[f64],
    dy: &[f64],
) -> f64 {
    let mut alpha_s = f64::INFINITY;
    let mut alpha_y = f64::INFINITY;
    for b in blocks {
        let sb = &s[b.start..b.start + b.dim];
        let yb = &y[b.start..b.start + b.dim];
        let dsb = &ds[b.start..b.start + b.dim];
        let dyb = &dy[b.start..b.start + b.dim];
        alpha_s = alpha_s.min(max_step(sb, dsb));
        alpha_y = alpha_y.min(max_step(yb, dyb));
    }
    let amax = alpha_s.min(alpha_y);
    if amax.is_infinite() {
        1.0
    } else {
        (0.5 * amax).min(1.0)
    }
}

/// Solves the dense KKT system for `(dx, dnu)`.
fn solve_kkt(
    mm: &[Vec<f64>],
    eq_rows: &[Vec<f64>],
    rhs: &[f64],
    n: usize,
    me: usize,
) -> Result<Vec<f64>, String> {
    let size = n + me;
    let mut k = vec![vec![0.0; size]; size];
    for i in 0..n {
        for j in 0..n {
            k[i][j] = mm[i][j];
        }
        k[i][i] += 1e-12;
    }
    for r in 0..me {
        for i in 0..n {
            k[i][n + r] = eq_rows[r][i];
            k[n + r][i] = eq_rows[r][i];
        }
    }
    let mat = Mat::from_rows(&k).map_err(|e| e.to_string())?;
    match dense::solve(&mat, rhs) {
        | Ok(sol) => Ok(sol),
        | Err(_) => {
            for i in 0..n {
                k[i][i] += 1e-8;
            }
            let mat2 = Mat::from_rows(&k).map_err(|e| e.to_string())?;
            dense::solve(&mat2, rhs).map_err(|_| "singular kkt system".to_string())
        },
    }
}

fn g_matvec(
    g: &[Vec<f64>],
    x: &[f64],
    m: usize,
    n: usize,
) -> Vec<f64> {
    let mut out = vec![0.0; m];
    for r in 0..m {
        let mut acc = 0.0;
        for k in 0..n {
            acc += g[r][k] * x[k];
        }
        out[r] = acc;
    }
    out
}

fn gt_matvec(
    g: &[Vec<f64>],
    y: &[f64],
    m: usize,
    n: usize,
) -> Vec<f64> {
    let mut out = vec![0.0; n];
    for r in 0..m {
        let yr = y[r];
        if yr != 0.0 {
            for k in 0..n {
                out[k] += yr * g[r][k];
            }
        }
    }
    out
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn dot(
    a: &[f64],
    b: &[f64],
) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// Validates shapes and finiteness, returning the number of variables.
fn validate(p: &SocpProblem) -> Result<usize, String> {
    let n = p.c.len();
    if p.lower.len() != n || p.upper.len() != n {
        return Err("bound vectors must match the number of variables".to_string());
    }
    if p.rel.len() != p.a.len() || p.b.len() != p.a.len() {
        return Err("linear rows, senses and right-hand sides must agree".to_string());
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
            return Err("linear rows must have one entry per variable".to_string());
        }
        if row.iter().any(|v| !v.is_finite()) {
            return Err("linear coefficients must be finite".to_string());
        }
    }
    if p.b.iter().any(|v| !v.is_finite()) {
        return Err("right-hand sides must be finite".to_string());
    }
    for cone in &p.cones {
        if cone.c.len() != n {
            return Err("cone right-hand sides must have one entry per variable".to_string());
        }
        if cone.b.len() != cone.a.len() {
            return Err("cone constant vector must match its rows".to_string());
        }
        if !cone.d.is_finite() {
            return Err("cone constants must be finite".to_string());
        }
        for row in &cone.a {
            if row.len() != n {
                return Err("cone rows must have one entry per variable".to_string());
            }
            if row.iter().any(|v| !v.is_finite()) {
                return Err("cone coefficients must be finite".to_string());
            }
        }
        if cone.b.iter().any(|v| !v.is_finite()) {
            return Err("cone constants must be finite".to_string());
        }
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_to_point() {
        // min t  s.t.  ||(x - 1, y - 2)|| <= t  ->  t = 0, x = 1, y = 2.
        let p = SocpProblem {
            c: vec![0.0, 0.0, 1.0],
            a: Vec::new(),
            rel: Vec::new(),
            b: Vec::new(),
            cones: vec![SocCone {
                a: vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]],
                b: vec![-1.0, -2.0],
                c: vec![0.0, 0.0, 1.0],
                d: 0.0,
            }],
            lower: vec![-10.0, -10.0, -10.0],
            upper: vec![None, None, None],
        };
        let s = solve_socp(&p, 400).unwrap();
        assert!(s.objective.abs() < 1e-4, "t = {}", s.objective);
        assert!((s.x[0] - 1.0).abs() < 1e-4, "x = {}", s.x[0]);
        assert!((s.x[1] - 2.0).abs() < 1e-4, "y = {}", s.x[1]);
    }

    #[test]
    fn robust_least_squares() {
        // min ||A x - b||_2 as an SOCP; the least-squares residual is 1/sqrt(3).
        let a_mat = vec![vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 1.0]];
        let b_vec = vec![1.0, 2.0, 4.0];
        // Rows of the cone map: [A, 0] acting on (x0, x1, t).
        let mut rows = Vec::new();
        for r in 0..3 {
            rows.push(vec![a_mat[r][0], a_mat[r][1], 0.0]);
        }
        let p = SocpProblem {
            c: vec![0.0, 0.0, 1.0],
            a: Vec::new(),
            rel: Vec::new(),
            b: Vec::new(),
            cones: vec![SocCone {
                a: rows,
                b: vec![-b_vec[0], -b_vec[1], -b_vec[2]],
                c: vec![0.0, 0.0, 1.0],
                d: 0.0,
            }],
            lower: vec![-10.0, -10.0, -10.0],
            upper: vec![None, None, None],
        };
        let s = solve_socp(&p, 400).unwrap();
        let expected = 1.0 / 3.0_f64.sqrt();
        assert!(
            (s.objective - expected).abs() < 1e-4,
            "t = {}, expected = {}",
            s.objective,
            expected
        );
        assert!((s.x[0] - 4.0 / 3.0).abs() < 1e-4, "x0 = {}", s.x[0]);
        assert!((s.x[1] - 7.0 / 3.0).abs() < 1e-4, "x1 = {}", s.x[1]);
    }

    #[test]
    fn dimension_mismatch_is_error() {
        let p = SocpProblem {
            c: vec![1.0, 1.0],
            a: Vec::new(),
            rel: Vec::new(),
            b: Vec::new(),
            cones: Vec::new(),
            lower: vec![-1.0],
            upper: vec![None],
        };
        assert!(solve_socp(&p, 10).is_err());
    }
}
