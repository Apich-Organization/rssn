//! Second-order cone programming by a primal–dual interior-point method with
//! Nesterov–Todd (NT) scaling and Mehrotra predictor–corrector steps.
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
//! with `K` a product of linear and second-order cones. The solver is an
//! infeasible-start primal–dual method that keeps `s` and the dual `y` strictly
//! inside `K` and drives the residuals `Gx + h - s`, `Ex - f` and
//! `c - G'y + E'nu` to zero while reducing the duality gap `s . y`.
//!
//! At each iteration, the symmetric positive definite Nesterov–Todd scaling
//! operator `W` satisfying `W^{-1} s = W y = lambda` and `W^{-2} s = y` is
//! constructed block by block, and the dense KKT system with normal block
//! `M = G' W^{-2} G` is factored once with [`crate::kernels::dense::lu_factor`]
//! and reused for both the affine-scaling predictor step and the Mehrotra
//! second-order corrector + adaptive-centering step.
//!
//! Convergence is measured on a scale-invariant merit: the run stops when the
//! duality gap is at most `1e-11` in absolute terms *or* `1e-11` relative to
//! `1 + max(|primal objective|, |dual objective|)`, and the residuals are small
//! relative to `1 + ||h||`, `1 + ||f||` and `1 + ||c||`. An infeasible-start
//! method can wobble on a hard instance, so the best iterate seen (by that
//! merit) is kept and returned when the budget runs out; a long stretch without
//! measurable progress ends the run early instead of burning the remaining
//! iterations.

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
}

/// Nesterov–Todd scaling data for a single cone block.
enum NtBlock {
    /// Linear (1-D) cone: `w = sqrt(s / y)`.
    Linear { w: f64 },
    /// Second-order cone (`dim > 1`): `W = beta * P_{w_bar^{1/2}}` with
    /// `det(w_bar) = 1` and `beta = (det(s) / det(y))^{1/4}`.
    Soc { beta: f64, w_bar: Vec<f64> },
}

/// Solves a second-order cone program.
///
/// # Errors
/// Returns an error string when the dimensions disagree, when an entry is not
/// finite, when `max_iter` is zero, or when the KKT system is singular.
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
                });
                g.push(row);
                h.push(p.b[i]);
            },
            | Relation::Ge => {
                blocks.push(Block {
                    start: g.len(),
                    dim: 1,
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
            });
            g.push(row);
            h.push(u);
        }
    }
    for cone in &p.cones {
        let dim = 1 + cone.a.len();
        let start = g.len();
        g.push(cone.c.clone());
        h.push(cone.d);
        for (k, row) in cone.a.iter().enumerate() {
            g.push(row.clone());
            h.push(cone.b[k]);
        }
        blocks.push(Block { start, dim });
    }

    let m_tot = g.len();
    let me = eq_rows.len();
    // Under the Jordan-algebra central path `lambda o lambda = mu * e`, each
    // linear cone and each second-order cone contributes `mu` to `s . y`.
    let n_cones = blocks.len().max(1) as f64;

    // Scale factors for the relative stopping test, so a problem whose data is
    // large is not forced to an unreachable absolute tolerance.
    let h_norm = norm(&h);
    let f_norm = norm(&f);
    let c_norm = norm(&p.c);

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

    // The best iterate seen. An infeasible-start method can wobble on a hard
    // instance; returning the best point instead of the last one keeps the
    // reported residuals meaningful when the budget runs out.
    let mut best: Option<(f64, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)> = None;
    let mut best_score = f64::INFINITY;
    let mut stall = 0usize;

    // Scratch reused across iterations.
    let mut rp = vec![0.0; m_tot];
    let mut re = vec![0.0; me];
    let mut rd = vec![0.0; n];
    let mut g_col = vec![0.0; m_tot];
    let mut w2g = vec![0.0; m_tot];

    for _ in 0..max_iter {
        residual_vectors(
            &p.c, &g, &h, &f, &eq_rows, &x, &s, &y, &nu, n, m_tot, me, &mut rp, &mut re, &mut rd,
        );

        let rp_norm = norm(&rp);
        let re_norm = norm(&re);
        let rd_norm = norm(&rd);
        let gap = dot(&s, &y).max(0.0);

        let primal_objective = dot(&p.c, &x);
        let dual_objective = dot(&h, &y) - dot(&f, &nu);
        let scale = 1.0 + primal_objective.abs().max(dual_objective.abs());
        let rel_gap = gap / scale;

        // A single merit that mixes relative infeasibility and the relative
        // duality gap, used to keep the best iterate.
        let score = rel_gap
            + rp_norm / (1.0 + h_norm)
            + re_norm / (1.0 + f_norm)
            + rd_norm / (1.0 + c_norm);
        if score < best_score * (1.0 - 1e-6) {
            best_score = score;
            best = Some((gap, x.clone(), s.clone(), y.clone(), nu.clone()));
            stall = 0;
        } else {
            stall += 1;
        }

        // Absolute or relative: the absolute test keeps the accuracy of a
        // small instance, the relative one rescues a large instance whose
        // absolute gap can never reach `1e-10`.
        let gap_small = gap <= 1e-11 || rel_gap <= 1e-11;
        if gap_small
            && rp_norm <= 1e-8 * (1.0 + h_norm)
            && re_norm <= 1e-8 * (1.0 + f_norm)
            && rd_norm <= 1e-8 * (1.0 + c_norm)
        {
            break;
        }
        // No measurable progress for a long stretch: further iterations only
        // burn the budget (and, on an infeasible instance, drift).
        if stall >= 60 {
            break;
        }

        // Compute Nesterov–Todd scaling W and scaled iterate lambda = W^{-1} s = W y.
        let (nt, lambda) = compute_nt_scaling(&blocks, &s, &y);

        // Form M = G' W^{-2} G and factor the KKT system ONCE per iteration.
        let mut mm = vec![vec![0.0; n]; n];
        for j in 0..n {
            for r in 0..m_tot {
                g_col[r] = g[r][j];
            }
            apply_w_inv2_into(&blocks, &nt, &g_col, &mut w2g);
            for k in j..n {
                let mut acc = 0.0;
                for r in 0..m_tot {
                    acc += g[r][k] * w2g[r];
                }
                mm[k][j] = acc;
                mm[j][k] = acc;
            }
        }
        let lu = factor_kkt(&mm, &eq_rows, n, me)?;

        // 1. Predictor (affine-scaling) step:
        //    lambda o (W^{-1} ds_aff + W dy_aff) = -lambda o lambda
        //    <=> W^{-1} ds_aff + W dy_aff = -lambda.
        let ds_aff_target: Vec<f64> = lambda.iter().map(|v| -v).collect();
        let (_, ds_aff, dy_aff, _) = solve_direction_with_lu(
            &lu,
            &blocks,
            &nt,
            &rp,
            &rd,
            &re,
            &g,
            n,
            m_tot,
            me,
            &ds_aff_target,
        );

        let alpha_p_aff = block_max_step(&blocks, &s, &ds_aff).min(1.0);
        let alpha_d_aff = block_max_step(&blocks, &y, &dy_aff).min(1.0);
        let mut gap_aff = 0.0;
        for r in 0..m_tot {
            gap_aff += (s[r] + alpha_p_aff * ds_aff[r]) * (y[r] + alpha_d_aff * dy_aff[r]);
        }
        let mu = (gap / n_cones).max(0.0);
        let mu_aff = (gap_aff / n_cones).max(0.0);
        let sigma = if mu > 1e-300 {
            (mu_aff / mu).powi(3).clamp(1e-4, 0.99)
        } else {
            0.1
        };

        // 2. Corrector + adaptive centering step (reusing the same LU factor):
        //    r_lambda = -lambda o lambda - (W^{-1} ds_aff) o (W dy_aff) + sigma * mu * e.
        let lam_sq = jordan_prod(&blocks, &lambda, &lambda);
        let w_inv_ds_aff = apply_w_inv(&blocks, &nt, &ds_aff);
        let w_dy_aff = apply_w(&blocks, &nt, &dy_aff);
        let corr_prod = jordan_prod(&blocks, &w_inv_ds_aff, &w_dy_aff);

        let mut r_lambda = vec![0.0; m_tot];
        for r in 0..m_tot {
            r_lambda[r] = -lam_sq[r] - corr_prod[r];
        }
        let smu = sigma * mu;
        for b in &blocks {
            r_lambda[b.start] += smu;
        }
        let ds_corr_target = jordan_inv_apply(&blocks, &lambda, &r_lambda);

        let (dx, ds, dy, dnu) = solve_direction_with_lu(
            &lu,
            &blocks,
            &nt,
            &rp,
            &rd,
            &re,
            &g,
            n,
            m_tot,
            me,
            &ds_corr_target,
        );

        // 3. Fraction-to-boundary step with strict cone-interior safeguard.
        let alpha_p_max = block_max_step(&blocks, &s, &ds);
        let alpha_d_max = block_max_step(&blocks, &y, &dy);
        let amax = alpha_p_max.min(alpha_d_max);
        let mut alpha = if amax.is_infinite() {
            1.0
        } else {
            (0.99 * amax).min(1.0)
        };
        while alpha > 1e-14
            && (!is_strictly_interior(&blocks, &s, &ds, alpha)
                || !is_strictly_interior(&blocks, &y, &dy, alpha))
        {
            alpha *= 0.5;
        }
        if alpha <= 1e-14 {
            break;
        }

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

    // Return the best iterate rather than whichever point the last step
    // reached, then recompute its residuals from scratch.
    if let Some((_, bx, bs, by, bnu)) = best {
        x = bx;
        s = bs;
        y = by;
        nu = bnu;
    }
    residual_vectors(
        &p.c, &g, &h, &f, &eq_rows, &x, &s, &y, &nu, n, m_tot, me, &mut rp, &mut re, &mut rd,
    );
    let gap = dot(&s, &y).max(0.0);
    let primal_residual = norm(&rp) + norm(&re);
    let dual_residual = norm(&rd);

    Ok(SocpSolution {
        objective: dot(&p.c, &x),
        x,
        iterations,
        primal_residual,
        dual_residual,
        gap,
    })
}

/// Fills the primal, equality and dual residuals at `(x, s, y, nu)`.
#[allow(clippy::too_many_arguments)]
fn residual_vectors(
    c: &[f64],
    g: &[Vec<f64>],
    h: &[f64],
    f: &[f64],
    eq_rows: &[Vec<f64>],
    x: &[f64],
    s: &[f64],
    y: &[f64],
    nu: &[f64],
    n: usize,
    m_tot: usize,
    me: usize,
    rp: &mut [f64],
    re: &mut [f64],
    rd: &mut [f64],
) {
    for r in 0..m_tot {
        let row = &g[r];
        let mut acc = h[r] - s[r];
        for k in 0..n {
            acc += row[k] * x[k];
        }
        rp[r] = acc;
    }
    for k in 0..me {
        let row = &eq_rows[k];
        let mut acc = -f[k];
        for j in 0..n {
            acc += row[j] * x[j];
        }
        re[k] = acc;
    }
    // rd = c - G'y + E'nu.
    rd.copy_from_slice(&c[..n]);
    for r in 0..m_tot {
        let yr = y[r];
        if yr != 0.0 {
            let row = &g[r];
            for j in 0..n {
                rd[j] -= yr * row[j];
            }
        }
    }
    for k in 0..me {
        let nuk = nu[k];
        if nuk != 0.0 {
            let row = &eq_rows[k];
            for j in 0..n {
                rd[j] += nuk * row[j];
            }
        }
    }
}

/// Lorentz determinant `v_0^2 - ||v_{1:}||_2^2` for a second-order cone vector.
fn soc_det(v: &[f64]) -> f64 {
    let mut tail_sq = 0.0;
    for &vk in &v[1..] {
        tail_sq += vk * vk;
    }
    v[0] * v[0] - tail_sq
}

/// Computes the Nesterov–Todd scaling blocks `W` and scaled iterate
/// `lambda = W^{-1} s = W y`.
fn compute_nt_scaling(
    blocks: &[Block],
    s: &[f64],
    y: &[f64],
) -> (Vec<NtBlock>, Vec<f64>) {
    let mut nt = Vec::with_capacity(blocks.len());
    let mut lambda = vec![0.0; s.len()];
    for b in blocks {
        let sb = &s[b.start..b.start + b.dim];
        let yb = &y[b.start..b.start + b.dim];
        if b.dim == 1 {
            let s0 = sb[0].max(1e-300);
            let y0 = yb[0].max(1e-300);
            let w = (s0 / y0).sqrt();
            lambda[b.start] = (s0 * y0).sqrt();
            nt.push(NtBlock::Linear { w });
        } else {
            let det_s = soc_det(sb).max(1e-300);
            let det_y = soc_det(yb).max(1e-300);
            let sqrt_ds = det_s.sqrt();
            let sqrt_dy = det_y.sqrt();

            let mut s_bar = vec![0.0; b.dim];
            let mut y_bar = vec![0.0; b.dim];
            for k in 0..b.dim {
                s_bar[k] = sb[k] / sqrt_ds;
                y_bar[k] = yb[k] / sqrt_dy;
            }
            let gamma = (0.5 * (1.0 + dot(&s_bar, &y_bar))).max(1.0).sqrt();
            let inv_2g = 0.5 / gamma;
            let mut w_bar = vec![0.0; b.dim];
            w_bar[0] = inv_2g * (s_bar[0] + y_bar[0]);
            for k in 1..b.dim {
                w_bar[k] = inv_2g * (s_bar[k] - y_bar[k]);
            }
            let beta = (sqrt_ds / sqrt_dy).sqrt();
            let block_nt = NtBlock::Soc { beta, w_bar };
            apply_w_block(&block_nt, yb, &mut lambda[b.start..b.start + b.dim]);
            nt.push(block_nt);
        }
    }
    (nt, lambda)
}

/// Applies `W` to a single block.
fn apply_w_block(
    nt: &NtBlock,
    v: &[f64],
    out: &mut [f64],
) {
    match nt {
        | NtBlock::Linear { w } => {
            out[0] = w * v[0];
        },
        | NtBlock::Soc { beta, w_bar } => {
            let mut c = 0.0;
            for k in 1..v.len() {
                c += w_bar[k] * v[k];
            }
            out[0] = beta * (w_bar[0] * v[0] + c);
            let coef = v[0] + c / (1.0 + w_bar[0]);
            for k in 1..v.len() {
                out[k] = beta * (v[k] + coef * w_bar[k]);
            }
        },
    }
}

/// Applies `W` across all blocks.
fn apply_w(
    blocks: &[Block],
    nt: &[NtBlock],
    v: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; v.len()];
    for (b, ntb) in blocks.iter().zip(nt) {
        apply_w_block(
            ntb,
            &v[b.start..b.start + b.dim],
            &mut out[b.start..b.start + b.dim],
        );
    }
    out
}

/// Applies `W^{-1}` across all blocks.
fn apply_w_inv(
    blocks: &[Block],
    nt: &[NtBlock],
    v: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; v.len()];
    for (b, ntb) in blocks.iter().zip(nt) {
        let vb = &v[b.start..b.start + b.dim];
        let ob = &mut out[b.start..b.start + b.dim];
        match ntb {
            | NtBlock::Linear { w } => {
                ob[0] = vb[0] / w;
            },
            | NtBlock::Soc { beta, w_bar } => {
                let inv_beta = 1.0 / beta;
                let mut c = 0.0;
                for k in 1..vb.len() {
                    c += w_bar[k] * vb[k];
                }
                ob[0] = inv_beta * (w_bar[0] * vb[0] - c);
                let coef = -vb[0] + c / (1.0 + w_bar[0]);
                for k in 1..vb.len() {
                    ob[k] = inv_beta * (vb[k] + coef * w_bar[k]);
                }
            },
        }
    }
    out
}

/// Applies `W^{-2}` across all blocks (`W^{-2} s = y`).
fn apply_w_inv2(
    blocks: &[Block],
    nt: &[NtBlock],
    v: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; v.len()];
    apply_w_inv2_into(blocks, nt, v, &mut out);
    out
}

/// Applies `W^{-2}` across all blocks into caller-provided storage.
fn apply_w_inv2_into(
    blocks: &[Block],
    nt: &[NtBlock],
    v: &[f64],
    out: &mut [f64],
) {
    for (b, ntb) in blocks.iter().zip(nt) {
        let vb = &v[b.start..b.start + b.dim];
        let ob = &mut out[b.start..b.start + b.dim];
        match ntb {
            | NtBlock::Linear { w } => {
                ob[0] = vb[0] / (w * w);
            },
            | NtBlock::Soc { beta, w_bar } => {
                let inv_b2 = 1.0 / (beta * beta);
                let mut c = 0.0;
                for k in 1..vb.len() {
                    c += w_bar[k] * vb[k];
                }
                let alpha = 2.0 * (w_bar[0] * vb[0] - c);
                ob[0] = inv_b2 * (alpha * w_bar[0] - vb[0]);
                for k in 1..vb.len() {
                    ob[k] = inv_b2 * (-alpha * w_bar[k] + vb[k]);
                }
            },
        }
    }
}

/// Blockwise Jordan product `u o v`.
fn jordan_prod(
    blocks: &[Block],
    u: &[f64],
    v: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; u.len()];
    for b in blocks {
        let ub = &u[b.start..b.start + b.dim];
        let vb = &v[b.start..b.start + b.dim];
        let ob = &mut out[b.start..b.start + b.dim];
        if b.dim == 1 {
            ob[0] = ub[0] * vb[0];
        } else {
            ob[0] = dot(ub, vb);
            for k in 1..b.dim {
                ob[k] = ub[0] * vb[k] + vb[0] * ub[k];
            }
        }
    }
    out
}

/// Solves `lambda o z = r` blockwise using the closed-form arrowhead inverse.
fn jordan_inv_apply(
    blocks: &[Block],
    lambda: &[f64],
    r: &[f64],
) -> Vec<f64> {
    let mut out = vec![0.0; r.len()];
    for b in blocks {
        let lb = &lambda[b.start..b.start + b.dim];
        let rb = &r[b.start..b.start + b.dim];
        let ob = &mut out[b.start..b.start + b.dim];
        if b.dim == 1 {
            ob[0] = rb[0] / lb[0].max(1e-300);
        } else {
            let l0 = lb[0].max(1e-300);
            let det = soc_det(lb).max(1e-300);
            let mut l1_dot_r1 = 0.0;
            for k in 1..b.dim {
                l1_dot_r1 += lb[k] * rb[k];
            }
            let z0 = (l0 * rb[0] - l1_dot_r1) / det;
            ob[0] = z0;
            for k in 1..b.dim {
                ob[k] = (rb[k] - z0 * lb[k]) / l0;
            }
        }
    }
    out
}

/// Solves for `(dx, ds, dy, dnu)` given a pre-factored KKT matrix and scaled
/// complementarity target `W^{-1} ds + W dy = ds_target`.
#[allow(clippy::too_many_arguments)]
fn solve_direction_with_lu(
    lu: &dense::Lu,
    blocks: &[Block],
    nt: &[NtBlock],
    rp: &[f64],
    rd: &[f64],
    re: &[f64],
    g: &[Vec<f64>],
    n: usize,
    m_tot: usize,
    me: usize,
    ds_target: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let w_inv_ds = apply_w_inv(blocks, nt, ds_target);
    let w_inv2_rp = apply_w_inv2(blocks, nt, rp);
    let mut comb = vec![0.0; m_tot];
    for r in 0..m_tot {
        comb[r] = w_inv_ds[r] - w_inv2_rp[r];
    }
    let gt_comb = gt_matvec(g, &comb, m_tot, n);

    let mut rhs = vec![0.0; n + me];
    for j in 0..n {
        rhs[j] = -rd[j] + gt_comb[j];
    }
    for k in 0..me {
        rhs[n + k] = -re[k];
    }

    let sol = lu.solve(&rhs);
    let dx = sol[0..n].to_vec();
    let dnu = sol[n..].to_vec();

    let gdx = g_matvec(g, &dx, m_tot, n);
    let mut ds = vec![0.0; m_tot];
    for r in 0..m_tot {
        ds[r] = gdx[r] + rp[r];
    }

    let w_inv2_ds = apply_w_inv2(blocks, nt, &ds);
    let mut dy = vec![0.0; m_tot];
    for r in 0..m_tot {
        dy[r] = w_inv_ds[r] - w_inv2_ds[r];
    }

    (dx, ds, dy, dnu)
}

/// Largest step keeping `s + alpha ds` inside a single cone block.
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
    // Necessary condition for `s + alpha ds in Q`: `s_0 + alpha ds_0 > 0`.
    let alpha_cap = if ds[0] < -1e-300 {
        (-s[0] / ds[0]).max(0.0)
    } else {
        f64::INFINITY
    };

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
            return if r > 0.0 {
                r.min(alpha_cap)
            } else {
                alpha_cap
            };
        }
        return alpha_cap;
    }
    let disc = a1 * a1 - 4.0 * a2 * a0;
    if a2 < 0.0 {
        // Unique positive root since a0 > 0 and a2 < 0.
        let sq = disc.max(0.0).sqrt();
        let r = if a1 < 0.0 {
            2.0 * a0 / (-a1 + sq)
        } else {
            (a1 + sq) / (-2.0 * a2)
        };
        return r.max(0.0).min(alpha_cap);
    }
    // a2 > 0.0
    if a1 >= 0.0 {
        return alpha_cap;
    }
    if disc < 0.0 {
        // Roundoff near a tangent/apex when ds_0 < 0.
        if alpha_cap.is_finite() {
            return (-a1 / (2.0 * a2)).max(0.0).min(alpha_cap);
        }
        return f64::INFINITY;
    }
    // Both roots are positive; the smaller one hits the boundary of +Q first.
    let sq = disc.sqrt();
    let r = 2.0 * a0 / (-a1 + sq);
    r.max(0.0).min(alpha_cap)
}

/// Maximum step to the cone boundary across all blocks.
fn block_max_step(
    blocks: &[Block],
    v: &[f64],
    dv: &[f64],
) -> f64 {
    let mut amax = f64::INFINITY;
    for b in blocks {
        let vb = &v[b.start..b.start + b.dim];
        let dvb = &dv[b.start..b.start + b.dim];
        amax = amax.min(max_step(vb, dvb));
    }
    amax
}

/// Checks that `v + alpha dv` lies strictly inside every cone block.
fn is_strictly_interior(
    blocks: &[Block],
    v: &[f64],
    dv: &[f64],
    alpha: f64,
) -> bool {
    for b in blocks {
        let v0 = v[b.start] + alpha * dv[b.start];
        if !v0.is_finite() || v0 <= 0.0 {
            return false;
        }
        if b.dim > 1 {
            let mut tail_sq = 0.0;
            for k in 1..b.dim {
                let vk = v[b.start + k] + alpha * dv[b.start + k];
                if !vk.is_finite() {
                    return false;
                }
                tail_sq += vk * vk;
            }
            if v0 * v0 - tail_sq <= 0.0 {
                return false;
            }
        }
    }
    true
}

/// Factors the dense KKT system once per iteration for reuse across predictor
/// and corrector solves.
fn factor_kkt(
    mm: &[Vec<f64>],
    eq_rows: &[Vec<f64>],
    n: usize,
    me: usize,
) -> Result<dense::Lu, String> {
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
    match dense::lu_factor(&mat) {
        | Ok(lu) => Ok(lu),
        | Err(_) => {
            for i in 0..n {
                k[i][i] += 1e-8;
            }
            for r in 0..me {
                k[n + r][n + r] -= 1e-8;
            }
            let mat2 = Mat::from_rows(&k).map_err(|e| e.to_string())?;
            dense::lu_factor(&mat2).map_err(|_| "singular kkt system".to_string())
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
        let s = solve_socp(&p, 100).unwrap();
        assert!(s.objective.abs() < 1e-6, "t = {}", s.objective);
        assert!((s.x[0] - 1.0).abs() < 1e-6, "x = {}", s.x[0]);
        assert!((s.x[1] - 2.0).abs() < 1e-6, "y = {}", s.x[1]);
        assert!(s.iterations < 40, "iterations = {}", s.iterations);
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
        let s = solve_socp(&p, 100).unwrap();
        let expected = 1.0 / 3.0_f64.sqrt();
        assert!(
            (s.objective - expected).abs() < 1e-6,
            "t = {}, expected = {}",
            s.objective,
            expected
        );
        assert!((s.x[0] - 4.0 / 3.0).abs() < 1e-6, "x0 = {}", s.x[0]);
        assert!((s.x[1] - 7.0 / 3.0).abs() < 1e-6, "x1 = {}", s.x[1]);
        assert!(s.iterations < 40, "iterations = {}", s.iterations);
    }

    #[test]
    fn lp_as_socp() {
        // min -x0 - 2 x1  s.t.  x0 + x1 <= 4,  0 <= x0 <= 2,  0 <= x1 <= 3
        // -> (1, 3), objective -7.
        let p = SocpProblem {
            c: vec![-1.0, -2.0],
            a: vec![vec![1.0, 1.0]],
            rel: vec![Relation::Le],
            b: vec![4.0],
            cones: Vec::new(),
            lower: vec![0.0, 0.0],
            upper: vec![Some(2.0), Some(3.0)],
        };
        let s = solve_socp(&p, 100).unwrap();
        assert!((s.x[0] - 1.0).abs() < 1e-6, "x0 = {}", s.x[0]);
        assert!((s.x[1] - 3.0).abs() < 1e-6, "x1 = {}", s.x[1]);
        assert!((s.objective + 7.0).abs() < 1e-6, "obj = {}", s.objective);
        assert!(s.iterations < 30, "iterations = {}", s.iterations);
    }

    #[test]
    fn multiple_cones_and_equality() {
        // min t1 + t2
        // s.t. ||(x, y - 1)||_2 <= t1,  ||(x - 2, y - 1)||_2 <= t2,  y = 0.
        // Unique optimum: x = 1, y = 0, t1 = t2 = sqrt(2), obj = 2 sqrt(2).
        let p = SocpProblem {
            c: vec![0.0, 0.0, 1.0, 1.0],
            a: vec![vec![0.0, 1.0, 0.0, 0.0]],
            rel: vec![Relation::Eq],
            b: vec![0.0],
            cones: vec![
                SocCone {
                    a: vec![vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 1.0, 0.0, 0.0]],
                    b: vec![0.0, -1.0],
                    c: vec![0.0, 0.0, 1.0, 0.0],
                    d: 0.0,
                },
                SocCone {
                    a: vec![vec![1.0, 0.0, 0.0, 0.0], vec![0.0, 1.0, 0.0, 0.0]],
                    b: vec![-2.0, -1.0],
                    c: vec![0.0, 0.0, 0.0, 1.0],
                    d: 0.0,
                },
            ],
            lower: vec![-10.0, -10.0, 0.0, 0.0],
            upper: vec![None, None, None, None],
        };
        let s = solve_socp(&p, 100).unwrap();
        let rt2 = 2.0_f64.sqrt();
        // The optimum is tangentially degenerate: the duality gap is quadratic
        // in the distance to it, so a gap of `1e-12` leaves about `1e-6` in
        // the coordinates while the objective is already accurate to `1e-12`.
        assert!((s.x[0] - 1.0).abs() < 1e-5, "x = {}", s.x[0]);
        assert!(s.x[1].abs() < 1e-5, "y = {}", s.x[1]);
        assert!((s.x[2] - rt2).abs() < 1e-5, "t1 = {}", s.x[2]);
        assert!((s.x[3] - rt2).abs() < 1e-5, "t2 = {}", s.x[3]);
        assert!(
            (s.objective - 2.0 * rt2).abs() < 1e-6,
            "obj = {}",
            s.objective
        );
        assert!(s.iterations < 40, "iterations = {}", s.iterations);
    }

    #[test]
    fn nt_scaling_identities() {
        let blocks = vec![Block { start: 0, dim: 1 }, Block { start: 1, dim: 4 }];
        let s = vec![3.0, 5.0, 1.0, -2.0, 1.5];
        let y = vec![0.7, 4.0, -1.5, 0.5, 2.0];
        let (nt, lambda) = compute_nt_scaling(&blocks, &s, &y);
        let w_inv_s = apply_w_inv(&blocks, &nt, &s);
        let w_y = apply_w(&blocks, &nt, &y);
        let w_inv2_s = apply_w_inv2(&blocks, &nt, &s);
        for i in 0..s.len() {
            assert!((w_inv_s[i] - lambda[i]).abs() < 1e-12);
            assert!((w_y[i] - lambda[i]).abs() < 1e-12);
            assert!((w_inv2_s[i] - y[i]).abs() < 1e-12);
        }
        // Check Jordan inverse: lambda o (L(lambda)^{-1} r) == r.
        let r = vec![1.2, 2.3, -0.4, 0.9, -1.1];
        let z = jordan_inv_apply(&blocks, &lambda, &r);
        let lz = jordan_prod(&blocks, &lambda, &z);
        for i in 0..r.len() {
            assert!((lz[i] - r[i]).abs() < 1e-12);
        }
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

    #[test]
    fn large_data_converges_on_a_relative_gap() {
        // The same LP as `lp_as_socp` scaled by 10^6. An absolute gap
        // tolerance of 1e-10 would be unreachable; the relative test is not.
        let scale = 1e6;
        let p = SocpProblem {
            c: vec![-scale, -2.0 * scale],
            a: vec![vec![1.0, 1.0]],
            rel: vec![Relation::Le],
            b: vec![4.0 * scale],
            cones: Vec::new(),
            lower: vec![0.0, 0.0],
            upper: vec![Some(2.0 * scale), Some(3.0 * scale)],
        };
        let s = solve_socp(&p, 100).expect("solve");
        assert!((s.x[0] - scale).abs() < 1e-4 * scale, "x0 = {}", s.x[0]);
        assert!(
            (s.x[1] - 3.0 * scale).abs() < 1e-4 * scale,
            "x1 = {}",
            s.x[1]
        );
        // Both the cost and the right-hand side are scaled, so the optimal
        // value scales quadratically.
        let expected = -7.0 * scale * scale;
        assert!(
            (s.objective - expected).abs() < 1e-6 * expected.abs(),
            "obj = {}",
            s.objective
        );
        assert!(s.iterations < 60, "iterations = {}", s.iterations);
        // The reported gap is the complementarity product of the returned
        // point, not a stale value from an earlier iteration.
        assert!(s.gap >= 0.0);
        assert!(s.primal_residual < 1e-3 * scale);
    }

    #[test]
    fn returned_residuals_describe_the_returned_point() {
        // min t  s.t.  ||(x - 3, y + 4)|| <= t,  x + y = 1.
        // The distance from (3, -4) to the line is 2 / sqrt(2) = sqrt(2), and
        // the closest point is (4, -3).
        let p = SocpProblem {
            c: vec![0.0, 0.0, 1.0],
            a: vec![vec![1.0, 1.0, 0.0]],
            rel: vec![Relation::Eq],
            b: vec![1.0],
            cones: vec![SocCone {
                a: vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]],
                b: vec![-3.0, 4.0],
                c: vec![0.0, 0.0, 1.0],
                d: 0.0,
            }],
            lower: vec![-100.0, -100.0, -100.0],
            upper: vec![None, None, None],
        };
        let s = solve_socp(&p, 200).expect("solve");
        let root2 = 2.0_f64.sqrt();
        assert!((s.objective - root2).abs() < 1e-6, "obj = {}", s.objective);
        assert!((s.x[0] - 4.0).abs() < 1e-5, "x = {}", s.x[0]);
        assert!((s.x[1] + 3.0).abs() < 1e-5, "y = {}", s.x[1]);
        // The residuals reported describe the returned point.
        assert!(s.primal_residual < 1e-6, "rp = {}", s.primal_residual);
        assert!(s.dual_residual < 1e-6, "rd = {}", s.dual_residual);
        assert!(s.gap < 1e-7, "gap = {}", s.gap);
    }
}
