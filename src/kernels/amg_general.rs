//! # Algebraic multigrid for non-symmetric systems
//!
//! [`AmgGeneral`] builds a hierarchy for a general sparse [`Csr`] matrix
//! (convection-diffusion, upwinded transport, M-matrices that are not
//! symmetric) with two coarsening strategies:
//!
//! * [`AmgKind::RugeStuben`]: classical Ruge-Stueben coarsening. The strength
//!   of connection of an M-matrix is `-a_ij >= theta max_k (-a_ik)`; the
//!   first pass picks the point with the largest measure `lambda_i` (number
//!   of points strongly depending on `i`) as coarse and makes the points that
//!   depend on it fine, the second pass promotes fine points until every pair
//!   of strongly connected fine points shares a coarse point. Direct
//!   interpolation uses the usual positive/negative splitting. **The
//!   restriction is interpolation built from `A^T`** (strength, interpolation
//!   weights, same `C/F` splitting), `R = P_{A^T}^T`, and the coarse operator
//!   is the Petrov-Galerkin product `R A P`;
//! * [`AmgKind::SmoothedAggregation`]: aggregation on the symmetrised
//!   strength graph `|a_ij| + |a_ji|`, tentative prolongator from the
//!   near-null space, and *Petrov-Galerkin* smoothing: `P = (I - omega
//!   D^{-1} A) T` for interpolation but `R = ((I - omega D^{-1} A^T) T)^T` for
//!   restriction, `A_c = R A P`.
//!
//! Smoothing is Gauss-Seidel (forward before, backward after the coarse
//! correction). The cycle is a good right preconditioner for GMRES and
//! BiCGSTAB ([`AmgGeneral::precondition`]) when used on non-symmetric
//! operators; plain stand-alone cycling ([`AmgGeneral::solve`]) converges for
//! diffusion-dominated and moderately convective problems.

use crate::kernels::dense::Lu;
use crate::kernels::dense::Mat;
use crate::kernels::dense::{
    self,
};
use crate::kernels::krylov::Csr;
use crate::kernels::multigrid::AmgError;
use crate::kernels::multigrid::Cycle;
use crate::kernels::multigrid::MgResult;
use crate::kernels::multigrid::aggregate;
use crate::kernels::multigrid::spgemm;
use crate::kernels::multigrid::transpose;

/// Coarsening strategy of [`AmgGeneral`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AmgKind {
    /// Classical Ruge-Stueben with `A^T`-based restriction.
    RugeStuben,
    /// Smoothed aggregation with Petrov-Galerkin restriction.
    SmoothedAggregation,
}

/// Parameters of [`AmgGeneral`].
#[derive(Debug, Clone, PartialEq)]
pub struct AmgGeneralOptions {
    /// Coarsening strategy.
    pub kind: AmgKind,
    /// Strength-of-connection threshold.
    pub theta: f64,
    /// Maximal size of the coarsest (directly solved) level.
    pub coarse_size: usize,
    /// Maximal number of levels.
    pub max_levels: usize,
    /// Gauss-Seidel sweeps before and after the coarse correction.
    pub sweeps: usize,
    /// Cycle shape.
    pub cycle: Cycle,
    /// Near-null-space vectors (smoothed aggregation); empty selects the
    /// constant vector.
    pub near_nullspace: Vec<Vec<f64>>,
}

impl Default for AmgGeneralOptions {
    fn default() -> Self {
        Self {
            kind: AmgKind::RugeStuben,
            theta: 0.25,
            coarse_size: 40,
            max_levels: 25,
            sweeps: 1,
            cycle: Cycle::V,
            near_nullspace: Vec::new(),
        }
    }
}

#[derive(Debug, Clone)]
struct Level {
    a: Csr,
    p: Csr,
    r: Csr,
    diag: Vec<f64>,
}

/// Algebraic multigrid hierarchy for general (non-symmetric) matrices.
#[derive(Debug)]
pub struct AmgGeneral {
    levels: Vec<Level>,
    coarse_a: Csr,
    coarse_lu: Option<Lu>,
    sweeps: usize,
    cycle: Cycle,
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Strong dependencies `S_i = {j != i : -a_ij >= theta max_k(-a_ik)}`.
fn strength(
    a: &Csr,
    theta: f64,
) -> Vec<Vec<usize>> {
    (0..a.rows)
        .map(|i| {
            let mut mx = 0.0_f64;
            for k in a.indptr[i]..a.indptr[i + 1] {
                if a.indices[k] != i {
                    mx = mx.max(-a.values[k]);
                }
            }
            if mx <= 0.0 {
                return Vec::new();
            }
            (a.indptr[i]..a.indptr[i + 1])
                .filter(|&k| a.indices[k] != i && -a.values[k] >= theta * mx)
                .map(|k| a.indices[k])
                .collect()
        })
        .collect()
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Pt {
    Undecided,
    Coarse,
    Fine,
}

/// Ruge-Stueben `C/F` splitting from the strength graph.
fn cf_splitting(s: &[Vec<usize>]) -> Vec<Pt> {
    let n = s.len();
    let mut st: Vec<Vec<usize>> = vec![Vec::new(); n]; // points strongly depending on i
    for (i, row) in s.iter().enumerate() {
        for &j in row {
            st[j].push(i);
        }
    }
    let mut state = vec![Pt::Undecided; n];
    let mut lambda: Vec<i64> = st.iter().map(|v| v.len() as i64).collect();
    let mut heap = std::collections::BinaryHeap::new();
    for i in 0..n {
        if s[i].is_empty() && st[i].is_empty() {
            state[i] = Pt::Fine; // isolated: handled by the smoother
        } else {
            heap.push((lambda[i], i));
        }
    }
    while let Some((l, i)) = heap.pop() {
        if state[i] != Pt::Undecided || l != lambda[i] {
            continue;
        }
        state[i] = Pt::Coarse;
        for &j in &st[i] {
            if state[j] == Pt::Undecided {
                state[j] = Pt::Fine;
                for &k in &s[j] {
                    if state[k] == Pt::Undecided {
                        lambda[k] += 1;
                        heap.push((lambda[k], k));
                    }
                }
            }
        }
        for &k in &s[i] {
            if state[k] == Pt::Undecided {
                lambda[k] -= 1;
                heap.push((lambda[k], k));
            }
        }
    }
    for p in &mut state {
        if *p == Pt::Undecided {
            *p = Pt::Fine;
        }
    }
    // second pass: strongly connected fine points must share a coarse point
    for i in 0..n {
        if state[i] != Pt::Fine {
            continue;
        }
        let ci: Vec<usize> = s[i]
            .iter()
            .copied()
            .filter(|&j| state[j] == Pt::Coarse)
            .collect();
        for &j in &s[i] {
            if state[j] != Pt::Fine {
                continue;
            }
            let shares = s[j].iter().any(|k| ci.contains(k));
            if !shares {
                state[j] = Pt::Coarse;
            }
        }
    }
    state
}

/// Direct interpolation for the matrix `a` with the given strength graph and
/// splitting; returns `P` (`n x nc`) and the coarse index of each `C` point.
fn direct_interpolation(
    a: &Csr,
    s: &[Vec<usize>],
    state: &[Pt],
    cidx: &[usize],
    nc: usize,
) -> Csr {
    let n = a.rows;
    let mut trip: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..n {
        if state[i] == Pt::Coarse {
            trip.push((i, cidx[i], 1.0));
            continue;
        }
        let ci: Vec<usize> = s[i]
            .iter()
            .copied()
            .filter(|&j| state[j] == Pt::Coarse)
            .collect();
        if ci.is_empty() {
            continue;
        }
        let mut aii = 0.0;
        let (mut sum_n, mut sum_p) = (0.0, 0.0);
        let (mut sum_nc, mut sum_pc) = (0.0, 0.0);
        for k in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[k];
            let v = a.values[k];
            if j == i {
                aii = v;
            } else if v < 0.0 {
                sum_n += v;
                if ci.contains(&j) {
                    sum_nc += v;
                }
            } else {
                sum_p += v;
                if ci.contains(&j) {
                    sum_pc += v;
                }
            }
        }
        let alpha = if sum_nc == 0.0 {
            0.0
        } else {
            sum_n / sum_nc
        };
        let (beta, dd) = if sum_pc == 0.0 {
            (0.0, aii + sum_p)
        } else {
            (sum_p / sum_pc, aii)
        };
        if dd == 0.0 {
            continue;
        }
        for k in a.indptr[i]..a.indptr[i + 1] {
            let j = a.indices[k];
            if j == i || !ci.contains(&j) {
                continue;
            }
            let v = a.values[k];
            let w = if v < 0.0 {
                -alpha * v / dd
            } else {
                -beta * v / dd
            };
            trip.push((i, cidx[j], w));
        }
    }
    Csr::from_triplets(n, nc, &trip)
}

/// Tentative prolongator by local QR of the near-null space on aggregates.
fn tentative(
    n: usize,
    agg: &[usize],
    nagg: usize,
    b: &[Vec<f64>],
) -> (Csr, Vec<Vec<f64>>) {
    let nb = b.len();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); nagg];
    for (i, &g) in agg.iter().enumerate() {
        members[g].push(i);
    }
    let mut trip: Vec<(usize, usize, f64)> = Vec::new();
    let mut bc_entries: Vec<(usize, usize, f64)> = Vec::new();
    let mut ncols = 0usize;
    for mem in &members {
        let mut q: Vec<Vec<f64>> = Vec::new();
        let mut rmat = vec![vec![0.0; nb]; nb];
        for c in 0..nb {
            let mut v: Vec<f64> = mem.iter().map(|&i| b[c][i]).collect();
            for (qi, qv) in q.iter().enumerate() {
                let dotp: f64 = v.iter().zip(qv).map(|(x, y)| x * y).sum();
                rmat[qi][c] = dotp;
                for (x, y) in v.iter_mut().zip(qv) {
                    *x -= dotp * y;
                }
            }
            let nv = norm2(&v);
            if nv > 1e-12 && q.len() < mem.len() {
                rmat[q.len()][c] = nv;
                q.push(v.iter().map(|x| x / nv).collect());
            }
        }
        for (qi, qv) in q.iter().enumerate() {
            for (pos, &i) in mem.iter().enumerate() {
                trip.push((i, ncols + qi, qv[pos]));
            }
        }
        for c in 0..nb {
            for r in 0..q.len() {
                bc_entries.push((c, ncols + r, rmat[r][c]));
            }
        }
        ncols += q.len();
    }
    let mut bc: Vec<Vec<f64>> = vec![vec![0.0; ncols]; nb];
    for (c, col, v) in bc_entries {
        bc[c][col] = v;
    }
    (Csr::from_triplets(n, ncols, &trip), bc)
}

/// `(I - omega D^{-1} A) T`.
fn smooth(
    a: &Csr,
    diag: &[f64],
    omega: f64,
    tent: &Csr,
) -> Csr {
    let mut tr: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..a.rows {
        tr.push((i, i, 1.0));
        for k in a.indptr[i]..a.indptr[i + 1] {
            tr.push((i, a.indices[k], -omega * a.values[k] / diag[i]));
        }
    }
    let s = Csr::from_triplets(a.rows, a.rows, &tr);
    spgemm(&s, tent)
}

impl AmgGeneral {
    /// Builds the hierarchy.
    ///
    /// # Errors
    /// [`AmgError::BadMatrix`] for a non-square matrix, a non-positive
    /// diagonal entry, or near-null-space vectors of the wrong length.
    pub fn build(
        a: &Csr,
        opts: &AmgGeneralOptions,
    ) -> Result<Self, AmgError> {
        if a.rows != a.cols || a.rows == 0 {
            return Err(AmgError::BadMatrix);
        }
        let mut b: Vec<Vec<f64>> = if opts.near_nullspace.is_empty() {
            vec![vec![1.0; a.rows]]
        } else {
            opts.near_nullspace.clone()
        };
        if b.iter().any(|v| v.len() != a.rows) {
            return Err(AmgError::BadMatrix);
        }
        let mut cur = a.clone();
        let mut levels: Vec<Level> = Vec::new();
        while cur.rows > opts.coarse_size && levels.len() + 1 < opts.max_levels {
            let diag = cur.diagonal();
            if diag.iter().any(|&d| d <= 0.0) {
                return Err(AmgError::BadMatrix);
            }
            let at = transpose(&cur);
            let (p, r) = match opts.kind {
                | AmgKind::RugeStuben => {
                    let s = strength(&cur, opts.theta);
                    let state = cf_splitting(&s);
                    let mut cidx = vec![usize::MAX; cur.rows];
                    let mut nc = 0;
                    for (i, st) in state.iter().enumerate() {
                        if *st == Pt::Coarse {
                            cidx[i] = nc;
                            nc += 1;
                        }
                    }
                    if nc == 0 || nc >= cur.rows {
                        break;
                    }
                    let p = direct_interpolation(&cur, &s, &state, &cidx, nc);
                    let st_t = strength(&at, opts.theta);
                    let pt = direct_interpolation(&at, &st_t, &state, &cidx, nc);
                    (p, transpose(&pt))
                },
                | AmgKind::SmoothedAggregation => {
                    // symmetrised strength graph
                    let mut trip: Vec<(usize, usize, f64)> = Vec::new();
                    for i in 0..cur.rows {
                        for k in cur.indptr[i]..cur.indptr[i + 1] {
                            let v = 0.5 * cur.values[k].abs();
                            trip.push((i, cur.indices[k], v));
                            trip.push((cur.indices[k], i, v));
                        }
                    }
                    let sym = Csr::from_triplets(cur.rows, cur.rows, &trip);
                    let sd = sym.diagonal();
                    let (agg, nagg) = aggregate(&sym, opts.theta.min(0.1), &sd);
                    let nb = b.len();
                    if nagg * nb >= cur.rows {
                        break;
                    }
                    let (tent, bc) = tentative(cur.rows, &agg, nagg, &b);
                    let mut rho = 0.0_f64;
                    for i in 0..cur.rows {
                        let sr: f64 = (cur.indptr[i]..cur.indptr[i + 1])
                            .map(|k| cur.values[k].abs())
                            .sum();
                        rho = rho.max(sr / diag[i]);
                    }
                    let omega = 4.0 / (3.0 * rho);
                    let p = smooth(&cur, &diag, omega, &tent);
                    let dt = at.diagonal();
                    let pt = smooth(&at, &dt, omega, &tent);
                    b = bc;
                    (p, transpose(&pt))
                },
            };
            if p.cols == 0 || p.cols >= cur.rows {
                break;
            }
            let ap = spgemm(&cur, &p);
            let coarse = spgemm(&r, &ap);
            levels.push(Level { a: cur, p, r, diag });
            cur = coarse;
            if opts.kind == AmgKind::RugeStuben {
                b = vec![vec![1.0; cur.rows]];
            }
        }
        let n = cur.rows;
        let mut dm = Mat::zeros(n, n);
        for i in 0..n {
            for k in cur.indptr[i]..cur.indptr[i + 1] {
                dm.set(i, cur.indices[k], cur.values[k]);
            }
        }
        let coarse_lu = dense::lu_factor(&dm).ok();
        Ok(Self {
            levels,
            coarse_a: cur,
            coarse_lu,
            sweeps: opts.sweeps.max(1),
            cycle: opts.cycle,
        })
    }

    /// Number of levels including the coarsest.
    #[must_use]
    pub const fn num_levels(&self) -> usize {
        self.levels.len() + 1
    }

    /// Sum of the nonzeros of all level operators divided by the fine one.
    #[must_use]
    pub fn operator_complexity(&self) -> f64 {
        let fine = self
            .levels
            .first()
            .map_or(self.coarse_a.values.len(), |l| l.a.values.len());
        let total: usize = self.levels.iter().map(|l| l.a.values.len()).sum::<usize>()
            + self.coarse_a.values.len();
        total as f64 / fine as f64
    }

    fn gauss_seidel(
        a: &Csr,
        x: &mut [f64],
        b: &[f64],
        diag: &[f64],
        sweeps: usize,
        forward: bool,
    ) {
        for _ in 0..sweeps {
            let order: Box<dyn Iterator<Item = usize>> = if forward {
                Box::new(0..a.rows)
            } else {
                Box::new((0..a.rows).rev())
            };
            for i in order {
                let mut s = b[i];
                for k in a.indptr[i]..a.indptr[i + 1] {
                    let j = a.indices[k];
                    if j != i {
                        s -= a.values[k] * x[j];
                    }
                }
                x[i] = s / diag[i];
            }
        }
    }

    fn cycle_at(
        &self,
        l: usize,
        x: &mut [f64],
        b: &[f64],
    ) {
        if l == self.levels.len() {
            if let Some(lu) = &self.coarse_lu {
                x.copy_from_slice(&lu.solve(b));
            } else {
                let d = self.coarse_a.diagonal();
                Self::gauss_seidel(&self.coarse_a, x, b, &d, 20, true);
            }
            return;
        }
        let lev = &self.levels[l];
        Self::gauss_seidel(&lev.a, x, b, &lev.diag, self.sweeps, true);
        let mut ax = vec![0.0; x.len()];
        lev.a.matvec(x, &mut ax);
        let res: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
        let mut rc = vec![0.0; lev.r.rows];
        lev.r.matvec(&res, &mut rc);
        let mut ec = vec![0.0; rc.len()];
        let gamma = if self.cycle == Cycle::W {
            2
        } else {
            1
        };
        for _ in 0..gamma {
            self.cycle_at(l + 1, &mut ec, &rc);
        }
        let mut corr = vec![0.0; x.len()];
        lev.p.matvec(&ec, &mut corr);
        for (xi, c) in x.iter_mut().zip(&corr) {
            *xi += c;
        }
        Self::gauss_seidel(&lev.a, x, b, &lev.diag, self.sweeps, false);
    }

    /// One cycle for `A x = b` starting from `x`.
    pub fn vcycle(
        &self,
        b: &[f64],
        x: &mut [f64],
    ) {
        self.cycle_at(0, x, b);
    }

    /// Preconditioner action `z = M^{-1} r` (one cycle from zero).
    pub fn precondition(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        z.fill(0.0);
        self.cycle_at(0, z, r);
    }

    /// Solves `A x = b` by repeated cycles.
    #[must_use]
    pub fn solve(
        &self,
        b: &[f64],
        x0: Option<&[f64]>,
        tol: f64,
        max_cycles: usize,
    ) -> MgResult {
        let fine = self.levels.first().map_or(&self.coarse_a, |l| &l.a);
        let mut x = x0.map_or_else(|| vec![0.0; b.len()], <[f64]>::to_vec);
        let bn = norm2(b).max(f64::MIN_POSITIVE);
        let res = |x: &[f64]| {
            let mut ax = vec![0.0; b.len()];
            fine.matvec(x, &mut ax);
            norm2(&b.iter().zip(&ax).map(|(p, q)| p - q).collect::<Vec<_>>()) / bn
        };
        let mut history = vec![res(&x)];
        let mut converged = history[0] <= tol;
        while !converged && history.len() <= max_cycles {
            self.cycle_at(0, &mut x, b);
            let r = res(&x);
            history.push(r);
            converged = r <= tol;
        }
        MgResult {
            x,
            history,
            converged,
        }
    }
}
