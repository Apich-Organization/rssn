//! Algebraic multigrid for non-symmetric systems (Ruge-Stueben, Petrov-Galerkin aggregation).
//!
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
//! `BiCGSTAB` ([`AmgGeneral::precondition`]) when used on non-symmetric
//! operators; plain stand-alone cycling ([`AmgGeneral::solve`]) converges for
//! diffusion-dominated and moderately convective problems. Cycles reuse sparse
//! level work vectors. For repeated calls, [`AmgGeneral::workspace_len`] gives
//! the size of a caller-owned workspace.

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
use crate::kernels::numeric_common::NumericError;
use crate::kernels::reuse::f64_scratch;

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
    coarse_diag: Vec<f64>,
    coarse_lu: Option<Lu>,
    workspace_len: usize,
    sweeps: usize,
    cycle: Cycle,
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn residual_into(
    a: &Csr,
    x: &[f64],
    b: &[f64],
    work: &mut [f64],
) {
    a.matvec(x, work);
    for (ax, rhs) in work.iter_mut().zip(b) {
        *ax = rhs - *ax;
    }
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
    pub(crate) fn heap_bytes(&self) -> usize {
        let levels = self.levels.iter().map(|level| {
            level.a.heap_bytes()
                + level.p.heap_bytes()
                + level.r.heap_bytes()
                + level.diag.capacity() * size_of::<f64>()
        });
        self.levels.capacity() * size_of::<Level>()
            + levels.sum::<usize>()
            + self.coarse_a.heap_bytes()
            + self.coarse_diag.capacity() * size_of::<f64>()
            + self.coarse_lu.as_ref().map_or(0, Lu::heap_bytes)
    }

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
        let coarse_diag = cur.diagonal();
        let workspace_len = levels.iter().fold(
            levels.first().map_or(0, |level| level.a.rows),
            |size, level| size.saturating_add(level.r.rows.saturating_mul(2)),
        );
        Ok(Self {
            levels,
            coarse_a: cur,
            coarse_diag,
            coarse_lu,
            workspace_len,
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

    #[expect(
        clippy::arithmetic_side_effects,
        clippy::indexing_slicing,
        reason = "row vectors match the operator and index ranges come from the CSR structure"
    )]
    fn gauss_seidel(
        a: &Csr,
        x: &mut [f64],
        b: &[f64],
        diag: &[f64],
        sweeps: usize,
        forward: bool,
    ) {
        let mut update_row = |i: usize| {
            let mut s = b[i];
            for k in a.indptr[i]..a.indptr[i + 1] {
                let j = a.indices[k];
                if j != i {
                    s -= a.values[k] * x[j];
                }
            }
            x[i] = s / diag[i];
        };
        for _ in 0..sweeps {
            if forward {
                (0..a.rows).for_each(&mut update_row);
            } else {
                (0..a.rows).rev().for_each(&mut update_row);
            }
        }
    }

    fn cycle_at(
        &self,
        levels: &[Level],
        x: &mut [f64],
        b: &[f64],
        fine_work: &mut [f64],
        workspace: &mut [f64],
    ) {
        let Some((lev, coarser)) = levels.split_first() else {
            if let Some(lu) = &self.coarse_lu
                && lu.solve_into(b, x).is_ok()
            {
                return;
            }
            Self::gauss_seidel(&self.coarse_a, x, b, &self.coarse_diag, 20, true);
            return;
        };
        let (rc, rest) = workspace.split_at_mut(lev.r.rows);
        let (ec, rest) = rest.split_at_mut(lev.r.rows);

        Self::gauss_seidel(&lev.a, x, b, &lev.diag, self.sweeps, true);
        let (fine, _) = fine_work.split_at_mut(lev.a.rows);
        residual_into(&lev.a, x, b, fine);
        lev.r.matvec(fine, rc);
        ec.fill(0.0);
        let gamma = if self.cycle == Cycle::W {
            2
        } else {
            1
        };
        // The fine residual has been restricted, so recursion may overwrite the fine buffer.
        for _ in 0..gamma {
            self.cycle_at(coarser, ec, rc, fine_work, rest);
        }
        let (fine, _) = fine_work.split_at_mut(lev.a.rows);
        lev.p.matvec(ec, fine);
        for (xi, c) in x.iter_mut().zip(fine) {
            *xi += *c;
        }
        Self::gauss_seidel(&lev.a, x, b, &lev.diag, self.sweeps, false);
    }

    fn fine(&self) -> &Csr {
        self.levels.first().map_or(&self.coarse_a, |l| &l.a)
    }

    fn apply_cycle(
        &self,
        b: &[f64],
        x: &mut [f64],
        workspace: &mut [f64],
    ) {
        let rows = self.levels.first().map_or(0, |level| level.a.rows);
        let (fine_work, coarse_work) = workspace.split_at_mut(rows);
        self.cycle_at(&self.levels, x, b, fine_work, coarse_work);
    }

    fn cycle_from_pool(
        &self,
        b: &[f64],
        x: &mut [f64],
    ) {
        let mut workspace = f64_scratch(self.workspace_len);
        workspace.resize(self.workspace_len, 0.0);
        self.apply_cycle(b, x, &mut workspace);
    }

    fn check_cycle_shapes(
        &self,
        b: &[f64],
        x: &[f64],
        workspace: &[f64],
    ) -> Result<(), NumericError> {
        let order = self.fine().rows;
        if b.len() != order || x.len() != order {
            return Err(NumericError::ShapeMismatch(format!(
                "AMG expects RHS and iterate of length {order}, got {} and {}",
                b.len(),
                x.len()
            )));
        }
        if workspace.len() < self.workspace_len {
            return Err(NumericError::ShapeMismatch(format!(
                "AMG requires {} workspace entries, got {}",
                self.workspace_len,
                workspace.len()
            )));
        }
        Ok(())
    }

    /// Number of `f64` entries needed by the caller-supplied cycle workspace.
    ///
    /// Storage is one finest-level vector shared by all levels and two vectors
    /// sized for each coarser level.
    /// A hierarchy consisting of only the terminal level needs no workspace.
    #[must_use]
    pub const fn workspace_len(&self) -> usize {
        self.workspace_len
    }

    /// One cycle for `A x = b` starting from `x`.
    ///
    /// # Panics
    /// Panics if either vector length differs from the order of the fine operator.
    pub fn vcycle(
        &self,
        b: &[f64],
        x: &mut [f64],
    ) {
        let order = self.fine().rows;
        assert!(
            b.len() == order && x.len() == order,
            "AMG vector length mismatch"
        );
        self.cycle_from_pool(b, x);
    }

    /// Preconditioner action `z = M^{-1} r` (one cycle from zero).
    ///
    /// # Panics
    /// Panics if either vector length differs from the order of the fine operator.
    pub fn precondition(
        &self,
        r: &[f64],
        z: &mut [f64],
    ) {
        let order = self.fine().rows;
        assert!(
            r.len() == order && z.len() == order,
            "AMG vector length mismatch"
        );
        z.fill(0.0);
        self.cycle_from_pool(r, z);
    }

    /// Applies one cycle with the current `x` and a caller-owned workspace.
    ///
    /// A slice of at least [`Self::workspace_len`] entries can be reused across calls.
    /// Its initial contents need not be initialized to zero.
    ///
    /// # Errors
    /// Returns [`NumericError::ShapeMismatch`] before modifying `x` or `workspace`
    /// if either vector has the wrong length or the workspace is too short.
    pub fn vcycle_with_workspace(
        &self,
        b: &[f64],
        x: &mut [f64],
        workspace: &mut [f64],
    ) -> Result<(), NumericError> {
        self.check_cycle_shapes(b, x, workspace)?;
        self.apply_cycle(b, x, workspace);
        Ok(())
    }

    /// Applies the preconditioner from zero with a caller-owned workspace.
    ///
    /// # Errors
    /// Returns [`NumericError::ShapeMismatch`] before modifying `z` or `workspace`
    /// if either vector has the wrong length or the workspace is too short.
    pub fn precondition_with_workspace(
        &self,
        r: &[f64],
        z: &mut [f64],
        workspace: &mut [f64],
    ) -> Result<(), NumericError> {
        self.check_cycle_shapes(r, z, workspace)?;
        z.fill(0.0);
        self.apply_cycle(r, z, workspace);
        Ok(())
    }

    /// Solves `A x = b` by repeated cycles.
    ///
    /// # Panics
    /// Panics if the right-hand side or supplied initial iterate has the wrong length.
    #[must_use]
    pub fn solve(
        &self,
        b: &[f64],
        x0: Option<&[f64]>,
        tol: f64,
        max_cycles: usize,
    ) -> MgResult {
        let fine = self.fine();
        assert!(
            b.len() == fine.rows && x0.is_none_or(|x| x.len() == fine.rows),
            "AMG vector length mismatch"
        );
        let mut x = x0.map_or_else(|| vec![0.0; b.len()], <[f64]>::to_vec);
        let work_len = self.workspace_len.max(fine.rows);
        let mut workspace = f64_scratch(work_len);
        workspace.resize(work_len, 0.0);
        let bn = norm2(b).max(f64::MIN_POSITIVE);
        let res = |x: &[f64], work: &mut [f64]| {
            let (ax, _) = work.split_at_mut(fine.rows);
            residual_into(fine, x, b, ax);
            norm2(ax) / bn
        };
        let r = res(&x, &mut workspace);
        let mut history = vec![r];
        let mut converged = r <= tol;
        while !converged && history.len() <= max_cycles {
            self.apply_cycle(b, &mut x, &mut workspace);
            let r = res(&x, &mut workspace);
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
