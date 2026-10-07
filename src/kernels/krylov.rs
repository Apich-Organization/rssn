//! # Sparse matrices and Krylov solvers
//!
//! A compressed-sparse-row matrix ([`Csr`]), Jacobi and ILU(0)
//! preconditioners, and the iterative solvers conjugate gradients (CG),
//! restarted GMRES and BiCGSTAB, all matrix-free: they only need a
//! matrix-vector product closure.
//!
//! A sparse direct solver complements them: [`Csr::sparse_lu`] performs a
//! left-looking Gilbert-Peierls LU factorisation with threshold pivoting
//! after a fill-reducing symmetric ordering (reverse Cuthill-McKee via
//! [`Csr::rcm_ordering`] or minimum degree via
//! [`Csr::min_degree_ordering`]).
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

/// Compressed sparse row matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct Csr {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Row pointers (`rows + 1` entries).
    pub indptr: Vec<usize>,
    /// Column indices, sorted within each row.
    pub indices: Vec<usize>,
    /// Stored values.
    pub values: Vec<f64>,
}

impl Csr {
    /// Builds a matrix from `(row, col, value)` triplets; duplicates are summed.
    #[must_use]
    pub fn from_triplets(rows: usize, cols: usize, trip: &[(usize, usize, f64)]) -> Self {
        let mut t: Vec<(usize, usize, f64)> =
            trip.iter().copied().filter(|&(r, c, _)| r < rows && c < cols).collect();
        t.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        let mut indptr = vec![0; rows + 1];
        let mut indices: Vec<usize> = Vec::new();
        let mut values: Vec<f64> = Vec::new();
        let mut last: Option<(usize, usize)> = None;
        for (r, c, v) in t {
            if last == Some((r, c)) {
                if let Some(x) = values.last_mut() {
                    *x += v;
                }
            } else {
                indices.push(c);
                values.push(v);
                indptr[r + 1] += 1;
                last = Some((r, c));
            }
        }
        for i in 0..rows {
            indptr[i + 1] += indptr[i];
        }
        Self { rows, cols, indptr, indices, values }
    }

    /// Matrix-vector product `y = A x`.
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        for i in 0..self.rows {
            let mut s = 0.0;
            for k in self.indptr[i]..self.indptr[i + 1] {
                s += self.values[k] * x[self.indices[k]];
            }
            y[i] = s;
        }
    }

    /// Diagonal of the matrix.
    #[must_use]
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.rows.min(self.cols))
            .map(|i| {
                (self.indptr[i]..self.indptr[i + 1])
                    .find(|&k| self.indices[k] == i)
                    .map_or(0.0, |k| self.values[k])
            })
            .collect()
    }

    /// Jacobi preconditioner apply function data: the inverse diagonal.
    #[must_use]
    pub fn jacobi_inverse(&self) -> Vec<f64> {
        self.diagonal().iter().map(|&d| if d == 0.0 { 1.0 } else { 1.0 / d }).collect()
    }

    /// ILU(0) factorisation (incomplete LU keeping the sparsity pattern).
    /// Returns `None` if a zero pivot is met.
    #[must_use]
    pub fn ilu0(&self) -> Option<Ilu0> {
        let n = self.rows;
        let mut lu = self.values.clone();
        let mut diag = vec![0usize; n];
        for i in 0..n {
            diag[i] = (self.indptr[i]..self.indptr[i + 1]).find(|&k| self.indices[k] == i)?;
        }
        let mut pos = vec![usize::MAX; n];
        for i in 0..n {
            for k in self.indptr[i]..self.indptr[i + 1] {
                pos[self.indices[k]] = k;
            }
            for k in self.indptr[i]..diag[i] {
                let j = self.indices[k];
                let piv = lu[diag[j]];
                if piv == 0.0 {
                    return None;
                }
                lu[k] /= piv;
                let l = lu[k];
                for kk in diag[j] + 1..self.indptr[j + 1] {
                    let p = pos[self.indices[kk]];
                    if p != usize::MAX {
                        lu[p] -= l * lu[kk];
                    }
                }
            }
            if lu[diag[i]] == 0.0 {
                return None;
            }
            for k in self.indptr[i]..self.indptr[i + 1] {
                pos[self.indices[k]] = usize::MAX;
            }
        }
        Some(Ilu0 { a: Self { values: lu, ..self.clone() }, diag })
    }
}

/// ILU(0) factors stored in the pattern of the original matrix.
#[derive(Debug, Clone)]
pub struct Ilu0 {
    a: Csr,
    diag: Vec<usize>,
}

impl Ilu0 {
    /// Solves `L U z = r` (applies the preconditioner `M^-1`).
    pub fn apply(&self, r: &[f64], z: &mut [f64]) {
        let n = self.a.rows;
        for i in 0..n {
            let mut s = r[i];
            for k in self.a.indptr[i]..self.diag[i] {
                s -= self.a.values[k] * z[self.a.indices[k]];
            }
            z[i] = s;
        }
        for i in (0..n).rev() {
            let mut s = z[i];
            for k in self.diag[i] + 1..self.a.indptr[i + 1] {
                s -= self.a.values[k] * z[self.a.indices[k]];
            }
            z[i] = s / self.a.values[self.diag[i]];
        }
    }
}

/// Outcome of an iterative solve.
#[derive(Debug, Clone, PartialEq)]
pub struct IterResult {
    /// Approximate solution.
    pub x: Vec<f64>,
    /// Final relative residual `|b - A x| / |b|`.
    pub residual: f64,
    /// Iterations performed (matrix-vector products for GMRES).
    pub iterations: usize,
    /// Whether the tolerance was reached.
    pub converged: bool,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn nrm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Identity preconditioner helper.
pub fn identity_precond(r: &[f64], z: &mut [f64]) {
    z.copy_from_slice(r);
}

/// Preconditioned conjugate gradients for SPD operators.
///
/// `apply(x, y)` computes `y = A x`; `precond(r, z)` computes `z = M^-1 r`.
pub fn cg<A, M>(apply: A, precond: M, b: &[f64], x0: Option<&[f64]>, tol: f64, max_iter: usize) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    let n = b.len();
    let mut x = x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    let mut ax = vec![0.0; n];
    apply(&x, &mut ax);
    let mut r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let mut z = vec![0.0; n];
    precond(&r, &mut z);
    let mut p = z.clone();
    let mut rz = dot(&r, &z);
    let mut ap = vec![0.0; n];
    let mut it = 0;
    while it < max_iter && nrm(&r) / bn > tol {
        apply(&p, &mut ap);
        let pap = dot(&p, &ap);
        if pap == 0.0 {
            break;
        }
        let alpha = rz / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        precond(&r, &mut z);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        it += 1;
    }
    let res = nrm(&r) / bn;
    IterResult { x, residual: res, iterations: it, converged: res <= tol }
}

/// Restarted GMRES(`restart`) with right preconditioning.
pub fn gmres<A, M>(
    apply: A,
    precond: M,
    b: &[f64],
    x0: Option<&[f64]>,
    restart: usize,
    tol: f64,
    max_iter: usize,
) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    let n = b.len();
    let m = restart.clamp(1, n.max(1));
    let mut x = x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let mut total = 0;
    let mut tmp = vec![0.0; n];
    let mut w = vec![0.0; n];
    let mut res;
    loop {
        apply(&x, &mut tmp);
        let r: Vec<f64> = (0..n).map(|i| b[i] - tmp[i]).collect();
        let beta = nrm(&r);
        res = beta / bn;
        if res <= tol || total >= max_iter {
            break;
        }
        let mut v: Vec<Vec<f64>> = vec![r.iter().map(|e| e / beta).collect()];
        let mut h = vec![vec![0.0; m]; m + 1];
        let (mut cs, mut sn) = (vec![0.0; m], vec![0.0; m]);
        let mut g = vec![0.0; m + 1];
        g[0] = beta;
        let mut k_used = 0;
        for j in 0..m {
            precond(&v[j], &mut tmp);
            apply(&tmp, &mut w);
            for i in 0..=j {
                h[i][j] = dot(&w, &v[i]);
                for q in 0..n {
                    w[q] -= h[i][j] * v[i][q];
                }
            }
            h[j + 1][j] = nrm(&w);
            let breakdown = h[j + 1][j] < 1e-300;
            if !breakdown {
                v.push(w.iter().map(|e| e / h[j + 1][j]).collect());
            }
            for i in 0..j {
                let t = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = t;
            }
            let d = h[j][j].hypot(h[j + 1][j]);
            if d == 0.0 {
                cs[j] = 1.0;
                sn[j] = 0.0;
            } else {
                cs[j] = h[j][j] / d;
                sn[j] = h[j + 1][j] / d;
            }
            h[j][j] = d;
            h[j + 1][j] = 0.0;
            g[j + 1] = -sn[j] * g[j];
            g[j] *= cs[j];
            total += 1;
            k_used = j + 1;
            if g[j + 1].abs() / bn <= tol || breakdown || total >= max_iter {
                break;
            }
        }
        let mut y = vec![0.0; k_used];
        for i in (0..k_used).rev() {
            let mut s = g[i];
            for j in i + 1..k_used {
                s -= h[i][j] * y[j];
            }
            y[i] = s / h[i][i];
        }
        let mut u = vec![0.0; n];
        for j in 0..k_used {
            for q in 0..n {
                u[q] += y[j] * v[j][q];
            }
        }
        precond(&u, &mut tmp);
        for q in 0..n {
            x[q] += tmp[q];
        }
    }
    IterResult { x, residual: res, iterations: total, converged: res <= tol }
}

/// Preconditioned `BiCGSTAB` for general nonsymmetric operators.
pub fn bicgstab<A, M>(apply: A, precond: M, b: &[f64], x0: Option<&[f64]>, tol: f64, max_iter: usize) -> IterResult
where
    A: Fn(&[f64], &mut [f64]),
    M: Fn(&[f64], &mut [f64]),
{
    let n = b.len();
    let mut x = x0.map_or_else(|| vec![0.0; n], <[f64]>::to_vec);
    let mut t = vec![0.0; n];
    apply(&x, &mut t);
    let mut r: Vec<f64> = (0..n).map(|i| b[i] - t[i]).collect();
    let r0 = r.clone();
    let bn = nrm(b).max(f64::MIN_POSITIVE);
    let (mut rho, mut alpha, mut omega) = (1.0, 1.0, 1.0);
    let mut v = vec![0.0; n];
    let mut p = vec![0.0; n];
    let mut ph = vec![0.0; n];
    let mut sh = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut it = 0;
    while it < max_iter && nrm(&r) / bn > tol {
        let rho_new = dot(&r0, &r);
        if rho_new == 0.0 {
            break;
        }
        let beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        precond(&p, &mut ph);
        apply(&ph, &mut v);
        let denom = dot(&r0, &v);
        if denom == 0.0 {
            break;
        }
        alpha = rho / denom;
        for i in 0..n {
            s[i] = r[i] - alpha * v[i];
        }
        it += 1;
        if nrm(&s) / bn <= tol {
            for i in 0..n {
                x[i] += alpha * ph[i];
            }
            r.copy_from_slice(&s);
            break;
        }
        precond(&s, &mut sh);
        apply(&sh, &mut t);
        let tt = dot(&t, &t);
        if tt == 0.0 {
            break;
        }
        omega = dot(&t, &s) / tt;
        for i in 0..n {
            x[i] += alpha * ph[i] + omega * sh[i];
            r[i] = s[i] - omega * t[i];
        }
        if omega == 0.0 {
            break;
        }
    }
    let res = nrm(&r) / bn;
    IterResult { x, residual: res, iterations: it, converged: res <= tol }
}

// ------------------------------------------------------------------
// Sparse direct solver: fill-reducing orderings and left-looking LU
// ------------------------------------------------------------------

/// Fill-reducing symmetric ordering used before sparse factorisation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SparseOrdering {
    /// Keep the given ordering.
    Natural,
    /// Reverse Cuthill-McKee (bandwidth / profile reduction).
    Rcm,
    /// Minimum degree on the elimination graph (fill reduction).
    MinDegree,
}

/// Errors of the sparse direct solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SparseLuError {
    /// The matrix is not square.
    NotSquare,
    /// A structurally or numerically zero pivot was met in this column.
    Singular(usize),
}

/// Breadth-first search from `start`; fills `dist` for the visited
/// vertices (which must be `usize::MAX` on entry) and returns them in
/// visiting order.
fn bfs_order(adj: &[Vec<usize>], start: usize, dist: &mut [usize]) -> Vec<usize> {
    let mut queue = vec![start];
    dist[start] = 0;
    let mut head = 0;
    while head < queue.len() {
        let v = queue[head];
        head += 1;
        for &w in &adj[v] {
            if dist[w] == usize::MAX {
                dist[w] = dist[v] + 1;
                queue.push(w);
            }
        }
    }
    queue
}

impl Csr {
    /// Adjacency lists of the symmetrised pattern `A + A^T` (no self loops).
    fn symmetric_adjacency(&self) -> Vec<Vec<usize>> {
        let n = self.rows.min(self.cols);
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            for k in self.indptr[i]..self.indptr[i + 1] {
                let j = self.indices[k];
                if j < n && j != i {
                    adj[i].push(j);
                    adj[j].push(i);
                }
            }
        }
        for a in &mut adj {
            a.sort_unstable();
            a.dedup();
        }
        adj
    }

    /// Half bandwidth `max |i - j|` over the stored entries.
    #[must_use]
    pub fn bandwidth(&self) -> usize {
        let mut bw = 0;
        for i in 0..self.rows {
            for k in self.indptr[i]..self.indptr[i + 1] {
                bw = bw.max(i.abs_diff(self.indices[k]));
            }
        }
        bw
    }

    /// Symmetric permutation `B[i][j] = A[perm[i]][perm[j]]`.
    #[must_use]
    pub fn permute_symmetric(&self, perm: &[usize]) -> Self {
        let n = self.rows;
        let mut inv = vec![0; n];
        for (new, &old) in perm.iter().enumerate() {
            inv[old] = new;
        }
        let mut trip = Vec::with_capacity(self.values.len());
        for r in 0..n {
            for k in self.indptr[r]..self.indptr[r + 1] {
                trip.push((inv[r], inv[self.indices[k]], self.values[k]));
            }
        }
        Self::from_triplets(n, self.cols, &trip)
    }

    /// Reverse Cuthill-McKee ordering of the symmetrised pattern, with a
    /// George-Liu pseudo-peripheral start vertex per connected component.
    /// Returns `perm` with `perm[new] = old`.
    #[must_use]
    pub fn rcm_ordering(&self) -> Vec<usize> {
        let adj = self.symmetric_adjacency();
        let n = adj.len();
        let deg: Vec<usize> = adj.iter().map(Vec::len).collect();
        let mut visited = vec![false; n];
        let mut order: Vec<usize> = Vec::with_capacity(n);
        let mut dist = vec![usize::MAX; n];
        for seed in 0..n {
            if visited[seed] {
                continue;
            }
            // pseudo-peripheral start vertex of the component
            let mut start = seed;
            let mut queue = bfs_order(&adj, start, &mut dist);
            let mut ecc = dist[*queue.last().unwrap_or(&start)];
            for _ in 0..16 {
                let cand = queue
                    .iter()
                    .copied()
                    .filter(|&v| dist[v] == ecc)
                    .min_by_key(|&v| deg[v])
                    .unwrap_or(start);
                for &v in &queue {
                    dist[v] = usize::MAX;
                }
                let q2 = bfs_order(&adj, cand, &mut dist);
                let e2 = dist[*q2.last().unwrap_or(&cand)];
                if e2 > ecc {
                    start = cand;
                    ecc = e2;
                    for &v in &q2 {
                        dist[v] = usize::MAX;
                    }
                    queue = bfs_order(&adj, start, &mut dist);
                } else {
                    for &v in &q2 {
                        dist[v] = usize::MAX;
                    }
                    queue = bfs_order(&adj, start, &mut dist);
                    break;
                }
            }
            for &v in &queue {
                dist[v] = usize::MAX;
            }
            // Cuthill-McKee numbering from `start`
            let begin = order.len();
            order.push(start);
            visited[start] = true;
            let mut head = begin;
            while head < order.len() {
                let v = order[head];
                head += 1;
                let mut nb: Vec<usize> = adj[v].iter().copied().filter(|&w| !visited[w]).collect();
                nb.sort_by_key(|&w| deg[w]);
                for w in nb {
                    visited[w] = true;
                    order.push(w);
                }
            }
        }
        order.reverse();
        order
    }

    /// Minimum-degree ordering on the elimination graph of the
    /// symmetrised pattern: repeatedly eliminates a vertex of smallest
    /// current degree, turning its neighbourhood into a clique (exact
    /// degrees, lazily maintained heap). Returns `perm[new] = old`.
    #[must_use]
    pub fn min_degree_ordering(&self) -> Vec<usize> {
        use std::cmp::Reverse;
        use std::collections::{BTreeSet, BinaryHeap};
        let adj0 = self.symmetric_adjacency();
        let n = adj0.len();
        let mut adj: Vec<BTreeSet<usize>> =
            adj0.into_iter().map(|v| v.into_iter().collect()).collect();
        let mut heap: BinaryHeap<Reverse<(usize, usize)>> =
            (0..n).map(|v| Reverse((adj[v].len(), v))).collect();
        let mut done = vec![false; n];
        let mut order = Vec::with_capacity(n);
        while let Some(Reverse((d, v))) = heap.pop() {
            if done[v] || d != adj[v].len() {
                continue;
            }
            done[v] = true;
            order.push(v);
            let nb: Vec<usize> = adj[v].iter().copied().collect();
            for &u in &nb {
                adj[u].remove(&v);
            }
            for (a, &u) in nb.iter().enumerate() {
                for &w in &nb[a + 1..] {
                    if adj[u].insert(w) {
                        adj[w].insert(u);
                    }
                }
            }
            for &u in &nb {
                heap.push(Reverse((adj[u].len(), u)));
            }
            adj[v].clear();
        }
        order
    }

    /// Sparse LU factorisation with a fill-reducing symmetric ordering
    /// and threshold partial pivoting (left-looking Gilbert-Peierls with
    /// a depth-first symbolic reach); see [`SparseLu`].
    ///
    /// # Errors
    /// [`SparseLuError::NotSquare`] or [`SparseLuError::Singular`].
    pub fn sparse_lu(&self, ordering: SparseOrdering) -> Result<SparseLu, SparseLuError> {
        self.sparse_lu_with(ordering, 0.1)
    }

    /// As [`Csr::sparse_lu`] with an explicit pivot threshold in `(0, 1]`:
    /// the diagonal entry is kept as pivot when it is at least
    /// `threshold` times the largest candidate (`1.0` is classical
    /// partial pivoting, small values favour sparsity).
    ///
    /// # Errors
    /// [`SparseLuError::NotSquare`] or [`SparseLuError::Singular`].
    pub fn sparse_lu_with(
        &self,
        ordering: SparseOrdering,
        threshold: f64,
    ) -> Result<SparseLu, SparseLuError> {
        if self.rows != self.cols {
            return Err(SparseLuError::NotSquare);
        }
        let n = self.rows;
        let perm: Vec<usize> = match ordering {
            SparseOrdering::Natural => (0..n).collect(),
            SparseOrdering::Rcm => self.rcm_ordering(),
            SparseOrdering::MinDegree => self.min_degree_ordering(),
        };
        let b = self.permute_symmetric(&perm);
        // CSC of B
        let mut cp = vec![0usize; n + 1];
        for &c in &b.indices {
            cp[c + 1] += 1;
        }
        for j in 0..n {
            cp[j + 1] += cp[j];
        }
        let mut ci = vec![0usize; b.indices.len()];
        let mut cx = vec![0.0; b.indices.len()];
        let mut next = cp.clone();
        for r in 0..n {
            for k in b.indptr[r]..b.indptr[r + 1] {
                let c = b.indices[k];
                ci[next[c]] = r;
                cx[next[c]] = b.values[k];
                next[c] += 1;
            }
        }
        let tol = threshold.clamp(1e-12, 1.0);
        let mut lp = vec![0usize; n + 1];
        let mut li: Vec<usize> = Vec::new();
        let mut lx: Vec<f64> = Vec::new();
        let mut up = vec![0usize; n + 1];
        let mut ui: Vec<usize> = Vec::new();
        let mut ux: Vec<f64> = Vec::new();
        let mut pinv = vec![usize::MAX; n];
        let mut x = vec![0.0; n];
        let mut xi = vec![0usize; n];
        let mut pstack = vec![0usize; n];
        let mut mark = vec![usize::MAX; n];
        for k in 0..n {
            lp[k] = li.len();
            up[k] = ui.len();
            // symbolic reach of column k of B in the graph of L
            let mut top = n;
            for p in cp[k]..cp[k + 1] {
                let start = ci[p];
                if mark[start] == k {
                    continue;
                }
                // iterative depth-first search; the stack lives in xi[0..=head]
                let mut head = 0usize;
                xi[0] = start;
                loop {
                    let j = xi[head];
                    let jnew = pinv[j];
                    if mark[j] != k {
                        mark[j] = k;
                        pstack[head] = if jnew == usize::MAX { 0 } else { lp[jnew] };
                    }
                    let p2 = if jnew == usize::MAX { 0 } else { lp[jnew + 1] };
                    let mut q = pstack[head];
                    let mut pushed = false;
                    while q < p2 {
                        let i = li[q];
                        q += 1;
                        if mark[i] == k {
                            continue;
                        }
                        pstack[head] = q;
                        head += 1;
                        xi[head] = i;
                        pushed = true;
                        break;
                    }
                    if !pushed {
                        top -= 1;
                        xi[top] = j;
                        if head == 0 {
                            break;
                        }
                        head -= 1;
                    }
                }
            }
            for p in top..n {
                x[xi[p]] = 0.0;
            }
            for p in cp[k]..cp[k + 1] {
                x[ci[p]] = cx[p];
            }
            // sparse triangular solve with the unit lower factor
            for p in top..n {
                let j = xi[p];
                let jn = pinv[j];
                if jn == usize::MAX {
                    continue;
                }
                let xj = x[j];
                for q in lp[jn] + 1..lp[jn + 1] {
                    x[li[q]] -= lx[q] * xj;
                }
            }
            // pivot search
            let mut ipiv = usize::MAX;
            let mut amax = -1.0;
            for p in top..n {
                let i = xi[p];
                if pinv[i] == usize::MAX {
                    let t = x[i].abs();
                    if t > amax {
                        amax = t;
                        ipiv = i;
                    }
                } else {
                    ui.push(pinv[i]);
                    ux.push(x[i]);
                }
            }
            if ipiv == usize::MAX || amax <= 0.0 || !amax.is_finite() {
                return Err(SparseLuError::Singular(k));
            }
            if pinv[k] == usize::MAX && x[k].abs() >= amax * tol {
                ipiv = k;
            }
            let pivot = x[ipiv];
            ui.push(k);
            ux.push(pivot);
            pinv[ipiv] = k;
            li.push(ipiv);
            lx.push(1.0);
            for p in top..n {
                let i = xi[p];
                if pinv[i] == usize::MAX {
                    li.push(i);
                    lx.push(x[i] / pivot);
                }
                x[i] = 0.0;
            }
        }
        lp[n] = li.len();
        up[n] = ui.len();
        for v in &mut li {
            *v = pinv[*v];
        }
        Ok(SparseLu { n, perm, pinv, lp, li, lx, up, ui, ux })
    }
}

/// Sparse LU factors `P B = L U` of a symmetrically permuted matrix
/// `B = A[perm][:, perm]` (CSC storage of `L` with unit diagonal and `U`).
#[derive(Debug, Clone)]
pub struct SparseLu {
    n: usize,
    perm: Vec<usize>,
    pinv: Vec<usize>,
    lp: Vec<usize>,
    li: Vec<usize>,
    lx: Vec<f64>,
    up: Vec<usize>,
    ui: Vec<usize>,
    ux: Vec<f64>,
}

impl SparseLu {
    /// Number of stored entries in `L` and `U` (a measure of fill-in).
    #[must_use]
    pub fn nnz(&self) -> usize {
        self.lx.len() + self.ux.len()
    }

    /// Dimension of the factored matrix.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.n
    }

    /// The symmetric fill-reducing permutation used (`perm[new] = old`).
    #[must_use]
    pub fn permutation(&self) -> &[usize] {
        &self.perm
    }

    /// Solves `A x = b`.
    #[must_use]
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut w = vec![0.0; n];
        for i in 0..n {
            w[self.pinv[i]] = b[self.perm[i]];
        }
        // forward: L w = w (unit diagonal stored first in each column)
        for j in 0..n {
            let wj = w[j];
            if wj != 0.0 {
                for q in self.lp[j] + 1..self.lp[j + 1] {
                    w[self.li[q]] -= self.lx[q] * wj;
                }
            }
        }
        // backward: U z = w (diagonal stored last in each column)
        for j in (0..n).rev() {
            let last = self.up[j + 1] - 1;
            w[j] /= self.ux[last];
            let wj = w[j];
            if wj != 0.0 {
                for q in self.up[j]..last {
                    w[self.ui[q]] -= self.ux[q] * wj;
                }
            }
        }
        let mut x = vec![0.0; n];
        for i in 0..n {
            x[self.perm[i]] = w[i];
        }
        x
    }
}
