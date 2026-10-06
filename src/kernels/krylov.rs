//! # Sparse matrices and Krylov solvers
//!
//! A compressed-sparse-row matrix ([`Csr`]), Jacobi and ILU(0)
//! preconditioners, and the iterative solvers conjugate gradients (CG),
//! restarted GMRES and BiCGSTAB, all matrix-free: they only need a
//! matrix-vector product closure.
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

/// Preconditioned BiCGSTAB for general nonsymmetric operators.
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
