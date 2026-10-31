//! # Pseudo-spectral solvers: Fourier (periodic) and Chebyshev (bounded), 1D–3D
//!
//! | Capability | 1D | 2D | 3D | Boundary conditions |
//! |---|---|---|---|---|
//! | N-D FFT ([`fft_nd`]), spectral derivatives ([`SpectralGrid`]) | yes | yes | yes | periodic |
//! | Spectral Helmholtz `-Δu + σu = f` ([`SpectralGrid::helmholtz`]) | yes | yes | yes | periodic |
//! | ETDRK4 for `u_t = L u + N(u, t)` ([`etdrk4`]) with 2/3 dealiasing | yes | yes | yes | periodic |
//! | Chebyshev collocation Helmholtz ([`ChebGrid::solve_helmholtz`]) | yes | yes | yes | Dirichlet, Neumann, Robin per face, time dependent data |
//! | Chebyshev heat equation θ-method ([`ChebGrid::heat`]) | yes | yes | yes | same |
//!
//! Fourier grids need power-of-two sizes (the FFT of
//! [`crate::kernels::transforms`] is radix-2). Chebyshev grids use
//! Gauss–Lobatto points, dense differentiation matrices and the sparse LU of
//! [`crate::kernels::krylov`] for the (dense-block) collocation systems, so
//! they are meant for modest sizes (up to a few thousand unknowns).

use crate::kernels::krylov::Csr;
use crate::kernels::krylov::SparseOrdering;
use crate::kernels::transforms::fft_slice;
use crate::kernels::transforms::ifft_slice;
use crate::sim::boundary::BcValue;
use crate::sim::boundary::BoundaryCondition;
use crate::sim::boundary::Face;
use crate::sim::boundary::FaceBcs;
use num_complex::Complex;
use std::f64::consts::PI;

/// In-place N-dimensional FFT of row-major (x fastest) data with extents `n`
/// (each a power of two). The inverse includes the `1/N` normalisation.
///
/// # Errors
/// If an extent is not a power of two or the length does not match.
pub fn fft_nd(
    data: &mut [Complex<f64>],
    n: &[usize],
    inverse: bool,
) -> Result<(), String> {
    if n.iter().any(|&m| !m.is_power_of_two()) || data.len() != n.iter().product::<usize>() {
        return Err("fft_nd needs power-of-two extents matching the data".into());
    }
    let mut stride = 1;
    for &len in n {
        let outer = data.len() / (stride * len);
        let mut line = vec![Complex::new(0.0, 0.0); len];
        for o in 0..outer {
            for s in 0..stride {
                let base = o * stride * len + s;
                for i in 0..len {
                    line[i] = data[base + i * stride];
                }
                if inverse {
                    ifft_slice(&mut line);
                } else {
                    fft_slice(&mut line);
                }
                for i in 0..len {
                    data[base + i * stride] = line[i];
                }
            }
        }
        stride *= len;
    }
    Ok(())
}

/// Periodic spectral grid (box with power-of-two points per axis).
#[derive(Clone, Debug)]
pub struct SpectralGrid {
    /// Active axes.
    pub dim: usize,
    /// Points per axis.
    pub n: [usize; 3],
    /// Box lengths.
    pub len: [f64; 3],
}

impl SpectralGrid {
    /// Creates a grid.
    ///
    /// # Errors
    /// For bad sizes.
    pub fn new(
        n: &[usize],
        len: &[f64],
    ) -> Result<Self, String> {
        let dim = n.len();
        if !(1..=3).contains(&dim)
            || len.len() != dim
            || n.iter().any(|&m| !m.is_power_of_two() || m < 4)
        {
            return Err("spectral grid needs 1-3 power-of-two axes with at least 4 points".into());
        }
        let mut g = Self {
            dim,
            n: [1; 3],
            len: [1.0; 3],
        };
        g.n[..dim].copy_from_slice(&n[..dim]);
        g.len[..dim].copy_from_slice(&len[..dim]);
        Ok(g)
    }

    /// Total number of points.
    #[must_use]
    pub fn size(&self) -> usize {
        self.n.iter().product()
    }

    /// Grid spacing along an axis.
    #[must_use]
    pub fn h(
        &self,
        a: usize,
    ) -> f64 {
        self.len[a] / self.n[a] as f64
    }

    /// Point coordinates (`x = i h`, origin at the corner).
    #[must_use]
    pub fn coord(
        &self,
        p: usize,
    ) -> [f64; 3] {
        let i = [
            p % self.n[0],
            (p / self.n[0]) % self.n[1],
            p / (self.n[0] * self.n[1]),
        ];
        let mut x = [0.0; 3];
        for a in 0..self.dim {
            x[a] = i[a] as f64 * self.h(a);
        }
        x
    }

    /// Wavenumber of mode index `j` along `axis` (negative frequencies wrapped).
    #[must_use]
    pub fn wavenumber(
        &self,
        axis: usize,
        j: usize,
    ) -> f64 {
        let n = self.n[axis];
        let m = if j <= n / 2 {
            j as f64
        } else {
            j as f64 - n as f64
        };
        2.0 * PI / self.len[axis] * m
    }

    /// Wavenumber vector of flat mode index `p`.
    #[must_use]
    pub fn kvec(
        &self,
        p: usize,
    ) -> [f64; 3] {
        let i = [
            p % self.n[0],
            (p / self.n[0]) % self.n[1],
            p / (self.n[0] * self.n[1]),
        ];
        let mut k = [0.0; 3];
        for a in 0..self.dim {
            k[a] = self.wavenumber(a, i[a]);
        }
        k
    }

    /// `|k|²` of flat mode index `p`.
    #[must_use]
    pub fn k2(
        &self,
        p: usize,
    ) -> f64 {
        self.kvec(p).iter().map(|k| k * k).sum()
    }

    /// Transform of real data to the spectral domain.
    #[must_use]
    pub fn to_spec(
        &self,
        u: &[f64],
    ) -> Vec<Complex<f64>> {
        let mut c: Vec<Complex<f64>> = u.iter().map(|&v| Complex::new(v, 0.0)).collect();
        let _ = fft_nd(&mut c, &self.n[..self.dim], false);
        c
    }

    /// Transform of complex data to the spectral domain.
    #[must_use]
    pub fn to_spec_c(
        &self,
        u: &[Complex<f64>],
    ) -> Vec<Complex<f64>> {
        let mut c = u.to_vec();
        let _ = fft_nd(&mut c, &self.n[..self.dim], false);
        c
    }

    /// Real part of the inverse transform.
    #[must_use]
    pub fn to_phys(
        &self,
        uh: &[Complex<f64>],
    ) -> Vec<f64> {
        let mut c = uh.to_vec();
        let _ = fft_nd(&mut c, &self.n[..self.dim], true);
        c.into_iter().map(|v| v.re).collect()
    }

    /// Complex inverse transform.
    #[must_use]
    pub fn to_phys_c(
        &self,
        uh: &[Complex<f64>],
    ) -> Vec<Complex<f64>> {
        let mut c = uh.to_vec();
        let _ = fft_nd(&mut c, &self.n[..self.dim], true);
        c
    }

    /// Spectral derivative `∂/∂x_axis` (the Nyquist mode is dropped).
    #[must_use]
    pub fn deriv(
        &self,
        uh: &[Complex<f64>],
        axis: usize,
    ) -> Vec<Complex<f64>> {
        uh.iter()
            .enumerate()
            .map(|(p, v)| {
                let i = [
                    p % self.n[0],
                    (p / self.n[0]) % self.n[1],
                    p / (self.n[0] * self.n[1]),
                ];
                if self.n[axis] > 1 && i[axis] == self.n[axis] / 2 {
                    Complex::new(0.0, 0.0)
                } else {
                    Complex::new(0.0, self.wavenumber(axis, i[axis])) * v
                }
            })
            .collect()
    }

    /// 2/3-rule dealiasing mask (`true` = keep).
    #[must_use]
    pub fn dealias_mask(&self) -> Vec<bool> {
        (0..self.size())
            .map(|p| {
                let i = [
                    p % self.n[0],
                    (p / self.n[0]) % self.n[1],
                    p / (self.n[0] * self.n[1]),
                ];
                (0..self.dim).all(|a| {
                    let m = if i[a] <= self.n[a] / 2 {
                        i[a]
                    } else {
                        self.n[a] - i[a]
                    };
                    3 * m < self.n[a]
                })
            })
            .collect()
    }

    /// Solves `-Δu + σ u = f` spectrally (periodic). For `σ = 0` the mean of
    /// `f` must vanish and the solution gets zero mean.
    ///
    /// # Errors
    /// For a singular mode (`σ = 0` with non-zero mean forcing).
    pub fn helmholtz(
        &self,
        sigma: f64,
        f: &[f64],
    ) -> Result<Vec<f64>, String> {
        let mut fh = self.to_spec(f);
        for (p, v) in fh.iter_mut().enumerate() {
            let d = self.k2(p) + sigma;
            if d.abs() < 1e-14 {
                if v.norm() > 1e-9 * (f.len() as f64) {
                    return Err("singular Helmholtz mode with non-zero forcing".into());
                }
                *v = Complex::new(0.0, 0.0);
            } else {
                *v /= d;
            }
        }
        Ok(self.to_phys(&fh))
    }
}

/// Nonlinear term `N(t, uh)` of a spectral PDE, returning its spectral transform.
pub type NonlinearTerm = dyn Fn(f64, &[Complex<f64>]) -> Vec<Complex<f64>>;

/// Fourth-order ETDRK4 time stepping of `u_t = L u + N(u, t)` in the spectral domain.
///
/// Integrates `u_t = L u + N(u, t)` in the spectral domain with the
/// fourth-order exponential time differencing Runge–Kutta scheme of
/// Cox and Matthews (contour-integral coefficients of Kassam and Trefethen).
///
/// `lin[p]` is the (diagonal) linear symbol; `nonlinear(t, uh)` returns the
/// spectral transform of `N`.
pub fn etdrk4(
    mut uh: Vec<Complex<f64>>,
    lin: &[Complex<f64>],
    nonlinear: &NonlinearTerm,
    t_end: f64,
    steps: usize,
) -> Vec<Complex<f64>> {
    let h = t_end / steps as f64;
    let n = uh.len();
    let m = 32;
    let mean = |f: &dyn Fn(Complex<f64>) -> Complex<f64>, hl: Complex<f64>| -> Complex<f64> {
        let mut s = Complex::new(0.0, 0.0);
        for j in 0..m {
            let r = Complex::new(0.0, PI * (f64::from(j) + 0.5) / f64::from(m)).exp();
            s += f(hl + r);
        }
        s / f64::from(m)
    };
    let mut e = vec![Complex::new(0.0, 0.0); n];
    let mut e2 = e.clone();
    let mut q = e.clone();
    let mut f1 = e.clone();
    let mut f2 = e.clone();
    let mut f3 = e.clone();
    for p in 0..n {
        let hl = lin[p] * h;
        e[p] = hl.exp();
        e2[p] = (hl * 0.5).exp();
        q[p] = mean(&|z| ((z * 0.5).exp() - 1.0) / z, hl) * h;
        f1[p] = mean(
            &|z| (-4.0 - z + z.exp() * (4.0 - 3.0 * z + z * z)) / (z * z * z),
            hl,
        ) * h;
        f2[p] = mean(&|z| (2.0 + z + z.exp() * (-2.0 + z)) / (z * z * z), hl) * h;
        f3[p] = mean(
            &|z| (-4.0 - 3.0 * z - z * z + z.exp() * (4.0 - z)) / (z * z * z),
            hl,
        ) * h;
    }
    let mut t = 0.0;
    for _ in 0..steps {
        let nv = nonlinear(t, &uh);
        let a: Vec<_> = (0..n).map(|p| e2[p] * uh[p] + q[p] * nv[p]).collect();
        let na = nonlinear(t + 0.5 * h, &a);
        let b: Vec<_> = (0..n).map(|p| e2[p] * uh[p] + q[p] * na[p]).collect();
        let nb = nonlinear(t + 0.5 * h, &b);
        let c: Vec<_> = (0..n)
            .map(|p| e2[p] * a[p] + q[p] * (nb[p] * 2.0 - nv[p]))
            .collect();
        let nc = nonlinear(t + h, &c);
        uh = (0..n)
            .map(|p| e[p] * uh[p] + f1[p] * nv[p] + f2[p] * (na[p] + nb[p]) * 2.0 + f3[p] * nc[p])
            .collect();
        t += h;
    }
    uh
}

/// Chebyshev–Gauss–Lobatto grid on a box (1D–3D, ascending nodes per axis).
#[derive(Clone, Debug)]
pub struct ChebGrid {
    /// Active axes.
    pub dim: usize,
    /// Polynomial degree per axis (`degree + 1` nodes).
    pub deg: [usize; 3],
    /// Lower corner.
    pub lo: [f64; 3],
    /// Upper corner.
    pub hi: [f64; 3],
}

/// Chebyshev differentiation matrix on `x_i = -cos(π i / n)` scaled to an interval of length `len`.
#[must_use]
pub fn cheb_diff_matrix(
    n: usize,
    len: f64,
) -> Vec<Vec<f64>> {
    let x: Vec<f64> = (0..=n).map(|i| -(PI * i as f64 / n as f64).cos()).collect();
    let c = |i: usize| if i == 0 || i == n { 2.0 } else { 1.0 };
    let mut d = vec![vec![0.0; n + 1]; n + 1];
    for i in 0..=n {
        for j in 0..=n {
            if i != j {
                let s = if (i + j) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                d[i][j] = c(i) / c(j) * s / (x[i] - x[j]);
            }
        }
    }
    for i in 0..=n {
        d[i][i] = -(0..=n).filter(|&j| j != i).map(|j| d[i][j]).sum::<f64>();
    }
    for r in &mut d {
        for v in r.iter_mut() {
            *v *= 2.0 / len;
        }
    }
    d
}

struct ChebOperator {
    lap: Csr,
    bc_rows: Csr,
    boundary: Vec<Option<(BcValue, [f64; 3])>>,
}

impl ChebGrid {
    /// Creates a grid with the given degrees and bounds.
    ///
    /// # Errors
    /// For bad dimensions or degrees below 2.
    pub fn new(
        deg: &[usize],
        lo: &[f64],
        hi: &[f64],
    ) -> Result<Self, String> {
        let dim = deg.len();
        if !(1..=3).contains(&dim)
            || lo.len() != dim
            || hi.len() != dim
            || deg.iter().any(|&d| d < 2)
        {
            return Err("invalid Chebyshev grid".into());
        }
        let mut g = Self {
            dim,
            deg: [0; 3],
            lo: [0.0; 3],
            hi: [1.0; 3],
        };
        g.deg[..dim].copy_from_slice(&deg[..dim]);
        g.lo[..dim].copy_from_slice(&lo[..dim]);
        g.hi[..dim].copy_from_slice(&hi[..dim]);
        Ok(g)
    }

    /// Number of unknowns.
    #[must_use]
    pub fn size(&self) -> usize {
        (0..3)
            .map(|a| {
                if a < self.dim {
                    self.deg[a] + 1
                } else {
                    1
                }
            })
            .product()
    }

    const fn ext(
        &self,
        a: usize,
    ) -> usize {
        if a < self.dim {
            self.deg[a] + 1
        } else {
            1
        }
    }

    const fn ijk(
        &self,
        p: usize,
    ) -> [usize; 3] {
        let (n0, n1) = (self.ext(0), self.ext(1));
        [p % n0, (p / n0) % n1, p / (n0 * n1)]
    }

    const fn index(
        &self,
        i: [usize; 3],
    ) -> usize {
        i[0] + self.ext(0) * (i[1] + self.ext(1) * i[2])
    }

    /// Node coordinate along an axis.
    #[must_use]
    pub fn node(
        &self,
        axis: usize,
        i: usize,
    ) -> f64 {
        let n = self.deg[axis] as f64;
        self.lo[axis] + 0.5 * (self.hi[axis] - self.lo[axis]) * (1.0 - (PI * i as f64 / n).cos())
    }

    /// Coordinates of flat index `p`.
    #[must_use]
    pub fn coord(
        &self,
        p: usize,
    ) -> [f64; 3] {
        let i = self.ijk(p);
        let mut x = [0.0; 3];
        for a in 0..self.dim {
            x[a] = self.node(a, i[a]);
        }
        x
    }

    /// Samples a function at all nodes.
    pub fn sample(
        &self,
        f: impl Fn([f64; 3]) -> f64,
    ) -> Vec<f64> {
        (0..self.size()).map(|p| f(self.coord(p))).collect()
    }

    /// Largest error against an exact function.
    pub fn max_error(
        &self,
        u: &[f64],
        exact: impl Fn([f64; 3]) -> f64,
    ) -> f64 {
        u.iter()
            .enumerate()
            .map(|(p, v)| (v - exact(self.coord(p))).abs())
            .fold(0.0, f64::max)
    }

    fn operator(
        &self,
        bcs: &FaceBcs,
    ) -> Result<ChebOperator, String> {
        let nn = self.size();
        let d1: Vec<Vec<Vec<f64>>> = (0..self.dim)
            .map(|a| cheb_diff_matrix(self.deg[a], self.hi[a] - self.lo[a]))
            .collect();
        let d2: Vec<Vec<Vec<f64>>> = d1
            .iter()
            .map(|d| {
                let n = d.len();
                (0..n)
                    .map(|i| {
                        (0..n)
                            .map(|j| (0..n).map(|k| d[i][k] * d[k][j]).sum())
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let mut lap = Vec::new();
        let mut bc = Vec::new();
        let mut boundary: Vec<Option<(BcValue, [f64; 3])>> = vec![None; nn];
        for p in 0..nn {
            let ijk = self.ijk(p);
            let x = self.coord(p);
            // boundary faces touched by this node
            let mut faces: Vec<Face> = Vec::new();
            for a in 0..self.dim {
                if ijk[a] == 0 {
                    faces.push(Face::of_axis(a, false));
                }
                if ijk[a] == self.deg[a] {
                    faces.push(Face::of_axis(a, true));
                }
            }
            if faces.is_empty() {
                for a in 0..self.dim {
                    for k in 0..=self.deg[a] {
                        let mut q = ijk;
                        q[a] = k;
                        lap.push((p, self.index(q), d2[a][ijk[a]][k]));
                    }
                }
                continue;
            }
            let face = faces
                .iter()
                .copied()
                .find(|f| bcs.get(*f).is_dirichlet())
                .unwrap_or(faces[0]);
            let a = face.axis();
            let sign = if face.is_max() { 1.0 } else { -1.0 };
            match bcs.get(face) {
                | BoundaryCondition::Dirichlet(g) => {
                    bc.push((p, p, 1.0));
                    boundary[p] = Some((g.clone(), x));
                },
                | BoundaryCondition::Neumann(g)
                | BoundaryCondition::Robin { g, alpha: _, beta: _ } => {
                    let (al, be) = match bcs.get(face) {
                        | BoundaryCondition::Robin { alpha, beta, .. } => (*alpha, *beta),
                        | _ => (0.0, 1.0),
                    };
                    bc.push((p, p, al));
                    for k in 0..=self.deg[a] {
                        let mut q = ijk;
                        q[a] = k;
                        bc.push((p, self.index(q), be * sign * d1[a][ijk[a]][k]));
                    }
                    boundary[p] = Some((g.clone(), x));
                },
                | BoundaryCondition::Symmetry => {
                    for k in 0..=self.deg[a] {
                        let mut q = ijk;
                        q[a] = k;
                        bc.push((p, self.index(q), sign * d1[a][ijk[a]][k]));
                    }
                    boundary[p] = Some((BcValue::constant(0.0), x));
                },
                | _ => {
                    return Err(
                        "Chebyshev grids support Dirichlet, Neumann, Robin and symmetry only"
                            .into(),
                    );
                },
            }
        }
        Ok(ChebOperator {
            lap: Csr::from_triplets(nn, nn, &lap),
            bc_rows: Csr::from_triplets(nn, nn, &bc),
            boundary,
        })
    }

    fn assemble(
        op: &ChebOperator,
        shift: f64,
        scale: f64,
    ) -> Csr {
        let nn = op.lap.rows;
        let mut trip = Vec::new();
        for p in 0..nn {
            if op.boundary[p].is_some() {
                for k in op.bc_rows.indptr[p]..op.bc_rows.indptr[p + 1] {
                    trip.push((p, op.bc_rows.indices[k], op.bc_rows.values[k]));
                }
            } else {
                trip.push((p, p, shift));
                for k in op.lap.indptr[p]..op.lap.indptr[p + 1] {
                    trip.push((p, op.lap.indices[k], scale * op.lap.values[k]));
                }
            }
        }
        Csr::from_triplets(nn, nn, &trip)
    }

    /// Solves `-Δu + σ u = f` with the given boundary conditions.
    ///
    /// # Errors
    /// For unsupported boundary types or a singular system.
    pub fn solve_helmholtz(
        &self,
        bcs: &FaceBcs,
        sigma: f64,
        f: impl Fn([f64; 3]) -> f64,
    ) -> Result<Vec<f64>, String> {
        let op = self.operator(bcs)?;
        let a = Self::assemble(&op, sigma, -1.0);
        let lu = a
            .sparse_lu(SparseOrdering::MinDegree)
            .map_err(|e| format!("LU failed: {e:?}"))?;
        let b: Vec<f64> = (0..self.size())
            .map(|p| {
                match &op.boundary[p] {
                    | Some((g, x)) => g.eval(0.0, *x),
                    | None => f(self.coord(p)),
                }
            })
            .collect();
        Ok(lu.solve(&b))
    }

    /// Integrates `u_t = α Δu + s(t, x)` with the θ-method (`θ` in `[1/2, 1]`).
    ///
    /// # Errors
    /// For unsupported boundary types or a singular system.
    #[allow(clippy::too_many_arguments)] // public solver signature (bcs, alpha, init, source, t_end, steps, theta) kept for API compatibility
    pub fn heat(
        &self,
        bcs: &FaceBcs,
        alpha: f64,
        init: impl Fn([f64; 3]) -> f64,
        source: Option<&dyn Fn(f64, [f64; 3]) -> f64>,
        t_end: f64,
        steps: usize,
        theta: f64,
    ) -> Result<Vec<f64>, String> {
        let op = self.operator(bcs)?;
        let dt = t_end / steps as f64;
        let a = Self::assemble(&op, 1.0, -dt * theta * alpha);
        let lu = a
            .sparse_lu(SparseOrdering::MinDegree)
            .map_err(|e| format!("LU failed: {e:?}"))?;
        let nn = self.size();
        let mut u = self.sample(init);
        let mut lu_u = vec![0.0; nn];
        for k in 0..steps {
            let t = k as f64 * dt;
            op.lap.matvec(&u, &mut lu_u);
            let b: Vec<f64> = (0..nn)
                .map(|p| {
                    if let Some((g, x)) = &op.boundary[p] {
                        g.eval(t + dt, *x)
                    } else {
                        let x = self.coord(p);
                        let s0 = source.map_or(0.0, |s| s(t, x));
                        let s1 = source.map_or(0.0, |s| s(t + dt, x));
                        u[p] + dt * (1.0 - theta) * (alpha * lu_u[p] + s0) + dt * theta * s1
                    }
                })
                .collect();
            u = lu.solve(&b);
        }
        Ok(u)
    }
}
