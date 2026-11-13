//! Three-dimensional finite-difference solvers and the dimensions the base
//! module lacked: 1D Poisson/Laplace and 2D/3D periodic advection.
//!
//! Everything here reuses the discrete [`Laplacian`](super::Laplacian) of the
//! 1D and 2D solvers, so the boundary treatment (Dirichlet nodes held fixed,
//! Neumann faces closed with a ghost node) and the time integrators
//! (Crank-Nicolson, leapfrog) are the same code in every dimension.
//!
//! Faces of a 3D grid are ordered `[x-, x+, y-, y+, z-, z+]`. The boundary
//! function of a face receives the two coordinates tangential to it, in axis
//! order: `(y, z)` on an `x` face, `(x, z)` on a `y` face and `(x, y)` on a
//! `z` face.

use super::AdvectionScheme;
use super::Boundary;
use super::Field1D;
use super::Grid1D;
use super::Laplacian;
use super::PdeError;
use super::PoissonSolver;
use super::axis_stencil;
use super::cfl_check;
use super::check_finite;
use super::check_steps;
use super::crank_nicolson;
use super::leapfrog;
use super::solve_poisson_system;
use crate::kernels::krylov::Csr;

/// Uniform 3D grid on `[x0, x1] x [y0, y1] x [z0, z1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grid3D {
    /// Lower `x` end.
    pub x0: f64,
    /// Upper `x` end.
    pub x1: f64,
    /// Lower `y` end.
    pub y0: f64,
    /// Upper `y` end.
    pub y1: f64,
    /// Lower `z` end.
    pub z0: f64,
    /// Upper `z` end.
    pub z1: f64,
    /// Nodes in `x`.
    pub nx: usize,
    /// Nodes in `y`.
    pub ny: usize,
    /// Nodes in `z`.
    pub nz: usize,
}

impl Grid3D {
    /// Unit cube with `n` nodes per axis.
    #[must_use]
    pub const fn unit_cube(n: usize) -> Self {
        Self {
            x0: 0.0,
            x1: 1.0,
            y0: 0.0,
            y1: 1.0,
            z0: 0.0,
            z1: 1.0,
            nx: n,
            ny: n,
            nz: n,
        }
    }

    /// Spacing in `x`.
    #[must_use]
    pub fn dx(&self) -> f64 {
        (self.x1 - self.x0) / (self.nx.max(2) - 1) as f64
    }

    /// Spacing in `y`.
    #[must_use]
    pub fn dy(&self) -> f64 {
        (self.y1 - self.y0) / (self.ny.max(2) - 1) as f64
    }

    /// Spacing in `z`.
    #[must_use]
    pub fn dz(&self) -> f64 {
        (self.z1 - self.z0) / (self.nz.max(2) - 1) as f64
    }

    /// Number of nodes.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Whether the grid has no nodes.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub(super) fn validate(&self) -> Result<(), PdeError> {
        let ordered = |a: f64, b: f64| a.is_finite() && b.is_finite() && b > a;
        if self.nx < 3
            || self.ny < 3
            || self.nz < 3
            || !ordered(self.x0, self.x1)
            || !ordered(self.y0, self.y1)
            || !ordered(self.z0, self.z1)
        {
            return Err(PdeError::InvalidInput(
                "3D grid needs >= 3 nodes per axis and positive extents",
            ));
        }
        Ok(())
    }
}

/// A scalar field on a 3D grid, stored with `x` fastest
/// (`u[(k * ny + j) * nx + i]`).
#[derive(Debug, Clone, PartialEq)]
pub struct Field3D {
    /// Nodes in `x`.
    pub nx: usize,
    /// Nodes in `y`.
    pub ny: usize,
    /// Nodes in `z`.
    pub nz: usize,
    /// Lower `x` end.
    pub x0: f64,
    /// Lower `y` end.
    pub y0: f64,
    /// Lower `z` end.
    pub z0: f64,
    /// Spacing in `x`.
    pub dx: f64,
    /// Spacing in `y`.
    pub dy: f64,
    /// Spacing in `z`.
    pub dz: f64,
    /// Values at the nodes.
    pub u: Vec<f64>,
    /// Time of the snapshot (`0` for stationary problems).
    pub t: f64,
}

impl Field3D {
    /// Value at node `(i, j, k)`.
    #[must_use]
    pub fn at(
        &self,
        i: usize,
        j: usize,
        k: usize,
    ) -> f64 {
        self.u[(k * self.ny + j) * self.nx + i]
    }

    /// Coordinates of node `(i, j, k)`.
    #[must_use]
    pub fn point(
        &self,
        i: usize,
        j: usize,
        k: usize,
    ) -> (f64, f64, f64) {
        (
            self.x0 + i as f64 * self.dx,
            self.y0 + j as f64 * self.dy,
            self.z0 + k as f64 * self.dz,
        )
    }

    /// Maximum absolute difference to `exact(x, y, z)` over all nodes.
    pub fn max_error(
        &self,
        exact: impl Fn(f64, f64, f64) -> f64,
    ) -> f64 {
        let mut m = 0.0_f64;
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let (x, y, z) = self.point(i, j, k);
                    m = m.max((self.at(i, j, k) - exact(x, y, z)).abs());
                }
            }
        }
        m
    }
}

fn field_3d(
    g: &Grid3D,
    u: Vec<f64>,
    t: f64,
) -> Field3D {
    Field3D {
        nx: g.nx,
        ny: g.ny,
        nz: g.nz,
        x0: g.x0,
        y0: g.y0,
        z0: g.z0,
        dx: g.dx(),
        dy: g.dy(),
        dz: g.dz(),
        u,
        t,
    }
}

/// Which face (if any) of a node carries a Dirichlet value, and its value.
fn dirichlet_value(
    g: &Grid3D,
    bc: &[Boundary; 6],
    (i, j, k): (usize, usize, usize),
    (x, y, z): (f64, f64, f64),
) -> Option<f64> {
    let faces = [
        (i == 0, 0, (y, z)),
        (i == g.nx - 1, 1, (y, z)),
        (j == 0, 2, (x, z)),
        (j == g.ny - 1, 3, (x, z)),
        (k == 0, 4, (x, y)),
        (k == g.nz - 1, 5, (x, y)),
    ];
    faces
        .into_iter()
        .find(|&(on, f, _)| on && bc[f].is_dirichlet())
        .map(|(_, f, (a, b))| bc[f].value(a, b))
}

/// The 7-point Laplacian of a 3D grid with the boundary treatment of the 1D
/// and 2D solvers.
fn build_laplacian_3d(
    g: &Grid3D,
    bc: &[Boundary; 6],
) -> Laplacian {
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let (hx, hy, hz) = (g.dx(), g.dy(), g.dz());
    let n = g.len();
    let mut trip = Vec::with_capacity(7 * n);
    let mut b = vec![0.0; n];
    let mut dirichlet = vec![None; n];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let id = (k * ny + j) * nx + i;
                let p = (
                    g.x0 + i as f64 * hx,
                    g.y0 + j as f64 * hy,
                    g.z0 + k as f64 * hz,
                );
                if let Some(v) = dirichlet_value(g, bc, (i, j, k), p) {
                    dirichlet[id] = Some(v);
                    continue;
                }
                let (x, y, z) = p;
                let mut centre = 0.0;
                // (index, count, spacing, stride, low face, tangential
                // coordinates of both faces of the axis)
                let axes = [
                    (i, nx, hx, 1, 0, (y, z)),
                    (j, ny, hy, nx, 2, (x, z)),
                    (k, nz, hz, nx * ny, 4, (x, y)),
                ];
                for (at, count, h, stride, face, (ta, tb)) in axes {
                    let glo = if at == 0 {
                        bc[face].value(ta, tb)
                    } else {
                        0.0
                    };
                    let ghi = if at == count - 1 {
                        bc[face + 1].value(ta, tb)
                    } else {
                        0.0
                    };
                    let (prev, c, next, r) =
                        axis_stencil(at, count, h, &bc[face], &bc[face + 1], glo, ghi);
                    centre += c;
                    if at > 0 {
                        trip.push((id, id - stride, prev));
                    }
                    if at + 1 < count {
                        trip.push((id, id + stride, next));
                    }
                    b[id] += r;
                }
                trip.push((id, id, centre));
            }
        }
    }
    Laplacian {
        a: Csr::from_triplets(n, n, &trip),
        b,
        dirichlet,
    }
}

fn sample_3d(
    g: &Grid3D,
    dirichlet: &[Option<f64>],
    f: &dyn Fn(f64, f64, f64) -> f64,
) -> Vec<f64> {
    let (hx, hy, hz) = (g.dx(), g.dy(), g.dz());
    let mut u = vec![0.0; g.len()];
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let id = (k * g.ny + j) * g.nx + i;
                u[id] = dirichlet[id].unwrap_or_else(|| {
                    f(
                        g.x0 + i as f64 * hx,
                        g.y0 + j as f64 * hy,
                        g.z0 + k as f64 * hz,
                    )
                });
            }
        }
    }
    u
}

/// Solves `u_t = alpha (u_xx + u_yy + u_zz)` with Crank-Nicolson;
/// `bc = [x-, x+, y-, y+, z-, z+]`.
///
/// # Errors
/// As [`heat_1d`](super::heat_1d).
pub fn heat_3d(
    alpha: f64,
    grid: &Grid3D,
    bc: &[Boundary; 6],
    u0: &dyn Fn(f64, f64, f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<Field3D, PdeError> {
    grid.validate()?;
    if alpha.is_nan() || alpha <= 0.0 {
        return Err(PdeError::InvalidInput("diffusivity must be positive"));
    }
    let dt = check_steps(t_end, steps)?;
    let lap = build_laplacian_3d(grid, bc);
    let u = crank_nicolson(&lap, alpha, dt, steps, sample_3d(grid, &lap.dirichlet, u0))?;
    Ok(field_3d(grid, u, t_end))
}

/// Solves `u_tt = c^2 (u_xx + u_yy + u_zz)` with leapfrog; stable for
/// `c dt sqrt(1/dx^2 + 1/dy^2 + 1/dz^2) <= 1`.
///
/// # Errors
/// As [`wave_1d`](super::wave_1d).
pub fn wave_3d(
    c: f64,
    grid: &Grid3D,
    bc: &[Boundary; 6],
    u0: &dyn Fn(f64, f64, f64) -> f64,
    v0: &dyn Fn(f64, f64, f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<Field3D, PdeError> {
    grid.validate()?;
    if c.is_nan() || c <= 0.0 {
        return Err(PdeError::InvalidInput("wave speed must be positive"));
    }
    let dt = check_steps(t_end, steps)?;
    let (hx, hy, hz) = (grid.dx(), grid.dy(), grid.dz());
    cfl_check(
        c * dt * (1.0 / (hx * hx) + 1.0 / (hy * hy) + 1.0 / (hz * hz)).sqrt(),
        1.0,
    )?;
    let lap = build_laplacian_3d(grid, bc);
    let vel = sample_3d(grid, &vec![None; grid.len()], v0);
    let start = sample_3d(grid, &lap.dirichlet, u0);
    let u = leapfrog(&lap, c, dt, steps, &start, &vel)?;
    Ok(field_3d(grid, u, t_end))
}

/// Solves `-(u_xx + u_yy + u_zz) = f` with the 7-point stencil;
/// `bc = [x-, x+, y-, y+, z-, z+]`.
///
/// # Errors
/// As [`poisson_2d`](super::poisson_2d).
pub fn poisson_3d(
    grid: &Grid3D,
    bc: &[Boundary; 6],
    f: &dyn Fn(f64, f64, f64) -> f64,
    solver: PoissonSolver,
) -> Result<Field3D, PdeError> {
    grid.validate()?;
    let n_dir = bc.iter().filter(|b| b.is_dirichlet()).count();
    let lap = build_laplacian_3d(grid, bc);
    let (hx, hy, hz) = (grid.dx(), grid.dy(), grid.dz());
    let (nx, ny) = (grid.nx, grid.ny);
    let rhs: Vec<f64> = (0..grid.len())
        .map(|id| {
            lap.dirichlet[id].unwrap_or_else(|| {
                let (i, j, k) = (id % nx, (id / nx) % ny, id / (nx * ny));
                f(
                    grid.x0 + i as f64 * hx,
                    grid.y0 + j as f64 * hy,
                    grid.z0 + k as f64 * hz,
                ) + lap.b[id]
            })
        })
        .collect();
    let u = solve_poisson_system(&lap, &rhs, n_dir, 6, solver)?;
    Ok(field_3d(grid, u, 0.0))
}

/// Solves `Laplace(u) = 0` on a 3D grid.
///
/// # Errors
/// As [`poisson_3d`].
pub fn laplace_3d(
    grid: &Grid3D,
    bc: &[Boundary; 6],
    solver: PoissonSolver,
) -> Result<Field3D, PdeError> {
    poisson_3d(grid, bc, &|_, _, _| 0.0, solver)
}

/// Solves `-u_xx = f` on a 1D grid; `bc = [left, right]`.
///
/// # Errors
/// As [`poisson_2d`](super::poisson_2d).
pub fn poisson_1d(
    grid: &Grid1D,
    bc: &[Boundary; 2],
    f: &dyn Fn(f64) -> f64,
    solver: PoissonSolver,
) -> Result<Field1D, PdeError> {
    grid.validate(3)?;
    let n_dir = bc.iter().filter(|b| b.is_dirichlet()).count();
    let h = grid.dx();
    let lap = super::build_laplacian(
        &super::LapGrid {
            nx: grid.n,
            ny: 1,
            x0: grid.x0,
            y0: 0.0,
            hx: h,
            hy: 1.0,
            two_d: false,
        },
        &super::bcs_1d(bc),
    );
    let rhs: Vec<f64> = (0..grid.n)
        .map(|i| lap.dirichlet[i].unwrap_or_else(|| f(grid.x0 + i as f64 * h) + lap.b[i]))
        .collect();
    let u = solve_poisson_system(&lap, &rhs, n_dir, 2, solver)?;
    Ok(Field1D {
        x: grid.nodes(),
        u,
        t: 0.0,
    })
}

/// Solves `u_xx = 0` on a 1D grid.
///
/// # Errors
/// As [`poisson_1d`].
pub fn laplace_1d(
    grid: &Grid1D,
    bc: &[Boundary; 2],
    solver: PoissonSolver,
) -> Result<Field1D, PdeError> {
    poisson_1d(grid, bc, &|_| 0.0, solver)
}

// ---------------------------------------------------------------------------
// Periodic advection in 2D and 3D
// ---------------------------------------------------------------------------

/// One periodic 1D sweep of `scheme` along `axis` of an `x`-fastest array
/// with extents `dims`.
fn sweep(
    u: &[f64],
    w: &mut [f64],
    dims: [usize; 3],
    axis: usize,
    nu: f64,
    scheme: AdvectionScheme,
) {
    let stride = [1, dims[0], dims[0] * dims[1]][axis];
    let n = dims[axis];
    for (id, out) in w.iter_mut().enumerate() {
        let at = (id / stride) % n;
        let base = id - at * stride;
        let um = u[base + ((at + n - 1) % n) * stride];
        let up = u[base + ((at + 1) % n) * stride];
        let c = u[id];
        *out = match scheme {
            | AdvectionScheme::Upwind => {
                if nu >= 0.0 {
                    c - nu * (c - um)
                } else {
                    c - nu * (up - c)
                }
            },
            | AdvectionScheme::LaxWendroff => {
                c - 0.5 * nu * (up - um) + 0.5 * nu * nu * (up - 2.0 * c + um)
            },
        };
    }
}

/// Strang-split periodic advection of `u_t + a . grad(u) = 0` on an
/// `x`-fastest array; every axis needs `|a_d dt / h_d| <= 1`.
fn advect(
    u: &mut Vec<f64>,
    dims: [usize; 3],
    h: [f64; 3],
    a: [f64; 3],
    scheme: AdvectionScheme,
    dt: f64,
    steps: usize,
) -> Result<(), PdeError> {
    let mut nus = [0.0; 3];
    for d in 0..3 {
        nus[d] = a[d] * dt / h[d];
        cfl_check(nus[d].abs(), 1.0)?;
    }
    let active: Vec<usize> = (0..3).filter(|&d| dims[d] > 1).collect();
    let mut w = vec![0.0; u.len()];
    for step in 0..steps {
        // Alternate the sweep order so the splitting error stays second order.
        let order: Vec<usize> = if step % 2 == 0 {
            active.clone()
        } else {
            active.iter().rev().copied().collect()
        };
        for axis in order {
            sweep(u, &mut w, dims, axis, nus[axis], scheme);
            std::mem::swap(u, &mut w);
        }
    }
    check_finite(u)
}

/// Solves `u_t + a_x u_x + a_y u_y = 0` on the periodic rectangle
/// `[x0, x1) x [y0, y1)` (the upper ends are identified with the lower ones).
///
/// # Errors
/// [`PdeError::Unstable`] when a directional Courant number exceeds one,
/// [`PdeError::InvalidInput`] or [`PdeError::NonFinite`].
pub fn advection_2d(
    a: [f64; 2],
    grid: &super::Grid2D,
    scheme: AdvectionScheme,
    u0: &dyn Fn(f64, f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<super::Field2D, PdeError> {
    grid.validate()?;
    let dt = check_steps(t_end, steps)?;
    let (nx, ny) = (grid.nx, grid.ny);
    let hx = (grid.x1 - grid.x0) / nx as f64;
    let hy = (grid.y1 - grid.y0) / ny as f64;
    let mut u = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            u.push(u0(grid.x0 + i as f64 * hx, grid.y0 + j as f64 * hy));
        }
    }
    advect(
        &mut u,
        [nx, ny, 1],
        [hx, hy, 1.0],
        [a[0], a[1], 0.0],
        scheme,
        dt,
        steps,
    )?;
    Ok(super::Field2D {
        nx,
        ny,
        x0: grid.x0,
        y0: grid.y0,
        dx: hx,
        dy: hy,
        u,
        t: t_end,
    })
}

/// Solves `u_t + a . grad(u) = 0` on the periodic box
/// `[x0, x1) x [y0, y1) x [z0, z1)`.
///
/// # Errors
/// As [`advection_2d`].
pub fn advection_3d(
    a: [f64; 3],
    grid: &Grid3D,
    scheme: AdvectionScheme,
    u0: &dyn Fn(f64, f64, f64) -> f64,
    t_end: f64,
    steps: usize,
) -> Result<Field3D, PdeError> {
    grid.validate()?;
    let dt = check_steps(t_end, steps)?;
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    let h = [
        (grid.x1 - grid.x0) / nx as f64,
        (grid.y1 - grid.y0) / ny as f64,
        (grid.z1 - grid.z0) / nz as f64,
    ];
    let mut u = Vec::with_capacity(nx * ny * nz);
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                u.push(u0(
                    grid.x0 + i as f64 * h[0],
                    grid.y0 + j as f64 * h[1],
                    grid.z0 + k as f64 * h[2],
                ));
            }
        }
    }
    advect(&mut u, [nx, ny, nz], h, a, scheme, dt, steps)?;
    Ok(Field3D {
        nx,
        ny,
        nz,
        x0: grid.x0,
        y0: grid.y0,
        z0: grid.z0,
        dx: h[0],
        dy: h[1],
        dz: h[2],
        u,
        t: t_end,
    })
}
