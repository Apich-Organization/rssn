//! Structured finite element, finite volume, DG and spectral simulation entry points.

use crate::sim::fdm_nd::Boundary1D;
use crate::sim::fdm_nd::Boundary2D;
use crate::sim::fdm_nd::Boundary3D;
use crate::sim::fdm_nd::BoundaryCondition;
use crate::sim::fdm_nd::FdmError;

/// Polynomial degree of Lagrange finite elements.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FemDegree {
    /// Piecewise linear (`P1`) Lagrange elements.
    P1,
    /// Piecewise quadratic (`P2`) Lagrange elements.
    P2,
}

impl FemDegree {
    const fn deg(self) -> usize {
        match self {
            | Self::P1 => 1,
            | Self::P2 => 2,
        }
    }
}

/// Solution of a 1D/2D/3D Finite Element Method (FEM) elliptic or parabolic problem.
#[derive(Clone, Debug, PartialEq)]
pub struct FemSolution {
    /// Coordinates `[x, y, z]` of every degree of freedom (vertex and edge midpoint nodes).
    pub dof_coords: Vec<[f64; 3]>,
    /// Nodal solution values `u_h(x_i)` at each degree of freedom.
    pub values: Vec<f64>,
    /// Topological dimension (`1`, `2`, or `3`).
    pub dim: usize,
    /// Element polynomial degree (`P1` or `P2`).
    pub degree: FemDegree,
}

impl FemSolution {
    /// Maximum nodal error `max_i |u_h(x_i) - exact(x_i)|`.
    #[must_use]
    pub fn max_nodal_error<E>(
        &self,
        mut exact: E,
    ) -> f64
    where
        E: FnMut([f64; 3]) -> f64,
    {
        self.dof_coords
            .iter()
            .zip(&self.values)
            .map(|(&pt, &val)| (val - exact(pt)).abs())
            .fold(0.0f64, f64::max)
    }
}

fn to_rssn_bc(bc: &BoundaryCondition) -> Result<crate::sim::boundary::BoundaryCondition, FdmError> {
    use crate::sim::boundary::BcValue;
    use crate::sim::boundary::BoundaryCondition as RssnBc;
    match bc {
        | BoundaryCondition::Dirichlet(v) => Ok(RssnBc::dirichlet(*v)),
        | BoundaryCondition::Neumann(g) => Ok(RssnBc::neumann(*g)),
        | BoundaryCondition::Robin { alpha, beta, gamma } => {
            Ok(RssnBc::Robin {
                alpha: *alpha,
                beta: *beta,
                g: BcValue::constant(*gamma),
            })
        },
        | BoundaryCondition::DirichletProfile(_) | BoundaryCondition::Periodic => {
            Err(FdmError::Domain(
                "Simplex FEM supports Dirichlet(v), Neumann(g), and Robin { alpha, beta, gamma } BCs"
                    .into(),
            ))
        },
    }
}

/// Solve a one-dimensional reaction-diffusion boundary problem by FEM.
///
/// Solves the 1D reaction-diffusion BVP `-(k u_x)_x + c u = f(x)` on `[x0, x1]`
/// with `nx` simplex intervals using `P1` or `P2` Lagrange finite elements
/// (`crate::sim::fem_simplex`) and boundary conditions `bc` (Dirichlet, Neumann, or Robin).
///
/// # Errors
///
/// Returns [`FdmError::Domain`] if the mesh, coefficients, or BCs are invalid.
pub fn solve_fem_1d<F>(
    bounds: (f64, f64),
    nx: usize,
    degree: FemDegree,
    k: f64,
    c: f64,
    source: F,
    bc: &Boundary1D,
) -> Result<FemSolution, FdmError>
where
    F: Fn(f64) -> f64 + Sync,
{
    use crate::sim::boundary::MarkerBcs;
    use crate::sim::fem_simplex::FemProblem;
    use crate::sim::fem_simplex::FemSpace;
    use crate::sim::fem_simplex::SimplexMesh;
    use crate::sim::fem_simplex::marker;
    use crate::sim::fem_simplex::solve_static;
    use crate::sim::structured::LinearSolver;

    if nx == 0 || bounds.1 <= bounds.0 || k <= 0.0 {
        return Err(FdmError::Domain("invalid 1D FEM domain or k <= 0".into()));
    }
    let mesh = SimplexMesh::interval(nx, bounds.0, bounds.1);
    let space = FemSpace::new(mesh, degree.deg()).map_err(FdmError::Domain)?;
    let bcs = MarkerBcs::new()
        .set(marker::XMIN, to_rssn_bc(&bc.left)?)
        .set(marker::XMAX, to_rssn_bc(&bc.right)?);
    let f_fn = |_t: f64, x: [f64; 3]| source(x[0]);
    let problem = FemProblem {
        k,
        c,
        f: &f_fn,
        bcs: &bcs,
    };
    let values =
        solve_static(&space, &problem, 0.0, LinearSolver::Auto).map_err(FdmError::Domain)?;
    Ok(FemSolution {
        dof_coords: space.dofs,
        values,
        dim: 1,
        degree,
    })
}

/// Solve a two-dimensional reaction-diffusion boundary problem by FEM.
///
/// Solves the 2D elliptic PDE `-∇·(k ∇u) + c u = f(x, y)` on `[x0, x1] × [y0, y1]`
/// triangulated into `2 * nx * ny` triangles using `P1` or `P2` Lagrange finite
/// elements (`crate::sim::fem_simplex`) and per-face boundary conditions `bc`
/// (Dirichlet, Neumann, or Robin).
///
/// # Errors
///
/// Returns [`FdmError::Domain`] if the mesh, coefficients, or BCs are invalid.
// Preserve the established solver API with explicit grid, rates, and callbacks.
#[allow(clippy::too_many_arguments)]
pub fn solve_fem_2d<F>(
    bounds: (f64, f64, f64, f64),
    nx: usize,
    ny: usize,
    degree: FemDegree,
    k: f64,
    c: f64,
    source: F,
    bc: &Boundary2D,
) -> Result<FemSolution, FdmError>
where
    F: Fn(f64, f64) -> f64 + Sync,
{
    use crate::sim::boundary::MarkerBcs;
    use crate::sim::fem_simplex::FemProblem;
    use crate::sim::fem_simplex::FemSpace;
    use crate::sim::fem_simplex::SimplexMesh;
    use crate::sim::fem_simplex::marker;
    use crate::sim::fem_simplex::solve_static;
    use crate::sim::structured::LinearSolver;

    if nx == 0 || ny == 0 || bounds.1 <= bounds.0 || bounds.3 <= bounds.2 || k <= 0.0 {
        return Err(FdmError::Domain("invalid 2D FEM domain or k <= 0".into()));
    }
    let mesh = SimplexMesh::rectangle(nx, ny, [bounds.0, bounds.2], [bounds.1, bounds.3]);
    let space = FemSpace::new(mesh, degree.deg()).map_err(FdmError::Domain)?;
    let bcs = MarkerBcs::new()
        .set(marker::XMIN, to_rssn_bc(&bc.left)?)
        .set(marker::XMAX, to_rssn_bc(&bc.right)?)
        .set(marker::YMIN, to_rssn_bc(&bc.bottom)?)
        .set(marker::YMAX, to_rssn_bc(&bc.top)?);
    let f_fn = |_t: f64, x: [f64; 3]| source(x[0], x[1]);
    let problem = FemProblem {
        k,
        c,
        f: &f_fn,
        bcs: &bcs,
    };
    let values =
        solve_static(&space, &problem, 0.0, LinearSolver::Auto).map_err(FdmError::Domain)?;
    Ok(FemSolution {
        dof_coords: space.dofs,
        values,
        dim: 2,
        degree,
    })
}

/// Solve a three-dimensional reaction-diffusion boundary problem by FEM.
///
/// Solves the 3D elliptic PDE `-∇·(k ∇u) + c u = f(x, y, z)` on a tetrahedral
/// box mesh `[lo, hi]` with `cells = [nx, ny, nz]` using `P1` or `P2` Lagrange
/// finite elements (`crate::sim::fem_simplex`) and per-face boundary conditions `bc`.
///
/// # Errors
///
/// Returns [`FdmError::Domain`] if the mesh, coefficients, or BCs are invalid.
// Preserve the established solver API with explicit grid, rates, and callbacks.
#[allow(clippy::too_many_arguments)]
pub fn solve_fem_3d<F>(
    lo: [f64; 3],
    hi: [f64; 3],
    cells: [usize; 3],
    degree: FemDegree,
    k: f64,
    c: f64,
    source: F,
    bc: &Boundary3D,
) -> Result<FemSolution, FdmError>
where
    F: Fn(f64, f64, f64) -> f64 + Sync,
{
    use crate::sim::boundary::MarkerBcs;
    use crate::sim::fem_simplex::FemProblem;
    use crate::sim::fem_simplex::FemSpace;
    use crate::sim::fem_simplex::SimplexMesh;
    use crate::sim::fem_simplex::marker;
    use crate::sim::fem_simplex::solve_static;
    use crate::sim::structured::LinearSolver;

    if cells.contains(&0)
        || (0..3).any(|a| hi[a] <= lo[a])
        || k <= 0.0
        || !k.is_finite()
        || !c.is_finite()
    {
        return Err(FdmError::Domain("invalid 3D FEM domain or k <= 0".into()));
    }
    let mesh = SimplexMesh::boxed(cells, lo, hi);
    let space = FemSpace::new(mesh, degree.deg()).map_err(FdmError::Domain)?;
    let bcs = MarkerBcs::new()
        .set(marker::XMIN, to_rssn_bc(&bc.x_min)?)
        .set(marker::XMAX, to_rssn_bc(&bc.x_max)?)
        .set(marker::YMIN, to_rssn_bc(&bc.y_min)?)
        .set(marker::YMAX, to_rssn_bc(&bc.y_max)?)
        .set(marker::ZMIN, to_rssn_bc(&bc.z_min)?)
        .set(marker::ZMAX, to_rssn_bc(&bc.z_max)?);
    let f_fn = |_t: f64, x: [f64; 3]| source(x[0], x[1], x[2]);
    let problem = FemProblem {
        k,
        c,
        f: &f_fn,
        bcs: &bcs,
    };
    let values =
        solve_static(&space, &problem, 0.0, LinearSolver::Auto).map_err(FdmError::Domain)?;
    Ok(FemSolution {
        dof_coords: space.dofs,
        values,
        dim: 3,
        degree,
    })
}

/// Solves a 2D elliptic PDE `-∇·(k ∇u) + c u = f(x, y)` on a triangle mesh with
/// adaptive Dörfler residual-based mesh refinement (`crate::sim::fem_triangles`).
///
/// # Errors
///
/// Returns [`FdmError::Domain`] if the mesh or parameters are invalid.
// Preserve the established solver API with explicit grid, rates, and callbacks.
#[allow(clippy::too_many_arguments)]
pub fn solve_fem_adaptive_2d<F, G>(
    bounds: (f64, f64, f64, f64),
    nx: usize,
    ny: usize,
    k: f64,
    c: f64,
    source: F,
    dirichlet_bc: G,
    refinement_levels: usize,
    dorfler_theta: f64,
) -> Result<FemSolution, FdmError>
where
    F: Fn(f64, f64) -> f64,
    G: Fn(f64, f64) -> f64,
{
    use crate::sim::fem_triangles::Problem;
    use crate::sim::fem_triangles::TriMesh;
    use crate::sim::fem_triangles::solve_adaptive;

    if nx == 0 || ny == 0 || bounds.1 <= bounds.0 || bounds.3 <= bounds.2 || k <= 0.0 {
        return Err(FdmError::Domain("invalid 2D adaptive FEM domain".into()));
    }
    let mesh = TriMesh::rectangle(bounds.0, bounds.1, bounds.2, bounds.3, nx, ny);
    let k_fn = |_: f64, _: f64| k;
    let c_fn = |_: f64, _: f64| c;
    let prob = Problem {
        k: &k_fn,
        c: &c_fn,
        f: &source,
        g: &dirichlet_bc,
    };
    let (_final_mesh, sol, _etas) = solve_adaptive(
        &mesh,
        &prob,
        refinement_levels,
        dorfler_theta.clamp(0.1, 0.95),
    )
    .map_err(FdmError::Domain)?;
    let dof_coords = sol.points.into_iter().map(|[x, y]| [x, y, 0.0]).collect();
    Ok(FemSolution {
        dof_coords,
        values: sol.values,
        dim: 2,
        degree: FemDegree::P1,
    })
}

// ---------------------------------------------------------------------------
// Finite Volume Method (FVM) and Discontinuous Galerkin (DG) Solvers
// ---------------------------------------------------------------------------

/// Slope limiter / numerical flux scheme for Finite Volume conservation laws.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FvmScheme {
    /// First-order upwind / Godunov flux.
    Upwind,
    /// First-order local Lax–Friedrichs (Rusanov) flux.
    LaxFriedrichs,
    /// Second-order MUSCL-Hancock reconstruction with TVD `minmod` slope limiter.
    MusclMinmod,
    /// Second-order MUSCL-Hancock reconstruction with TVD `van_leer` slope limiter.
    MusclVanLeer,
}

/// Solve a scalar conservation law with finite volumes.
///
/// Solves a 1D scalar conservation law `u_t + f(u)_x = 0` on `num_cells` finite
/// volumes covering `[0, domain_size]` using `scheme` (`Upwind`, `LaxFriedrichs`,
/// `MusclMinmod`, or `MusclVanLeer`, leveraging `crate::sim::physics_fvm` limiters
/// and fluxes) with periodic or transmissive boundary conditions.
///
/// Returns `(cell_centers, cell_averages)` at `t = steps * dt`.
///
/// # Errors
///
/// Returns [`FdmError::Domain`] for invalid grid, time step, or non-finite fluxes.
// Preserve the established solver API with explicit grid, rates, and callbacks.
#[allow(clippy::too_many_arguments)]
pub fn solve_fvm_conservation_1d<IC, FL, DFL>(
    num_cells: usize,
    domain_size: f64,
    dt: f64,
    steps: usize,
    periodic: bool,
    scheme: FvmScheme,
    initial: IC,
    flux: FL,
    max_wave_speed: DFL,
) -> Result<(Vec<f64>, Vec<f64>), FdmError>
where
    IC: Fn(f64) -> f64,
    FL: Fn(f64) -> f64,
    DFL: Fn(f64) -> f64,
{
    use crate::sim::physics_fvm::minmod;
    use crate::sim::physics_fvm::van_leer;

    if num_cells < 3
        || domain_size <= 0.0
        || dt <= 0.0
        || !domain_size.is_finite()
        || !dt.is_finite()
    {
        return Err(FdmError::Domain("invalid 1D FVM parameters".into()));
    }
    let dx = domain_size / num_cells as f64;
    let centers: Vec<f64> = (0..num_cells).map(|i| (i as f64 + 0.5) * dx).collect();
    let mut u: Vec<f64> = centers.iter().map(|&x| initial(x)).collect();
    if !u.iter().all(|v| v.is_finite()) {
        return Err(FdmError::Domain("non-finite initial FVM state".into()));
    }

    let get_cell = |state: &[f64], idx: isize| -> f64 {
        let n = state.len() as isize;
        if periodic {
            state[((idx % n + n) % n) as usize]
        } else {
            state[idx.clamp(0, n - 1) as usize]
        }
    };

    let num_flux = |ul: f64, ur: f64| -> f64 {
        let fl = flux(ul);
        let fr = flux(ur);
        let alpha = max_wave_speed(ul).abs().max(max_wave_speed(ur).abs());
        match scheme {
            | FvmScheme::LaxFriedrichs => f64::midpoint(fl, fr) - 0.5 * (dx / dt) * (ur - ul),
            | FvmScheme::Upwind | FvmScheme::MusclMinmod | FvmScheme::MusclVanLeer => {
                f64::midpoint(fl, fr) - 0.5 * alpha * (ur - ul)
            },
        }
    };

    for _ in 0..steps {
        // Reconstruct left and right cell-interface states at every cell i
        let mut u_left_edge = vec![0.0; num_cells];
        let mut u_right_edge = vec![0.0; num_cells];
        for i in 0..num_cells {
            let idx = i as isize;
            let um1 = get_cell(&u, idx - 1);
            let u0 = u[i];
            let up1 = get_cell(&u, idx + 1);
            let slope = match scheme {
                | FvmScheme::Upwind | FvmScheme::LaxFriedrichs => 0.0,
                | FvmScheme::MusclMinmod => minmod(u0 - um1, up1 - u0),
                | FvmScheme::MusclVanLeer => van_leer(u0 - um1, up1 - u0),
            };
            let ul = u0 - 0.5 * slope;
            let ur = u0 + 0.5 * slope;
            // MUSCL-Hancock half-step evolution within cell
            let dflux = if matches!(scheme, FvmScheme::MusclMinmod | FvmScheme::MusclVanLeer) {
                0.5 * (dt / dx) * (flux(ur) - flux(ul))
            } else {
                0.0
            };
            u_left_edge[i] = ul - dflux;
            u_right_edge[i] = ur - dflux;
        }

        // Compute interface fluxes F_{i+1/2} for i = -1..num_cells-1
        let mut next = vec![0.0; num_cells];
        for i in 0..num_cells {
            let im1 = if i == 0 {
                if periodic { num_cells - 1 } else { 0 }
            } else {
                i - 1
            };
            let ip1 = if i + 1 == num_cells {
                if periodic { 0 } else { num_cells - 1 }
            } else {
                i + 1
            };
            let f_minus = num_flux(u_right_edge[im1], u_left_edge[i]);
            let f_plus = num_flux(u_right_edge[i], u_left_edge[ip1]);
            next[i] = u[i] - (dt / dx) * (f_plus - f_minus);
        }
        if !next.iter().all(|v| v.is_finite()) {
            return Err(FdmError::NotConverged);
        }
        u = next;
    }

    Ok((centers, u))
}

/// Solves the 1D inviscid Burgers' equation `u_t + (u^2 / 2)_x = 0` using
/// `crate::sim::physics_fvm::solve_burgers_1d` or high-resolution MUSCL-Hancock FVM.
///
/// # Errors
///
/// Returns [`FdmError::Domain`] for invalid grid or step parameters.
pub fn solve_fvm_burgers_1d<IC>(
    num_cells: usize,
    domain_size: f64,
    dt: f64,
    steps: usize,
    periodic: bool,
    scheme: FvmScheme,
    initial: IC,
) -> Result<(Vec<f64>, Vec<f64>), FdmError>
where
    IC: Fn(f64) -> f64,
{
    if !periodic && scheme == FvmScheme::LaxFriedrichs {
        let mut mesh = crate::sim::physics_fvm::Mesh::new(num_cells, domain_size, &initial);
        let dx = mesh.dx;
        let centers: Vec<f64> = (0..num_cells).map(|i| (i as f64 + 0.5) * dx).collect();
        let vals = crate::sim::physics_fvm::solve_burgers_1d(&mut mesh, dt, steps);
        return Ok((centers, vals));
    }
    solve_fvm_conservation_1d(
        num_cells,
        domain_size,
        dt,
        steps,
        periodic,
        scheme,
        initial,
        |u| 0.5 * u * u,
        f64::abs,
    )
}

/// Solves the 1D Shallow Water Equations (`h_t + (hu)_x = 0`, `(hu)_t + (hu^2 + g h^2 / 2)_x = 0`)
/// using `crate::sim::physics_fvm::solve_shallow_water_1d`.
///
/// Returns `(height_h, momentum_hu)` at `t = steps * dt`.
///
/// # Errors
///
/// [`FdmError::ShapeMismatch`] if `initial_h.len() != initial_hu.len()`,
/// [`FdmError::Domain`] for non-positive parameters.
pub fn solve_fvm_shallow_water_1d(
    initial_h: &[f64],
    initial_hu: &[f64],
    dx: f64,
    dt: f64,
    steps: usize,
    gravity: f64,
) -> Result<(Vec<f64>, Vec<f64>), FdmError> {
    if initial_h.is_empty() || initial_h.len() != initial_hu.len() {
        return Err(FdmError::ShapeMismatch(
            "initial_h and initial_hu must have the same non-zero length".into(),
        ));
    }
    if dx <= 0.0 || dt <= 0.0 || gravity <= 0.0 {
        return Err(FdmError::Domain(
            "shallow water dx, dt, and gravity must be positive".into(),
        ));
    }
    let states = crate::sim::physics_fvm::solve_shallow_water_1d(
        initial_h.to_vec(),
        initial_hu.to_vec(),
        dx,
        dt,
        steps,
        gravity,
    );
    let h_out: Vec<f64> = states.iter().map(|s| s.h).collect();
    let hu_out: Vec<f64> = states.iter().map(|s| s.hu).collect();
    Ok((h_out, hu_out))
}

/// Solve periodic linear advection with discontinuous Galerkin elements.
///
/// Solves 1D linear advection `u_t + a u_x = 0` using Nodal Discontinuous Galerkin
/// (`crate::sim::dg::solve_dg_1d`) of polynomial degree `poly_degree >= 1` on `n_elem`
/// elements over `[x0, x1]` with periodic boundary conditions.
///
/// Returns `(lgl_nodes_x, nodal_values)` at `t = t_end`.
///
/// # Errors
///
/// Returns [`FdmError::Domain`] if parameters are invalid.
pub fn solve_dg_advection_1d<IC>(
    bounds: (f64, f64),
    n_elem: usize,
    poly_degree: usize,
    velocity: f64,
    t_end: f64,
    cfl: f64,
    initial: IC,
) -> Result<(Vec<f64>, Vec<f64>), FdmError>
where
    IC: Fn(f64) -> f64,
{
    use crate::sim::dg::DgBc;
    use crate::sim::dg::DgMesh;
    use crate::sim::dg::LinearAdvection;
    use crate::sim::dg::solve_dg;

    if n_elem == 0 || poly_degree == 0 || bounds.1 <= bounds.0 || t_end < 0.0 || cfl <= 0.0 {
        return Err(FdmError::Domain("invalid 1D DG parameters".into()));
    }
    let mesh = DgMesh::line(n_elem, bounds.0, bounds.1, poly_degree);
    let law = LinearAdvection { a: [velocity, 0.0] };
    let bcs = [
        DgBc::Periodic,
        DgBc::Periodic,
        DgBc::Periodic,
        DgBc::Periodic,
    ];
    let init_fn = |x: [f64; 2]| vec![initial(x[0])];
    let sol = solve_dg(&mesh, &law, &bcs, &init_fn, t_end, cfl.min(0.9), None)
        .map_err(FdmError::Domain)?;
    let np = mesh.nodes_per_element();
    let mut xs = Vec::with_capacity(n_elem * np);
    for e in 0..n_elem {
        for i in 0..np {
            xs.push(sol.coord(e, i)[0]);
        }
    }
    Ok((xs, sol.q))
}

// ---------------------------------------------------------------------------
// Spectral Solvers (Fourier Pseudospectral & Chebyshev Collocation)
// ---------------------------------------------------------------------------

/// Solve a periodic Helmholtz problem by Fourier inversion.
///
/// Solves the periodic Helmholtz / screened Poisson equation `-Δu + σ u = f`
/// in 1D, 2D, or 3D on a periodic box of lengths `lengths` and power-of-two grid
/// sizes `extents` using Fourier pseudospectral inversion (`crate::sim::spectral::SpectralGrid`).
///
/// For `σ == 0.0` (pure periodic Poisson `-Δu = f`), the source `f` must have zero mean
/// and the returned solution has zero mean.
///
/// # Errors
///
/// [`FdmError::Domain`] if extents are not powers of two (`>= 4`) or if `σ == 0`
/// with non-zero mean source, [`FdmError::ShapeMismatch`] if `source.len()` does not match.
pub fn solve_spectral_helmholtz(
    extents: &[usize],
    lengths: &[f64],
    sigma: f64,
    source: &[f64],
) -> Result<Vec<f64>, FdmError> {
    use crate::sim::spectral::SpectralGrid;

    let grid = SpectralGrid::new(extents, lengths).map_err(FdmError::Domain)?;
    if source.len() != grid.size() {
        return Err(FdmError::ShapeMismatch(format!(
            "spectral grid of size {} given {} source values",
            grid.size(),
            source.len()
        )));
    }
    grid.helmholtz(sigma, source).map_err(FdmError::Domain)
}

/// Integrate a periodic reaction-diffusion equation by ETDRK4.
///
/// Integrates a periodic reaction-diffusion or heat equation
/// `u_t = \alpha Δu - \sigma u + N(t, u)` in 1D, 2D, or 3D on a periodic Fourier
/// pseudospectral grid using 4th-order Exponential Time Differencing Runge–Kutta
/// (`crate::sim::spectral::etdrk4`) with 2/3 dealiasing on the nonlinear term.
///
/// # Errors
///
/// [`FdmError::Domain`] for invalid grid extents or time span,
/// [`FdmError::ShapeMismatch`] if `initial.len()` does not match the grid size.
// Preserve the established solver API with explicit grid, rates, and callbacks.
#[allow(clippy::too_many_arguments)]
pub fn solve_spectral_etdrk4<NL>(
    extents: &[usize],
    lengths: &[f64],
    alpha: f64,
    sigma: f64,
    t_end: f64,
    steps: usize,
    initial: &[f64],
    nonlinear_phys: NL,
) -> Result<Vec<f64>, FdmError>
where
    NL: Fn(f64, &[f64]) -> Vec<f64> + 'static,
{
    use crate::sim::spectral::SpectralGrid;
    use crate::sim::spectral::etdrk4;

    let grid = SpectralGrid::new(extents, lengths).map_err(FdmError::Domain)?;
    let n_total = grid.size();
    if initial.len() != n_total {
        return Err(FdmError::ShapeMismatch(format!(
            "spectral grid of size {n_total} given {} initial values",
            initial.len()
        )));
    }
    if !alpha.is_finite() || alpha < 0.0 || !t_end.is_finite() || t_end < 0.0 || steps == 0 {
        return Err(FdmError::Domain("invalid ETDRK4 parameters".into()));
    }
    let uh0 = grid.to_spec(initial);
    let zero_c = uh0[0] * 0.0;
    let lin: Vec<_> = (0..n_total)
        .map(|p| zero_c + (-alpha * grid.k2(p) - sigma))
        .collect();
    let mask = grid.dealias_mask();
    let grid_out = grid.clone();
    let nl_spec = move |t: f64, uh: &[_]| {
        let mut uh_dealiased = uh.to_vec();
        for (v, &keep) in uh_dealiased.iter_mut().zip(&mask) {
            if !keep {
                *v = zero_c;
            }
        }
        let u_phys = grid.to_phys(&uh_dealiased);
        let n_phys = nonlinear_phys(t, &u_phys);
        let mut nh = grid.to_spec(&n_phys);
        for (v, &keep) in nh.iter_mut().zip(&mask) {
            if !keep {
                *v = zero_c;
            }
        }
        nh
    };
    let uh_end = etdrk4(uh0, &lin, &nl_spec, t_end, steps);
    Ok(grid_out.to_phys(&uh_end))
}

/// Solve a bounded Helmholtz problem by Chebyshev collocation.
///
/// Solves the 1D/2D/3D Helmholtz / Poisson equation `-Δu + σ u = f(x)` on a bounded
/// box `[lo, hi]` using Chebyshev–Gauss–Lobatto spectral collocation
/// (`crate::sim::spectral::ChebGrid::solve_helmholtz`) with Dirichlet, Neumann, or
/// Robin boundary conditions.
///
/// Returns `(chebyshev_node_coords, nodal_values)`.
///
/// # Errors
///
/// Returns [`FdmError::Domain`] for invalid degrees (`< 2`), bounds, or BCs.
pub fn solve_chebyshev_helmholtz<F>(
    degrees: &[usize],
    lo: &[f64],
    hi: &[f64],
    sigma: f64,
    source: F,
    bc: &BoundaryCondition,
) -> Result<(Vec<[f64; 3]>, Vec<f64>), FdmError>
where
    F: Fn([f64; 3]) -> f64,
{
    use crate::sim::boundary::FaceBcs;
    use crate::sim::spectral::ChebGrid;

    let grid = ChebGrid::new(degrees, lo, hi).map_err(FdmError::Domain)?;
    let rssn_bc = to_rssn_bc(bc)?;
    let face_bcs = FaceBcs::uniform(&rssn_bc);
    let values = grid
        .solve_helmholtz(&face_bcs, sigma, &source)
        .map_err(FdmError::Domain)?;
    let coords = (0..grid.size()).map(|p| grid.coord(p)).collect();
    Ok((coords, values))
}
