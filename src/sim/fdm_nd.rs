//! 1D, 2D, and 3D finite-difference PDE solvers (Poisson, Heat, Wave, Advection)
//! with Dirichlet, Neumann, Robin, and periodic boundary conditions, plus
//! explicit, Crank-Nicolson, and ADI time-stepping.

/// Error returned by 1D/2D/3D finite-difference PDE solvers.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FdmError {
    /// Empty input or zero grid dimension.
    Empty,
    /// Mismatched dimensions or field lengths.
    ShapeMismatch(String),
    /// Invalid domain, spacing, CFL condition, or BC parameter.
    Domain(String),
    /// Linear system was singular.
    Singular,
    /// Iterative solver did not converge.
    NotConverged,
}

impl std::fmt::Display for FdmError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Empty => f.write_str("empty input"),
            | Self::ShapeMismatch(s) => write!(f, "shape mismatch: {s}"),
            | Self::Domain(s) => write!(f, "domain error: {s}"),
            | Self::Singular => f.write_str("singular matrix"),
            | Self::NotConverged => f.write_str("did not converge"),
        }
    }
}

impl std::error::Error for FdmError {}

/// 2D uniform Cartesian grid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    /// Points along x.
    pub nx: usize,
    /// Points along y.
    pub ny: usize,
    /// The spacing.
    pub h: f64,
}

impl Grid {
    fn check(&self) -> Result<(), FdmError> {
        if self.nx < 3 || self.ny < 3 {
            return Err(FdmError::Domain(
                "a PDE grid needs at least three points per axis".to_string(),
            ));
        }
        if self.h <= 0.0 || !self.h.is_finite() {
            return Err(FdmError::Domain(
                "a PDE grid needs a positive spacing".to_string(),
            ));
        }
        Ok(())
    }

    /// The number of points.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.nx * self.ny
    }

    /// True when the grid has no points.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// A field on a grid: `values[y * nx + x]`.
#[derive(Clone, Debug, PartialEq)]
pub struct Field {
    /// The grid it lives on.
    pub grid: Grid,
    /// The samples, row-major.
    pub values: Vec<f64>,
}

/// The Laplacian of a field with zero boundary data, five-point stencil.
fn laplacian(
    grid: Grid,
    values: &[f64],
) -> Vec<f64> {
    let (nx, ny, h) = (grid.nx, grid.ny, grid.h);
    let mut out = vec![0.0; values.len()];
    for y in 1..ny - 1 {
        for x in 1..nx - 1 {
            let index = y * nx + x;
            out[index] =
                (values[index + 1] + values[index - 1] + values[index + nx] + values[index - nx]
                    - 4.0 * values[index])
                    / (h * h);
        }
    }
    out
}

/// A Poisson problem `∇²u = f` with zero boundary data.

pub struct PoissonSpec {
    /// The grid.
    pub grid: Grid,
    /// The source term `f`.
    pub source: Vec<f64>,
    /// Maximum SOR sweeps.
    pub iterations: usize,
    /// The residual tolerance.
    pub tolerance: f64,
}

/// Solves `∇²u = f` by successive over-relaxation.
///
/// # Errors
///
/// [`FdmError::ShapeMismatch`] when the source does not fit the grid, and
/// [`FdmError::NotConverged`] when the sweep limit is reached.
pub fn solve_poisson(spec: &PoissonSpec) -> Result<Field, FdmError> {
    spec.grid.check()?;
    if spec.source.len() != spec.grid.len() {
        return Err(FdmError::ShapeMismatch(format!(
            "a {}x{} grid needs {} source values, found {}",
            spec.grid.ny,
            spec.grid.nx,
            spec.grid.len(),
            spec.source.len()
        )));
    }
    let (nx, ny, h) = (spec.grid.nx, spec.grid.ny, spec.grid.h);
    let mut u = vec![0.0; spec.grid.len()];
    let omega = 1.5;
    let tolerance = if spec.tolerance > 0.0 {
        spec.tolerance
    } else {
        1e-10
    };
    for _ in 0..spec.iterations.max(1) {
        let mut change = 0.0f64;
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let index = y * nx + x;
                let neighbour = (u[index + 1] + u[index - 1] + u[index + nx] + u[index - nx]) / 4.0;
                let target = neighbour - h * h * spec.source[index] / 4.0;
                let updated = u[index] + omega * (target - u[index]);
                change = change.max((updated - u[index]).abs());
                u[index] = updated;
            }
        }
        if change < tolerance {
            return Ok(Field {
                grid: spec.grid,
                values: u,
            });
        }
    }
    Err(FdmError::NotConverged)
}

/// A heat equation `u_t = α ∇²u`.
#[derive(Clone, Debug, PartialEq)]
pub struct HeatSpec {
    /// The grid.
    pub grid: Grid,
    /// The diffusivity.
    pub alpha: f64,
    /// How many time steps to take.
    pub steps: usize,
    /// The time step.
    pub dt: f64,
    /// The initial field.
    pub initial: Vec<f64>,
}

/// Solves the heat equation with an explicit scheme, returning every snapshot
/// including the initial one.
///
/// # Errors
///
/// [`FdmError::ShapeMismatch`] when the initial field does not fit, and
/// [`FdmError::Domain`] when `α dt / h²` exceeds the stability limit.
pub fn solve_heat(spec: &HeatSpec) -> Result<Vec<Field>, FdmError> {
    spec.grid.check()?;
    if spec.initial.len() != spec.grid.len() {
        return Err(FdmError::ShapeMismatch(
            "the initial field does not fit the grid".to_string(),
        ));
    }
    let spacing_squared = spec.grid.h * spec.grid.h;
    let ratio = spec.alpha * spec.dt / spacing_squared;
    if !(ratio.is_finite() && ratio <= 0.25) {
        return Err(FdmError::Domain(format!(
            "the explicit heat scheme needs `alpha*dt/h^2 <= 0.25`, found {ratio}"
        )));
    }
    let mut field = spec.initial.clone();
    let mut out = vec![Field {
        grid: spec.grid,
        values: field.clone(),
    }];
    for _ in 0..spec.steps {
        let lap = laplacian(spec.grid, &field);
        for index in 0..field.len() {
            field[index] += spec.alpha * spec.dt * lap[index];
        }
        out.push(Field {
            grid: spec.grid,
            values: field.clone(),
        });
    }
    Ok(out)
}

/// A wave equation `u_tt = c² ∇²u`.
#[derive(Clone, Debug, PartialEq)]
pub struct WaveSpec {
    /// The grid.
    pub grid: Grid,
    /// The wave speed.
    pub c: f64,
    /// How many time steps to take.
    pub steps: usize,
    /// The time step.
    pub dt: f64,
    /// The initial displacement.
    pub initial: Vec<f64>,
    /// The initial velocity.
    pub velocity: Vec<f64>,
}

/// Solves the wave equation with the leapfrog scheme.
///
/// # Errors
///
/// [`FdmError::ShapeMismatch`] when a field does not fit, and
/// [`FdmError::Domain`] when the CFL condition fails.
pub fn solve_wave(spec: &WaveSpec) -> Result<Vec<Field>, FdmError> {
    spec.grid.check()?;
    if spec.initial.len() != spec.grid.len() || spec.velocity.len() != spec.grid.len() {
        return Err(FdmError::ShapeMismatch(
            "the initial fields do not fit the grid".to_string(),
        ));
    }
    let courant = spec.c * spec.dt / spec.grid.h;
    if !(courant.is_finite() && courant.abs() <= std::f64::consts::FRAC_1_SQRT_2 + 1e-12) {
        return Err(FdmError::Domain(format!(
            "the leapfrog scheme needs `c*dt/h <= 1/sqrt(2)`, found {courant}"
        )));
    }
    let squared = courant * courant;
    let mut previous = spec.initial.clone();
    let mut current = spec.initial.clone();
    // The first step needs the velocity: u^1 = u^0 + dt v^0 + dt²/2 c² ∇²u^0.
    let lap = laplacian(spec.grid, &current);
    for index in 0..current.len() {
        current[index] += spec.dt * spec.velocity[index] + 0.5 * squared * lap[index];
    }
    let mut out = vec![
        Field {
            grid: spec.grid,
            values: previous.clone(),
        },
        Field {
            grid: spec.grid,
            values: current.clone(),
        },
    ];
    for _ in 1..spec.steps {
        let lap = laplacian(spec.grid, &current);
        let mut next = vec![0.0; current.len()];
        for index in 0..current.len() {
            next[index] = 2.0 * current[index] - previous[index] + squared * lap[index];
        }
        previous = current;
        current = next;
        out.push(Field {
            grid: spec.grid,
            values: current.clone(),
        });
    }
    Ok(out)
}

/// An advection equation `u_t + v·∇u = 0`.
#[derive(Clone, Debug, PartialEq)]
pub struct AdvectionSpec {
    /// The grid.
    pub grid: Grid,
    /// The velocity along x.
    pub vx: f64,
    /// The velocity along y.
    pub vy: f64,
    /// How many time steps to take.
    pub steps: usize,
    /// The time step.
    pub dt: f64,
    /// The initial field.
    pub initial: Vec<f64>,
}

/// Solves advection with the first-order upwind scheme.
///
/// # Errors
///
/// [`FdmError::ShapeMismatch`] when the field does not fit, and
/// [`FdmError::Domain`] when the CFL condition fails.
pub fn solve_advection(spec: &AdvectionSpec) -> Result<Vec<Field>, FdmError> {
    spec.grid.check()?;
    if spec.initial.len() != spec.grid.len() {
        return Err(FdmError::ShapeMismatch(
            "the initial field does not fit the grid".to_string(),
        ));
    }
    let (nx, ny, h) = (spec.grid.nx, spec.grid.ny, spec.grid.h);
    let cx = spec.vx * spec.dt / h;
    let cy = spec.vy * spec.dt / h;
    if !(cx.is_finite() && cy.is_finite()) || cx.abs() + cy.abs() > 1.0 {
        return Err(FdmError::Domain(format!(
            "the upwind scheme needs `(|vx|+|vy|)*dt/h <= 1`, found {}",
            cx.abs() + cy.abs()
        )));
    }
    let mut field = spec.initial.clone();
    let mut out = vec![Field {
        grid: spec.grid,
        values: field.clone(),
    }];
    for _ in 0..spec.steps {
        let mut next = field.clone();
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let index = y * nx + x;
                let dx = if cx >= 0.0 {
                    field[index] - field[index - 1]
                } else {
                    field[index + 1] - field[index]
                };
                let dy = if cy >= 0.0 {
                    field[index] - field[index - nx]
                } else {
                    field[index + nx] - field[index]
                };
                next[index] = field[index] - cx * dx - cy * dy;
            }
        }
        field = next;
        out.push(Field {
            grid: spec.grid,
            values: field.clone(),
        });
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Boundary Conditions (Dirichlet, Neumann, Robin, Periodic)
// ---------------------------------------------------------------------------

/// Boundary condition on a domain face.
///
/// The normal derivative `\partial u / \partial n` is always taken along the
/// **outward** unit normal of the boundary face.

#[derive(Clone, Debug, PartialEq)]
pub enum BoundaryCondition {
    /// Constant Dirichlet condition `u = value` on the boundary face
    /// (`Dirichlet(0.0)` is homogeneous zero Dirichlet).
    Dirichlet(f64),
    /// Spatially varying Dirichlet profile `u = values[idx]` along the boundary face.
    DirichletProfile(Vec<f64>),
    /// Constant outward Neumann flux condition `\partial u / \partial n = flux`.
    Neumann(f64),
    /// Constant Robin (mixed/convective) condition
    /// `\alpha u + \beta (\partial u / \partial n) = \gamma`.
    Robin {
        /// Coefficient of `u`.
        alpha: f64,
        /// Coefficient of the outward normal derivative `\partial u / \partial n`.
        beta: f64,
        /// Right-hand-side datum `\gamma`.
        gamma: f64,
    },
    /// Periodic wrap-around across opposite faces of the axis.
    Periodic,
}

impl Default for BoundaryCondition {
    fn default() -> Self {
        Self::Dirichlet(0.0)
    }
}

impl BoundaryCondition {
    /// Homogeneous Dirichlet boundary condition (`u = 0`).
    #[must_use]
    pub const fn zero_dirichlet() -> Self {
        Self::Dirichlet(0.0)
    }

    /// Homogeneous Neumann (insulating / reflective) boundary condition (`\partial u / \partial n = 0`).
    #[must_use]
    pub const fn zero_neumann() -> Self {
        Self::Neumann(0.0)
    }

    fn dirichlet_at(
        &self,
        face_idx: usize,
    ) -> f64 {
        match self {
            | Self::Dirichlet(v) => *v,
            | Self::DirichletProfile(vals) => vals.get(face_idx).copied().unwrap_or(0.0),
            | _ => 0.0,
        }
    }
}

/// Boundary conditions for a 1D interval `[x_min, x_max]`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Boundary1D {
    /// Condition at the left endpoint (`i = 0`).
    pub left: BoundaryCondition,
    /// Condition at the right endpoint (`i = nx - 1`).
    pub right: BoundaryCondition,
}

impl Boundary1D {
    /// Uniform boundary condition on both ends.
    #[must_use]
    pub fn uniform(bc: BoundaryCondition) -> Self {
        Self {
            left: bc.clone(),
            right: bc,
        }
    }

    /// Periodic boundary conditions on `[0, L)`.
    #[must_use]
    pub const fn periodic() -> Self {
        Self {
            left: BoundaryCondition::Periodic,
            right: BoundaryCondition::Periodic,
        }
    }

    fn validate(&self) -> Result<(), FdmError> {
        let lp = matches!(self.left, BoundaryCondition::Periodic);
        let rp = matches!(self.right, BoundaryCondition::Periodic);
        if lp != rp {
            return Err(FdmError::Domain(
                "Periodic boundary condition must be set on both left and right faces".into(),
            ));
        }
        for bc in [&self.left, &self.right] {
            if let BoundaryCondition::Robin { alpha, beta, gamma } = bc
                && (!alpha.is_finite()
                    || !beta.is_finite()
                    || !gamma.is_finite()
                    || (alpha.abs() < 1e-15 && beta.abs() < 1e-15))
            {
                return Err(FdmError::Domain(
                    "invalid Robin boundary condition coefficients".into(),
                ));
            }
        }
        Ok(())
    }
}

/// Boundary conditions for a 2D rectangle `[x_min, x_max] × [y_min, y_max]`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Boundary2D {
    /// Condition at `x = x_min` (`x = 0`).
    pub left: BoundaryCondition,
    /// Condition at `x = x_max` (`x = nx - 1`).
    pub right: BoundaryCondition,
    /// Condition at `y = y_min` (`y = 0`).
    pub bottom: BoundaryCondition,
    /// Condition at `y = y_max` (`y = ny - 1`).
    pub top: BoundaryCondition,
}

impl Boundary2D {
    /// Same boundary condition on all four faces.
    #[must_use]
    pub fn uniform(bc: BoundaryCondition) -> Self {
        Self {
            left: bc.clone(),
            right: bc.clone(),
            bottom: bc.clone(),
            top: bc,
        }
    }

    /// Doubly periodic boundary conditions.
    #[must_use]
    pub const fn periodic() -> Self {
        Self {
            left: BoundaryCondition::Periodic,
            right: BoundaryCondition::Periodic,
            bottom: BoundaryCondition::Periodic,
            top: BoundaryCondition::Periodic,
        }
    }

    fn validate(&self) -> Result<(), FdmError> {
        let xp = matches!(self.left, BoundaryCondition::Periodic);
        if xp != matches!(self.right, BoundaryCondition::Periodic) {
            return Err(FdmError::Domain(
                "Periodic BC in x requires both left and right faces to be Periodic".into(),
            ));
        }
        let yp = matches!(self.bottom, BoundaryCondition::Periodic);
        if yp != matches!(self.top, BoundaryCondition::Periodic) {
            return Err(FdmError::Domain(
                "Periodic BC in y requires both bottom and top faces to be Periodic".into(),
            ));
        }
        Ok(())
    }
}

/// Boundary conditions for a 3D box `[x_min, x_max] × [y_min, y_max] × [z_min, z_max]`.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Boundary3D {
    /// Condition at `x = 0`.
    pub x_min: BoundaryCondition,
    /// Condition at `x = nx - 1`.
    pub x_max: BoundaryCondition,
    /// Condition at `y = 0`.
    pub y_min: BoundaryCondition,
    /// Condition at `y = ny - 1`.
    pub y_max: BoundaryCondition,
    /// Condition at `z = 0`.
    pub z_min: BoundaryCondition,
    /// Condition at `z = nz - 1`.
    pub z_max: BoundaryCondition,
}

impl Boundary3D {
    /// Same boundary condition on all six faces.
    #[must_use]
    pub fn uniform(bc: BoundaryCondition) -> Self {
        Self {
            x_min: bc.clone(),
            x_max: bc.clone(),
            y_min: bc.clone(),
            y_max: bc.clone(),
            z_min: bc.clone(),
            z_max: bc,
        }
    }

    /// Triply periodic boundary conditions.
    #[must_use]
    pub const fn periodic() -> Self {
        Self {
            x_min: BoundaryCondition::Periodic,
            x_max: BoundaryCondition::Periodic,
            y_min: BoundaryCondition::Periodic,
            y_max: BoundaryCondition::Periodic,
            z_min: BoundaryCondition::Periodic,
            z_max: BoundaryCondition::Periodic,
        }
    }

    fn validate(&self) -> Result<(), FdmError> {
        for face in [
            &self.x_min,
            &self.x_max,
            &self.y_min,
            &self.y_max,
            &self.z_min,
            &self.z_max,
        ] {
            Boundary1D::uniform(face.clone()).validate()?;
        }

        if matches!(self.x_min, BoundaryCondition::Periodic)
            != matches!(self.x_max, BoundaryCondition::Periodic)
            || matches!(self.y_min, BoundaryCondition::Periodic)
                != matches!(self.y_max, BoundaryCondition::Periodic)
            || matches!(self.z_min, BoundaryCondition::Periodic)
                != matches!(self.z_max, BoundaryCondition::Periodic)
        {
            return Err(FdmError::Domain(
                "Periodic BC must be paired on opposite faces of each axis".into(),
            ));
        }
        Ok(())
    }
}

/// Time-stepping schemes for parabolic / diffusion PDEs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum TimeStepping {
    /// Forward Euler explicit time-stepping (requires CFL `\alpha \Delta t / h^2 \le 1/(2d)`).
    #[default]
    Explicit,
    /// First-order L-stable Backward Euler implicit time-stepping (unconditionally stable).
    BackwardEuler,
    /// Second-order A-stable Crank–Nicolson trapezoidal implicit time-stepping (unconditionally stable).
    CrankNicolson,
    /// Alternating Direction Implicit (Peaceman–Rachford in 2D, Douglas–Gunn in 3D;
    /// equivalent to Crank–Nicolson in 1D). Unconditionally stable with `O(N)` tridiagonal
    /// line solves per step.
    Adi,
}

// ---------------------------------------------------------------------------
// 1D and 3D Grids and Fields
// ---------------------------------------------------------------------------

/// A regular 1D grid with `nx` points and spacing `h`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid1D {
    /// Number of grid points.
    pub nx: usize,
    /// Grid spacing `h > 0`.
    pub h: f64,
}

impl Grid1D {
    fn check(&self) -> Result<(), FdmError> {
        if self.nx < 3 {
            return Err(FdmError::Domain(
                "a 1D PDE grid needs at least three points".into(),
            ));
        }
        if self.h <= 0.0 || !self.h.is_finite() {
            return Err(FdmError::Domain(
                "a 1D PDE grid needs a positive finite spacing".into(),
            ));
        }
        Ok(())
    }

    /// Number of grid points.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.nx
    }

    /// True if `nx == 0`.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.nx == 0
    }
}

/// A scalar field on a [`Grid1D`].
#[derive(Clone, Debug, PartialEq)]
pub struct Field1D {
    /// The 1D grid.
    pub grid: Grid1D,
    /// Nodal values `values[x]`.
    pub values: Vec<f64>,
}

/// A regular 3D grid with `nx × ny × nz` points and spacing `h`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid3D {
    /// Points along x.
    pub nx: usize,
    /// Points along y.
    pub ny: usize,
    /// Points along z.
    pub nz: usize,
    /// Grid spacing `h > 0`.
    pub h: f64,
}

impl Grid3D {
    fn check(&self) -> Result<(), FdmError> {
        if self.nx < 3 || self.ny < 3 || self.nz < 3 {
            return Err(FdmError::Domain(
                "a 3D PDE grid needs at least three points per axis".into(),
            ));
        }
        if self.h <= 0.0 || !self.h.is_finite() {
            return Err(FdmError::Domain(
                "a 3D PDE grid needs a positive finite spacing".into(),
            ));
        }
        Ok(())
    }

    /// Total number of grid points `nx * ny * nz`.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// True if empty.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Flat index `(z * ny + y) * nx + x`.
    #[must_use]
    pub const fn index(
        &self,
        x: usize,
        y: usize,
        z: usize,
    ) -> usize {
        (z * self.ny + y) * self.nx + x
    }
}

/// A scalar field on a [`Grid3D`]: `values[(z * ny + y) * nx + x]`.
#[derive(Clone, Debug, PartialEq)]
pub struct Field3D {
    /// The 3D grid.
    pub grid: Grid3D,
    /// Nodal values in `z`-major, `y`-middle, `x`-fastest order.
    pub values: Vec<f64>,
}

impl Field3D {
    /// Returns the value at grid node `(x, y, z)`.
    #[must_use]
    pub fn at(
        &self,
        x: usize,
        y: usize,
        z: usize,
    ) -> f64 {
        self.values[self.grid.index(x, y, z)]
    }
}

// ---------------------------------------------------------------------------
// Tridiagonal, Cyclic-Tridiagonal, and Ghost-Node Helpers
// ---------------------------------------------------------------------------

/// Thomas algorithm for a tridiagonal system with lower diagonal `a[1..n]`,
/// main diagonal `b[0..n]`, upper diagonal `c[0..n-1]`, and right-hand side `d[0..n]`.
fn solve_tridiagonal(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64],
) -> Result<Vec<f64>, FdmError> {
    let n = b.len();
    let mut c_star = vec![0.0; n];
    let mut d_star = vec![0.0; n];
    if b[0].abs() < 1e-300 || !b[0].is_finite() {
        return Err(FdmError::Singular);
    }
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];
    for i in 1..n {
        let denom = b[i] - a[i] * c_star[i - 1];
        if denom.abs() < 1e-300 || !denom.is_finite() {
            return Err(FdmError::Singular);
        }
        if i < n - 1 {
            c_star[i] = c[i] / denom;
        }
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / denom;
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d_star[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_star[i] - c_star[i] * x[i + 1];
    }
    if !x.iter().all(|v| v.is_finite()) {
        return Err(FdmError::Singular);
    }
    Ok(x)
}

/// Solves a periodic cyclic-tridiagonal system of size `n` with uniform diagonal
/// `diag` and off-diagonal `off` (`A_{i, i-1} = A_{i, i+1} = off`, `A_{0, n-1} = A_{n-1, 0} = off`)
/// via the Sherman–Morrison formula.
fn solve_cyclic_tridiagonal(
    diag: f64,
    off: f64,
    rhs: &[f64],
) -> Result<Vec<f64>, FdmError> {
    let n = rhs.len();
    let gamma = -diag;
    let mut b = vec![diag; n];
    b[0] -= gamma;
    b[n - 1] -= (off * off) / gamma;
    let a = vec![off; n];
    let mut c = vec![off; n];
    c[n - 1] = 0.0;

    let y = solve_tridiagonal(&a, &b, &c, rhs)?;
    let mut u = vec![0.0; n];
    u[0] = gamma;
    u[n - 1] = off;
    let q = solve_tridiagonal(&a, &b, &c, &u)?;

    let vx = y[0] + (off / gamma) * y[n - 1];
    let vq = q[0] + (off / gamma) * q[n - 1];
    let denom = 1.0 + vq;
    if denom.abs() < 1e-300 {
        return Err(FdmError::Singular);
    }
    let factor = vx / denom;
    Ok((0..n).map(|i| y[i] - factor * q[i]).collect())
}

/// Evaluates the 1D second-derivative stencil `(u_{i-1} - 2 u_i + u_{i+1}) / h^2`
/// along a line `u[0..n]` respecting `left` and `right` boundary conditions.
///
/// For Neumann `\partial u / \partial n = g` at `i = 0`, the outward normal points
/// left (`-x`), so `(u_{-1} - u_1) / (2h) = g => u_{-1} = u_1 + 2 h g`, giving
/// `u_{xx}(0) = 2 (u_1 - u_0 + h g) / h^2`.
/// For Robin `\alpha u + \beta (\partial u / \partial n) = \gamma` (`\beta \ne 0`),
/// `\partial u / \partial n = (\gamma - \alpha u_b) / \beta`.
fn line_second_diff(
    u: &[f64],
    h: f64,
    left: &BoundaryCondition,
    right: &BoundaryCondition,
    face_idx: usize,
) -> Vec<f64> {
    let n = u.len();
    let inv_h2 = 1.0 / (h * h);
    let mut out = vec![0.0; n];

    if matches!(left, BoundaryCondition::Periodic) {
        for i in 0..n {
            let im1 = if i == 0 { n - 1 } else { i - 1 };
            let ip1 = if i + 1 == n { 0 } else { i + 1 };
            out[i] = (u[im1] - 2.0 * u[i] + u[ip1]) * inv_h2;
        }
        return out;
    }

    for i in 1..n - 1 {
        out[i] = (u[i - 1] - 2.0 * u[i] + u[i + 1]) * inv_h2;
    }

    // Left endpoint i = 0
    out[0] = match left {
        | BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_) => 0.0,
        | BoundaryCondition::Neumann(g) => 2.0 * (u[1] - u[0] + h * g) * inv_h2,
        | BoundaryCondition::Robin { alpha, beta, gamma } => {
            if beta.abs() < 1e-15 {
                0.0
            } else {
                let dudn = (gamma - alpha * u[0]) / beta;
                2.0 * (u[1] - u[0] + h * dudn) * inv_h2
            }
        },
        | BoundaryCondition::Periodic => unreachable!(),
    };

    // Right endpoint i = n - 1
    out[n - 1] = match right {
        | BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_) => 0.0,
        | BoundaryCondition::Neumann(g) => 2.0 * (u[n - 2] - u[n - 1] + h * g) * inv_h2,
        | BoundaryCondition::Robin { alpha, beta, gamma } => {
            if beta.abs() < 1e-15 {
                0.0
            } else {
                let dudn = (gamma - alpha * u[n - 1]) / beta;
                2.0 * (u[n - 2] - u[n - 1] + h * dudn) * inv_h2
            }
        },
        | BoundaryCondition::Periodic => unreachable!(),
    };

    let _ = face_idx;
    out
}

/// Solves the 1D implicit diffusion operator `(I - r h^2 \partial_{xx}) u = rhs`
/// along a line of length `n` with `r = \theta \alpha \Delta t / h^2`, enforcing
/// `left` and `right` boundary conditions.
fn solve_line_implicit(
    rhs: &[f64],
    r: f64,
    h: f64,
    left: &BoundaryCondition,
    right: &BoundaryCondition,
    face_idx: usize,
) -> Result<Vec<f64>, FdmError> {
    let n = rhs.len();
    if matches!(left, BoundaryCondition::Periodic) {
        return solve_cyclic_tridiagonal(1.0 + 2.0 * r, -r, rhs);
    }

    let mut a = vec![-r; n];
    let mut b = vec![1.0 + 2.0 * r; n];
    let mut c = vec![-r; n];
    let mut d = rhs.to_vec();
    a[0] = 0.0;
    c[n - 1] = 0.0;

    // Left boundary row i = 0
    match left {
        | BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_) => {
            b[0] = 1.0;
            c[0] = 0.0;
            d[0] = left.dirichlet_at(face_idx);
        },
        | BoundaryCondition::Neumann(g) => {
            // u_0 - r * 2 * (u_1 - u_0 + h g) = rhs_0
            b[0] = 1.0 + 2.0 * r;
            c[0] = -2.0 * r;
            d[0] = rhs[0] + 2.0 * r * h * g;
        },
        | BoundaryCondition::Robin { alpha, beta, gamma } => {
            if beta.abs() < 1e-15 {
                b[0] = 1.0;
                c[0] = 0.0;
                d[0] = gamma / alpha;
            } else {
                // dudn = (gamma - alpha u_0) / beta
                b[0] = 1.0 + 2.0 * r * (1.0 + h * alpha / beta);
                c[0] = -2.0 * r;
                d[0] = rhs[0] + 2.0 * r * h * (gamma / beta);
            }
        },
        | BoundaryCondition::Periodic => unreachable!(),
    }

    // Right boundary row i = n - 1
    match right {
        | BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_) => {
            a[n - 1] = 0.0;
            b[n - 1] = 1.0;
            d[n - 1] = right.dirichlet_at(face_idx);
        },
        | BoundaryCondition::Neumann(g) => {
            a[n - 1] = -2.0 * r;
            b[n - 1] = 1.0 + 2.0 * r;
            d[n - 1] = rhs[n - 1] + 2.0 * r * h * g;
        },
        | BoundaryCondition::Robin { alpha, beta, gamma } => {
            if beta.abs() < 1e-15 {
                a[n - 1] = 0.0;
                b[n - 1] = 1.0;
                d[n - 1] = gamma / alpha;
            } else {
                a[n - 1] = -2.0 * r;
                b[n - 1] = 1.0 + 2.0 * r * (1.0 + h * alpha / beta);
                d[n - 1] = rhs[n - 1] + 2.0 * r * h * (gamma / beta);
            }
        },
        | BoundaryCondition::Periodic => unreachable!(),
    }

    solve_tridiagonal(&a, &b, &c, &d)
}

/// Applies Dirichlet boundary values directly onto a 1D vector.
fn enforce_dirichlet_1d(
    u: &mut [f64],
    bc: &Boundary1D,
) {
    let n = u.len();
    if matches!(
        bc.left,
        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
    ) {
        u[0] = bc.left.dirichlet_at(0);
    }
    if matches!(
        bc.right,
        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
    ) {
        u[n - 1] = bc.right.dirichlet_at(0);
    }
}

/// Applies Dirichlet boundary values directly onto a 2D field.
fn enforce_dirichlet_2d(
    u: &mut [f64],
    grid: Grid,
    bc: &Boundary2D,
) {
    let (nx, ny) = (grid.nx, grid.ny);
    if matches!(
        bc.left,
        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
    ) {
        for y in 0..ny {
            u[y * nx] = bc.left.dirichlet_at(y);
        }
    }
    if matches!(
        bc.right,
        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
    ) {
        for y in 0..ny {
            u[y * nx + (nx - 1)] = bc.right.dirichlet_at(y);
        }
    }
    if matches!(
        bc.bottom,
        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
    ) {
        for x in 0..nx {
            u[x] = bc.bottom.dirichlet_at(x);
        }
    }
    if matches!(
        bc.top,
        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
    ) {
        for x in 0..nx {
            u[(ny - 1) * nx + x] = bc.top.dirichlet_at(x);
        }
    }
}

// ---------------------------------------------------------------------------
// 1D, 2D (with BCs & Implicit/CN/ADI), and 3D Finite-Difference Solvers
// ---------------------------------------------------------------------------

/// Solves the 1D Poisson equation `u_{xx} = source` on `grid` with boundary
/// conditions `bc` using direct `O(N)` tridiagonal / cyclic-tridiagonal elimination.
///
/// Note: Pure Neumann or Periodic problems with zero reaction have a constant
/// nullspace; in that case the first node is pinned to zero mean.
///
/// # Errors
///
/// [`FdmError::Domain`] for invalid grid or boundary conditions,
/// [`FdmError::ShapeMismatch`] if `source.len() != grid.nx`.
pub fn solve_poisson_1d(
    grid: Grid1D,
    source: &[f64],
    bc: &Boundary1D,
) -> Result<Field1D, FdmError> {
    grid.check()?;
    bc.validate()?;
    let n = grid.nx;
    let h = grid.h;
    if source.len() != n {
        return Err(FdmError::ShapeMismatch(format!(
            "1D grid of {n} points given {} source values",
            source.len()
        )));
    }
    let h2 = h * h;

    if matches!(bc.left, BoundaryCondition::Periodic) {
        // Pin mode by small shift if singular, then subtract mean
        let mean_src: f64 = source.iter().sum::<f64>() / n as f64;
        let rhs: Vec<f64> = source.iter().map(|&f| -h2 * (f - mean_src)).collect();
        let mut u = solve_cyclic_tridiagonal(2.0 + 1e-12, -1.0, &rhs)?;
        let mean_u: f64 = u.iter().sum::<f64>() / n as f64;
        for v in &mut u {
            *v -= mean_u;
        }
        return Ok(Field1D { grid, values: u });
    }

    let mut a = vec![1.0; n];
    let mut b = vec![-2.0; n];
    let mut c = vec![1.0; n];
    let mut d: Vec<f64> = source.iter().map(|&f| h2 * f).collect();
    a[0] = 0.0;
    c[n - 1] = 0.0;

    match &bc.left {
        | BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_) => {
            b[0] = 1.0;
            c[0] = 0.0;
            d[0] = bc.left.dirichlet_at(0);
        },
        | BoundaryCondition::Neumann(g) => {
            b[0] = -2.0;
            c[0] = 2.0;
            d[0] = h2 * source[0] - 2.0 * h * g;
        },
        | BoundaryCondition::Robin { alpha, beta, gamma } => {
            if beta.abs() < 1e-15 {
                b[0] = 1.0;
                c[0] = 0.0;
                d[0] = gamma / alpha;
            } else {
                b[0] = -2.0 * (1.0 + h * alpha / beta);
                c[0] = 2.0;
                d[0] = h2 * source[0] - 2.0 * h * (gamma / beta);
            }
        },
        | BoundaryCondition::Periodic => unreachable!(),
    }

    match &bc.right {
        | BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_) => {
            // If pure Neumann on left and Dirichlet on right, system is nonsingular.
            a[n - 1] = 0.0;
            b[n - 1] = 1.0;
            d[n - 1] = bc.right.dirichlet_at(0);
        },
        | BoundaryCondition::Neumann(g) => {
            if matches!(bc.left, BoundaryCondition::Neumann(_)) {
                // Pure Neumann on both ends: pin u[0] = 0 to remove constant nullspace
                b[0] = 1.0;
                c[0] = 0.0;
                d[0] = 0.0;
            }
            a[n - 1] = 2.0;
            b[n - 1] = -2.0;
            d[n - 1] = h2 * source[n - 1] - 2.0 * h * g;
        },
        | BoundaryCondition::Robin { alpha, beta, gamma } => {
            if beta.abs() < 1e-15 {
                a[n - 1] = 0.0;
                b[n - 1] = 1.0;
                d[n - 1] = gamma / alpha;
            } else {
                a[n - 1] = 2.0;
                b[n - 1] = -2.0 * (1.0 + h * alpha / beta);
                d[n - 1] = h2 * source[n - 1] - 2.0 * h * (gamma / beta);
            }
        },
        | BoundaryCondition::Periodic => unreachable!(),
    }

    let values = solve_tridiagonal(&a, &b, &c, &d)?;
    Ok(Field1D { grid, values })
}

/// Solves the 2D Poisson equation `∇²u = spec.source` with general boundary
/// conditions `bc` (Dirichlet, non-zero Dirichlet, Neumann, Robin, or Periodic).
///
/// # Errors
///
/// [`FdmError::Domain`] for invalid grid or BCs,
/// [`FdmError::ShapeMismatch`] if `spec.source.len() != spec.grid.len()`, and
/// [`FdmError::NotConverged`] if SOR does not converge within `spec.iterations`.
pub fn solve_poisson_with_bc(
    spec: &PoissonSpec,
    bc: &Boundary2D,
) -> Result<Field, FdmError> {
    spec.grid.check()?;
    bc.validate()?;
    if spec.source.len() != spec.grid.len() {
        return Err(FdmError::ShapeMismatch(format!(
            "a {}x{} grid needs {} source values, found {}",
            spec.grid.ny,
            spec.grid.nx,
            spec.grid.len(),
            spec.source.len()
        )));
    }
    let (nx, ny, h) = (spec.grid.nx, spec.grid.ny, spec.grid.h);
    let h2 = h * h;
    let mut u = vec![0.0; spec.grid.len()];
    enforce_dirichlet_2d(&mut u, spec.grid, bc);

    let omega = 1.5;
    let tolerance = if spec.tolerance > 0.0 {
        spec.tolerance
    } else {
        1e-10
    };
    let xp = matches!(bc.left, BoundaryCondition::Periodic);
    let yp = matches!(bc.bottom, BoundaryCondition::Periodic);
    let pure_singular = (xp
        || (matches!(bc.left, BoundaryCondition::Neumann(_))
            && matches!(bc.right, BoundaryCondition::Neumann(_))))
        && (yp
            || (matches!(bc.bottom, BoundaryCondition::Neumann(_))
                && matches!(bc.top, BoundaryCondition::Neumann(_))));

    for _ in 0..spec.iterations.max(1) {
        let mut change = 0.0f64;
        for y in 0..ny {
            for x in 0..nx {
                let idx = y * nx + x;
                // Check if this boundary node is fixed by Dirichlet
                if !xp
                    && x == 0
                    && matches!(
                        bc.left,
                        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
                    )
                {
                    u[idx] = bc.left.dirichlet_at(y);
                    continue;
                }
                if !xp
                    && x == nx - 1
                    && matches!(
                        bc.right,
                        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
                    )
                {
                    u[idx] = bc.right.dirichlet_at(y);
                    continue;
                }
                if !yp
                    && y == 0
                    && matches!(
                        bc.bottom,
                        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
                    )
                {
                    u[idx] = bc.bottom.dirichlet_at(x);
                    continue;
                }
                if !yp
                    && y == ny - 1
                    && matches!(
                        bc.top,
                        BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_)
                    )
                {
                    u[idx] = bc.top.dirichlet_at(x);
                    continue;
                }

                // Neighbor contributions along x
                let (sum_x, diag_x) = if xp {
                    let xm = if x == 0 { nx - 1 } else { x - 1 };
                    let xp_i = if x + 1 == nx { 0 } else { x + 1 };
                    (u[y * nx + xm] + u[y * nx + xp_i], 2.0)
                } else if x == 0 {
                    match &bc.left {
                        | BoundaryCondition::Neumann(g) => (2.0 * u[y * nx + 1] + 2.0 * h * g, 2.0),
                        | BoundaryCondition::Robin { alpha, beta, gamma } => {
                            (
                                2.0 * u[y * nx + 1] + 2.0 * h * gamma / beta,
                                2.0 * (1.0 + h * alpha / beta),
                            )
                        },
                        | _ => (0.0, 2.0),
                    }
                } else if x == nx - 1 {
                    match &bc.right {
                        | BoundaryCondition::Neumann(g) => {
                            (2.0 * u[y * nx + nx - 2] + 2.0 * h * g, 2.0)
                        },
                        | BoundaryCondition::Robin { alpha, beta, gamma } => {
                            (
                                2.0 * u[y * nx + nx - 2] + 2.0 * h * gamma / beta,
                                2.0 * (1.0 + h * alpha / beta),
                            )
                        },
                        | _ => (0.0, 2.0),
                    }
                } else {
                    (u[idx - 1] + u[idx + 1], 2.0)
                };

                // Neighbor contributions along y
                let (sum_y, diag_y) = if yp {
                    let ym = if y == 0 { ny - 1 } else { y - 1 };
                    let yp_i = if y + 1 == ny { 0 } else { y + 1 };
                    (u[ym * nx + x] + u[yp_i * nx + x], 2.0)
                } else if y == 0 {
                    match &bc.bottom {
                        | BoundaryCondition::Neumann(g) => (2.0 * u[nx + x] + 2.0 * h * g, 2.0),
                        | BoundaryCondition::Robin { alpha, beta, gamma } => {
                            (
                                2.0 * u[nx + x] + 2.0 * h * gamma / beta,
                                2.0 * (1.0 + h * alpha / beta),
                            )
                        },
                        | _ => (0.0, 2.0),
                    }
                } else if y == ny - 1 {
                    match &bc.top {
                        | BoundaryCondition::Neumann(g) => {
                            (2.0 * u[(ny - 2) * nx + x] + 2.0 * h * g, 2.0)
                        },
                        | BoundaryCondition::Robin { alpha, beta, gamma } => {
                            (
                                2.0 * u[(ny - 2) * nx + x] + 2.0 * h * gamma / beta,
                                2.0 * (1.0 + h * alpha / beta),
                            )
                        },
                        | _ => (0.0, 2.0),
                    }
                } else {
                    (u[idx - nx] + u[idx + nx], 2.0)
                };

                let target = (sum_x + sum_y - h2 * spec.source[idx]) / (diag_x + diag_y);
                let updated = u[idx] + omega * (target - u[idx]);
                change = change.max((updated - u[idx]).abs());
                u[idx] = updated;
            }
        }
        if pure_singular {
            let mean: f64 = u.iter().sum::<f64>() / u.len() as f64;
            for v in &mut u {
                *v -= mean;
            }
        }
        if change < tolerance {
            return Ok(Field {
                grid: spec.grid,
                values: u,
            });
        }
    }
    Err(FdmError::NotConverged)
}

fn fixed_value_3d(
    grid: Grid3D,
    bc: &Boundary3D,
    x: usize,
    y: usize,
    z: usize,
) -> Option<f64> {
    let faces = [
        (x == 0, &bc.x_min, z * grid.ny + y),
        (x + 1 == grid.nx, &bc.x_max, z * grid.ny + y),
        (y == 0, &bc.y_min, z * grid.nx + x),
        (y + 1 == grid.ny, &bc.y_max, z * grid.nx + x),
        (z == 0, &bc.z_min, y * grid.nx + x),
        (z + 1 == grid.nz, &bc.z_max, y * grid.nx + x),
    ];
    for (active, face, index) in faces {
        if !active {
            continue;
        }
        match face {
            | BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_) => {
                return Some(face.dirichlet_at(index));
            },
            | BoundaryCondition::Robin { alpha, beta, gamma } if beta.abs() < 1e-15 => {
                return Some(gamma / alpha);
            },
            | _ => {},
        }
    }
    None
}

fn enforce_dirichlet_3d(
    u: &mut [f64],
    grid: Grid3D,
    bc: &Boundary3D,
) {
    for z in 0..grid.nz {
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                if let Some(value) = fixed_value_3d(grid, bc, x, y, z) {
                    u[grid.index(x, y, z)] = value;
                }
            }
        }
    }
}

fn axis_stencil_3d(
    u: &[f64],
    grid: Grid3D,
    bc: &Boundary3D,
    xyz: [usize; 3],
    axis: usize,
) -> (f64, f64) {
    let dims = [grid.nx, grid.ny, grid.nz];
    let position = xyz[axis];
    let n = dims[axis];
    let (left, right) = match axis {
        | 0 => (&bc.x_min, &bc.x_max),
        | 1 => (&bc.y_min, &bc.y_max),
        | _ => (&bc.z_min, &bc.z_max),
    };
    let at = |p| {
        let mut point = xyz;
        point[axis] = p;
        u[grid.index(point[0], point[1], point[2])]
    };
    if matches!(left, BoundaryCondition::Periodic) {
        return (at((position + n - 1) % n) + at((position + 1) % n), 2.0);
    }
    if position > 0 && position + 1 < n {
        return (at(position - 1) + at(position + 1), 2.0);
    }
    let face = if position == 0 { left } else { right };
    let neighbor = at(if position == 0 { 1 } else { n - 2 });
    match face {
        | BoundaryCondition::Neumann(g) => (2.0 * neighbor + 2.0 * grid.h * g, 2.0),
        | BoundaryCondition::Robin { alpha, beta, gamma } => {
            (
                2.0 * neighbor + 2.0 * grid.h * gamma / beta,
                2.0 * (1.0 + grid.h * alpha / beta),
            )
        },
        | _ => (0.0, 0.0),
    }
}

fn laplacian_3d(
    u: &[f64],
    grid: Grid3D,
    bc: &Boundary3D,
    axis: Option<usize>,
) -> Vec<f64> {
    let mut output = vec![0.0; u.len()];
    for z in 0..grid.nz {
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                if fixed_value_3d(grid, bc, x, y, z).is_some() {
                    continue;
                }
                let index = grid.index(x, y, z);
                let mut sum = 0.0;
                for a in 0..3 {
                    if axis.is_none() || axis == Some(a) {
                        let (neighbors, diagonal) = axis_stencil_3d(u, grid, bc, [x, y, z], a);
                        sum += neighbors - diagonal * u[index];
                    }
                }
                output[index] = sum / (grid.h * grid.h);
            }
        }
    }
    output
}

fn homogeneous_3d(bc: &Boundary3D) -> Boundary3D {
    let face = |condition: &BoundaryCondition| {
        match condition {
            | BoundaryCondition::Dirichlet(_) | BoundaryCondition::DirichletProfile(_) => {
                BoundaryCondition::Dirichlet(0.0)
            },
            | BoundaryCondition::Neumann(_) => BoundaryCondition::Neumann(0.0),
            | BoundaryCondition::Robin { alpha, beta, .. } => {
                BoundaryCondition::Robin {
                    alpha: *alpha,
                    beta: *beta,
                    gamma: 0.0,
                }
            },
            | BoundaryCondition::Periodic => BoundaryCondition::Periodic,
        }
    };
    Boundary3D {
        x_min: face(&bc.x_min),
        x_max: face(&bc.x_max),
        y_min: face(&bc.y_min),
        y_max: face(&bc.y_max),
        z_min: face(&bc.z_min),
        z_max: face(&bc.z_max),
    }
}

/// Solves the 3D Poisson equation `∇²u = source` on `grid` with boundary conditions
/// `bc` using Successive Over-Relaxation (SOR).
///
/// # Errors
///
/// [`FdmError::Domain`] for invalid grid or BCs,
/// [`FdmError::ShapeMismatch`] if `source.len() != grid.len()`, and
/// [`FdmError::NotConverged`] if `iterations` is exceeded.
pub fn solve_poisson_3d(
    grid: Grid3D,
    source: &[f64],
    bc: &Boundary3D,
    iterations: usize,
    tolerance: f64,
) -> Result<Field3D, FdmError> {
    grid.check()?;
    bc.validate()?;
    if source.len() != grid.len() {
        return Err(FdmError::ShapeMismatch(format!(
            "3D grid of {} points given {} source values",
            grid.len(),
            source.len()
        )));
    }
    if source.iter().any(|v| !v.is_finite())
        || !tolerance.is_finite()
        || tolerance <= 0.0
        || iterations == 0
    {
        return Err(FdmError::Domain(
            "invalid Poisson source or iteration controls".into(),
        ));
    }
    let mut u = vec![0.0; grid.len()];
    enforce_dirichlet_3d(&mut u, grid, bc);
    let pure_singular = [
        &bc.x_min, &bc.x_max, &bc.y_min, &bc.y_max, &bc.z_min, &bc.z_max,
    ]
    .iter()
    .all(|v| {
        matches!(
            v,
            BoundaryCondition::Periodic | BoundaryCondition::Neumann(_)
        )
    });
    for _ in 0..iterations {
        let mut change = 0.0_f64;
        for z in 0..grid.nz {
            for y in 0..grid.ny {
                for x in 0..grid.nx {
                    if fixed_value_3d(grid, bc, x, y, z).is_some() {
                        continue;
                    }
                    let index = grid.index(x, y, z);
                    let mut sum = 0.0;
                    let mut diagonal = 0.0;
                    for axis in 0..3 {
                        let (v, d) = axis_stencil_3d(&u, grid, bc, [x, y, z], axis);
                        sum += v;
                        diagonal += d;
                    }
                    if diagonal.abs() < 1e-15 {
                        return Err(FdmError::Singular);
                    }
                    let next = u[index]
                        + 1.3 * ((sum - grid.h * grid.h * source[index]) / diagonal - u[index]);
                    change = change.max((next - u[index]).abs());
                    u[index] = next;
                }
            }
        }
        if pure_singular {
            let mean = u.iter().sum::<f64>() / u.len() as f64;
            for v in &mut u {
                *v -= mean;
            }
        }
        if u.iter().any(|v| !v.is_finite()) {
            return Err(FdmError::NotConverged);
        }
        if change < tolerance {
            return Ok(Field3D { grid, values: u });
        }
    }
    Err(FdmError::NotConverged)
}

/// Solves the 1D heat equation `u_t = \alpha u_{xx}` on `grid` with boundary
/// conditions `bc` and time-stepping `stepping` (`Explicit`, `BackwardEuler`,
/// `CrankNicolson`, or `Adi`).
///
/// Implicit schemes (`BackwardEuler`, `CrankNicolson`, `Adi`) are unconditionally
/// stable for any `\Delta t > 0`.
///
/// # Errors
///
/// [`FdmError::Domain`] for invalid parameters or unstable explicit step,
/// [`FdmError::ShapeMismatch`] if `initial.len() != grid.nx`.
pub fn solve_heat_1d(
    grid: Grid1D,
    alpha: f64,
    dt: f64,
    steps: usize,
    initial: &[f64],
    bc: &Boundary1D,
    stepping: TimeStepping,
) -> Result<Vec<Field1D>, FdmError> {
    grid.check()?;
    bc.validate()?;
    if initial.len() != grid.nx {
        return Err(FdmError::ShapeMismatch(
            "1D initial field does not fit grid".into(),
        ));
    }
    if !alpha.is_finite() || alpha < 0.0 || !dt.is_finite() || dt <= 0.0 {
        return Err(FdmError::Domain(
            "alpha >= 0 and dt > 0 must be finite".into(),
        ));
    }
    let h = grid.h;
    let spacing_squared = h * h;
    let r = alpha * dt / spacing_squared;
    if stepping == TimeStepping::Explicit && r > 0.5 {
        return Err(FdmError::Domain(format!(
            "1D explicit heat scheme requires alpha*dt/h^2 <= 0.5, found {r}"
        )));
    }

    let mut u = initial.to_vec();
    enforce_dirichlet_1d(&mut u, bc);
    let mut out = Vec::with_capacity(steps + 1);
    out.push(Field1D {
        grid,
        values: u.clone(),
    });

    for _ in 0..steps {
        u = match stepping {
            | TimeStepping::Explicit => {
                let d2 = line_second_diff(&u, h, &bc.left, &bc.right, 0);
                let mut next: Vec<f64> = u
                    .iter()
                    .zip(&d2)
                    .map(|(&ui, &l)| ui + alpha * dt * l)
                    .collect();
                enforce_dirichlet_1d(&mut next, bc);
                next
            },
            | TimeStepping::BackwardEuler => solve_line_implicit(&u, r, h, &bc.left, &bc.right, 0)?,
            | TimeStepping::CrankNicolson | TimeStepping::Adi => {
                let d2 = line_second_diff(&u, h, &bc.left, &bc.right, 0);
                let rhs: Vec<f64> = u
                    .iter()
                    .zip(&d2)
                    .map(|(&ui, &l)| ui + 0.5 * alpha * dt * l)
                    .collect();
                solve_line_implicit(&rhs, 0.5 * r, h, &bc.left, &bc.right, 0)?
            },
        };
        out.push(Field1D {
            grid,
            values: u.clone(),
        });
    }
    Ok(out)
}

/// Solve the two-dimensional heat equation with prescribed face data.
///
/// Solves the 2D heat equation `u_t = \alpha ∇²u` with general boundary conditions
/// `bc` and time-stepping `stepping` (`Explicit`, `BackwardEuler`, `CrankNicolson`,
/// or `Adi` Peaceman–Rachford Alternating Direction Implicit).
///
/// `BackwardEuler`, `CrankNicolson`, and `Adi` are unconditionally stable for
/// arbitrarily large diffusion Courant numbers `\alpha \Delta t / h^2`.
///
/// # Errors
///
/// [`FdmError::Domain`] for invalid grid, BCs, or unstable explicit CFL,
/// [`FdmError::ShapeMismatch`] if `spec.initial.len() != spec.grid.len()`.
pub fn solve_heat_with_bc(
    spec: &HeatSpec,
    bc: &Boundary2D,
    stepping: TimeStepping,
) -> Result<Vec<Field>, FdmError> {
    spec.grid.check()?;
    bc.validate()?;
    if spec.initial.len() != spec.grid.len() {
        return Err(FdmError::ShapeMismatch(
            "the initial field does not fit the grid".into(),
        ));
    }
    let (nx, ny, h) = (spec.grid.nx, spec.grid.ny, spec.grid.h);
    let r = spec.alpha * spec.dt / (h * h);
    if !r.is_finite() || spec.alpha < 0.0 || spec.dt <= 0.0 {
        return Err(FdmError::Domain("invalid heat parameters".into()));
    }
    if stepping == TimeStepping::Explicit && r > 0.25 {
        return Err(FdmError::Domain(format!(
            "the explicit heat scheme needs `alpha*dt/h^2 <= 0.25`, found {r}"
        )));
    }

    let mut u = spec.initial.clone();
    enforce_dirichlet_2d(&mut u, spec.grid, bc);
    let mut out = Vec::with_capacity(spec.steps + 1);
    out.push(Field {
        grid: spec.grid,
        values: u.clone(),
    });

    // Helper: compute \partial_{xx} u for every row y
    let diff_xx = |field: &[f64]| -> Vec<f64> {
        let mut dxx = vec![0.0; nx * ny];
        for y in 0..ny {
            let row = &field[y * nx..(y + 1) * nx];
            let d2 = line_second_diff(row, h, &bc.left, &bc.right, y);
            dxx[y * nx..(y + 1) * nx].copy_from_slice(&d2);
        }
        dxx
    };

    // Helper: compute \partial_{yy} u for every column x
    let diff_yy = |field: &[f64]| -> Vec<f64> {
        let mut dyy = vec![0.0; nx * ny];
        let mut col = vec![0.0; ny];
        for x in 0..nx {
            for y in 0..ny {
                col[y] = field[y * nx + x];
            }
            let d2 = line_second_diff(&col, h, &bc.bottom, &bc.top, x);
            for y in 0..ny {
                dyy[y * nx + x] = d2[y];
            }
        }
        dyy
    };

    // Helper: implicit line solves along x: (I - rx h^2 \partial_{xx}) u = rhs
    let implicit_x = |rhs_field: &[f64], rx: f64| -> Result<Vec<f64>, FdmError> {
        let mut res = vec![0.0; nx * ny];
        for y in 0..ny {
            let row = &rhs_field[y * nx..(y + 1) * nx];
            let sol = solve_line_implicit(row, rx, h, &bc.left, &bc.right, y)?;
            res[y * nx..(y + 1) * nx].copy_from_slice(&sol);
        }
        Ok(res)
    };

    // Helper: implicit line solves along y: (I - ry h^2 \partial_{yy}) u = rhs
    let implicit_y = |rhs_field: &[f64], ry: f64| -> Result<Vec<f64>, FdmError> {
        let mut res = vec![0.0; nx * ny];
        let mut col = vec![0.0; ny];
        for x in 0..nx {
            for y in 0..ny {
                col[y] = rhs_field[y * nx + x];
            }
            let sol = solve_line_implicit(&col, ry, h, &bc.bottom, &bc.top, x)?;
            for y in 0..ny {
                res[y * nx + x] = sol[y];
            }
        }
        Ok(res)
    };

    for _ in 0..spec.steps {
        u = match stepping {
            | TimeStepping::Explicit => {
                let dxx = diff_xx(&u);
                let dyy = diff_yy(&u);
                let mut next: Vec<f64> = (0..u.len())
                    .map(|i| u[i] + spec.alpha * spec.dt * (dxx[i] + dyy[i]))
                    .collect();
                enforce_dirichlet_2d(&mut next, spec.grid, bc);
                next
            },
            | TimeStepping::BackwardEuler => {
                // Fractional-step / LOD Backward Euler: (I - dt alpha d_xx)(I - dt alpha d_yy) u^{n+1} = u^n
                let u_half = implicit_x(&u, r)?;
                let mut u_next = implicit_y(&u_half, r)?;
                enforce_dirichlet_2d(&mut u_next, spec.grid, bc);
                u_next
            },
            | TimeStepping::CrankNicolson | TimeStepping::Adi => {
                // Peaceman-Rachford ADI (2nd-order in time and space, unconditionally stable):
                // Step 1 (half step): (I - r/2 h^2 d_xx) u* = (I + r/2 h^2 d_yy) u^n
                let dyy = diff_yy(&u);
                let rhs1: Vec<f64> = (0..u.len())
                    .map(|i| u[i] + 0.5 * spec.alpha * spec.dt * dyy[i])
                    .collect();
                let u_star = implicit_x(&rhs1, 0.5 * r)?;

                // Step 2 (half step): (I - r/2 h^2 d_yy) u^{n+1} = (I + r/2 h^2 d_xx) u*
                let dxx_star = diff_xx(&u_star);
                let rhs2: Vec<f64> = (0..u.len())
                    .map(|i| u_star[i] + 0.5 * spec.alpha * spec.dt * dxx_star[i])
                    .collect();
                let mut u_next = implicit_y(&rhs2, 0.5 * r)?;
                enforce_dirichlet_2d(&mut u_next, spec.grid, bc);
                u_next
            },
        };
        out.push(Field {
            grid: spec.grid,
            values: u.clone(),
        });
    }

    Ok(out)
}

/// Solve the three-dimensional heat equation with prescribed face data.
///
/// Solves the 3D heat equation `u_t = \alpha ∇²u` on `grid` with boundary conditions
/// `bc` and time-stepping `stepping` (`Explicit`, `BackwardEuler`, `CrankNicolson`,
/// or `Adi` Douglas–Gunn 3D Alternating Direction Implicit).
///
/// # Errors
///
/// [`FdmError::Domain`] for invalid parameters or unstable explicit step (`r > 1/6`),
/// [`FdmError::ShapeMismatch`] if `initial.len() != grid.len()`.
pub fn solve_heat_3d(
    grid: Grid3D,
    alpha: f64,
    dt: f64,
    steps: usize,
    initial: &[f64],
    bc: &Boundary3D,
    stepping: TimeStepping,
) -> Result<Vec<Field3D>, FdmError> {
    grid.check()?;
    bc.validate()?;
    if initial.len() != grid.len() {
        return Err(FdmError::ShapeMismatch(
            "3D initial field does not fit grid".into(),
        ));
    }
    if initial.iter().any(|v| !v.is_finite())
        || !alpha.is_finite()
        || alpha < 0.0
        || !dt.is_finite()
        || dt <= 0.0
    {
        return Err(FdmError::Domain("invalid 3D heat inputs".into()));
    }
    let spacing_squared = grid.h * grid.h;
    let r = alpha * dt / spacing_squared;
    if !r.is_finite() || (stepping == TimeStepping::Explicit && r > 1.0 / 6.0) {
        return Err(FdmError::Domain(
            "invalid or unstable 3D diffusion Courant number".into(),
        ));
    }
    let mut u = initial.to_vec();
    enforce_dirichlet_3d(&mut u, grid, bc);
    let mut output = vec![Field3D {
        grid,
        values: u.clone(),
    }];
    let homogeneous = homogeneous_3d(bc);
    let boundary_source = laplacian_3d(&vec![0.0; u.len()], grid, bc, None);
    for _ in 0..steps {
        let lap = laplacian_3d(&u, grid, bc, None);
        match stepping {
            | TimeStepping::Explicit => {
                for i in 0..u.len() {
                    u[i] += alpha * dt * lap[i];
                }
            },
            | TimeStepping::BackwardEuler | TimeStepping::CrankNicolson => {
                let theta = if stepping == TimeStepping::BackwardEuler {
                    1.0
                } else {
                    0.5
                };
                let rhs: Vec<_> = (0..u.len())
                    .map(|i| {
                        u[i] + alpha * dt * ((1.0 - theta) * lap[i] + theta * boundary_source[i])
                    })
                    .collect();
                let apply = |state: &[f64], out: &mut [f64]| {
                    let linear = laplacian_3d(state, grid, &homogeneous, None);
                    for i in 0..state.len() {
                        out[i] = state[i] - theta * alpha * dt * linear[i];
                    }
                };
                let result = crate::kernels::krylov::gmres(
                    apply,
                    |x, y| y.copy_from_slice(x),
                    &rhs,
                    Some(&u),
                    40.min(u.len()),
                    1e-11,
                    20 * u.len(),
                );
                if !result.converged {
                    return Err(FdmError::NotConverged);
                }
                u = result.x;
            },
            | TimeStepping::Adi => {
                // Douglas theta=1/2 splitting: each sweep subtracts its old
                // directional derivative before solving the implicit factor.
                let mut stage: Vec<_> = (0..u.len()).map(|i| u[i] + alpha * dt * lap[i]).collect();
                for axis in 0..3 {
                    let old_axis = laplacian_3d(&u, grid, bc, Some(axis));
                    for i in 0..stage.len() {
                        stage[i] -= 0.5 * alpha * dt * old_axis[i];
                    }
                    let dims = [grid.nx, grid.ny, grid.nz];
                    let n = dims[axis];
                    let (left, right) = match axis {
                        | 0 => (&bc.x_min, &bc.x_max),
                        | 1 => (&bc.y_min, &bc.y_max),
                        | _ => (&bc.z_min, &bc.z_max),
                    };
                    for a in 0..dims[(axis + 1) % 3] {
                        for b in 0..dims[(axis + 2) % 3] {
                            let mut point = [0; 3];
                            point[(axis + 1) % 3] = a;
                            point[(axis + 2) % 3] = b;
                            let line: Vec<_> = (0..n)
                                .map(|j| {
                                    point[axis] = j;
                                    stage[grid.index(point[0], point[1], point[2])]
                                })
                                .collect();
                            let face_index = match axis {
                                | 0 => point[2] * grid.ny + point[1],
                                | 1 => point[2] * grid.nx + point[0],
                                | _ => point[1] * grid.nx + point[0],
                            };
                            let solved = solve_line_implicit(
                                &line,
                                0.5 * r,
                                grid.h,
                                left,
                                right,
                                face_index,
                            )?;
                            for j in 0..n {
                                point[axis] = j;
                                stage[grid.index(point[0], point[1], point[2])] = solved[j];
                            }
                        }
                    }
                    enforce_dirichlet_3d(&mut stage, grid, bc);
                }
                u = stage;
            },
        }
        enforce_dirichlet_3d(&mut u, grid, bc);
        if u.iter().any(|v| !v.is_finite()) {
            return Err(FdmError::NotConverged);
        }
        output.push(Field3D {
            grid,
            values: u.clone(),
        });
    }
    Ok(output)
}

/// Solves the 1D wave equation `u_{tt} = c^2 u_{xx}` using the leapfrog scheme
/// with boundary conditions `bc`.
///
/// # Errors
///
/// [`FdmError::Domain`] if CFL `|c| dt / h > 1`, [`FdmError::ShapeMismatch`]
/// if `initial` or `velocity` does not match `grid.nx`.
pub fn solve_wave_1d(
    grid: Grid1D,
    c: f64,
    dt: f64,
    steps: usize,
    initial: &[f64],
    velocity: &[f64],
    bc: &Boundary1D,
) -> Result<Vec<Field1D>, FdmError> {
    grid.check()?;
    bc.validate()?;
    let n = grid.nx;
    if initial.len() != n || velocity.len() != n {
        return Err(FdmError::ShapeMismatch(
            "1D wave initial or velocity field does not fit grid".into(),
        ));
    }
    let courant = c * dt / grid.h;
    if !courant.is_finite() || courant.abs() > 1.0 + 1e-12 || dt <= 0.0 {
        return Err(FdmError::Domain(format!(
            "1D wave leapfrog requires |c|*dt/h <= 1, found {courant}"
        )));
    }
    let c2_dt2 = (c * dt) * (c * dt);
    let mut prev = initial.to_vec();
    enforce_dirichlet_1d(&mut prev, bc);
    let mut out = Vec::with_capacity(steps + 1);
    out.push(Field1D {
        grid,
        values: prev.clone(),
    });
    if steps == 0 {
        return Ok(out);
    }

    let d2_0 = line_second_diff(&prev, grid.h, &bc.left, &bc.right, 0);
    let mut curr: Vec<f64> = (0..n)
        .map(|i| prev[i] + dt * velocity[i] + 0.5 * c2_dt2 * d2_0[i])
        .collect();
    enforce_dirichlet_1d(&mut curr, bc);
    out.push(Field1D {
        grid,
        values: curr.clone(),
    });

    for _ in 1..steps {
        let d2 = line_second_diff(&curr, grid.h, &bc.left, &bc.right, 0);
        let mut next: Vec<f64> = (0..n)
            .map(|i| 2.0 * curr[i] - prev[i] + c2_dt2 * d2[i])
            .collect();
        enforce_dirichlet_1d(&mut next, bc);
        prev = curr;
        curr = next;
        out.push(Field1D {
            grid,
            values: curr.clone(),
        });
    }
    Ok(out)
}

/// Solves the 2D wave equation `u_{tt} = c^2 ∇²u` with general boundary conditions `bc`.
///
/// # Errors
///
/// [`FdmError::Domain`] if CFL `|c| dt / h > 1/\sqrt{2}`, [`FdmError::ShapeMismatch`]
/// if `initial` or `velocity` does not match `spec.grid.len()`.
pub fn solve_wave_with_bc(
    spec: &WaveSpec,
    bc: &Boundary2D,
) -> Result<Vec<Field>, FdmError> {
    spec.grid.check()?;
    bc.validate()?;
    let (nx, ny, h) = (spec.grid.nx, spec.grid.ny, spec.grid.h);
    if spec.initial.len() != nx * ny || spec.velocity.len() != nx * ny {
        return Err(FdmError::ShapeMismatch(
            "the initial fields do not fit the grid".into(),
        ));
    }
    let courant = spec.c * spec.dt / h;
    if !(courant.is_finite() && courant.abs() <= std::f64::consts::FRAC_1_SQRT_2 + 1e-12) {
        return Err(FdmError::Domain(format!(
            "the leapfrog scheme needs `c*dt/h <= 1/sqrt(2)`, found {courant}"
        )));
    }
    let c2_dt2 = (spec.c * spec.dt) * (spec.c * spec.dt);

    let lap_bc = |field: &[f64]| -> Vec<f64> {
        let mut lap = vec![0.0; nx * ny];
        for y in 0..ny {
            let row = &field[y * nx..(y + 1) * nx];
            let d2x = line_second_diff(row, h, &bc.left, &bc.right, y);
            for x in 0..nx {
                lap[y * nx + x] += d2x[x];
            }
        }
        let mut col = vec![0.0; ny];
        for x in 0..nx {
            for y in 0..ny {
                col[y] = field[y * nx + x];
            }
            let d2y = line_second_diff(&col, h, &bc.bottom, &bc.top, x);
            for y in 0..ny {
                lap[y * nx + x] += d2y[y];
            }
        }
        lap
    };

    let mut prev = spec.initial.clone();
    enforce_dirichlet_2d(&mut prev, spec.grid, bc);
    let lap0 = lap_bc(&prev);
    let mut curr: Vec<f64> = (0..prev.len())
        .map(|i| prev[i] + spec.dt * spec.velocity[i] + 0.5 * c2_dt2 * lap0[i])
        .collect();
    enforce_dirichlet_2d(&mut curr, spec.grid, bc);

    let mut out = vec![
        Field {
            grid: spec.grid,
            values: prev.clone(),
        },
        Field {
            grid: spec.grid,
            values: curr.clone(),
        },
    ];
    for _ in 1..spec.steps {
        let lap = lap_bc(&curr);
        let mut next: Vec<f64> = (0..curr.len())
            .map(|i| 2.0 * curr[i] - prev[i] + c2_dt2 * lap[i])
            .collect();
        enforce_dirichlet_2d(&mut next, spec.grid, bc);
        prev = curr;
        curr = next;
        out.push(Field {
            grid: spec.grid,
            values: curr.clone(),
        });
    }
    Ok(out)
}

/// Solves the 3D wave equation `u_{tt} = c^2 ∇²u` on `grid` using the leapfrog scheme.
///
/// # Errors
///
/// [`FdmError::Domain`] if CFL `|c| dt / h > 1/\sqrt{3}`, [`FdmError::ShapeMismatch`]
/// if `initial` or `velocity` does not match `grid.len()`.
pub fn solve_wave_3d(
    grid: Grid3D,
    c: f64,
    dt: f64,
    steps: usize,
    initial: &[f64],
    velocity: &[f64],
) -> Result<Vec<Field3D>, FdmError> {
    solve_wave_3d_with_bc(
        grid,
        c,
        dt,
        steps,
        initial,
        velocity,
        &Boundary3D::default(),
    )
}

/// Three-dimensional leapfrog wave solver with Dirichlet, Neumann, Robin, or periodic faces.
///
/// # Errors
/// Invalid finite inputs, incompatible field dimensions, or unstable Courant number.
pub fn solve_wave_3d_with_bc(
    grid: Grid3D,
    c: f64,
    dt: f64,
    steps: usize,
    initial: &[f64],
    velocity: &[f64],
    bc: &Boundary3D,
) -> Result<Vec<Field3D>, FdmError> {
    grid.check()?;
    bc.validate()?;
    if initial.len() != grid.len() || velocity.len() != grid.len() {
        return Err(FdmError::ShapeMismatch(
            "3D wave fields must match the grid".into(),
        ));
    }
    if !c.is_finite()
        || !dt.is_finite()
        || dt <= 0.0
        || (c * dt / grid.h).abs() > 1.0 / 3.0_f64.sqrt()
        || initial.iter().chain(velocity).any(|v| !v.is_finite())
    {
        return Err(FdmError::Domain(
            "invalid or unstable 3D wave inputs".into(),
        ));
    }
    let mut prev = initial.to_vec();
    enforce_dirichlet_3d(&mut prev, grid, bc);
    let mut output = vec![Field3D {
        grid,
        values: prev.clone(),
    }];
    if steps == 0 {
        return Ok(output);
    }
    let lap0 = laplacian_3d(&prev, grid, bc, None);
    let coefficient = c * c * dt * dt;
    let mut current: Vec<_> = (0..grid.len())
        .map(|i| prev[i] + dt * velocity[i] + 0.5 * coefficient * lap0[i])
        .collect();
    enforce_dirichlet_3d(&mut current, grid, bc);
    output.push(Field3D {
        grid,
        values: current.clone(),
    });
    for _ in 1..steps {
        let lap = laplacian_3d(&current, grid, bc, None);
        let mut next: Vec<_> = (0..grid.len())
            .map(|i| 2.0 * current[i] - prev[i] + coefficient * lap[i])
            .collect();
        enforce_dirichlet_3d(&mut next, grid, bc);
        prev = current;
        current = next;
        if current.iter().any(|v| !v.is_finite()) {
            return Err(FdmError::NotConverged);
        }
        output.push(Field3D {
            grid,
            values: current.clone(),
        });
    }
    Ok(output)
}

fn upwind_delta_3d(
    u: &[f64],
    grid: Grid3D,
    bc: &Boundary3D,
    point: [usize; 3],
    axis: usize,
    courant: f64,
) -> f64 {
    let n = [grid.nx, grid.ny, grid.nz][axis];
    let position = point[axis];
    let (left, right) = match axis {
        | 0 => (&bc.x_min, &bc.x_max),
        | 1 => (&bc.y_min, &bc.y_max),
        | _ => (&bc.z_min, &bc.z_max),
    };
    let at = |p| {
        let mut q = point;
        q[axis] = p;
        u[grid.index(q[0], q[1], q[2])]
    };
    if courant >= 0.0 {
        if position > 0 {
            return at(position) - at(position - 1);
        }
        if matches!(left, BoundaryCondition::Periodic) {
            return at(0) - at(n - 1);
        }
    } else {
        if position + 1 < n {
            return at(position + 1) - at(position);
        }
        if matches!(right, BoundaryCondition::Periodic) {
            return at(0) - at(n - 1);
        }
    }
    let face = if position == 0 { left } else { right };
    let flux = match face {
        | BoundaryCondition::Neumann(g) => *g,
        | BoundaryCondition::Robin { alpha, beta, gamma } => (gamma - alpha * at(position)) / beta,
        | _ => 0.0,
    };
    if position == 0 {
        -grid.h * flux
    } else {
        grid.h * flux
    }
}

/// Solves 1D linear advection `u_t + v u_x = 0` on `grid` with boundary conditions `bc`.
///
/// # Errors
///
/// [`FdmError::Domain`] if `|v| dt / h > 1`, [`FdmError::ShapeMismatch`] if
/// `initial.len() != grid.nx`.
pub fn solve_advection_1d(
    grid: Grid1D,
    v: f64,
    dt: f64,
    steps: usize,
    initial: &[f64],
    bc: &Boundary1D,
) -> Result<Vec<Field1D>, FdmError> {
    grid.check()?;
    bc.validate()?;
    let n = grid.nx;
    if initial.len() != n {
        return Err(FdmError::ShapeMismatch(
            "1D advection initial field does not fit grid".into(),
        ));
    }
    let cfl = v * dt / grid.h;
    if !cfl.is_finite() || cfl.abs() > 1.0 + 1e-12 || dt <= 0.0 {
        return Err(FdmError::Domain(format!(
            "1D advection upwind requires |v|*dt/h <= 1, found {cfl}"
        )));
    }
    let periodic = matches!(bc.left, BoundaryCondition::Periodic);
    let mut u = initial.to_vec();
    if !periodic {
        enforce_dirichlet_1d(&mut u, bc);
    }
    let mut out = Vec::with_capacity(steps + 1);
    out.push(Field1D {
        grid,
        values: u.clone(),
    });

    for _ in 0..steps {
        let mut next = u.clone();
        if periodic {
            for i in 0..n {
                let im1 = if i == 0 { n - 1 } else { i - 1 };
                let ip1 = if i + 1 == n { 0 } else { i + 1 };
                let du = if cfl >= 0.0 {
                    u[i] - u[im1]
                } else {
                    u[ip1] - u[i]
                };
                next[i] = u[i] - cfl * du;
            }
        } else {
            for i in 1..n - 1 {
                let du = if cfl >= 0.0 {
                    u[i] - u[i - 1]
                } else {
                    u[i + 1] - u[i]
                };
                next[i] = u[i] - cfl * du;
            }
            if cfl >= 0.0 && matches!(bc.right, BoundaryCondition::Neumann(_)) {
                next[n - 1] = u[n - 1] - cfl * (u[n - 1] - u[n - 2]);
            }
            if cfl < 0.0 && matches!(bc.left, BoundaryCondition::Neumann(_)) {
                next[0] = u[0] - cfl * (u[1] - u[0]);
            }
            enforce_dirichlet_1d(&mut next, bc);
        }
        u = next;
        out.push(Field1D {
            grid,
            values: u.clone(),
        });
    }
    Ok(out)
}

/// Solves 2D linear advection `u_t + v·∇u = 0` with general boundary conditions `bc`
/// (including `Periodic`, `Dirichlet`, and outflow `Neumann`).
///
/// # Errors
///
/// [`FdmError::Domain`] if `(|vx| + |vy|) dt / h > 1`, [`FdmError::ShapeMismatch`]
/// if `spec.initial.len() != spec.grid.len()`.
pub fn solve_advection_with_bc(
    spec: &AdvectionSpec,
    bc: &Boundary2D,
) -> Result<Vec<Field>, FdmError> {
    spec.grid.check()?;
    bc.validate()?;
    if spec.initial.len() != spec.grid.len() {
        return Err(FdmError::ShapeMismatch(
            "the initial field does not fit the grid".into(),
        ));
    }
    let grid = Grid3D {
        nx: spec.grid.nx,
        ny: spec.grid.ny,
        nz: 1,
        h: spec.grid.h,
    };
    let b3 = Boundary3D {
        x_min: bc.left.clone(),
        x_max: bc.right.clone(),
        y_min: bc.bottom.clone(),
        y_max: bc.top.clone(),
        z_min: BoundaryCondition::Periodic,
        z_max: BoundaryCondition::Periodic,
    };
    let c = [spec.vx * spec.dt / grid.h, spec.vy * spec.dt / grid.h];
    if c.iter().any(|v| !v.is_finite())
        || c[0].abs() + c[1].abs() > 1.0
        || spec.dt <= 0.0
        || spec.initial.iter().any(|v| !v.is_finite())
    {
        return Err(FdmError::Domain(
            "invalid or unstable advection inputs".into(),
        ));
    }
    let mut u = spec.initial.clone();
    enforce_dirichlet_3d(&mut u, grid, &b3);
    let mut output = vec![Field {
        grid: spec.grid,
        values: u.clone(),
    }];
    for _ in 0..spec.steps {
        let mut next = u.clone();
        for y in 0..grid.ny {
            for x in 0..grid.nx {
                if fixed_value_3d(grid, &b3, x, y, 0).is_some() {
                    continue;
                }
                let i = grid.index(x, y, 0);
                for a in 0..2 {
                    next[i] -= c[a] * upwind_delta_3d(&u, grid, &b3, [x, y, 0], a, c[a]);
                }
            }
        }
        enforce_dirichlet_3d(&mut next, grid, &b3);
        u = next;
        output.push(Field {
            grid: spec.grid,
            values: u.clone(),
        });
    }
    Ok(output)
}

/// Solves 3D linear advection `u_t + v·∇u = 0` on `grid` with velocity `(vx, vy, vz)`
/// and boundary conditions `bc`.
///
/// # Errors
///
/// [`FdmError::Domain`] if `(|vx| + |vy| + |vz|) dt / h > 1`,
/// [`FdmError::ShapeMismatch`] if `initial.len() != grid.len()`.
pub fn solve_advection_3d(
    grid: Grid3D,
    velocity: (f64, f64, f64),
    dt: f64,
    steps: usize,
    initial: &[f64],
    bc: &Boundary3D,
) -> Result<Vec<Field3D>, FdmError> {
    grid.check()?;
    bc.validate()?;
    if initial.len() != grid.len() {
        return Err(FdmError::ShapeMismatch(
            "3D advection initial field does not fit grid".into(),
        ));
    }
    let c = [
        velocity.0 * dt / grid.h,
        velocity.1 * dt / grid.h,
        velocity.2 * dt / grid.h,
    ];
    if c.iter().any(|v| !v.is_finite())
        || c.iter().map(|v| v.abs()).sum::<f64>() > 1.0
        || dt <= 0.0
        || initial.iter().any(|v| !v.is_finite())
    {
        return Err(FdmError::Domain(
            "invalid or unstable 3D advection inputs".into(),
        ));
    }
    let mut u = initial.to_vec();
    enforce_dirichlet_3d(&mut u, grid, bc);
    let mut output = vec![Field3D {
        grid,
        values: u.clone(),
    }];
    for _ in 0..steps {
        let mut next = u.clone();
        for z in 0..grid.nz {
            for y in 0..grid.ny {
                for x in 0..grid.nx {
                    if fixed_value_3d(grid, bc, x, y, z).is_some() {
                        continue;
                    }
                    let i = grid.index(x, y, z);
                    for a in 0..3 {
                        next[i] -= c[a] * upwind_delta_3d(&u, grid, bc, [x, y, z], a, c[a]);
                    }
                }
            }
        }
        enforce_dirichlet_3d(&mut next, grid, bc);
        u = next;
        output.push(Field3D {
            grid,
            values: u.clone(),
        });
    }
    Ok(output)
}

// ---------------------------------------------------------------------------
// Finite Element Method (FEM) Solvers (1D, 2D, 3D via rssn::sim::fem_*)
// ---------------------------------------------------------------------------
