//! # Common boundary-condition description
//!
//! One vocabulary shared by every continuum solver of [`crate::sim`]:
//! structured finite differences, multigrid, spectral Chebyshev, finite
//! elements, boundary-value fluids and elasticity.
//!
//! | Type | Meaning |
//! |------|---------|
//! | [`BcValue`] | boundary datum `g(t, x)`: constant, `g(t)`, `g(x)` or a space-time closure |
//! | [`BoundaryCondition`] | scalar condition: Dirichlet `u = g`, Neumann `∂u/∂n = g`, Robin `αu + β∂u/∂n = g`, periodic, symmetry, absorbing |
//! | [`Face`] / [`FaceBcs`] | the six faces of a box, one condition per face (mixed per face; 1D uses X faces, 2D X and Y faces) |
//! | [`MarkerBcs`] | per-marker assignment for unstructured meshes (marker id → condition) |
//! | [`FluidBoundary`] / [`FluidFaceBcs`] | no-slip, slip, inflow, outflow, periodic for incompressible flow |
//! | [`ElasticBoundary`] | displacement (per component), traction, free for elasticity |
//!
//! `∂/∂n` is always the derivative along the **outward** unit normal.


use std::collections::HashMap;
use std::fmt;
use std::sync::Arc;

/// Closure type behind a [`BcValue`]: `g(t, [x, y, z])`.
pub type SpaceTimeFn = dyn Fn(f64, [f64; 3]) -> f64 + Send + Sync;

/// Boundary datum `g(t, x)`, possibly time and space dependent.
#[derive(Clone)]
pub struct BcValue(Arc<SpaceTimeFn>);

impl fmt::Debug for BcValue {
    fn fmt(
        &self,
        f: &mut fmt::Formatter<'_>,
    ) -> fmt::Result {
        f.write_str("BcValue(..)")
    }
}

impl BcValue {
    /// Constant datum.
    #[must_use]
    pub fn constant(v: f64) -> Self {
        Self(Arc::new(move |_, _| v))
    }

    /// Time-dependent datum `g(t)`.
    pub fn time(g: impl Fn(f64) -> f64 + Send + Sync + 'static) -> Self {
        Self(Arc::new(move |t, _| g(t)))
    }

    /// Space-dependent datum `g(x)`.
    pub fn space(g: impl Fn([f64; 3]) -> f64 + Send + Sync + 'static) -> Self {
        Self(Arc::new(move |_, x| g(x)))
    }

    /// Space-time datum `g(t, x)`.
    pub fn spacetime(g: impl Fn(f64, [f64; 3]) -> f64 + Send + Sync + 'static) -> Self {
        Self(Arc::new(g))
    }

    /// Evaluates the datum.
    #[must_use]
    pub fn eval(
        &self,
        t: f64,
        x: [f64; 3],
    ) -> f64 {
        (self.0)(t, x)
    }
}

/// Scalar boundary condition on one face or marker.
#[derive(Clone, Debug)]
pub enum BoundaryCondition {
    /// `u = g`.
    Dirichlet(BcValue),
    /// `∂u/∂n = g` (outward normal).
    Neumann(BcValue),
    /// `alpha u + beta ∂u/∂n = g`.
    Robin {
        /// Coefficient of `u`.
        alpha: f64,
        /// Coefficient of the normal derivative.
        beta: f64,
        /// Datum.
        g: BcValue,
    },
    /// Periodic wrap-around; must be set on both faces of the axis.
    Periodic,
    /// Mirror symmetry, equal to homogeneous Neumann.
    Symmetry,
    /// First-order radiation condition `u_t + c ∂u/∂n = 0` for waves
    /// (only honoured by the wave solvers).
    Absorbing,
}

impl BoundaryCondition {
    /// Homogeneous or constant Dirichlet condition.
    #[must_use]
    pub fn dirichlet(v: f64) -> Self {
        Self::Dirichlet(BcValue::constant(v))
    }

    /// Homogeneous or constant Neumann condition.
    #[must_use]
    pub fn neumann(v: f64) -> Self {
        Self::Neumann(BcValue::constant(v))
    }

    /// Robin condition with constant data.
    #[must_use]
    pub fn robin(
        alpha: f64,
        beta: f64,
        g: f64,
    ) -> Self {
        Self::Robin {
            alpha,
            beta,
            g: BcValue::constant(g),
        }
    }

    /// Whether this is a Dirichlet condition.
    #[must_use]
    pub const fn is_dirichlet(&self) -> bool {
        matches!(self, Self::Dirichlet(_))
    }

    /// Whether this is periodic.
    #[must_use]
    pub const fn is_periodic(&self) -> bool {
        matches!(self, Self::Periodic)
    }
}

/// Face of an axis-aligned box.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Face {
    /// `x = x_min`.
    XMin,
    /// `x = x_max`.
    XMax,
    /// `y = y_min`.
    YMin,
    /// `y = y_max`.
    YMax,
    /// `z = z_min`.
    ZMin,
    /// `z = z_max`.
    ZMax,
}

impl Face {
    /// All faces in storage order.
    pub const ALL: [Self; 6] = [
        Self::XMin,
        Self::XMax,
        Self::YMin,
        Self::YMax,
        Self::ZMin,
        Self::ZMax,
    ];

    /// Axis index 0, 1 or 2.
    #[must_use]
    pub const fn axis(self) -> usize {
        match self {
            | Self::XMin | Self::XMax => 0,
            | Self::YMin | Self::YMax => 1,
            | Self::ZMin | Self::ZMax => 2,
        }
    }

    /// Whether this is the upper face of its axis.
    #[must_use]
    pub const fn is_max(self) -> bool {
        matches!(self, Self::XMax | Self::YMax | Self::ZMax)
    }

    /// Storage index `2 * axis + is_max`.
    #[must_use]
    pub const fn index(self) -> usize {
        2 * self.axis() + self.is_max() as usize
    }

    /// Face of `axis` (lower or upper).
    #[must_use]
    pub const fn of_axis(
        axis: usize,
        max: bool,
    ) -> Self {
        Self::ALL[2 * axis + max as usize]
    }
}

/// One scalar condition per face of a box.
#[derive(Clone, Debug)]
pub struct FaceBcs {
    faces: [BoundaryCondition; 6],
}

impl FaceBcs {
    /// The same condition on every face (periodic is allowed).
    #[must_use]
    pub fn uniform(bc: &BoundaryCondition) -> Self {
        Self {
            faces: std::array::from_fn(|_| bc.clone()),
        }
    }

    /// Homogeneous Dirichlet on every face.
    #[must_use]
    pub fn dirichlet_zero() -> Self {
        Self::uniform(&BoundaryCondition::dirichlet(0.0))
    }

    /// Sets the condition on one face.
    #[must_use]
    pub fn with(
        mut self,
        face: Face,
        bc: BoundaryCondition,
    ) -> Self {
        self.faces[face.index()] = bc;
        self
    }

    /// Sets the condition on both faces of an axis.
    #[must_use]
    pub fn with_axis(
        mut self,
        axis: usize,
        bc: &BoundaryCondition,
    ) -> Self {
        self.faces[2 * axis.min(2)] = bc.clone();
        self.faces[2 * axis.min(2) + 1] = bc.clone();
        self
    }

    /// The same description with all boundary data set to zero (used for
    /// error / residual equations in multigrid).
    #[must_use]
    pub fn homogeneous(&self) -> Self {
        let z = || BcValue::constant(0.0);
        Self {
            faces: std::array::from_fn(|i| {
                match &self.faces[i] {
                    | BoundaryCondition::Dirichlet(_) => BoundaryCondition::Dirichlet(z()),
                    | BoundaryCondition::Neumann(_) => BoundaryCondition::Neumann(z()),
                    | BoundaryCondition::Robin { alpha, beta, .. } => {
                        BoundaryCondition::Robin {
                            alpha: *alpha,
                            beta: *beta,
                            g: z(),
                        }
                    },
                    | other => other.clone(),
                }
            }),
        }
    }

    /// Condition on a face.
    #[must_use]
    pub const fn get(
        &self,
        face: Face,
    ) -> &BoundaryCondition {
        &self.faces[face.index()]
    }

    /// Whether `axis` is periodic (both faces must agree).
    ///
    /// # Errors
    /// If only one face of the axis is periodic.
    pub fn axis_periodic(
        &self,
        axis: usize,
    ) -> Result<bool, String> {
        let a = self.get(Face::of_axis(axis, false)).is_periodic();
        let b = self.get(Face::of_axis(axis, true)).is_periodic();
        if a == b {
            Ok(a)
        } else {
            Err(format!("axis {axis}: periodic must be set on both faces"))
        }
    }
}

/// Per-marker assignment for unstructured meshes.
#[derive(Clone, Debug, Default)]
pub struct MarkerBcs<T> {
    map: HashMap<i32, T>,
    /// Condition used for markers that are not listed.
    pub default: Option<T>,
}

impl<T: Clone> MarkerBcs<T> {
    /// Empty assignment without default.
    #[must_use]
    pub fn new() -> Self {
        Self {
            map: HashMap::new(),
            default: None,
        }
    }

    /// Assignment where unlisted markers use `default`.
    #[must_use]
    pub fn with_default(default: T) -> Self {
        Self {
            map: HashMap::new(),
            default: Some(default),
        }
    }

    /// Assigns a condition to a marker.
    #[must_use]
    pub fn set(
        mut self,
        marker: i32,
        bc: T,
    ) -> Self {
        self.map.insert(marker, bc);
        self
    }

    /// Condition of a marker, or the default.
    #[must_use]
    pub fn get(
        &self,
        marker: i32,
    ) -> Option<&T> {
        self.map.get(&marker).or(self.default.as_ref())
    }
}

/// Fluid boundary type on one face.
#[derive(Clone, Debug)]
pub enum FluidBoundary {
    /// Solid wall moving with the given velocity (zero for a fixed wall).
    NoSlip([BcValue; 3]),
    /// Free-slip wall / symmetry plane: `u·n = 0`, no tangential stress.
    Slip,
    /// Prescribed velocity inflow (Dirichlet velocity, Neumann pressure).
    Inflow([BcValue; 3]),
    /// Outflow: zero normal-derivative of velocity, prescribed pressure.
    Outflow(BcValue),
    /// Periodic wrap-around.
    Periodic,
}

impl FluidBoundary {
    /// Fixed no-slip wall.
    #[must_use]
    pub fn wall() -> Self {
        Self::NoSlip([
            BcValue::constant(0.0),
            BcValue::constant(0.0),
            BcValue::constant(0.0),
        ])
    }

    /// Wall moving with a constant velocity.
    #[must_use]
    pub fn moving_wall(v: [f64; 3]) -> Self {
        Self::NoSlip(v.map(BcValue::constant))
    }

    /// Inflow with constant velocity.
    #[must_use]
    pub fn inflow(v: [f64; 3]) -> Self {
        Self::Inflow(v.map(BcValue::constant))
    }

    /// Outflow at zero pressure.
    #[must_use]
    pub fn outflow() -> Self {
        Self::Outflow(BcValue::constant(0.0))
    }
}

/// One fluid boundary per face.
#[derive(Clone, Debug)]
pub struct FluidFaceBcs {
    faces: [FluidBoundary; 6],
}

impl FluidFaceBcs {
    /// The same condition on every face.
    #[must_use]
    pub fn uniform(b: &FluidBoundary) -> Self {
        Self {
            faces: std::array::from_fn(|_| b.clone()),
        }
    }

    /// Sets one face.
    #[must_use]
    pub fn with(
        mut self,
        face: Face,
        b: FluidBoundary,
    ) -> Self {
        self.faces[face.index()] = b;
        self
    }

    /// Sets both faces of an axis.
    #[must_use]
    pub fn with_axis(
        mut self,
        axis: usize,
        b: &FluidBoundary,
    ) -> Self {
        self.faces[2 * axis.min(2)] = b.clone();
        self.faces[2 * axis.min(2) + 1] = b.clone();
        self
    }

    /// Condition of a face.
    #[must_use]
    pub const fn get(
        &self,
        face: Face,
    ) -> &FluidBoundary {
        &self.faces[face.index()]
    }
}

/// Elasticity boundary condition on a face or marker.
#[derive(Clone, Debug)]
pub enum ElasticBoundary {
    /// Prescribed displacement per component (`None` leaves the component free,
    /// e.g. a roller).
    Displacement([Option<BcValue>; 3]),
    /// Prescribed traction vector `σ n = t`.
    Traction([BcValue; 3]),
    /// Traction-free surface.
    Free,
}

impl ElasticBoundary {
    /// Fully clamped boundary.
    #[must_use]
    pub fn clamped() -> Self {
        Self::Displacement([
            Some(BcValue::constant(0.0)),
            Some(BcValue::constant(0.0)),
            Some(BcValue::constant(0.0)),
        ])
    }

    /// Constant traction.
    #[must_use]
    pub fn traction(t: [f64; 3]) -> Self {
        Self::Traction(t.map(BcValue::constant))
    }

    /// Constant normal-pressure traction `σ n = -p n` is built by the caller from `n`;
    /// this helper fixes only one displacement component to zero (a roller).
    #[must_use]
    pub fn roller(component: usize) -> Self {
        let mut d = [None, None, None];
        if let Some(slot) = d.get_mut(component) {
            *slot = Some(BcValue::constant(0.0));
        }
        Self::Displacement(d)
    }
}
