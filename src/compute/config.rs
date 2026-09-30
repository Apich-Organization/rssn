//! # RSSN Compute Configuration
//!
//! Provides `ComputeConfig` and `TargetRepresentation` to configure the unified
//! identity transformation pipeline.

use std::collections::HashMap;
use crate::symbolic::egraph::heuristics::HeuristicBudget;

/// The target representation domain for the identity transformation.
#[derive(Clone, Debug, PartialEq)]
pub enum TargetRepresentation {
    /// Canonical closed-form symbolic expression.
    Symbolic,
    /// Numerical scalar or tensor evaluation with specified relative tolerance.
    Numerical {
        /// Absolute/relative tolerance for numerical convergence.
        tolerance: f64,
        /// Maximum iterations for iterative solvers.
        max_iterations: usize,
    },
    /// Exact discrete/algebraic representation (e.g. rationals, polynomial rings).
    Exact,
    /// Automatically infer target: numerical if all variables are bound, symbolic otherwise.
    Auto,
}

impl Default for TargetRepresentation {
    fn default() -> Self {
        Self::Auto
    }
}

/// Configuration for the unified `compute` pipeline.
#[derive(Clone, Debug)]
pub struct ComputeConfig {
    /// Desired target representation kind.
    pub target: TargetRepresentation,
    /// Execution and memory budget for E-Graph equality saturation.
    pub budget: HeuristicBudget,
    /// Whether calculus rules (differentiation, integration, limits, series, vector calculus) are enabled.
    pub enable_calculus: bool,
    /// Whether algebraic rules (constant folding, polynomial algebra, Gröbner) are enabled.
    pub enable_algebra: bool,
    /// Whether ordinary differential equation (ODE) solving rules are enabled.
    pub enable_ode: bool,
    /// Whether partial differential equation (PDE) solving rules are enabled.
    pub enable_pde: bool,
    /// Whether linear algebra and matrix rules are enabled.
    pub enable_matrix: bool,
    /// Whether integral transform rules (Laplace, Fourier, Z-transform) are enabled.
    pub enable_transforms: bool,
    /// Whether special function rules (Bessel, Gamma, Zeta, Erf) are enabled.
    pub enable_special_functions: bool,
    /// Whether optimization rules are enabled.
    pub enable_optimization: bool,
    /// Whether statistics and probability rules are enabled.
    pub enable_stats: bool,
    /// Physics simulation and quantum operator rules enabled.
    pub enable_physics: bool,
    /// Concrete variable bindings for evaluation (e.g., "x" -> 3.0).
    pub bindings: HashMap<String, f64>,
    /// Initial condition mappings for ODE solving (e.g., "y" -> 1.0).
    pub ode_initial_conditions: HashMap<String, f64>,
    /// Integration range (x0, x_end) for numerical ODE solving.
    pub ode_range: Option<(f64, f64)>,
    /// Number of integration steps for numerical ODE stepper.
    pub ode_steps: Option<usize>,
}

impl Default for ComputeConfig {
    fn default() -> Self {
        Self::new()
    }
}

impl ComputeConfig {
    /// Creates a default compute configuration.
    #[must_use]
    pub fn new() -> Self {
        Self {
            target: TargetRepresentation::Auto,
            budget: HeuristicBudget::default(),
            enable_calculus: false,
            enable_algebra: true, // Basic algebraic simplification on by default
            enable_ode: false,
            enable_pde: false,
            enable_matrix: false,
            enable_transforms: false,
            enable_special_functions: false,
            enable_optimization: false,
            enable_stats: false,
            enable_physics: false,
            bindings: HashMap::new(),
            ode_initial_conditions: HashMap::new(),
            ode_range: None,
            ode_steps: None,
        }
    }

    /// Enables calculus rules (differentiation, integration, limits, series, vector calculus).
    #[must_use]
    pub fn with_calculus(mut self) -> Self {
        self.enable_calculus = true;
        self
    }

    /// Enables ODE solving rules (symbolic separation/linear + numerical Runge-Kutta).
    #[must_use]
    pub fn with_ode(mut self) -> Self {
        self.enable_ode = true;
        self
    }

    /// Enables PDE solving rules.
    #[must_use]
    pub fn with_pde(mut self) -> Self {
        self.enable_pde = true;
        self
    }

    /// Enables advanced algebraic rules (Gröbner bases, polynomial factorization, radicals).
    #[must_use]
    pub fn with_algebra(mut self) -> Self {
        self.enable_algebra = true;
        self
    }

    /// Enables linear algebra and matrix rules (symbolic + numerical SVD, QR, Eigen, Inversion).
    #[must_use]
    pub fn with_matrix(mut self) -> Self {
        self.enable_matrix = true;
        self
    }

    /// Enables integral transforms (Laplace, Fourier, Z-transform).
    #[must_use]
    pub fn with_transforms(mut self) -> Self {
        self.enable_transforms = true;
        self
    }

    /// Enables special functions rules (Bessel, Gamma, Zeta, Erf).
    #[must_use]
    pub fn with_special_functions(mut self) -> Self {
        self.enable_special_functions = true;
        self
    }

    /// Enables optimization rules (gradient descent, Newton, BFGS).
    #[must_use]
    pub fn with_optimization(mut self) -> Self {
        self.enable_optimization = true;
        self
    }

    /// Enables probability and statistics rules.
    #[must_use]
    pub fn with_stats(mut self) -> Self {
        self.enable_stats = true;
        self
    }

    /// Enables physics simulation and quantum operator rules.
    #[must_use]
    pub fn with_physics(mut self) -> Self {
        self.enable_physics = true;
        self
    }

    /// Injects all available domain transformation rules.
    #[must_use]
    pub fn with_all(mut self) -> Self {
        self.enable_calculus = true;
        self.enable_algebra = true;
        self.enable_ode = true;
        self.enable_pde = true;
        self.enable_matrix = true;
        self.enable_transforms = true;
        self.enable_special_functions = true;
        self.enable_optimization = true;
        self.enable_stats = true;
        self.enable_physics = true;
        self
    }

    /// Explicitly targets a symbolic closed-form result.
    #[must_use]
    pub fn target_symbolic(mut self) -> Self {
        self.target = TargetRepresentation::Symbolic;
        self
    }

    /// Explicitly targets a numerical result with specified tolerance.
    #[must_use]
    pub fn target_numerical(mut self, tolerance: f64) -> Self {
        self.target = TargetRepresentation::Numerical {
            tolerance,
            max_iterations: 1000,
        };
        self
    }

    /// Binds a variable name to a concrete numerical float value.
    #[must_use]
    pub fn bind(mut self, var: impl Into<String>, val: f64) -> Self {
        self.bindings.insert(var.into(), val);
        self
    }

    /// Sets the heuristic saturation budget.
    #[must_use]
    pub fn with_budget(mut self, budget: HeuristicBudget) -> Self {
        self.budget = budget;
        self
    }

    /// Sets an initial condition for ODE solving (e.g. y(0) = 1.0).
    #[must_use]
    pub fn with_initial_condition(mut self, func: impl Into<String>, val: f64) -> Self {
        self.ode_initial_conditions.insert(func.into(), val);
        self
    }

    /// Sets the integration range `(x0, x_end)` for numerical ODE solvers.
    #[must_use]
    pub fn with_ode_range(mut self, x0: f64, x_end: f64) -> Self {
        self.ode_range = Some((x0, x_end));
        self
    }

    /// Sets the number of integration steps for numerical ODE solvers.
    #[must_use]
    pub fn with_ode_steps(mut self, steps: usize) -> Self {
        self.ode_steps = Some(steps);
        self
    }
}
