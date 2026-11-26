//! Numerical optimization on the `basin` framework.
//!
//! [`minimize`] runs one of the solvers listed by [`SolverKind`] on a cost
//! function with an optional analytic gradient and optional box bounds, and
//! returns a plain [`MinimizeResult`]. A problem is described by [`FnProblem`],
//! so a caller does not have to implement any trait.
//!
//! The wrapper covers the solver families `basin` provides:
//!
//! - **First order**: [`SolverKind::GradientDescent`],
//!   [`SolverKind::NonlinearCg`], [`SolverKind::TrustRegion`],
//!   [`SolverKind::ProjectedGradientDescent`] and [`SolverKind::Tnc`].
//! - **Quasi-Newton**: [`SolverKind::Bfgs`] and [`SolverKind::Lbfgs`].
//! - **Derivative free**: [`SolverKind::NelderMead`],
//!   [`SolverKind::NelderMeadProjected`], [`SolverKind::Bobyqa`],
//!   [`SolverKind::Newuoa`], [`SolverKind::Lincoa`], [`SolverKind::Cobyla`],
//!   [`SolverKind::Direct`], [`SolverKind::Mads`] and [`SolverKind::Slsqp`].
//! - **Global and stochastic**: [`SolverKind::CmaEs`], [`SolverKind::De`],
//!   [`SolverKind::GlobalBestPso`], [`SolverKind::RandomSearch`],
//!   [`SolverKind::SolisWets`], [`SolverKind::Ssga`] and
//!   [`SolverKind::Gbnm`].
//! - **Scalar**: [`SolverKind::Brent`] and [`SolverKind::GoldenSection`] for a
//!   one-dimensional parameter.
//!
//! Box bounds belong to the problem, not the solver: a bounded solver reads
//! them from the problem, an unbounded solver ignores them, and a solver that
//! needs bounds but receives none fails with a clear error.
//!
//! `basin` is re-exported so a caller can reach a facility this wrapper does
//! not wrap (least-squares solvers, robust losses, constrained blocks,
//! observers, checkpointing).

use std::sync::Arc;

use basin::Bfgs;
use basin::BisectionRoot;
use basin::Bobyqa;
use basin::BoundedCmaEs;
use basin::BoxBounds;
use basin::Brent;
use basin::BrentRoot;
use basin::CmaEs;
use basin::Cobyla;
use basin::ConstraintDescription;
use basin::ConstraintSet;
use basin::Constraints;
use basin::CostFunction;
use basin::De;
use basin::DenseMatrix;
use basin::Direct;
use basin::Executor;
use basin::Gbnm;
use basin::GlobalBestPso;
use basin::GoldenSection;
use basin::Gradient;
use basin::GradientDescent;
use basin::HalleyRoot;
use basin::IncumbentState;
use basin::Lbfgs;
use basin::line_search::MoreThuente;
use basin::Lbfgsb;
use basin::Lincoa;
use basin::Mads;
use basin::NelderMead;
use basin::NewtonRoot;
use basin::OptimizationResult;
use basin::Newuoa;
use basin::NonlinearCg;
use basin::PointState;
use basin::PopulationProgress;
use basin::ProjectedGradientDescent;
use basin::RandomSearch;
use basin::SecantRoot;
use basin::Slsqp;
use basin::SolisWets;
use basin::Ssga;
use basin::Tnc;
use basin::Toms748Root;

pub use basin;

/// Bound used for a side a caller leaves free.
const UNBOUNDED: f64 = f64::MAX;

/// A cost function with an optional analytic gradient and box bounds.
///
/// Clones share their closures, so one problem can feed several runs.
#[derive(Clone)]
pub struct FnProblem {
    cost: Arc<dyn Fn(&[f64]) -> Result<f64, String> + Send + Sync>,
    gradient: Option<Arc<dyn Fn(&[f64]) -> Result<Vec<f64>, String> + Send + Sync>>,
    set: ConstraintSet<BoxBounds<Vec<f64>>>,
}

impl FnProblem {
    /// Builds a problem from a cost closure and an optional gradient closure.
    ///
    /// `bounds` is `(lower, upper)`; a free side uses [`f64::INFINITY`]. When
    /// `bounds` is `None` every variable is free.
    #[must_use]
    pub fn new(
        cost: impl Fn(&[f64]) -> f64 + Send + Sync + 'static,
        gradient: Option<impl Fn(&[f64]) -> Vec<f64> + Send + Sync + 'static>,
        dimension: usize,
        bounds: Option<(Vec<f64>, Vec<f64>)>,
    ) -> Self {
        Self::try_new(
            move |x| Ok(cost(x)),
            gradient
                .map(|gradient| move |x: &[f64]| -> Result<Vec<f64>, String> { Ok(gradient(x)) }),
            dimension,
            bounds,
        )
    }

    /// Builds a problem from fallible closures.
    #[must_use]
    pub fn try_new(
        cost: impl Fn(&[f64]) -> Result<f64, String> + Send + Sync + 'static,
        gradient: Option<impl Fn(&[f64]) -> Result<Vec<f64>, String> + Send + Sync + 'static>,
        dimension: usize,
        bounds: Option<(Vec<f64>, Vec<f64>)>,
    ) -> Self {
        let (lower, upper) =
            bounds.unwrap_or_else(|| (vec![-UNBOUNDED; dimension], vec![UNBOUNDED; dimension]));
        Self {
            cost: Arc::new(cost),
            gradient: gradient.map(|gradient| {
                let closure: Arc<dyn Fn(&[f64]) -> Result<Vec<f64>, String> + Send + Sync> =
                    Arc::new(gradient);
                closure
            }),
            set: ConstraintSet::default().with_bounds(BoxBounds::new(Some(lower), Some(upper))),
        }
    }

    /// Whether the problem carries an analytic gradient.
    #[must_use]
    pub const fn has_gradient(&self) -> bool {
        self.gradient.is_some()
    }

    /// The same closures over a different box.
    #[must_use]
    pub fn with_bounds(
        &self,
        lower: Vec<f64>,
        upper: Vec<f64>,
    ) -> Self {
        Self {
            cost: Arc::clone(&self.cost),
            gradient: self.gradient.as_ref().map(Arc::clone),
            set: ConstraintSet::default().with_bounds(BoxBounds::new(Some(lower), Some(upper))),
        }
    }

    /// Cost and gradient at `x`.
    ///
    /// # Errors
    /// Returns an error when the cost fails or the problem has no gradient.
    pub fn value_and_gradient(
        &self,
        x: &[f64],
    ) -> Result<(f64, Vec<f64>), String> {
        let value = (self.cost)(x)?;
        match &self.gradient {
            | Some(gradient) => Ok((value, gradient(x)?)),
            | None => Err("the problem has no analytic gradient".to_string()),
        }
    }
}

impl CostFunction for FnProblem {
    type Error = String;
    type Output = f64;
    type Param = Vec<f64>;

    fn cost(
        &self,
        param: &Vec<f64>,
    ) -> Result<f64, String> {
        (self.cost)(param)
    }
}

impl Gradient for FnProblem {
    type Gradient = Vec<f64>;

    fn gradient(
        &self,
        param: &Vec<f64>,
    ) -> Result<Vec<f64>, String> {
        match &self.gradient {
            | Some(gradient) => gradient(param),
            | None => Err("this solver needs an analytic gradient".to_string()),
        }
    }
}

impl Constraints for FnProblem {
    type Set = ConstraintSet<BoxBounds<Vec<f64>>>;

    fn constraints(&self) -> <Self::Set as ConstraintDescription>::Borrowed<'_> {
        self.set.borrowed()
    }
}

/// The optimization problems the bundled test functions describe.
#[derive(Debug, Clone, Copy)]
pub enum ProblemType {
    /// Rosenbrock's banana function.
    Rosenbrock,
    /// Sum of squares.
    Sphere,
    /// Rastrigin's multimodal function.
    Rastrigin,
    /// Ackley's multimodal function.
    Ackley,
    /// A function the caller supplies.
    Custom,
}

/// Iteration budget and convergence target for a run.
#[derive(Debug, Clone)]
pub struct OptimizationConfig {
    /// Maximum number of iterations.
    pub max_iters: u64,
    /// Cost at or below which the run counts as converged.
    pub tolerance: f64,
    /// Which test problem the configuration describes.
    pub problem_type: ProblemType,
    /// Number of variables.
    pub dimension: usize,
}

impl Default for OptimizationConfig {
    fn default() -> Self {
        Self {
            max_iters: 1000,
            tolerance: 1e-6,
            problem_type: ProblemType::Rosenbrock,
            dimension: 2,
        }
    }
}

/// One `basin` solver, with the parameters that change its behaviour.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SolverKind {
    /// Steepest descent with a constant step.
    GradientDescent {
        /// Step length.
        step: f64,
    },
    /// Dense BFGS quasi-Newton.
    Bfgs,
    /// Limited-memory BFGS.
    Lbfgs,
    /// Limited-memory BFGS clamped to the box bounds (needs bounds).
    Lbfgsb,
    /// Nonlinear conjugate gradients (Hager–Zhang line search).
    NonlinearCg,
    /// Nelder–Mead simplex.
    NelderMead,
    /// Nelder–Mead clamped to the box bounds.
    NelderMeadProjected,
    /// BOBYQA model-based trust region (needs bounds).
    Bobyqa,
    /// NEWUOA model-based trust region.
    Newuoa,
    /// LINCOA with linear constraints (needs bounds).
    Lincoa,
    /// COBYLA linear approximation (needs bounds).
    Cobyla,
    /// DIRECT deterministic global search (needs bounds).
    Direct,
    /// Mesh-adaptive direct search (needs bounds).
    Mads,
    /// Sequential least squares programming (needs bounds).
    Slsqp,
    /// Projected gradient descent (needs bounds).
    ProjectedGradientDescent {
        /// Step length.
        alpha: f64,
    },
    /// Truncated Newton with bounds.
    Tnc,
    /// CMA-ES, a stochastic evolution strategy.
    CmaEs {
        /// Initial step size.
        sigma: f64,
    },
    /// CMA-ES clamped to the box bounds (needs bounds).
    BoundedCmaEs {
        /// Initial step size.
        sigma: f64,
    },
    /// Differential evolution.
    De,
    /// Global-best particle swarm (needs bounds).
    GlobalBestPso,
    /// Random search with `lambda` samples per iteration.
    RandomSearch {
        /// Samples per iteration.
        lambda: usize,
    },
    /// Solis–Wets adaptive random search.
    SolisWets,
    /// Steady-state genetic algorithm.
    Ssga,
    /// Globalized and bounded Nelder–Mead.
    Gbnm,
    /// Brent's scalar minimiser (one variable).
    Brent,
    /// Golden-section scalar minimiser (one variable).
    GoldenSection,
}

impl SolverKind {
    /// Every solver this wrapper exposes.
    pub const ALL: [Self; 26] = [
        Self::GradientDescent { step: 0.1 },
        Self::Bfgs,
        Self::Lbfgs,
        Self::Lbfgsb,
        Self::NonlinearCg,
        Self::NelderMead,
        Self::NelderMeadProjected,
        Self::Bobyqa,
        Self::Newuoa,
        Self::Lincoa,
        Self::Cobyla,
        Self::Direct,
        Self::Mads,
        Self::Slsqp,
        Self::ProjectedGradientDescent { alpha: 0.1 },
        Self::Tnc,
        Self::CmaEs { sigma: 0.3 },
        Self::BoundedCmaEs { sigma: 0.3 },
        Self::De,
        Self::GlobalBestPso,
        Self::RandomSearch { lambda: 10 },
        Self::SolisWets,
        Self::Ssga,
        Self::Gbnm,
        Self::Brent,
        Self::GoldenSection,
    ];

    /// Whether the solver requires finite box bounds.
    #[must_use]
    pub const fn needs_bounds(self) -> bool {
        matches!(
            self,
            Self::Lbfgsb
                | Self::Bobyqa
                | Self::Lincoa
                | Self::Cobyla
                | Self::Direct
                | Self::Mads
                | Self::Slsqp
                | Self::ProjectedGradientDescent { .. }
                | Self::GlobalBestPso
                | Self::BoundedCmaEs { .. }
        )
    }

    /// Whether the solver evaluates the analytic gradient.
    #[must_use]
    pub const fn needs_gradient(self) -> bool {
        matches!(
            self,
            Self::GradientDescent { .. }
                | Self::Bfgs
                | Self::Lbfgs
                | Self::NonlinearCg
                | Self::ProjectedGradientDescent { .. }
                | Self::Tnc
        )
    }
}

/// Execution options for [`minimize`].
#[derive(Debug, Clone)]
pub struct MinimizeOptions {
    /// Iteration budget.
    pub max_iter: u64,
    /// Cost at or below which the result is reported as converged. Passed to
    /// the run as the target objective, so a solver stops as soon as it reaches
    /// it and reports `TargetCost` — except `Cobyla`, `Slsqp` and `Tnc`, whose
    /// seeded state carries no incumbent record for the framework to test.
    pub tolerance: f64,
    /// Euclidean gradient norm (or, for a projected solver, the projected
    /// gradient) at or below which the solver stops. Only the solvers that
    /// report a gradient use it; the default is the square root of the
    /// machine epsilon.
    pub gradient_tolerance: f64,
    /// Seed for the stochastic solvers, so runs are reproducible.
    pub seed: u64,
    /// Box bounds `(lower, upper)`; `None` leaves every variable free.
    pub bounds: Option<(Vec<f64>, Vec<f64>)>,
}

impl Default for MinimizeOptions {
    fn default() -> Self {
        Self {
            max_iter: 1_000,
            tolerance: 1e-6,
            gradient_tolerance: f64::EPSILON.sqrt(),
            seed: 0,
            bounds: None,
        }
    }
}

impl From<&OptimizationConfig> for MinimizeOptions {
    fn from(config: &OptimizationConfig) -> Self {
        Self {
            max_iter: config.max_iters,
            tolerance: config.tolerance,
            ..Self::default()
        }
    }
}

/// What a run produced.
///
/// The fingerprint derives let a result cross the native ABI unchanged, the
/// same way every other record RSSN publishes does.
#[derive(
    Debug,
    Clone,
    PartialEq,
    bincode_next::Encode,
    bincode_next::Decode,
    bincode_next::Fingerprint,
)]
pub struct MinimizeResult {
    /// Best parameter vector found.
    pub x: Vec<f64>,
    /// Cost at `x`.
    pub cost: f64,
    /// Iterations performed.
    pub iterations: u64,
    /// Cost evaluations performed.
    pub cost_evals: u64,
    /// Whether `cost` reached the requested tolerance.
    pub converged: bool,
    /// Why the solver stopped, named as the framework names it
    /// (`TargetCost`, `MaxIter`, `GradientTolerance`, …).
    pub termination: String,
}

/// The best iterate a run observed, with its cost.
///
/// A run whose costs are all infinite never establishes an incumbent, so the
/// best-record accessors would panic; the final iterate is the honest answer
/// there.
fn best_point<S>(
    result: &OptimizationResult<S>,
) -> (S::Param, S::Float)
where
    S: IncumbentState,
    S::Param: Clone,
{
    match result.state.incumbent_record() {
        | Some(best) => (best.param.clone(), best.cost),
        | None => (result.param().clone(), result.cost()),
    }
}

/// Runs `solver` on `problem` and converts the framework result.
///
/// The target objective is only set when the tolerance is finite, because the
/// framework asserts on a non-finite target; an unreachable target is the way
/// a caller asks for "run until something else stops you".
macro_rules! solve_with {
    ($problem:expr, $solver:expr, $start:expr, $options:expr) => {{
        let start: Vec<f64> = $start;
        let executor = Executor::from_start($problem, $solver, start).max_iter($options.max_iter);
        let executor = match $options.tolerance.is_finite() {
            | true => executor.target_objective($options.tolerance),
            | false => executor,
        };
        let result = executor.run().map_err(|error| error.to_string())?;
        let (x, cost) = best_point(&result);
        Ok(MinimizeResult {
            x,
            cost,
            iterations: result.iter(),
            cost_evals: result.cost_evals(),
            converged: cost <= $options.tolerance,
            termination: format!("{:?}", result.report.code()),
        })
    }};
}

/// A one-variable view of a [`FnProblem`], for the scalar solvers.
struct ScalarView<'a> {
    problem: &'a FnProblem,
    set: ConstraintSet<BoxBounds<f64>>,
}

impl<'a> ScalarView<'a> {
    fn new(
        problem: &'a FnProblem,
        start: f64,
    ) -> Self {
        let bounds = problem.set.bounds();
        let lower = bounds
            .lower()
            .and_then(|values| values.first().copied())
            .unwrap_or(-UNBOUNDED);
        let upper = bounds
            .upper()
            .and_then(|values| values.first().copied())
            .unwrap_or(UNBOUNDED);
        // A scalar solver brackets inside its box, so a free variable gets a
        // finite bracket around the start instead of the whole line.
        let (lower, upper) = if lower <= -UNBOUNDED && upper >= UNBOUNDED {
            let reach = 10.0 * start.abs().max(1.0);
            (start - reach, start + reach)
        } else {
            (lower, upper)
        };
        Self {
            problem,
            set: ConstraintSet::default().with_bounds(BoxBounds::new(Some(lower), Some(upper))),
        }
    }
}

impl CostFunction for ScalarView<'_> {
    type Error = String;
    type Output = f64;
    type Param = f64;

    fn cost(
        &self,
        x: &f64,
    ) -> Result<f64, String> {
        (self.problem.cost)(std::slice::from_ref(x))
    }
}

impl Constraints for ScalarView<'_> {
    type Set = ConstraintSet<BoxBounds<f64>>;

    fn constraints(&self) -> <Self::Set as ConstraintDescription>::Borrowed<'_> {
        self.set.borrowed()
    }
}

/// As [`solve_with`], for a solver whose seeded state carries no incumbent
/// record: the framework cannot express a target objective for it, so the run
/// stops on the solver's own convergence test or the iteration budget.
macro_rules! solve_without_target {
    ($problem:expr, $solver:expr, $start:expr, $options:expr) => {{
        let start: Vec<f64> = $start;
        let result = Executor::from_start($problem, $solver, start)
            .max_iter($options.max_iter)
            .run()
            .map_err(|error| error.to_string())?;
        let (x, cost) = best_point(&result);
        Ok(MinimizeResult {
            x,
            cost,
            iterations: result.iter(),
            cost_evals: result.cost_evals(),
            converged: cost <= $options.tolerance,
            termination: format!("{:?}", result.report.code()),
        })
    }};
}

/// Runs a solver whose initial state is a single point.
macro_rules! solve_point {
    ($problem:expr, $solver:expr, $start:expr, $options:expr) => {{
        let start: Vec<f64> = $start;
        let executor =
            Executor::new($problem, $solver, PointState::new(start)).max_iter($options.max_iter);
        let executor = match $options.tolerance.is_finite() {
            | true => executor.target_objective($options.tolerance),
            | false => executor,
        };
        let result = executor.run().map_err(|error| error.to_string())?;
        let (x, cost) = best_point(&result);
        Ok(MinimizeResult {
            x,
            cost,
            iterations: result.iter(),
            cost_evals: result.cost_evals(),
            converged: cost <= $options.tolerance,
            termination: format!("{:?}", result.report.code()),
        })
    }};
}

/// Runs a scalar solver whose initial state is a single number.
macro_rules! solve_scalar {
    ($problem:expr, $solver:expr, $start:expr, $options:expr) => {{
        let start: f64 = $start;
        let executor =
            Executor::new($problem, $solver, PointState::new(start)).max_iter($options.max_iter);
        let executor = match $options.tolerance.is_finite() {
            | true => executor.target_objective($options.tolerance),
            | false => executor,
        };
        let result = executor.run().map_err(|error| error.to_string())?;
        let (x, cost) = best_point(&result);
        Ok(MinimizeResult {
            x: vec![x],
            cost,
            iterations: result.iter(),
            cost_evals: result.cost_evals(),
            converged: cost <= $options.tolerance,
            termination: format!("{:?}", result.report.code()),
        })
    }};
}

/// Runs a population solver seeded with one member.
macro_rules! solve_population {
    ($problem:expr, $solver:expr, $start:expr, $options:expr) => {{
        let _ = $start;
        let state = PopulationProgress::empty();
        let executor = Executor::new($problem, $solver, state)
            .require_evaluated_state()
            .max_iter($options.max_iter);
        let executor = match $options.tolerance.is_finite() {
            | true => executor.target_objective($options.tolerance),
            | false => executor,
        };
        let result = executor
            .run_with_solver()
            .map_err(|error| error.to_string())?
            .into_result();
        let (x, cost) = best_point(&result);
        Ok(MinimizeResult {
            x,
            cost,
            iterations: result.iter(),
            cost_evals: result.cost_evals(),
            converged: cost <= $options.tolerance,
            termination: format!("{:?}", result.report.code()),
        })
    }};
}

/// The stationarity test, for a solver that reports a gradient.
///
/// The framework's cost-change test is deliberately not used: it stops when
/// two consecutive costs are within the tolerance, and a simplex or population
/// method reports the same best cost on many consecutive iterations, so that
/// test would stop such a solver on its first repeat. Those families keep
/// their own tolerance test — simplex size, mesh size, search-distribution
/// collapse — which is the test their algorithm actually defines.
macro_rules! stationary {
    ($solver:expr, $options:expr) => {
        $solver.with_absolute_gradient_tolerance($options.gradient_tolerance)
    };
}

/// The stationarity test for a solver whose bounds define the projection.
macro_rules! projected {
    ($solver:expr, $options:expr) => {
        $solver.with_absolute_projected_gradient_tolerance($options.gradient_tolerance)
    };
}

/// Minimises `problem` with `kind`.
///
/// `start` is the initial point; a scalar solver ([`SolverKind::Brent`],
/// [`SolverKind::GoldenSection`]) requires exactly one variable.
///
/// # Errors
/// Returns an error when the problem has no gradient but the solver needs one,
/// when the solver needs bounds and none were supplied, when a scalar solver
/// receives more than one variable, or when the solver aborts.
pub fn minimize(
    problem: &FnProblem,
    start: &[f64],
    kind: SolverKind,
    options: &MinimizeOptions,
) -> Result<MinimizeResult, String> {
    if start.is_empty() {
        return Err("the initial point is empty".to_string());
    }
    if matches!(kind, SolverKind::Brent | SolverKind::GoldenSection) && start.len() != 1 {
        return Err(format!(
            "{kind:?} minimises one variable, but the start has {} entries",
            start.len()
        ));
    }
    if kind.needs_bounds() && options.bounds.is_none() {
        return Err(format!("{kind:?} requires box bounds"));
    }
    if kind.needs_gradient() && !problem.has_gradient() {
        return Err(format!("{kind:?} requires an analytic gradient"));
    }
    let problem = match &options.bounds {
        | Some((lower, upper)) => problem.with_bounds(lower.clone(), upper.clone()),
        | None => problem.clone(),
    };
    let x0 = start.to_vec();
    match kind {
        | SolverKind::GradientDescent { step } => {
            solve_with!(problem, stationary!(GradientDescent::new(step), options), x0, options)
        },
        // A shallow objective needs a long step: the More–Thuente search is
        // left free to take one, which the Wolfe search's default bound rules
        // out.
        | SolverKind::Bfgs => solve_with!(
            problem,
            stationary!(
                Bfgs::with_line_search(MoreThuente::new().stpmax(f64::MAX).maxfev(100)),
                options
            ),
            x0,
            options
        ),
        | SolverKind::Lbfgs => solve_with!(
            problem,
            stationary!(Lbfgs::new().unbounded(), options),
            x0,
            options
        ),
        // The bounded variant exposes no gradient test: its iterate is
        // projected onto the box, so a plain gradient norm says nothing.
        | SolverKind::Lbfgsb => solve_with!(problem, Lbfgsb::new(), x0, options),
        | SolverKind::NonlinearCg => {
            solve_with!(problem, stationary!(NonlinearCg::new(), options), x0, options)
        },
        // Derivative-free, global and scalar families: no framework check, for
        // the reason recorded at `stationary!` — each keeps the tolerance test
        // its own algorithm defines.
        | SolverKind::NelderMead => solve_with!(problem, NelderMead::new(), x0, options),
        | SolverKind::NelderMeadProjected => {
            solve_with!(
                problem,
                NelderMead::new().projected(),
                x0,
                options
            )
        },
        | SolverKind::Bobyqa => solve_with!(problem, Bobyqa::new(), x0, options),
        | SolverKind::Newuoa => solve_with!(problem, Newuoa::new(), x0, options),
        | SolverKind::Lincoa => solve_with!(problem, Lincoa::new(), x0, options),
        // `Cobyla`, `Slsqp` and `Tnc` seed a state with no incumbent record, so
        // the framework cannot express a target objective for them; they stop
        // on their own internal test or the iteration budget.
        | SolverKind::Cobyla => solve_without_target!(problem, Cobyla::new(), x0, options),
        | SolverKind::Direct => solve_point!(problem, Direct::new(), x0, options),
        | SolverKind::Mads => solve_with!(problem, Mads::new(), x0, options),
        | SolverKind::Slsqp => solve_without_target!(problem, Slsqp::new(), x0, options),
        | SolverKind::ProjectedGradientDescent { alpha } => {
            solve_with!(
                problem,
                projected!(ProjectedGradientDescent::new(alpha), options),
                x0,
                options
            )
        },
        | SolverKind::Tnc => solve_without_target!(problem, Tnc::new(), x0, options),
        | SolverKind::CmaEs { sigma } => {
            solve_with!(
                problem,
                CmaEs::<Vec<f64>, DenseMatrix>::new(options.seed, sigma),
                x0,
                options
            )
        },
        | SolverKind::BoundedCmaEs { sigma } => {
            solve_with!(
                problem,
                BoundedCmaEs::<Vec<f64>, DenseMatrix>::new(options.seed, sigma),
                x0,
                options
            )
        },
        | SolverKind::De => {
            solve_population!(problem, De::new(options.seed), x0, options)
        },
        | SolverKind::GlobalBestPso => {
            solve_population!(
                problem,
                GlobalBestPso::new(options.seed),
                x0,
                options
            )
        },
        | SolverKind::RandomSearch { lambda } => {
            solve_population!(
                problem,
                RandomSearch::new(lambda, options.seed),
                x0,
                options
            )
        },
        | SolverKind::SolisWets => {
            solve_with!(
                problem,
                SolisWets::new(options.seed),
                x0,
                options
            )
        },
        | SolverKind::Ssga => {
            solve_population!(problem, Ssga::new(options.seed), x0, options)
        },
        | SolverKind::Gbnm => {
            solve_with!(problem, Gbnm::new(options.seed), x0, options)
        },
        | SolverKind::Brent => {
            let scalar = ScalarView::new(&problem, x0[0]);
            solve_scalar!(scalar, Brent::new(), x0[0], options)
        },
        | SolverKind::GoldenSection => {
            let scalar = ScalarView::new(&problem, x0[0]);
            solve_scalar!(scalar, GoldenSection::new(), x0[0], options)
        },
    }
}

/// A scalar root finder.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootSolver {
    /// Bisection on a sign-changing bracket.
    Bisection,
    /// Brent's method on a sign-changing bracket.
    Brent,
    /// Safeguarded secant iteration.
    Secant,
    /// TOMS Algorithm 748 (`k = 2`).
    Toms748,
    /// Safeguarded Newton iteration; uses a central difference when no
    /// derivative is supplied.
    Newton,
    /// Safeguarded Halley iteration; uses central differences when no
    /// derivatives are supplied.
    Halley,
}

/// Options for [`find_root`].
#[derive(Debug, Clone, Copy)]
pub struct RootOptions {
    /// Iteration budget.
    pub max_iter: u64,
    /// Absolute position tolerance.
    pub tolerance: f64,
}

impl Default for RootOptions {
    fn default() -> Self {
        Self {
            max_iter: 200,
            tolerance: 1e-12,
        }
    }
}

/// What a root search produced.
#[derive(Debug, Clone, PartialEq)]
pub struct RootReport {
    /// The root the search settled on.
    pub root: f64,
    /// The function value at `root`.
    pub residual: f64,
    /// Iterations performed.
    pub iterations: u64,
    /// Function evaluations performed.
    pub function_evals: u64,
    /// Whether the solver's own convergence test passed.
    pub converged: bool,
}

/// Central difference of a scalar function.
fn central_difference(
    f: &impl Fn(f64) -> f64,
    x: f64,
) -> f64 {
    let h = 1e-6 * x.abs().max(1.0);
    (f(x + h) - f(x - h)) / (2.0 * h)
}

/// Finds a root of `f` inside `bracket`.
///
/// `derivative` is used by [`RootSolver::Newton`] and
/// [`RootSolver::Halley`]; when it is `None` a central difference stands in,
/// and Halley's second derivative is differenced from the same source.
///
/// # Errors
/// Returns an error when the bracket is not finite and ordered, when the
/// endpoints do not bracket a sign change, or when the solver reports a
/// non-convergent failure.
pub fn find_root(
    f: impl Fn(f64) -> f64,
    derivative: Option<impl Fn(f64) -> f64>,
    bracket: (f64, f64),
    solver: RootSolver,
    options: &RootOptions,
) -> Result<RootReport, String> {
    let (lower, upper) = bracket;
    let call = |x: f64| Ok::<f64, String>(f(x));
    let derivative_of = |x: f64| {
        match &derivative {
            | Some(derivative) => derivative(x),
            | None => central_difference(&f, x),
        }
    };
    let result = match solver {
        | RootSolver::Bisection => {
            BisectionRoot::new(lower, upper)
                .with_max_iter(options.max_iter)
                .with_absolute_position_tolerance(options.tolerance)
                .solve(call)
                .map_err(|error| format!("{error:?}"))?
        },
        | RootSolver::Brent => {
            BrentRoot::new(lower, upper)
                .with_max_iter(options.max_iter)
                .with_absolute_position_tolerance(options.tolerance)
                .solve(call)
                .map_err(|error| format!("{error:?}"))?
        },
        | RootSolver::Secant => {
            SecantRoot::new(lower, upper)
                .with_max_iter(options.max_iter)
                .with_absolute_position_tolerance(options.tolerance)
                .solve(call)
                .map_err(|error| format!("{error:?}"))?
        },
        | RootSolver::Toms748 => {
            Toms748Root::new(lower, upper)
                .with_max_iter(options.max_iter)
                .with_absolute_position_tolerance(options.tolerance)
                .solve(call)
                .map_err(|error| format!("{error:?}"))?
        },
        | RootSolver::Newton => {
            NewtonRoot::new(lower, upper)
                .with_max_iter(options.max_iter)
                .with_absolute_position_tolerance(options.tolerance)
                .solve(
                    |x: f64| Ok::<f64, String>(f(x)),
                    |x: f64| Ok::<f64, String>(derivative_of(x)),
                )
                .map_err(|error| format!("{error:?}"))?
        },
        | RootSolver::Halley => {
            let second = |x: f64| {
                let h = 1e-4 * x.abs().max(1.0);
                (derivative_of(x + h) - derivative_of(x - h)) / (2.0 * h)
            };
            HalleyRoot::new(lower, upper)
                .with_max_iter(options.max_iter)
                .with_absolute_position_tolerance(options.tolerance)
                .solve(
                    |x: f64| Ok::<f64, String>(f(x)),
                    |x: f64| Ok::<f64, String>(derivative_of(x)),
                    |x: f64| Ok::<f64, String>(second(x)),
                )
                .map_err(|error| format!("{error:?}"))?
        },
    };
    Ok(RootReport {
        root: result.root(),
        residual: result.value(),
        iterations: result.iterations(),
        function_evals: result.function_evals(),
        converged: result.converged(),
    })
}

/// Rosenbrock's banana function.

/// Rosenbrock's banana function `sum 100 (x[i+1] - x[i]^2)^2 + (1 - x[i])^2`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Rosenbrock;

impl Rosenbrock {
    /// The problem over `dimension` variables.
    #[must_use]
    pub fn problem(dimension: usize) -> FnProblem {
        FnProblem::new(
            |x| {
                let mut total = 0.0;
                for pair in x.windows(2) {
                    let (a, b) = (pair[0], pair[1]);
                    total += 100.0 * (b - a * a).powi(2) + (1.0 - a).powi(2);
                }
                total
            },
            Some(|x: &[f64]| {
                let n = x.len();
                let mut g = vec![0.0; n];
                for i in 0..n.saturating_sub(1) {
                    let (a, b) = (x[i], x[i + 1]);
                    g[i] += -400.0 * a * (b - a * a) - 2.0 * (1.0 - a);
                    g[i + 1] += 200.0 * (b - a * a);
                }
                g
            }),
            dimension.max(2),
            None,
        )
    }

    /// The function value at `x`.
    #[must_use]
    pub fn cost(x: &[f64]) -> f64 {
        let mut total = 0.0;
        for pair in x.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            total += 100.0 * (b - a * a).powi(2) + (1.0 - a).powi(2);
        }
        total
    }
}

/// The sphere function `sum x[i]^2`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Sphere;

impl Sphere {
    /// The problem over `dimension` variables.
    #[must_use]
    pub fn problem(dimension: usize) -> FnProblem {
        FnProblem::new(
            |x| x.iter().map(|value| value * value).sum(),
            Some(|x: &[f64]| x.iter().map(|value| 2.0 * value).collect()),
            dimension,
            None,
        )
    }

    /// The function value at `x`.
    #[must_use]
    pub fn cost(x: &[f64]) -> f64 {
        x.iter().map(|value| value * value).sum()
    }
}

/// Rastrigin's function `10 n + sum x[i]^2 - 10 cos(2 pi x[i])`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rastrigin {
    /// Number of variables the function is defined over.
    pub dimension: usize,
}

impl Default for Rastrigin {
    fn default() -> Self {
        Self { dimension: 2 }
    }
}

impl Rastrigin {
    /// The problem over this instance's dimension.
    #[must_use]
    pub fn problem(&self) -> FnProblem {
        FnProblem::new(
            |x| {
                let n = f64::from(u32::try_from(x.len()).unwrap_or(0));
                10.0 * n
                    + x.iter()
                        .map(|value| {
                            value * value - 10.0 * (2.0 * std::f64::consts::PI * value).cos()
                        })
                        .sum::<f64>()
            },
            Some(|x: &[f64]| {
                x.iter()
                    .map(|value| {
                        2.0 * value
                            + 20.0
                                * std::f64::consts::PI
                                * (2.0 * std::f64::consts::PI * value).sin()
                    })
                    .collect()
            }),
            self.dimension,
            None,
        )
    }

    /// The function value at `x`.
    #[must_use]
    pub fn cost(x: &[f64]) -> f64 {
        let n = f64::from(u32::try_from(x.len()).unwrap_or(0));
        10.0 * n
            + x.iter()
                .map(|value| value * value - 10.0 * (2.0 * std::f64::consts::PI * value).cos())
                .sum::<f64>()
    }
}

/// Half the mean square of `residuals`, scaled when the raw sum overflows or
/// underflows.
fn half_mean_square(residuals: &[f64]) -> f64 {
    let count = residuals.len() as f64;
    let total: f64 = residuals.iter().map(|value| value * value).sum();
    if total.is_infinite() || (total > 0.0 && total < f64::MIN_POSITIVE * count) {
        let scale = residuals
            .iter()
            .map(|value| value.abs())
            .fold(0.0, f64::max);
        if scale.is_finite() && scale > 0.0 {
            let scaled: f64 = residuals.iter().map(|value| (value / scale).powi(2)).sum();
            return (scale * (scaled / (2.0 * count))) * scale;
        }
    }
    total / (2.0 * count)
}

/// Ordinary least squares with an implied intercept, measured as half the
/// mean squared residual (so its gradient is the mean residual times the
/// design column).
#[derive(Debug, Clone, PartialEq)]
pub struct LinearRegression {
    /// Predictors, one row per observation.
    pub x: Vec<Vec<f64>>,
    /// Responses, one per observation.
    pub y: Vec<f64>,
}

impl LinearRegression {
    /// Builds a least-squares problem from predictors and responses.
    ///
    /// # Errors
    /// Returns an error when the row counts differ, the design is empty, or a
    /// value is not finite.
    pub fn new(
        x: Vec<Vec<f64>>,
        y: Vec<f64>,
    ) -> Result<Self, String> {
        if x.is_empty() || x.len() != y.len() {
            return Err("predictors and responses must be non-empty and aligned".to_string());
        }
        let width = x[0].len();
        if width == 0 || x.iter().any(|row| row.len() != width) {
            return Err("predictor rows must share a non-zero width".to_string());
        }
        if x.iter().flatten().any(|value| !value.is_finite())
            || y.iter().any(|value| !value.is_finite())
        {
            return Err("least-squares data must be finite".to_string());
        }
        Ok(Self { x, y })
    }

    /// Number of coefficients, including the intercept.
    #[must_use]
    pub fn coefficient_count(&self) -> usize {
        self.x.first().map_or(0, |row| row.len() + 1)
    }

    /// Residual vector at `beta`.
    #[must_use]
    pub fn residuals(
        &self,
        beta: &[f64],
    ) -> Vec<f64> {
        self.x
            .iter()
            .zip(&self.y)
            .map(|(row, response)| {
                let mut fitted = beta.first().copied().unwrap_or(0.0);
                for (value, coefficient) in row.iter().zip(beta.iter().skip(1)) {
                    fitted += value * coefficient;
                }
                fitted - response
            })
            .collect()
    }

    /// Half the mean squared residual and its coefficient gradient.
    ///
    /// The intercept is the first coefficient, followed by one slope per
    /// column of the design. Both outputs come from the same residual pass.
    ///
    /// # Errors
    /// Returns an error when the design is empty, the row counts differ, or
    /// the coefficient count is not `columns + 1`.
    pub fn value_and_gradient(
        &self,
        beta: &[f64],
    ) -> Result<(f64, Vec<f64>), String> {
        if self.x.is_empty() || self.x.len() != self.y.len() {
            return Err("invalid regression data shape".to_string());
        }
        if beta.len() != self.coefficient_count() {
            return Err("expected the intercept and one slope per feature column".to_string());
        }
        let residuals = self.residuals(beta);
        let cost = half_mean_square(&residuals);
        let count = residuals.len() as f64;
        let columns = self.x.first().map_or(0, Vec::len);
        let mut gradient = vec![0.0; columns + 1];
        gradient[0] = residuals.iter().sum();
        for (row, &residual) in self.x.iter().zip(&residuals) {
            for (entry, &predictor) in gradient.iter_mut().skip(1).zip(row) {
                *entry += predictor * residual;
            }
        }
        for entry in &mut gradient {
            *entry /= count;
        }
        // Normalise term by term when the raw sum overflows but its mean does
        // not.
        if gradient.iter().any(|value| !value.is_finite()) {
            if !gradient[0].is_finite() {
                gradient[0] = residuals.iter().map(|value| value / count).sum();
            }
            for index in 0..columns {
                if !gradient[index + 1].is_finite() {
                    gradient[index + 1] = self
                        .x
                        .iter()
                        .zip(&residuals)
                        .map(|(row, value)| row[index] * (value / count))
                        .sum();
                }
            }
        }
        Ok((cost, gradient))
    }

    /// The problem form used by [`minimize`].
    #[must_use]
    pub fn problem(&self) -> FnProblem {
        let coefficients = self.coefficient_count();
        let cost_model = self.clone();
        let gradient_model = self.clone();
        FnProblem::new(
            move |beta| {
                cost_model
                    .value_and_gradient(beta)
                    .map_or(f64::INFINITY, |(cost, _)| cost)
            },
            Some(move |beta: &[f64]| {
                gradient_model
                    .value_and_gradient(beta)
                    .map_or_else(|_| vec![0.0; coefficients], |(_, gradient)| gradient)
            }),
            coefficients,
            None,
        )
    }
}

/// The solver entry points the kernel ABI binds to.
pub struct EquationOptimizer;

impl EquationOptimizer {
    /// Minimises `problem` with steepest descent.
    ///
    /// # Errors
    /// Propagates [`minimize`]'s errors.
    pub fn solve_with_gradient_descent(
        problem: &FnProblem,
        initial: &[f64],
        config: &OptimizationConfig,
    ) -> Result<MinimizeResult, String> {
        minimize(
            problem,
            initial,
            SolverKind::GradientDescent { step: 0.1 },
            &MinimizeOptions::from(config),
        )
    }

    /// Minimises `problem` with nonlinear conjugate gradients.
    ///
    /// # Errors
    /// Propagates [`minimize`]'s errors.
    pub fn auto_solve_conjugate_gradient(
        problem: &FnProblem,
        initial: &[f64],
        config: &OptimizationConfig,
    ) -> Result<MinimizeResult, String> {
        minimize(
            problem,
            initial,
            SolverKind::NonlinearCg,
            &MinimizeOptions::from(config),
        )
    }

    /// Minimises `problem` with BFGS.
    ///
    /// # Errors
    /// Propagates [`minimize`]'s errors.
    pub fn solve_with_bfgs(
        problem: &FnProblem,
        initial: &[f64],
        config: &OptimizationConfig,
    ) -> Result<MinimizeResult, String> {
        minimize(
            problem,
            initial,
            SolverKind::Bfgs,
            &MinimizeOptions::from(config),
        )
    }

    /// Minimises `problem` with a global-best particle swarm.
    ///
    /// `bounds` is `(lower, upper)`; the swarm needs a finite box.
    ///
    /// # Errors
    /// Propagates [`minimize`]'s errors.
    pub fn solve_with_pso(
        problem: &FnProblem,
        bounds: (Vec<f64>, Vec<f64>),
        config: &OptimizationConfig,
    ) -> Result<MinimizeResult, String> {
        let start = bounds.0.clone();
        let mut options = MinimizeOptions::from(config);
        options.bounds = Some(bounds);
        minimize(problem, &start, SolverKind::GlobalBestPso, &options)
    }

    /// Runs any [`SolverKind`] on `problem`.
    ///
    /// # Errors
    /// Propagates [`minimize`]'s errors.
    pub fn auto_solve(
        problem: &FnProblem,
        initial: &[f64],
        kind: SolverKind,
        config: &OptimizationConfig,
    ) -> Result<MinimizeResult, String> {
        minimize(problem, initial, kind, &MinimizeOptions::from(config))
    }
}

/// Labels a finished run from its cost.
pub struct ResultAnalyzer;

impl ResultAnalyzer {
    /// Prints a short report for `result`.
    pub fn print_optimization_result(result: &MinimizeResult) {
        println!("Optimization Results:");
        println!("  Converged: {}", result.converged);
        println!("  Best solution: {:?}", result.x);
        println!("  Best value: {:.6}", result.cost);
        println!("  Iterations: {}", result.iterations);
        println!("  Function evaluations: {}", result.cost_evals);
    }

    /// Labels the convergence quality of `result`.
    #[must_use]
    pub fn analyze_convergence(result: &MinimizeResult) -> String {
        if result.cost < 1e-6 {
            "Excellent convergence".to_string()
        } else if result.cost < 1e-3 {
            "Good convergence".to_string()
        } else if result.cost < 1e-1 {
            "Moderate convergence".to_string()
        } else {
            "Poor convergence".to_string()
        }
    }
}
