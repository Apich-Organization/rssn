//! A comprehensive optimization module for solving various types of equations
//! using multiple optimization algorithms from the argmin 0.11 library.

use std::f64::consts::PI;

use argmin::core::ArgminFloat;
use argmin::core::CostFunction;
use argmin::core::Error;
use argmin::core::Executor;
use argmin::core::Gradient;
use argmin::core::IterState;
use argmin::core::OptimizationResult;
use argmin::core::PopulationState;
use argmin::core::Solver;
use argmin::core::State;
use argmin::solver::conjugategradient::NonlinearConjugateGradient;
use argmin::solver::conjugategradient::beta::PolakRibierePlus;
use argmin::solver::gradientdescent::SteepestDescent;
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::particleswarm::ParticleSwarm;
use argmin::solver::quasinewton::BFGS;
use argmin_math::ArgminConj;
use argmin_math::ArgminDot;
use argmin_math::ArgminL2Norm;
use argmin_math::ArgminMul;
use argmin_math::ArgminScaledAdd;
use argmin_math::ArgminSub;
use ndarray::Array1;
use ndarray::Array2;
use ndarray::ArrayView1;
use ndarray::ArrayViewMut1;
use ndarray::linalg::general_mat_vec_mul;
use rand_v09::rngs::StdRng as ParticleSwarmRng;

use crate::kernels::numeric_common::scale_product;
use crate::kernels::reuse::f64_scratch;

type P = Array1<f64>;

type G = Array1<f64>;

type F = f64;

#[allow(dead_code)]
type H = Array2<f64>;

#[allow(dead_code)]
type BFGSState = IterState<P, G, (), H, (), F>;

type MThLineSearch = MoreThuenteLineSearch<P, G, F>;

/// Types of optimization problems supported
#[derive(Debug, Clone, Copy)]
pub enum ProblemType {
    /// Rosenbrock function.
    Rosenbrock,
    /// Sphere function.
    Sphere,
    /// Rastrigin function.
    Rastrigin,
    /// Ackley function.
    Ackley,
    /// Custom function provided by the user.
    Custom,
}

/// Configuration for optimization algorithms
#[derive(Debug, Clone)]
pub struct OptimizationConfig {
    /// Maximum number of iterations.
    pub max_iters: u64,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Type of problem being solved.
    pub problem_type: ProblemType,
    /// Number of dimensions in the parameter space.
    pub dimension: usize,
}

impl Default for OptimizationConfig {
    /// Provides a default configuration for optimization algorithms.
    fn default() -> Self {
        Self {
            max_iters: 1000,
            tolerance: 1e-6,
            problem_type: ProblemType::Rosenbrock,
            dimension: 2,
        }
    }
}

/// The chained Rosenbrock function, summed over adjacent coordinates.
///
/// Each pair `(x, y)` contributes `(a - x)^2 + b (y - x^2)^2`.
pub struct Rosenbrock {
    /// Target of the linear terms (usually 1.0).
    pub a: f64,
    /// Weight of the squared coupling residuals (usually 100.0).
    pub b: f64,
}

impl Default for Rosenbrock {
    fn default() -> Self {
        Self { a: 1.0, b: 100.0 }
    }
}

// Returns a mantissa for the coupling residual and a factor whose square
// restores its magnitude. Scaling before the fused difference preserves
// finite weighted terms even when the unweighted residual is out of range.
fn rescale_rosenbrock_residual(
    x: f64,
    y: f64,
    residual: f64,
) -> (f64, f64) {
    const UP: f64 = f64::from_bits(1535_u64 << 52); // 2^512
    const DOWN: f64 = f64::from_bits(511_u64 << 52); // 2^-512
    const SMALL_X: f64 = f64::from_bits(573_u64 << 52); // 2^-450
    const SMALL_Y: f64 = f64::from_bits(123_u64 << 52); // 2^-900
    if residual.is_infinite() && x.is_finite() && y.is_finite() {
        let scaled = x * DOWN;
        (scaled.mul_add(-scaled, (y * DOWN) * DOWN), UP)
    } else if residual.abs() < f64::MIN_POSITIVE
        && x.abs() < SMALL_X
        && y.abs() < SMALL_Y
        && !crate::kernels::exact_eq(x, 0.0)
    {
        let scaled = x * UP;
        (scaled.mul_add(-scaled, (y * UP) * UP), DOWN)
    } else {
        (residual, 1.0)
    }
}

fn rosenbrock_scaled_product(
    value: f64,
    factors: &[f64],
) -> f64 {
    // The caller handles zero/non-finite factors before requesting scaling.
    scale_product(value, factors, &[]).unwrap_or_else(|_| f64::INFINITY.copysign(value))
}

const ROSENBROCK_BOUNDED_SCALE: u64 = 1123_u64 << 52; // 2^100
const ROSENBROCK_COST_THRESHOLD: f64 = f64::from_bits(173_u64 << 52); // 2^-850

fn ordinary_rosenbrock_scale(value: f64) -> bool {
    const LOWER: u64 = 923_u64 << 52; // 2^-100
    let magnitude = value.to_bits() & 0x7fff_ffff_ffff_ffff;
    magnitude == 0 || (LOWER..=ROSENBROCK_BOUNDED_SCALE).contains(&magnitude)
}

impl Rosenbrock {
    fn cost_pairs<const SCALED: bool>(
        &self,
        param: &Array1<f64>,
    ) -> f64 {
        let term = |x: f64, y: f64| {
            let linear = self.a - x;
            let penalty = if SCALED {
                self.pair_penalty(x, y)
            } else {
                self.b * x.mul_add(-x, y).powi(2)
            };
            linear.mul_add(linear, penalty)
        };
        param.as_slice().map_or_else(
            || {
                (0..param.len() - 1)
                    .map(|i| term(param[i], param[i + 1]))
                    .sum()
            },
            |values| values.windows(2).map(|pair| term(pair[0], pair[1])).sum(),
        )
    }

    fn gradient_pairs<const SCALED: bool>(
        &self,
        param: &Array1<f64>,
        gradient: &mut Vec<f64>,
    ) {
        let (a, b) = (self.a, self.b);
        let weighted_residual = |x: f64, y: f64| b * x.mul_add(-x, y);
        let left = |x: f64, weighted: f64| (-2.0_f64).mul_add(a - x, -4.0 * (x * weighted));
        if !SCALED && let Some(values) = param.as_slice() {
            // Each interior output is independent. Three shifted slices
            // expose contiguous reads and a single disjoint output stream.
            let end = values.len() - 1;
            gradient.resize(values.len(), 0.0);
            gradient[0] = left(values[0], weighted_residual(values[0], values[1]));
            for (((entry, &x), &y), &z) in gradient[1..end]
                .iter_mut()
                .zip(&values[..end - 1])
                .zip(&values[1..end])
                .zip(&values[2..])
            {
                *entry = 2.0 * weighted_residual(x, y) + left(y, weighted_residual(y, z));
            }
            gradient[end] = 2.0 * weighted_residual(values[end - 1], values[end]);
            return;
        }
        // Each interior coordinate receives the right derivative of the
        // previous pair and the left derivative of the current one.
        let mut previous = 0.0;
        let mut derivative = |x, y| {
            let (left, right) = if SCALED {
                self.pair_derivatives(x, y)
            } else {
                let weighted = weighted_residual(x, y);
                (left(x, weighted), 2.0 * weighted)
            };
            let value = previous + left;
            previous = right;
            value
        };
        if let Some(values) = param.as_slice() {
            gradient.extend(values.windows(2).map(|pair| derivative(pair[0], pair[1])));
        } else {
            gradient.extend((0..param.len() - 1).map(|i| derivative(param[i], param[i + 1])));
        }
        gradient.push(previous);
    }

    #[inline]
    fn pair_penalty(
        &self,
        x: f64,
        y: f64,
    ) -> f64 {
        let residual = x.mul_add(-x, y);
        let square = residual * residual;
        if square.is_normal() || crate::kernels::exact_eq(residual, 0.0) {
            return self.b * square;
        }
        self.scaled_penalty(x, y, residual)
    }

    #[cold]
    fn scaled_penalty(
        &self,
        x: f64,
        y: f64,
        residual: f64,
    ) -> f64 {
        let (residual, scale) = rescale_rosenbrock_residual(x, y, residual);
        if crate::kernels::exact_eq(scale, 1.0)
            || !self.b.is_finite()
            || crate::kernels::exact_eq(residual, 0.0)
        {
            (self.b * residual) * residual
        } else {
            let magnitude = residual.abs();
            rosenbrock_scaled_product(self.b, &[magnitude, magnitude, scale, scale, scale, scale])
        }
    }

    #[inline]
    fn pair_derivatives(
        &self,
        x: f64,
        y: f64,
    ) -> (f64, f64) {
        let residual = x.mul_add(-x, y);
        let weighted = self.b * residual;
        let mixed = x * weighted;
        let (coupling, next) =
            if weighted.is_normal() && (mixed.is_normal() || crate::kernels::exact_eq(x, 0.0)) {
                (-4.0 * mixed, 2.0 * weighted)
            } else {
                self.scaled_derivatives(x, y, residual)
            };
        ((-2.0_f64).mul_add(self.a - x, coupling), next)
    }

    #[cold]
    fn scaled_derivatives(
        &self,
        x: f64,
        y: f64,
        residual: f64,
    ) -> (f64, f64) {
        let (residual, scale) = rescale_rosenbrock_residual(x, y, residual);
        if !residual.is_finite() || !x.is_finite() || !self.b.is_finite() {
            let weighted = self.b * residual;
            return (-4.0 * (x * weighted), 2.0 * weighted);
        }
        if crate::kernels::exact_eq(residual, 0.0) {
            return (0.0, 0.0);
        }
        let signed = if residual.is_sign_negative() {
            -self.b
        } else {
            self.b
        };
        let magnitude = residual.abs();
        let right = rosenbrock_scaled_product(signed, &[2.0, magnitude, scale, scale]);
        let left = if crate::kernels::exact_eq(x, 0.0) {
            0.0
        } else {
            let weight = if x.is_sign_negative() {
                signed
            } else {
                -signed
            };
            rosenbrock_scaled_product(weight, &[4.0, magnitude, x.abs(), scale, scale])
        };
        (left, right)
    }
}

impl CostFunction for Rosenbrock {
    type Output = f64;
    type Param = Array1<f64>;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, Error> {
        if param.len() < 2 {
            return Err(Error::msg("Rosenbrock needs at least two coordinates"));
        }
        if crate::kernels::exact_eq(self.b, 0.0) {
            return Ok(param
                .iter()
                .take(param.len() - 1)
                .map(|x| (self.a - x).powi(2))
                .sum());
        }
        // With this weight, even 2^64 subnormal square roundings
        // contribute less than 2^-911 in total. At this cost cutoff they
        // are below the precision of the ordinary result.
        if self.b.abs() <= f64::from_bits(ROSENBROCK_BOUNDED_SCALE) {
            let ordinary = self.cost_pairs::<false>(param);
            if ordinary.is_finite() && ordinary.abs() >= ROSENBROCK_COST_THRESHOLD {
                return Ok(ordinary);
            }
        }
        Ok(self.cost_pairs::<true>(param))
    }
}

impl Gradient for Rosenbrock {
    type Gradient = Array1<f64>;
    type Param = Array1<f64>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, Error> {
        let n = param.len();
        if n < 2 {
            return Err(Error::msg("Rosenbrock needs at least two coordinates"));
        }
        let mut gradient = Vec::with_capacity(n);
        if crate::kernels::exact_eq(self.b, 0.0) {
            gradient.extend(param.iter().take(n - 1).map(|x| -2.0 * (self.a - x)));
            gradient.push(0.0);
        } else {
            // Nonzero inputs in range are multiples of 2^-152. Residuals
            // are multiples of 2^-304, so the bounded nonzero products used
            // for both gradient contributions stay in the normal range.
            let ordinary = n >= 16
                && ordinary_rosenbrock_scale(self.b)
                && param.as_slice().map_or_else(
                    || (0..n).fold(true, |ok, i| ok & ordinary_rosenbrock_scale(param[i])),
                    |values| {
                        values
                            .iter()
                            .fold(true, |ok, &x| ok & ordinary_rosenbrock_scale(x))
                    },
                );
            if ordinary {
                self.gradient_pairs::<false>(param, &mut gradient);
            } else {
                self.gradient_pairs::<true>(param, &mut gradient);
            }
        }
        Ok(Array1::from_vec(gradient))
    }
}

/// Sphere function optimization (convex function)
pub struct Sphere;

impl CostFunction for Sphere {
    type Output = f64;
    type Param = Array1<f64>;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, Error> {
        Ok(param.iter().map(|&x| x * x).sum())
    }
}

impl Gradient for Sphere {
    type Gradient = Array1<f64>;
    type Param = Array1<f64>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, Error> {
        Ok(param.iter().map(|&x| 2.0 * x).collect())
    }
}

/// Rastrigin function optimization (multimodal function)
pub struct Rastrigin {
    /// Parameter a (usually 10.0).
    pub a: f64,
}

impl Default for Rastrigin {
    fn default() -> Self {
        Self { a: 10.0 }
    }
}

impl CostFunction for Rastrigin {
    type Output = f64;
    type Param = Array1<f64>;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, Error> {
        let n = param.len() as f64;

        let sum: f64 = param
            .iter()
            .map(|&x| x.mul_add(x, -(self.a * (2.0 * PI * x).cos())))
            .sum();

        Ok(self.a.mul_add(n, sum))
    }
}

/// Linear regression problem optimization
pub struct LinearRegression {
    /// Feature matrix X.
    pub x: Array2<f64>,
    /// Target vector y.
    pub y: Array1<f64>,
}

impl LinearRegression {
    /// Creates a new `LinearRegression` problem instance.
    ///
    /// # Errors
    /// Returns an error if the input data dimensions are mismatched or empty.
    #[allow(clippy::suspicious_operation_groupings)]
    pub fn new(
        x: Array2<f64>,
        y: Array1<f64>,
    ) -> Result<Self, Error> {
        if x.is_empty() || y.is_empty() || x.nrows() != y.len() {
            return Err(Error::msg(
                "Input data dimension \
                 mismatch",
            ));
        }

        Ok(Self { x, y })
    }

    /// Evaluates the half mean-squared residual and its coefficient gradient.
    ///
    /// The intercept is the first coefficient, followed by one slope per
    /// column of `x`. Both outputs use the same residual evaluation.
    ///
    /// # Errors
    /// Returns an error if the data are empty, the row and target counts
    /// differ, or the coefficient count differs from `x.ncols() + 1`.
    pub fn value_and_gradient(
        &self,
        param: &Array1<f64>,
    ) -> Result<(f64, Array1<f64>), Error> {
        self.with_residuals(param, |residuals| {
            let cost = half_mean_square(residuals.iter().copied(), residuals.len());
            (cost, self.gradient_from_residuals(residuals))
        })
    }

    fn parameters<'p>(
        &self,
        param: &'p Array1<f64>,
    ) -> Result<(f64, ArrayView1<'p, f64>), Error> {
        if self.x.is_empty() || self.x.nrows() != self.y.len() {
            return Err(Error::msg("invalid regression data shape"));
        }
        if param.len().checked_sub(1) != Some(self.x.ncols()) {
            return Err(Error::msg(
                "expected the intercept and one slope per feature column",
            ));
        }
        Ok((param[0], param.slice(ndarray::s![1..])))
    }

    fn single_column_residuals(
        &self,
        intercept: f64,
        slope: f64,
    ) -> impl Iterator<Item = f64> + Clone + '_ {
        self.x
            .column(0)
            .into_iter()
            .zip(&self.y)
            .map(move |(&x, &target)| regression_residual(x * slope, intercept, target))
    }

    fn with_row_buffer<T>(
        &self,
        evaluate: impl FnOnce(ArrayViewMut1<'_, f64>) -> T,
    ) -> T {
        const STACK_ROWS: usize = 64;
        let n = self.y.len();
        if n <= STACK_ROWS {
            let mut storage = [0.0; STACK_ROWS];
            evaluate(ArrayViewMut1::from(&mut storage[..n]))
        } else {
            let mut scratch = f64_scratch(n);
            scratch.resize(n, 0.0);
            evaluate(ArrayViewMut1::from(&mut scratch[..]))
        }
    }

    fn with_residuals<T>(
        &self,
        param: &Array1<f64>,
        evaluate: impl FnOnce(ArrayView1<'_, f64>) -> T,
    ) -> Result<T, Error> {
        let (intercept, slopes) = self.parameters(param)?;
        Ok(self.with_row_buffer(|mut residuals| {
            if slopes.len() == 1 {
                for (entry, value) in residuals
                    .iter_mut()
                    .zip(self.single_column_residuals(intercept, slopes[0]))
                {
                    *entry = value;
                }
            } else {
                general_mat_vec_mul(1.0, &self.x, &slopes, 0.0, &mut residuals);
                for (prediction, &target) in residuals.iter_mut().zip(&self.y) {
                    *prediction = regression_residual(*prediction, intercept, target);
                }
            }
            evaluate(residuals.view())
        }))
    }

    fn gradient_from_residuals(
        &self,
        residuals: ArrayView1<'_, f64>,
    ) -> Array1<f64> {
        let count = residuals.len() as f64;
        let columns = self.x.ncols();
        let mut gradient = Array1::zeros(columns + 1);
        gradient[0] = residuals.sum();
        if columns > 1
            && residuals.len() >= 64
            && let Some(matrix) = self.x.as_slice()
            && let Some(coefficients) = gradient.as_slice_mut()
        {
            // Stream contiguous rows while keeping the order of terms in each
            // column's dot product.
            let slopes = &mut coefficients[1..];
            for (row, &residual) in matrix.chunks_exact(columns).zip(&residuals) {
                for (entry, &predictor) in slopes.iter_mut().zip(row) {
                    *entry += predictor * residual;
                }
            }
        } else {
            general_mat_vec_mul(
                1.0,
                &self.x.t(),
                &residuals,
                0.0,
                &mut gradient.slice_mut(ndarray::s![1..]),
            );
        }
        gradient /= count;
        if gradient.iter().any(|v| !v.is_finite()) {
            // Normalize the terms before summing when the unscaled sum
            // overflows, even though its mean can be finite.
            if !gradient[0].is_finite() {
                gradient[0] = residuals.iter().map(|r| r / count).sum();
            }
            for (entry, column) in gradient.iter_mut().skip(1).zip(self.x.columns()) {
                if !entry.is_finite() {
                    *entry = column
                        .iter()
                        .zip(&residuals)
                        .map(|(x, r)| x * (r / count))
                        .sum();
                }
            }
        }
        gradient
    }
}

fn regression_residual(
    prediction: f64,
    intercept: f64,
    target: f64,
) -> f64 {
    prediction + intercept - target
}

fn half_mean_square(
    residuals: impl Iterator<Item = f64> + Clone,
    count: usize,
) -> f64 {
    let count = count as f64;
    let total: f64 = residuals.clone().map(|r| r * r).sum();
    if total.is_infinite() || (total > 0.0 && total < f64::MIN_POSITIVE * count) {
        let scale = residuals.clone().map(f64::abs).fold(0.0, f64::max);
        if scale.is_finite() && scale > 0.0 {
            let scaled: f64 = residuals.map(|r| (r / scale).powi(2)).sum();
            return (scale * (scaled / (2.0 * count))) * scale;
        }
    }
    total / (2.0 * count)
}

impl CostFunction for LinearRegression {
    type Output = f64;
    type Param = Array1<f64>;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, Error> {
        let (intercept, slopes) = self.parameters(param)?;
        if slopes.len() == 1 {
            let residuals = self.single_column_residuals(intercept, slopes[0]);
            return Ok(half_mean_square(residuals, self.y.len()));
        }
        Ok(self.with_row_buffer(|mut predictions| {
            general_mat_vec_mul(1.0, &self.x, &slopes, 0.0, &mut predictions);
            let residuals = predictions
                .iter()
                .zip(&self.y)
                .map(|(&prediction, &target)| regression_residual(prediction, intercept, target));
            half_mean_square(residuals, self.y.len())
        }))
    }
}

impl Gradient for LinearRegression {
    type Gradient = Array1<f64>;
    type Param = Array1<f64>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, Error> {
        self.with_residuals(param, |residuals| self.gradient_from_residuals(residuals))
    }
}

/// Main optimization solver
pub struct EquationOptimizer;

// MoreThuenteLineSearch<Vector, ScaledVector, Scalar>
type LineSearch = MoreThuenteLineSearch<Array1<f64>, Array1<f64>, f64>;

// IterState<Vector, ScaledVector, Alpha, Direction, LineSearchState, Scalar>
type GDIterState = IterState<Array1<f64>, Array1<f64>, (), (), (), f64>;

// SteepestDescent<LineSearch>
type GDOptimizer = SteepestDescent<LineSearch>;

// OptimizationResult<CostFunction, Optimizer, IterState>
type GDResult<C> = OptimizationResult<C, GDOptimizer, GDIterState>;

// Result<OptimizationResult<...>, Error>
type SolveResult<C, Error> = Result<GDResult<C>, Error>;

type CGLineSearch<P, G, F> = MoreThuenteLineSearch<P, G, F>;

type CGIterState<P, G, F> = IterState<P, G, (), (), (), F>;

// Genuine nonlinear conjugate gradient with a nonnegative Polak–Ribière update.
type CGOptimizer<P, G, F> =
    NonlinearConjugateGradient<P, CGLineSearch<P, G, F>, PolakRibierePlus, F>;

// OptimizationResult<CostFunction, Optimizer, IterState>
type CGResult<C, P, G, F> = OptimizationResult<C, CGOptimizer<P, G, F>, CGIterState<P, G, F>>;

type AutoSolveResult<C, P, G, F, Error> = Result<CGResult<C, P, G, F>, Error>;

type Vector = Array1<f64>;

type Scalar = f64;

type HessianApprox = Array2<f64>;

// 1. Line Search Type
// MoreThuenteLineSearch<Vector, ScaledVector, Scalar>
type BFGSLineSearch = MoreThuenteLineSearch<Vector, Vector, Scalar>;

// IterState<Vector, ScaledVector, Alpha, Direction (HessianApprox), LineSearchState, Scalar>
type BFGSIterState = IterState<Vector, Vector, (), HessianApprox, (), Scalar>;

// 3. Optimizer Type
// BFGS<LineSearch, Scalar>
type BFGSOptimizer = BFGS<BFGSLineSearch, Scalar>;

// OptimizationResult<CostFunction, Optimizer, IterState>
type BFGSResult<C> = OptimizationResult<C, BFGSOptimizer, BFGSIterState>;

type BFGSSolveResult<C, E> = Result<BFGSResult<C>, E>;

type PsoVector = Array1<f64>;

type PsoScalar = f64;

type PsoRng = ParticleSwarmRng;

// Particle<Vector, Scalar>
type PsoParticle = argmin::solver::particleswarm::Particle<PsoVector, PsoScalar>;

// PopulationState<ParticleType, Scalar>
type PsoState = PopulationState<PsoParticle, PsoScalar>;

// 3. Optimizer Type
// ParticleSwarm<Vector, Scalar, Rng>
type PsoOptimizer = ParticleSwarm<PsoVector, PsoScalar, PsoRng>;

// OptimizationResult<CostFunction, Optimizer, IterState>
type PsoResult<C> = OptimizationResult<C, PsoOptimizer, PsoState>;

type PsoSolveResult<C, E> = Result<PsoResult<C>, E>;

impl EquationOptimizer {
    /// Solve using gradient descent
    ///
    /// # Errors
    /// Returns an error if the optimization process fails to initialize or execute.
    pub fn solve_with_gradient_descent<C>(
        cost_function: C,
        initial_param: Array1<f64>,
        config: &OptimizationConfig,
    ) -> SolveResult<C, Error>
    where
        C: CostFunction<Param = Array1<f64>, Output = f64>
            + Gradient<Param = Array1<f64>, Gradient = Array1<f64>>,
    {
        let linesearch = MoreThuenteLineSearch::new();

        let solver = SteepestDescent::new(linesearch);

        let res = Executor::new(cost_function, solver)
            .configure(|state| {
                state
                    .param(initial_param)
                    .max_iters(config.max_iters)
                    .target_cost(config.tolerance)
            })
            .run()?;

        Ok(res)
    }

    /// Automatically configures and solves a problem using the Conjugate Gradient method.
    ///
    /// # Errors
    /// Returns an error if the optimization process fails.
    pub fn auto_solve_conjugate_gradient<C>(
        cost_function: C,
        init_param: P,
        options: &OptimizationConfig,
    ) -> AutoSolveResult<C, P, G, F, Error>
    where
        C: CostFunction<Param = P, Output = F> + Gradient<Param = P, Gradient = G>,
        P: ArgminDot<P, F>
            + ArgminScaledAdd<P, F, P>
            + ArgminSub<P, P>
            + ArgminConj
            + ArgminMul<F, P>
            + Clone
            + std::fmt::Debug,
        G: ArgminL2Norm<F>,
        F: ArgminFloat + ArgminL2Norm<F>,
    {
        // Nonlinear CG only produces descent directions when the line search
        // satisfies the strong Wolfe conditions with c2 < 1/2; the default
        // c2 = 0.9 (meant for quasi-Newton methods) let Polak-Ribiere+ step
        // into an ascent direction on curved valleys such as Rosenbrock's.
        // Powell's orthogonality restart resets to steepest descent when
        // successive gradients stop being nearly orthogonal.
        let constant = |value: f64| {
            <F as num_traits::FromPrimitive>::from_f64(value)
                .ok_or_else(|| Error::msg("conjugate gradient: constant not representable"))
        };
        let linesearch: MThLineSearch =
            MoreThuenteLineSearch::new().with_c(constant(1e-4)?, constant(0.1)?)?;

        let solver = NonlinearConjugateGradient::new(linesearch, PolakRibierePlus::new())
            .restart_orthogonality(constant(0.1)?);

        let res = Executor::new(cost_function, solver)
            .configure(|state| {
                state
                    .param(init_param)
                    .max_iters(options.max_iters)
                    .target_cost(options.tolerance)
            })
            .run()?;

        Ok(res)
    }

    /// Solve using BFGS quasi-Newton method
    ///
    /// # Errors
    /// Returns an error if the BFGS optimization fails.
    pub fn solve_with_bfgs<C>(
        cost_function: C,
        initial_param: Array1<f64>,
        config: &OptimizationConfig,
    ) -> BFGSSolveResult<C, Error>
    where
        C: CostFunction<Param = Array1<f64>, Output = f64>
            + Gradient<Param = Array1<f64>, Gradient = Array1<f64>>,
    {
        let linesearch = MoreThuenteLineSearch::new();

        let solver = BFGS::new(linesearch);

        let res = Executor::new(cost_function, solver)
            .configure(|state| {
                state
                    .inv_hessian(Array2::eye(initial_param.len()))
                    .param(initial_param)
                    .max_iters(config.max_iters)
                    .target_cost(config.tolerance)
            })
            .run()?;

        Ok(res)
    }

    /// Solve using particle swarm optimization (for non-differentiable functions)
    ///
    /// # Errors
    /// Returns an error if the PSO optimization fails.
    pub fn solve_with_pso<C>(
        cost_function: C,
        bounds: (Array1<f64>, Array1<f64>),
        config: &OptimizationConfig,
    ) -> PsoSolveResult<C, Error>
    where
        C: CostFunction<Param = Array1<f64>, Output = f64>,
    {
        let solver = ParticleSwarm::new(bounds, 40);

        let res = Executor::new(cost_function, solver)
            .configure(|state| {
                state
                    .max_iters(config.max_iters)
                    .target_cost(config.tolerance)
            })
            .run()?;

        Ok(res)
    }

    /// Runs `solver` on `problem` with an `argmin` executor.
    ///
    /// `problem` is any argmin problem type (the cost function / gradient
    /// provider), not the parameter vector. `configure` receives the solver's
    /// fresh state and returns it set up, e.g. `|s| s.param(x0).max_iters(100)`
    /// for iterative solvers or `|s| s.max_iters(200)` for population solvers;
    /// it is required because argmin's default iteration limit is unbounded.
    ///
    /// # Errors
    /// Returns an error if the optimization process fails.
    pub fn auto_solve<C, S, I, Cfg>(
        problem: C,
        solver: S,
        configure: Cfg,
    ) -> Result<OptimizationResult<C, S, I>, Error>
    where
        S: Solver<C, I>,
        I: State,
        Cfg: FnOnce(I) -> I,
    {
        Executor::new(problem, solver).configure(configure).run()
    }
}

/// Result analysis tools
pub struct ResultAnalyzer;

impl ResultAnalyzer {
    /// Prints the optimization results to the console.
    pub fn print_optimization_result<S: State<Param = Array1<f64>, Float = f64>>(state: &S) {
        println!("Optimization Results:");

        println!("  Converged: {}", state.get_best_cost() < 1e-4);

        if let Some(best_param) = state.get_best_param() {
            println!(
                "  Best solution: \
                 {best_param:?}"
            );
        } else {
            println!(
                "  Best solution: Not \
                 available"
            );
        }

        println!("  Best value: {:.6}", state.get_best_cost());

        println!("  Iterations: {}", state.get_iter());

        let func_counts = state.get_func_counts();

        println!(
            "  Function evaluations: \
             {}",
            func_counts.get("cost").unwrap_or(&0)
        );

        if let Some(grad_counts) = func_counts.get("gradient")
            && *grad_counts > 0
        {
            println!(
                "  Gradient \
                     evaluations: \
                     {grad_counts}"
            );
        }
    }

    /// Analyzes the convergence state and returns a summary string.
    pub fn analyze_convergence<S: State<Float = f64>>(state: &S) -> String {
        let cost = state.get_best_cost();

        if cost < 1e-6 {
            "Excellent convergence".to_string()
        } else if cost < 1e-3 {
            "Good convergence".to_string()
        } else if cost < 1e-1 {
            "Moderate convergence".to_string()
        } else {
            "Poor convergence".to_string()
        }
    }
}
