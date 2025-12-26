//! A comprehensive optimization module for solving various types of equations
//! using multiple optimization algorithms from the argmin 0.11 library.

use argmin::core::{
    ArgminFloat, CostFunction, Error, Executor, Gradient, IterState, Operator, OptimizationResult,
    PopulationState, Solver, State,
};
use argmin::solver::gradientdescent::SteepestDescent;
use argmin::solver::linesearch::condition::ArmijoCondition;
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::particleswarm::ParticleSwarm;
use argmin::solver::quasinewton::BFGS;
use argmin_math::ArgminConj;
use argmin_math::ArgminDot;
use argmin_math::ArgminL2Norm;
use argmin_math::ArgminMul;
use argmin_math::ArgminScaledAdd;
use argmin_math::ArgminSub;
use ndarray::{Array1, Array2};
use rand_v09::rngs::StdRng as ParticleSwarmRng;
use std::f64::consts::PI;

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
    Rosenbrock,
    Sphere,
    Rastrigin,
    Ackley,
    Custom,
}

/// Configuration for optimization algorithms
#[derive(Debug, Clone)]

pub struct OptimizationConfig {
    pub max_iters: u64,
    pub tolerance: f64,
    pub problem_type: ProblemType,
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

/// Rosenbrock function optimization (classical test function)

pub struct Rosenbrock {
    pub a: f64,
    pub b: f64,
}

impl Default for Rosenbrock {
    fn default() -> Self {

        Self { a: 1.0, b: 100.0 }
    }
}

impl CostFunction for Rosenbrock {
    type Param = Array1<f64>;

    type Output = f64;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, Error> {

        if param.len() < 2 {

            return Err(Error::msg("Parameter dimension must be at least 2"));
        }

        let mut sum = 0.0;

        for i in 0..param.len() - 1 {

            let x = param[i];

            let y = param[i + 1];

            sum += (self.a - x).powi(2)
                + self.b
                    * x.mul_add(-x, y)
                        .powi(2);
        }

        Ok(sum)
    }
}

impl Gradient for Rosenbrock {
    type Param = Array1<f64>;

    type Gradient = Array1<f64>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, Error> {

        let n = param.len();

        if n < 2 {

            return Err(Error::msg("Parameter dimension must be at least 2"));
        }

        let mut grad = Array1::zeros(n);

        for i in 0..n - 1 {

            let x = param[i];

            let y = param[i + 1];

            if i == 0 {

                grad[i] = -2.0 * (self.a - x) - 4.0 * self.b * x * x.mul_add(-x, y);
            } else {

                grad[i] += 2.0 * self.b * param[i - 1].mul_add(-param[i - 1], param[i]);
            }

            grad[i + 1] = 2.0 * self.b * x.mul_add(-x, y);
        }

        Ok(grad)
    }
}

/// Sphere function optimization (convex function)

pub struct Sphere;

impl CostFunction for Sphere {
    type Param = Array1<f64>;

    type Output = f64;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, Error> {

        Ok(param
            .iter()
            .map(|&x| x * x)
            .sum())
    }
}

impl Gradient for Sphere {
    type Param = Array1<f64>;

    type Gradient = Array1<f64>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, Error> {

        Ok(param
            .iter()
            .map(|&x| 2.0 * x)
            .collect())
    }
}

/// Rastrigin function optimization (multimodal function)

pub struct Rastrigin {
    pub a: f64,
}

impl Default for Rastrigin {
    fn default() -> Self {

        Self { a: 10.0 }
    }
}

impl CostFunction for Rastrigin {
    type Param = Array1<f64>;

    type Output = f64;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, Error> {

        let n = param.len() as f64;

        let sum: f64 = param
            .iter()
            .map(|&x| x.mul_add(x, -(self.a * (2.0 * PI * x).cos())))
            .sum();

        Ok(self
            .a
            .mul_add(n, sum))
    }
}

/// Linear regression problem optimization

pub struct LinearRegression {
    pub x: Array2<f64>,
    pub y: Array1<f64>,
}

impl LinearRegression {
    pub fn new(
        x: Array2<f64>,
        y: Array1<f64>,
    ) -> Result<Self, Error> {

        if x.is_empty() || y.is_empty() || x.nrows() != y.len() {

            return Err(Error::msg("Input data dimension mismatch"));
        }

        Ok(Self { x, y })
    }
}

impl CostFunction for LinearRegression {
    type Param = Array1<f64>;

    type Output = f64;

    fn cost(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Output, Error> {

        let mut total_error = 0.0;

        let predictions = self
            .x
            .dot(&param.slice(ndarray::s![1..]));

        for (prediction, &target) in predictions
            .iter()
            .zip(self.y.iter())
        {

            total_error += (prediction + param[0] - target).powi(2);
        }

        Ok(total_error / (2.0 * self.y.len() as f64))
    }
}

impl Gradient for LinearRegression {
    type Param = Array1<f64>;

    type Gradient = Array1<f64>;

    fn gradient(
        &self,
        param: &Self::Param,
    ) -> Result<Self::Gradient, Error> {

        let m = self.y.len() as f64;

        let predictions = self
            .x
            .dot(&param.slice(ndarray::s![1..]))
            + param[0];

        let errors = predictions - &self.y;

        let mut grad = Array1::zeros(param.len());

        grad[0] = errors.sum() / m;

        let x_t = self.x.t();

        let grad_rest = x_t.dot(&errors) / m;

        grad.slice_mut(ndarray::s![1..])
            .assign(&grad_rest);

        Ok(grad)
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

// SteepestDescent<LineSearch>
type CGOptimizer<P, G, F> = SteepestDescent<CGLineSearch<P, G, F>>;

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

    pub fn auto_solve_conjugate_gradient<C>(
        cost_function: C,
        init_param: P,
        options: &OptimizationConfig,
    ) -> AutoSolveResult<C, P, G, F, Error>
    where
        C: CostFunction<Param = P, Output = F>
            + Gradient<Param = P, Gradient = G>
            + Operator<Param = P, Output = P>,
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

        let _linesearch_condition: ArmijoCondition<F> =
            ArmijoCondition::new(0.0001).expect("Failed to create Armijo condition");

        let linesearch: MThLineSearch = MoreThuenteLineSearch::new();

        let solver: SteepestDescent<MThLineSearch> = SteepestDescent::new(linesearch);

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
                    .param(initial_param.clone())
                    .inv_hessian(Array2::eye(initial_param.len()))
                    .max_iters(config.max_iters)
                    .target_cost(config.tolerance)
            })
            .run()?;

        Ok(res)
    }

    /// Solve using particle swarm optimization (for non-differentiable functions)

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

    /// Automatically select solver and solve

    pub fn auto_solve<S, I>(
        problem: P,
        solver: S,
    ) -> Result<OptimizationResult<P, S, I>, Error>
    where
        S: Solver<P, I>,
        I: State<Param = Array1<f64>>,
    {

        Executor::new(problem, solver).run()
    }
}

/// Result analysis tools

pub struct ResultAnalyzer;

impl ResultAnalyzer {
    pub fn print_optimization_result<S: State<Param = Array1<f64>, Float = f64>>(state: &S) {

        println!("Optimization Results:");

        println!("  Converged: {}", state.get_best_cost() < 1e-4);

        if let Some(best_param) = state.get_best_param() {

            println!("  Best solution: {best_param:?}");
        } else {

            println!("  Best solution: Not available");
        }

        println!("  Best value: {:.6}", state.get_best_cost());

        println!("  Iterations: {}", state.get_iter());

        let func_counts = state.get_func_counts();

        println!(
            "  Function evaluations: {}",
            func_counts
                .get("cost")
                .unwrap_or(&0)
        );

        if let Some(grad_counts) = func_counts.get("gradient") {

            if *grad_counts > 0 {

                println!("  Gradient evaluations: {grad_counts}");
            }
        }
    }

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
