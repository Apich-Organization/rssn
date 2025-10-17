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

// Type aliases for clarity (P=Parameter, G=Gradient, F=Float)
type P = Array1<f64>; // Assuming Array1<f64> from ndarray
type G = Array1<f64>;
type F = f64;
//type StateType = IterState<P, G, (), (), (), F>;
//type SolverType = ConjugateGradient<P, F>;
// Type alias for the LineSearch type
//type BacktrackingLS = BacktrackingLineSearch<P, G, ArmijoCondition<F>, F>;
// Used in tests.
#[allow(dead_code)]
type H = Array2<f64>;
// Used in tests.
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

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        if param.len() < 2 {
            return Err(Error::msg("Parameter dimension must be at least 2"));
        }

        let mut sum = 0.0;
        for i in 0..param.len() - 1 {
            let x = param[i];
            let y = param[i + 1];
            sum += (self.a - x).powi(2) + self.b * (y - x.powi(2)).powi(2);
        }
        Ok(sum)
    }
}

impl Gradient for Rosenbrock {
    type Param = Array1<f64>;
    type Gradient = Array1<f64>;

    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, Error> {
        let n = param.len();
        if n < 2 {
            return Err(Error::msg("Parameter dimension must be at least 2"));
        }

        let mut grad = Array1::zeros(n);

        for i in 0..n - 1 {
            let x = param[i];
            let y = param[i + 1];

            if i == 0 {
                grad[i] = -2.0 * (self.a - x) - 4.0 * self.b * x * (y - x.powi(2));
            } else {
                grad[i] += 2.0 * self.b * (param[i] - param[i - 1].powi(2));
            }

            grad[i + 1] = 2.0 * self.b * (y - x.powi(2));
        }

        Ok(grad)
    }
}

/// Sphere function optimization (convex function)
pub struct Sphere;

impl CostFunction for Sphere {
    type Param = Array1<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        Ok(param.iter().map(|&x| x * x).sum())
    }
}

impl Gradient for Sphere {
    type Param = Array1<f64>;
    type Gradient = Array1<f64>;

    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, Error> {
        Ok(param.iter().map(|&x| 2.0 * x).collect())
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

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        let n = param.len() as f64;
        let sum: f64 = param
            .iter()
            .map(|&x| x * x - self.a * (2.0 * PI * x).cos())
            .sum();
        Ok(self.a * n + sum)
    }
}

/// Linear regression problem optimization
pub struct LinearRegression {
    pub x: Array2<f64>,
    pub y: Array1<f64>,
}

impl LinearRegression {
    pub fn new(x: Array2<f64>, y: Array1<f64>) -> Result<Self, Error> {
        if x.is_empty() || y.is_empty() || x.nrows() != y.len() {
            return Err(Error::msg("Input data dimension mismatch"));
        }
        Ok(Self { x, y })
    }
}

impl CostFunction for LinearRegression {
    type Param = Array1<f64>;
    type Output = f64;

    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        let mut total_error = 0.0;
        let predictions = self.x.dot(&param.slice(ndarray::s![1..]));
        for (prediction, &target) in predictions.iter().zip(self.y.iter()) {
            total_error += (prediction + param[0] - target).powi(2);
        }
        Ok(total_error / (2.0 * self.y.len() as f64))
    }
}

impl Gradient for LinearRegression {
    type Param = Array1<f64>;
    type Gradient = Array1<f64>;

    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, Error> {
        let m = self.y.len() as f64;
        let predictions = self.x.dot(&param.slice(ndarray::s![1..])) + param[0];
        let errors = predictions - &self.y;
        let mut grad = Array1::zeros(param.len());
        grad[0] = errors.sum() / m;
        let x_t = self.x.t();
        let grad_rest = x_t.dot(&errors) / m;
        grad.slice_mut(ndarray::s![1..]).assign(&grad_rest);
        Ok(grad)
    }
}

/// Main optimization solver
pub struct EquationOptimizer;

impl EquationOptimizer {
    /// Solve using gradient descent
    pub fn solve_with_gradient_descent<C>(
        cost_function: C,
        initial_param: Array1<f64>,
        config: &OptimizationConfig,
    ) -> Result<
        OptimizationResult<
            C,
            SteepestDescent<MoreThuenteLineSearch<Array1<f64>, Array1<f64>, f64>>,
            IterState<Array1<f64>, Array1<f64>, (), (), (), f64>,
        >,
        Error,
    >
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
    ) -> Result<
        OptimizationResult<
            C,
            SteepestDescent<MoreThuenteLineSearch<P, G, F>>,
            IterState<P, G, (), (), (), F>,
        >,
        Error,
    >
    where
        // Problem Trait Bounds
        C: CostFunction<Param = P, Output = F>
            + Gradient<Param = P, Gradient = G>
            // FIX: E0220 - Operator has NO 'Float' associated type. Removed it.
            + Operator<Param = P, Output = P>,

        // Math Trait Bounds (F = f64 is usually ArgminFloat)
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
        // 1. Initialize ArmijoCondition (Returns Result, MUST unwrap)
        // The explicit type ensures we get the concrete struct
        let _linesearch_condition: ArmijoCondition<F> =
            ArmijoCondition::new(0.0001).expect("Failed to create Armijo condition");

        let linesearch: MThLineSearch = MoreThuenteLineSearch::new();

        let solver: SteepestDescent<MThLineSearch> = SteepestDescent::new(linesearch);

        // 4. Run the Executor
        let res = Executor::new(cost_function, solver)
            .configure(|state| {
                state
                    .param(init_param)
                    .max_iters(options.max_iters)
                    .target_cost(options.tolerance)
            })
            .run()?; // This '?' handles the final optimization Result

        // 5. Return the success result
        Ok(res)
    }

    /*
        /// Solve using conjugate gradient method
        pub fn run_conjugate_gradient<C>(
            cost_function: C,
            init_param: Array1<f64>,
            config: &OptimizationConfig,
        ) -> Result<OptimizationResult<C, ConjugateGradient<ArrayBase<OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>, f64>, IterState<ArrayBase<OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>, (), (), (), ArrayBase<OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>, f64>>, argmin_math::Error>
        where
            C: CostFunction<Param = P, Output = F>
                + Gradient<Param = P, Gradient = G>
                + Operator<Param = P, Output = P>, // FIX: E0220 - Added Float to the Operator bound for Hessian-vector product
            P: ArgminDot<P, F> + ArgminScaledAdd<P, F, P>, // Added necessary math bounds for P (Array1<f64>)
        {

        let linesearch_condition = ArmijoCondition::new(0.0001)
            .expect("Failed to create Armijo condition");

        let linesearch = BacktrackingLineSearch::new(linesearch_condition)
            .initial_step_length(1.0);

        // Conjugate gradient outer solver (takes the parameter type P)
        // let solver: ConjugateGradient<P, F> = ConjugateGradient::new(linesearch_condition);

        // Executor uses the now correctly-typed solver
        let res = Executor::new(cost_function, solver)
            .configure(|state| {
                state
                    .param(init_param.clone())
                    .max_iters(config.max_iters)
                    .target_cost(config.tolerance)
            })
            .run()?;

        Ok(res)
        }
    */

    /// Solve using BFGS quasi-Newton method
    pub fn solve_with_bfgs<C>(
        cost_function: C,
        initial_param: Array1<f64>,
        config: &OptimizationConfig,
    ) -> Result<
        OptimizationResult<
            C,
            BFGS<MoreThuenteLineSearch<Array1<f64>, Array1<f64>, f64>, f64>,
            IterState<Array1<f64>, Array1<f64>, (), Array2<f64>, (), f64>,
        >,
        Error,
    >
    where
        C: CostFunction<Param = Array1<f64>, Output = f64>
            + Gradient<Param = Array1<f64>, Gradient = Array1<f64>>,
    {
        let linesearch = MoreThuenteLineSearch::new();
        let solver = BFGS::new(linesearch);

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

    /// Solve using particle swarm optimization (for non-differentiable functions)
    pub fn solve_with_pso<C>(
        cost_function: C,
        bounds: (Array1<f64>, Array1<f64>),
        config: &OptimizationConfig,
    ) -> Result<
        OptimizationResult<
            C,
            ParticleSwarm<Array1<f64>, f64, ParticleSwarmRng>,
            PopulationState<argmin::solver::particleswarm::Particle<Array1<f64>, f64>, f64>,
        >,
        Error,
    >
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
    pub fn auto_solve<S, I>(problem: P, solver: S) -> Result<OptimizationResult<P, S, I>, Error>
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
            println!("  Best solution: {:?}", best_param);
        } else {
            println!("  Best solution: Not available");
        }
        println!("  Best value: {:.6}", state.get_best_cost());
        println!("  Iterations: {}", state.get_iter());

        let func_counts = state.get_func_counts();
        println!(
            "  Function evaluations: {}",
            func_counts.get("cost").unwrap_or(&0)
        );
        if let Some(grad_counts) = func_counts.get("gradient") {
            if *grad_counts > 0 {
                println!("  Gradient evaluations: {}", grad_counts);
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

#[cfg(test)]
mod tests {
    use super::*;
    use argmin::core::Executor;

    #[test]
    #[allow(unused_variables)]
    pub(crate) fn test_rosenbrock_optimization() {
        // ... (config and initial_guess setup)
        let config = OptimizationConfig {
            problem_type: ProblemType::Rosenbrock,
            max_iters: 1000,
            tolerance: 1e-8,
            dimension: 2,
        };
        let problem = Rosenbrock::default();
        let initial_guess = Array1::from(vec![-1.2, 1.0]);

        let linesearch = MoreThuenteLineSearch::new();
        let solver = BFGS::new(linesearch);

        let dimension = initial_guess.len();
        let initial_inverse_hessian = Array2::eye(dimension);

        // IterState::new(initial_param, initial_cost, initial_grad, initial_hessian, initial_inverse_hessian, initial_jacobian)
        let initial_state: BFGSState = IterState::new()
            .param(initial_guess)
            .inv_hessian(initial_inverse_hessian)
            .max_iters(config.max_iters)
            .target_cost(config.tolerance);

        let result = Executor::new(problem, solver)
            .configure(|state: BFGSState| initial_state)
            .run()
            .unwrap();

        let best_param = result.state.get_best_param().unwrap();
        let best_cost = result.state.get_best_cost();

        // Rosenbrock function has global minimum at (1,1) with value 0
        assert!(best_cost < 1e-4);
        assert!((best_param[0] - 1.0).abs() < 0.1);
        assert!((best_param[1] - 1.0).abs() < 0.1);
    }

    #[test]
    pub(crate) fn test_linear_regression() {
        // Generate test data: y = 2 + 3x
        let x = Array2::from_shape_vec((5, 1), vec![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
        let y = Array1::from(vec![5.0, 8.0, 11.0, 14.0, 17.0]);

        let problem = match LinearRegression::new(x, y) {
            Ok(p) => p,
            Err(e) => panic!("Failed to create LinearRegression problem: {}", e),
        };
        let config = OptimizationConfig {
            problem_type: ProblemType::Custom,
            max_iters: 1000,
            tolerance: 1e-6,
            dimension: 2,
        };

        let result = match EquationOptimizer::solve_with_gradient_descent(
            problem,
            Array1::from(vec![0.0, 0.0]), // Initial parameters [intercept, slope]
            &config,
        ) {
            Ok(r) => r,
            Err(e) => panic!("Solver failed for linear regression: {}", e),
        };

        let best_param = match result.state.get_best_param() {
            Some(p) => p,
            None => panic!("Best param should not be None after successful optimization"),
        };

        // Check if close to true parameters [2, 3]
        assert!((best_param[0] - 2.0).abs() < 0.5);
        assert!((best_param[1] - 3.0).abs() < 0.5);
    }

    #[test]
    pub(crate) fn test_sphere_function() {
        let config = OptimizationConfig {
            problem_type: ProblemType::Sphere,
            max_iters: 500,
            tolerance: 1e-8,
            dimension: 3,
        };
        let problem = Sphere;
        let initial_guess = Array1::from(vec![2.0, -1.5, 3.0]);

        // FIX: Replaced EquationOptimizer::auto_solve with the correctly-signed solve_with_gradient_descent
        let result =
            EquationOptimizer::solve_with_gradient_descent(problem, initial_guess, &config)
                .unwrap();

        let best_cost = result.state.get_best_cost();

        // Sphere function has minimum 0 at origin
        assert!(best_cost < 1e-6);
    }
}
