//! A comprehensive optimization module for solving various types of equations
//! using multiple optimization algorithms from the argmin 0.11 library.
//! This module is temporarily disabled for now because of conflicts between
//! the 'num-bigint' (v0.4.x) and 'statrs' (v0.18.0) dependencies,
//! which require 'rand v0.8', and 'argmin' (v0.11.0), which requires 'rand v0.9'.
//! This module will be re-enabled once the math libraries upgrade their 'rand' dependency.

use argmin::core::{CostFunction, Error, Gradient};
/*
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::gradientdescent::SteepestDescent;
use argmin::solver::conjugategradient::ConjugateGradient;
use argmin::solver::quasinewton::BFGS;
use argmin::solver::particleswarm::ParticleSwarm;
use argmin::core::Executor;
use argmin::core::State;
*/
use std::f64::consts::PI;

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
    ///
    /// Default values are: `max_iters = 1000`, `tolerance = 1e-6`,
    /// `problem_type = Rosenbrock`, `dimension = 2`.
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
    /// Provides default parameters for the Rosenbrock function.
    ///
    /// Default values are: `a = 1.0`, `b = 100.0`.
    fn default() -> Self {
        Self { a: 1.0, b: 100.0 }
    }
}

impl CostFunction for Rosenbrock {
    type Param = Vec<f64>;
    type Output = f64;

    /// Computes the cost (function value) of the Rosenbrock function.
    ///
    /// The Rosenbrock function is a non-convex function used as a performance test problem
    /// for optimization algorithms. It is defined as `f(x,y) = (a-x)² + b(y-x²)²`.
    ///
    /// # Arguments
    /// * `param` - The input parameters `[x, y, ...]`.
    ///
    /// # Returns
    /// A `Result` containing the function value, or an `Error` if parameter dimension is too small.
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
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    /// Computes the gradient of the Rosenbrock function.
    ///
    /// # Arguments
    /// * `param` - The input parameters `[x, y, ...]`.
    ///
    /// # Returns
    /// A `Result` containing the gradient vector, or an `Error` if parameter dimension is too small.
    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, Error> {
        let n = param.len();
        if n < 2 {
            return Err(Error::msg("Parameter dimension must be at least 2"));
        }

        let mut grad = vec![0.0; n];

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
    type Param = Vec<f64>;
    type Output = f64;

    /// Computes the cost (function value) of the Sphere function.
    ///
    /// The Sphere function is a simple convex function often used to test optimization algorithms.
    /// It is defined as `f(x) = Σ x_i²`.
    ///
    /// # Arguments
    /// * `param` - The input parameters `[x_1, x_2, ...]`.
    ///
    /// # Returns
    /// A `Result` containing the function value.
    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        Ok(param.iter().map(|&x| x * x).sum())
    }
}

impl Gradient for Sphere {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    /// Computes the gradient of the Sphere function.
    ///
    /// # Arguments
    /// * `param` - The input parameters `[x_1, x_2, ...]`.
    ///
    /// # Returns
    /// A `Result` containing the gradient vector.
    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, Error> {
        Ok(param.iter().map(|&x| 2.0 * x).collect())
    }
}

/// Rastrigin function optimization (multimodal function)
pub struct Rastrigin {
    pub a: f64,
}

impl Default for Rastrigin {
    /// Provides default parameters for the Rastrigin function.
    ///
    /// Default value for `a` is `10.0`.
    fn default() -> Self {
        Self { a: 10.0 }
    }
}

impl CostFunction for Rastrigin {
    type Param = Vec<f64>;
    type Output = f64;

    /// Computes the cost (function value) of the Rastrigin function.
    ///
    /// The Rastrigin function is a non-convex, multimodal function used as a performance
    /// test problem for optimization algorithms. It is defined as `f(x) = A*n + Σ (x_i² - A*cos(2πx_i))`.
    ///
    /// # Arguments
    /// * `param` - The input parameters `[x_1, x_2, ...]`.
    ///
    /// # Returns
    /// A `Result` containing the function value.
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
    pub x: Vec<Vec<f64>>, // Feature matrix
    pub y: Vec<f64>,      // Target values
}

impl LinearRegression {
    /// Creates a new `LinearRegression` problem instance.
    ///
    /// # Arguments
    /// * `x` - The feature matrix, where each inner `Vec<f64>` is a data point's features.
    /// * `y` - The target values corresponding to each data point.
    ///
    /// # Returns
    /// A `Result` containing the `LinearRegression` instance or an `Error` if data dimensions mismatch.
    pub fn new(x: Vec<Vec<f64>>, y: Vec<f64>) -> Result<Self, Error> {
        if x.is_empty() || y.is_empty() || x.len() != y.len() {
            return Err(Error::msg("Input data dimension mismatch"));
        }
        Ok(Self { x, y })
    }
}

impl CostFunction for LinearRegression {
    type Param = Vec<f64>; // Parameters: [intercept, coefficients...]
    type Output = f64;

    /// Computes the cost (sum of squared errors) for the linear regression problem.
    ///
    /// The cost function is typically `(1 / 2m) * Σ (h(x_i) - y_i)²`, where `h(x_i)` is the prediction.
    ///
    /// # Arguments
    /// * `param` - The current regression parameters `[intercept, slope1, slope2, ...]`.
    ///
    /// # Returns
    /// A `Result` containing the computed cost.
    fn cost(&self, param: &Self::Param) -> Result<Self::Output, Error> {
        let mut total_error = 0.0;

        for (features, &target) in self.x.iter().zip(self.y.iter()) {
            let mut prediction = param[0]; // Intercept
            for (i, &feature) in features.iter().enumerate() {
                prediction += param[i + 1] * feature;
            }

            total_error += (prediction - target).powi(2);
        }

        Ok(total_error / (2.0 * self.y.len() as f64))
    }
}

impl Gradient for LinearRegression {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    /// Computes the gradient of the cost function for the linear regression problem.
    ///
    /// The gradient indicates the direction of steepest ascent of the cost function.
    /// For linear regression, the gradient of the squared error cost function is used
    /// to update the parameters in gradient descent algorithms.
    ///
    /// # Arguments
    /// * `param` - The current regression parameters.
    ///
    /// # Returns
    /// A `Result` containing the gradient vector.
    fn gradient(&self, param: &Self::Param) -> Result<Self::Gradient, Error> {
        let m = self.y.len() as f64;
        let n = param.len();
        let mut grad = vec![0.0; n];

        for (features, &target) in self.x.iter().zip(self.y.iter()) {
            let mut prediction = param[0];
            for (i, &feature) in features.iter().enumerate() {
                prediction += param[i + 1] * feature;
            }

            let error = prediction - target;

            grad[0] += error;
            for (i, &feature) in features.iter().enumerate() {
                grad[i + 1] += error * feature;
            }
        }

        for g in grad.iter_mut() {
            *g /= m;
        }

        Ok(grad)
    }
}

/*
/// Main optimization solver
pub struct EquationOptimizer;

impl EquationOptimizer {
    /// Solve using gradient descent
    pub fn solve_with_gradient_descent<C>(
        cost_function: C,
        initial_param: Vec<f64>,
        config: &OptimizationConfig,
    ) -> Result<Result<Vec<f64>, (), f64>, Error>
    where
        C: CostFunction<Param = Vec<f64>, Output = f64> + Gradient<Param = Vec<f64>, Gradient = Vec<f64>>,
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

    /// Solve using conjugate gradient method
    pub fn solve_with_conjugate_gradient<C>(
        cost_function: C,
        initial_param: Vec<f64>,
        config: &OptimizationConfig,
    ) -> Result<Result<Vec<f64>, (), f64>, Error>
    where
        C: CostFunction<Param = Vec<f64>, Output = f64> + Gradient<Param = Vec<f64>, Gradient = Vec<f64>>,
    {
        let linesearch = MoreThuenteLineSearch::new();
        let solver = ConjugateGradient::new(linesearch);

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

    /// Solve using BFGS quasi-Newton method
    pub fn solve_with_bfgs<C>(
        cost_function: C,
        initial_param: Vec<f64>,
        config: &OptimizationConfig,
    ) -> Result<Result<Vec<f64>, (), f64>, Error>
    where
        C: CostFunction<Param = Vec<f64>, Output = f64> + Gradient<Param = Vec<f64>, Gradient = Vec<f64>>,
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
        bounds: (Vec<f64>, Vec<f64>),
        config: &OptimizationConfig,
    ) -> Result<Result<Vec<f64>, (), f64>, Error>
    where
        C: CostFunction<Param = Vec<f64>, Output = f64>,
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
    pub fn auto_solve(
        problem_type: ProblemType,
        initial_param: Option<Vec<f64>>,
        config: &OptimizationConfig,
    ) -> Result<Box<dyn State<Param = Vec<f64>, Float = f64>>, Error> {
        let dim = config.dimension;

        match problem_type {
            ProblemType::Rosenbrock => {
                let problem = Rosenbrock::default();
                let init_param = initial_param.unwrap_or_else(|| vec![-1.2, 1.0]);
                let result = Self::solve_with_bfgs(problem, init_param, config)?;
                Ok(Box::new(result.state))
            }
            ProblemType::Sphere => {
                let problem = Sphere;
                let init_param = initial_param.unwrap_or_else(|| vec![2.0; dim]);
                let result = Self::solve_with_conjugate_gradient(problem, init_param, config)?;
                Ok(Box::new(result.state))
            }
            ProblemType::Rastrigin => {
                let problem = Rastrigin::default();
                let bounds = (vec![-5.12; dim], vec![5.12; dim]);
                let result = Self::solve_with_pso(problem, bounds, config)?;
                Ok(Box::new(result.state))
            }
            _ => {
                let problem = Rosenbrock::default();
                let init_param = initial_param.unwrap_or_else(|| vec![-1.2, 1.0]);
                let result = Self::solve_with_gradient_descent(problem, init_param, config)?;
                Ok(Box::new(result.state))
            }
        }
    }
}

/// Result analysis tools
pub struct ResultAnalyzer;

impl ResultAnalyzer {
    pub fn print_optimization_result(state: &dyn State<Param = Vec<f64>, Float = f64>) {
        println!("Optimization Results:");
        println!("  Converged: {}", state.get_best_cost() < 1e-4);
        println!("  Best solution: {:?}", state.get_best_param().unwrap());
        println!("  Best value: {:.6}", state.get_best_cost());
        println!("  Iterations: {}", state.get_iter());

        let func_counts = state.get_func_counts();
        println!("  Function evaluations: {}", func_counts.0);
        if func_counts.1 > 0 {
            println!("  Gradient evaluations: {}", func_counts.1);
        }
    }

    pub fn analyze_convergence(state: &dyn State<Float = f64>) -> String {
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

    #[test]
    pub(crate) fn test_rosenbrock_optimization() {
        let config = OptimizationConfig {
            problem_type: ProblemType::Rosenbrock,
            max_iters: 1000,
            tolerance: 1e-8,
            dimension: 2,
        };

        let state = EquationOptimizer::auto_solve(
            ProblemType::Rosenbrock,
            Some(vec![-1.2, 1.0]),
            &config,
        ).unwrap();

        let best_param = state.get_best_param().unwrap();
        let best_cost = state.get_best_cost();

        // Rosenbrock function has global minimum at (1,1) with value 0
        assert!(best_cost < 1e-4);
        assert!((best_param[0] - 1.0).abs() < 0.1);
        assert!((best_param[1] - 1.0).abs() < 0.1);
    }

    #[test]
    pub(crate) fn test_linear_regression() {
        // Generate test data: y = 2 + 3x
        let x = vec![
            vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]
        ];
        let y = vec![5.0, 8.0, 11.0, 14.0, 17.0];

        let problem = LinearRegression::new(x, y).unwrap();
        let config = OptimizationConfig {
            problem_type: ProblemType::Custom,
            max_iters: 1000,
            tolerance: 1e-6,
            dimension: 2,
        };

        let result = EquationOptimizer::solve_with_gradient_descent(
            problem,
            vec![0.0, 0.0],  // Initial parameters [intercept, slope]
            &config,
        ).unwrap();

        let best_param = result.state.get_best_param().unwrap();

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

        let state = EquationOptimizer::auto_solve(
            ProblemType::Sphere,
            Some(vec![2.0, -1.5, 3.0]),
            &config,
        ).unwrap();

        let best_cost = state.get_best_cost();

        // Sphere function has minimum 0 at origin
        assert!(best_cost < 1e-6);
    }
}

/*
// Example usage in main function
pub(crate) fn main() -> Result<(), Error> {
    println!("=== Multiple Equation Types Optimization Module Example ===");

    // Example 1: Rosenbrock function optimization
    println!("\n1. Rosenbrock Function Optimization:");
    let config = OptimizationConfig {
        problem_type: ProblemType::Rosenbrock,
        max_iters: 1000,
        tolerance: 1e-8,
        dimension: 2,
    };

    let state = EquationOptimizer::auto_solve(
        ProblemType::Rosenbrock,
        Some(vec![-1.2, 1.0]),
        &config,
    )?;

    ResultAnalyzer::print_optimization_result(&*state);
    println!("  Convergence analysis: {}", ResultAnalyzer::analyze_convergence(&*state));

    // Example 2: Linear regression problem
    println!("\n2. Linear Regression Problem Optimization:");
    let x = vec![
        vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]
    ];
    let y = vec![5.1, 7.9, 10.8, 14.2, 16.9];  // y ≈ 2 + 3x

    let problem = LinearRegression::new(x, y).unwrap();
    let config = OptimizationConfig {
        problem_type: ProblemType::Custom,
        max_iters: 500,
        tolerance: 1e-6,
        dimension: 2,
    };

    let result = EquationOptimizer::solve_with_gradient_descent(
        problem,
        vec![0.0, 0.0],
        &config,
    )?;

    ResultAnalyzer::print_optimization_result(&result.state);

    // Example 3: Multimodal function optimization
    println!("\n3. Rastrigin Function Optimization (Multimodal):");
    let config = OptimizationConfig {
        problem_type: ProblemType::Rastrigin,
        max_iters: 1000,
        tolerance: 1e-4,
        dimension: 2,
    };

    let state = EquationOptimizer::auto_solve(
        ProblemType::Rastrigin,
        None,
        &config,
    )?;

    ResultAnalyzer::print_optimization_result(&*state);

    Ok(())
}
*/
*/
