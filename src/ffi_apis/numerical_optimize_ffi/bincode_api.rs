use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::numerical::optimize::*;
use argmin::core::State;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]

struct OptimizeRequest {
    problem_type: String,
    init_param: Vec<f64>,
    max_iters: u64,
    tolerance: f64,
    rosenbrock_a: Option<f64>,
    rosenbrock_b: Option<f64>,
}

#[derive(Serialize, Deserialize)]

struct OptimizeResponse {
    success: bool,
    best_param: Option<Vec<f64>>,
    best_cost: Option<f64>,
    iterations: Option<u64>,
    error: Option<String>,
}

#[no_mangle]

pub unsafe extern "C" fn numerical_optimize_solve_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let request: OptimizeRequest = match from_bincode_buffer(&buffer) {
        Some(req) => req,
        None => {

            let response = OptimizeResponse {
                success: false,
                best_param: None,
                best_cost: None,
                iterations: None,
                error: Some("Invalid Bincode Input".to_string()),
            };

            return to_bincode_buffer(&response);
        }
    };

    let init_param = Array1::from(
        request
            .init_param
            .clone(),
    );

    let config = OptimizationConfig {
        max_iters: request.max_iters,
        tolerance: request.tolerance,
        problem_type: ProblemType::Custom,
        dimension: request
            .init_param
            .len(),
    };

    let response = match request
        .problem_type
        .as_str()
    {
        "Rosenbrock" => {

            let a = request
                .rosenbrock_a
                .unwrap_or(1.0);

            let b = request
                .rosenbrock_b
                .unwrap_or(100.0);

            let problem = Rosenbrock { a, b };

            match EquationOptimizer::solve_with_gradient_descent(problem, init_param, &config) {
                Ok(res) => OptimizeResponse {
                    success: true,
                    best_param: Some(
                        res.state
                            .get_best_param()
                            .unwrap()
                            .to_vec(),
                    ),
                    best_cost: Some(
                        res.state
                            .get_best_cost(),
                    ),
                    iterations: Some(res.state.get_iter()),
                    error: None,
                },
                Err(e) => OptimizeResponse {
                    success: false,
                    best_param: None,
                    best_cost: None,
                    iterations: None,
                    error: Some(e.to_string()),
                },
            }
        }
        "Sphere" => {

            let problem = Sphere;

            match EquationOptimizer::solve_with_gradient_descent(problem, init_param, &config) {
                Ok(res) => OptimizeResponse {
                    success: true,
                    best_param: Some(
                        res.state
                            .get_best_param()
                            .unwrap()
                            .to_vec(),
                    ),
                    best_cost: Some(
                        res.state
                            .get_best_cost(),
                    ),
                    iterations: Some(res.state.get_iter()),
                    error: None,
                },
                Err(e) => OptimizeResponse {
                    success: false,
                    best_param: None,
                    best_cost: None,
                    iterations: None,
                    error: Some(e.to_string()),
                },
            }
        }
        _ => OptimizeResponse {
            success: false,
            best_param: None,
            best_cost: None,
            iterations: None,
            error: Some(format!("Unknown problem type: {}", request.problem_type)),
        },
    };

    to_bincode_buffer(&response)
}
