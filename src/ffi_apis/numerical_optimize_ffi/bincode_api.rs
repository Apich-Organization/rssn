use argmin::core::State;
use ndarray::Array1;
use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::numerical::optimize::EquationOptimizer;
use crate::numerical::optimize::OptimizationConfig;
use crate::numerical::optimize::ProblemType;
use crate::numerical::optimize::Rosenbrock;
use crate::numerical::optimize::Sphere;

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

/// Solves a numerical optimization problem using bincode serialization.
///
/// This function performs gradient descent optimization on well-known test problems
/// (Rosenbrock, Sphere) via FFI. The optimization configuration and problem parameters
/// are deserialized from the input buffer.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `OptimizeRequest` with:
///   - `problem_type`: Name of the optimization problem ("Rosenbrock" or "Sphere")
///   - `init_param`: Initial parameter vector for optimization
///   - `max_iters`: Maximum number of iterations allowed
///   - `tolerance`: Convergence tolerance threshold
///   - `rosenbrock_a`: Optional parameter `a` for Rosenbrock function (default: 1.0)
///   - `rosenbrock_b`: Optional parameter `b` for Rosenbrock function (default: 100.0)
///
/// # Returns
///
/// A bincode-encoded buffer containing `OptimizeResponse` with:
/// - `success`: Whether optimization succeeded
/// - `best_param`: Optimal parameter vector found
/// - `best_cost`: Minimum cost achieved
/// - `iterations`: Number of iterations performed
/// - `error`: Error message if optimization failed
///
/// # Safety
///
/// This function is unsafe because it receives raw pointers through FFI.
/// The caller must ensure the input buffer contains valid bincode data.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn numerical_optimize_solve_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let request: OptimizeRequest =
        match from_bincode_buffer(
            &buffer,
        ) {
            | Some(req) => req,
            | None => {

                let response = OptimizeResponse {
                success : false,
                best_param : None,
                best_cost : None,
                iterations : None,
                error : Some("Invalid Bincode Input".to_string()),
            };

                return to_bincode_buffer(&response);
            },
        };

    let init_param = Array1::from(
        request
            .init_param
            .clone(),
    );

    let config = OptimizationConfig {
        max_iters: request.max_iters,
        tolerance: request.tolerance,
        problem_type:
            ProblemType::Custom,
        dimension: request
            .init_param
            .len(),
    };

    let response = match request
        .problem_type
        .as_str()
    {
        | "Rosenbrock" => {

            let a = request
                .rosenbrock_a
                .unwrap_or(1.0);

            let b = request
                .rosenbrock_b
                .unwrap_or(100.0);

            let problem = Rosenbrock {
                a,
                b,
            };

            match EquationOptimizer::solve_with_gradient_descent(
                problem,
                init_param,
                &config,
            ) {
                | Ok(res) => {
                    OptimizeResponse {
                        success : true,
                        best_param : Some(
                            res.state
                                .get_best_param()
                                .unwrap()
                                .to_vec(),
                        ),
                        best_cost : Some(
                            res.state
                                .get_best_cost(),
                        ),
                        iterations : Some(res.state.get_iter()),
                        error : None,
                    }
                },
                | Err(e) => {
                    OptimizeResponse {
                        success : false,
                        best_param : None,
                        best_cost : None,
                        iterations : None,
                        error : Some(e.to_string()),
                    }
                },
            }
        },
        | "Sphere" => {

            let problem = Sphere;

            match EquationOptimizer::solve_with_gradient_descent(
                problem,
                init_param,
                &config,
            ) {
                | Ok(res) => {
                    OptimizeResponse {
                        success : true,
                        best_param : Some(
                            res.state
                                .get_best_param()
                                .unwrap()
                                .to_vec(),
                        ),
                        best_cost : Some(
                            res.state
                                .get_best_cost(),
                        ),
                        iterations : Some(res.state.get_iter()),
                        error : None,
                    }
                },
                | Err(e) => {
                    OptimizeResponse {
                        success : false,
                        best_param : None,
                        best_cost : None,
                        iterations : None,
                        error : Some(e.to_string()),
                    }
                },
            }
        },
        | _ => {
            OptimizeResponse {
                success: false,
                best_param: None,
                best_cost: None,
                iterations: None,
                error: Some(format!(
                    "Unknown problem \
                     type: {}",
                    request
                        .problem_type
                )),
            }
        },
    };

    to_bincode_buffer(&response)
}
