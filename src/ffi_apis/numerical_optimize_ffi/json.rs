use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use argmin::core::State;
use ndarray::Array1;
use serde::Deserialize;
use serde::Serialize;

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

/// Solves a numerical optimization problem using JSON serialization.
///
/// This function performs gradient descent optimization on well-known test problems
/// (Rosenbrock, Sphere) via FFI. The optimization configuration and problem parameters
/// are deserialized from a JSON string.
///
/// # Arguments
///
/// * `json_ptr` - A null-terminated C string containing JSON with:
///   - `problem_type`: Name of the optimization problem ("Rosenbrock" or "Sphere")
///   - `init_param`: Initial parameter vector for optimization
///   - `max_iters`: Maximum number of iterations allowed
///   - `tolerance`: Convergence tolerance threshold
///   - `rosenbrock_a`: Optional parameter `a` for Rosenbrock function (default: 1.0)
///   - `rosenbrock_b`: Optional parameter `b` for Rosenbrock function (default: 100.0)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `OptimizeResponse` with:
/// - `success`: Whether optimization succeeded
/// - `best_param`: Optimal parameter vector found
/// - `best_cost`: Minimum cost achieved
/// - `iterations`: Number of iterations performed
/// - `error`: Error message if optimization failed
///
/// The caller must free the returned string using `numerical_optimize_free_json`.
/// Returns null if the input pointer is invalid.
#[unsafe(no_mangle)]

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub extern "C" fn numerical_optimize_solve_json(
    json_ptr: *const c_char
) -> *mut c_char {

    if json_ptr.is_null() {

        return std::ptr::null_mut();
    }

    let c_str = unsafe {

        CStr::from_ptr(json_ptr)
    };

    let json_str = match c_str.to_str()
    {
        | Ok(s) => s,
        | Err(_) => {
            return std::ptr::null_mut()
        },
    };

    let request: OptimizeRequest =
        match serde_json::from_str(
            json_str,
        ) {
            | Ok(req) => req,
            | Err(e) => {

                let response =
                    OptimizeResponse {
                        success: false,
                        best_param:
                            None,
                        best_cost: None,
                        iterations:
                            None,
                        error: Some(
                            format!(
                    "Invalid JSON: {e}"
                ),
                        ),
                    };

                let json_resp = serde_json::to_string(&response).unwrap();

                return CString::new(
                    json_resp,
                )
                .unwrap()
                .into_raw();
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
            ProblemType::Custom, /* Placeholder, not used by logic below effectively */
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

    let json_resp =
        serde_json::to_string(
            &response,
        )
        .unwrap();

    CString::new(json_resp)
        .unwrap()
        .into_raw()
}

/// Frees a JSON string previously allocated by `numerical_optimize_solve_json`.
///
/// # Arguments
///
/// * `ptr` - Pointer to a C string to deallocate
///
/// # Safety
///
/// The caller must ensure the pointer was previously returned by
/// `numerical_optimize_solve_json` and has not already been freed.
/// Passing a null pointer is safe but has no effect.
#[unsafe(no_mangle)]

pub extern "C" fn numerical_optimize_free_json(
    ptr: *mut c_char
) {

    if !ptr.is_null() {

        unsafe {

            let _ =
                CString::from_raw(ptr);
        }
    }
}
