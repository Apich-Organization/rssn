use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;

use argmin::core::State;
use ndarray::Array1;
use serde::Deserialize;
use serde::Serialize;

use crate::numerical::optimize::*;

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
                    "Invalid JSON: {}",
                    e
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

#[no_mangle]

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
