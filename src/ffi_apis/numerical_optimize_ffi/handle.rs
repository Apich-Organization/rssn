use std::ptr;
use std::slice;

use argmin::core::State;
use ndarray::Array1;

use crate::numerical::optimize::*;

pub struct FfiOptimizationResult {
    pub best_param: Vec<f64>,
    pub best_cost: f64,
    pub iterations: u64,
}

#[no_mangle]

pub extern "C" fn numerical_optimize_rosenbrock_gd_handle(
    a: f64,
    b: f64,
    init_param_ptr: *const f64,
    init_param_len: usize,
    max_iters: u64,
    tolerance: f64,
) -> *mut FfiOptimizationResult {

    if init_param_ptr.is_null() {

        return ptr::null_mut();
    }

    let init_param_slice = unsafe {

        slice::from_raw_parts(
            init_param_ptr,
            init_param_len,
        )
    };

    let init_param = Array1::from(
        init_param_slice.to_vec(),
    );

    let problem = Rosenbrock { a, b };

    let config = OptimizationConfig {
        max_iters,
        tolerance,
        problem_type:
            ProblemType::Rosenbrock,
        dimension: init_param_len,
    };

    match EquationOptimizer::solve_with_gradient_descent(
        problem, init_param, &config,
    ) {
        | Ok(res) => {

            // Use explicit trait methods if needed, but normally method syntax works
            let best_param = res
                .state
                .get_best_param()
                .unwrap()
                .to_vec();

            let best_cost = res
                .state
                .get_best_cost();

            let iterations = res.state.get_iter();

            let result = FfiOptimizationResult {
                best_param,
                best_cost,
                iterations,
            };

            Box::into_raw(Box::new(result))
        },
        | Err(_) => ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn numerical_optimize_rosenbrock_bfgs_handle(
    a: f64,
    b: f64,
    init_param_ptr: *const f64,
    init_param_len: usize,
    max_iters: u64,
    tolerance: f64,
) -> *mut FfiOptimizationResult {

    if init_param_ptr.is_null() {

        return ptr::null_mut();
    }

    let init_param_slice = unsafe {

        slice::from_raw_parts(
            init_param_ptr,
            init_param_len,
        )
    };

    let init_param = Array1::from(
        init_param_slice.to_vec(),
    );

    let problem = Rosenbrock { a, b };

    let config = OptimizationConfig {
        max_iters,
        tolerance,
        problem_type:
            ProblemType::Rosenbrock,
        dimension: init_param_len,
    };

    match EquationOptimizer::solve_with_bfgs(
        problem, init_param, &config,
    ) {
        | Ok(res) => {

            let best_param = res
                .state
                .get_best_param()
                .unwrap()
                .to_vec();

            let best_cost = res
                .state
                .get_best_cost();

            let iterations = res.state.get_iter();

            let result = FfiOptimizationResult {
                best_param,
                best_cost,
                iterations,
            };

            Box::into_raw(Box::new(result))
        },
        | Err(_) => ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn numerical_optimize_sphere_gd_handle(
    init_param_ptr: *const f64,
    init_param_len: usize,
    max_iters: u64,
    tolerance: f64,
) -> *mut FfiOptimizationResult {

    if init_param_ptr.is_null() {

        return ptr::null_mut();
    }

    let init_param_slice = unsafe {

        slice::from_raw_parts(
            init_param_ptr,
            init_param_len,
        )
    };

    let init_param = Array1::from(
        init_param_slice.to_vec(),
    );

    let problem = Sphere;

    let config = OptimizationConfig {
        max_iters,
        tolerance,
        problem_type:
            ProblemType::Sphere,
        dimension: init_param_len,
    };

    match EquationOptimizer::solve_with_gradient_descent(
        problem, init_param, &config,
    ) {
        | Ok(res) => {

            let best_param = res
                .state
                .get_best_param()
                .unwrap()
                .to_vec();

            let best_cost = res
                .state
                .get_best_cost();

            let iterations = res.state.get_iter();

            let result = FfiOptimizationResult {
                best_param,
                best_cost,
                iterations,
            };

            Box::into_raw(Box::new(result))
        },
        | Err(_) => ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn numerical_optimize_get_result_cost_handle(
    handle: *const FfiOptimizationResult
) -> f64 {

    if handle.is_null() {

        return f64::NAN;
    }

    unsafe {

        (*handle).best_cost
    }
}

#[no_mangle]

pub extern "C" fn numerical_optimize_get_result_iterations_handle(
    handle: *const FfiOptimizationResult
) -> u64 {

    if handle.is_null() {

        return 0;
    }

    unsafe {

        (*handle).iterations
    }
}

#[no_mangle]

pub extern "C" fn numerical_optimize_get_result_param_len_handle(
    handle: *const FfiOptimizationResult
) -> usize {

    if handle.is_null() {

        return 0;
    }

    unsafe {

        (*handle)
            .best_param
            .len()
    }
}

#[no_mangle]

pub extern "C" fn numerical_optimize_get_result_param_handle(
    handle: *const FfiOptimizationResult,
    buffer: *mut f64,
) -> bool {

    if handle.is_null()
        || buffer.is_null()
    {

        return false;
    }

    unsafe {

        let res = &*handle;

        ptr::copy_nonoverlapping(
            res.best_param
                .as_ptr(),
            buffer,
            res.best_param.len(),
        );
    }

    true
}

#[no_mangle]

pub extern "C" fn numerical_optimize_drop_result_handle(
    handle: *mut FfiOptimizationResult
) {

    if !handle.is_null() {

        unsafe {

            let _ =
                Box::from_raw(handle);
        }
    }
}
