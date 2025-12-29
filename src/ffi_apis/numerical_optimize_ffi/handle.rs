use std::ptr;
use std::slice;

use argmin::core::State;
use ndarray::Array1;

use crate::numerical::optimize::EquationOptimizer;
use crate::numerical::optimize::OptimizationConfig;
use crate::numerical::optimize::ProblemType;
use crate::numerical::optimize::Rosenbrock;
use crate::numerical::optimize::Sphere;

/// FFI-compatible optimization result containing the solution and convergence information.
///
/// This structure holds the outcome of a numerical optimization procedure and is designed
/// to be safely passed across FFI boundaries.

pub struct FfiOptimizationResult {
    /// Optimal parameter vector that minimizes the objective function.
    pub best_param: Vec<f64>,
    /// Minimum cost (objective function value) achieved at the optimal parameters.
    pub best_cost: f64,
    /// Number of iterations performed during optimization.
    pub iterations: u64,
}

/// Optimizes the Rosenbrock function using gradient descent via handle-based FFI.
///
/// The Rosenbrock function is a non-convex test function commonly used to evaluate
/// optimization algorithms: f(x,y) = (a-x)² + b(y-x²)²
///
/// # Arguments
///
/// * `a` - First parameter of the Rosenbrock function (typically 1.0)
/// * `b` - Second parameter of the Rosenbrock function (typically 100.0)
/// * `init_param_ptr` - Pointer to the initial parameter array
/// * `init_param_len` - Length of the initial parameter array
/// * `max_iters` - Maximum number of optimization iterations
/// * `tolerance` - Convergence tolerance for termination
///
/// # Returns
///
/// A raw pointer to `FfiOptimizationResult` containing the optimization outcome,
/// or null if the input is invalid or optimization fails. The caller must free
/// the result using `numerical_optimize_drop_result_handle`.
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

    let problem = Rosenbrock {
        a,
        b,
    };

    let config = OptimizationConfig {
        max_iters,
        tolerance,
        problem_type:
            ProblemType::Rosenbrock,
        dimension: init_param_len,
    };

    match EquationOptimizer::solve_with_gradient_descent(
        problem,
        init_param,
        &config,
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

/// Optimizes the Rosenbrock function using BFGS quasi-Newton method via handle-based FFI.
///
/// The BFGS (Broyden-Fletcher-Goldfarb-Shanno) algorithm is a quasi-Newton method that
/// approximates the Hessian matrix for faster convergence on smooth problems.
///
/// # Arguments
///
/// * `a` - First parameter of the Rosenbrock function (typically 1.0)
/// * `b` - Second parameter of the Rosenbrock function (typically 100.0)
/// * `init_param_ptr` - Pointer to the initial parameter array
/// * `init_param_len` - Length of the initial parameter array
/// * `max_iters` - Maximum number of optimization iterations
/// * `tolerance` - Convergence tolerance for termination
///
/// # Returns
///
/// A raw pointer to `FfiOptimizationResult` containing the optimization outcome,
/// or null if the input is invalid or optimization fails. The caller must free
/// the result using `numerical_optimize_drop_result_handle`.
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

    let problem = Rosenbrock {
        a,
        b,
    };

    let config = OptimizationConfig {
        max_iters,
        tolerance,
        problem_type:
            ProblemType::Rosenbrock,
        dimension: init_param_len,
    };

    match EquationOptimizer::solve_with_bfgs(
        problem,
        init_param,
        &config,
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

/// Optimizes the sphere function using gradient descent via handle-based FFI.
///
/// The sphere function is a convex test function: f(x) = `Σx_i²`, commonly used
/// to verify that optimization algorithms can find the global minimum at the origin.
///
/// # Arguments
///
/// * `init_param_ptr` - Pointer to the initial parameter array
/// * `init_param_len` - Length of the initial parameter array
/// * `max_iters` - Maximum number of optimization iterations
/// * `tolerance` - Convergence tolerance for termination
///
/// # Returns
///
/// A raw pointer to `FfiOptimizationResult` containing the optimization outcome,
/// or null if the input is invalid or optimization fails. The caller must free
/// the result using `numerical_optimize_drop_result_handle`.
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
        problem,
        init_param,
        &config,
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

/// Retrieves the optimal cost from an optimization result handle.
///
/// # Arguments
///
/// * `handle` - Pointer to an `FfiOptimizationResult` instance
///
/// # Returns
///
/// The minimum cost achieved, or `NaN` if the handle is null.
#[no_mangle]

pub const extern "C" fn numerical_optimize_get_result_cost_handle(
    handle : *const FfiOptimizationResult
) -> f64 {

    if handle.is_null() {

        return f64::NAN;
    }

    unsafe {

        (*handle).best_cost
    }
}

/// Retrieves the iteration count from an optimization result handle.
///
/// # Arguments
///
/// * `handle` - Pointer to an `FfiOptimizationResult` instance
///
/// # Returns
///
/// The number of iterations performed, or 0 if the handle is null.
#[no_mangle]

pub const extern "C" fn numerical_optimize_get_result_iterations_handle(
    handle : *const FfiOptimizationResult
) -> u64 {

    if handle.is_null() {

        return 0;
    }

    unsafe {

        (*handle).iterations
    }
}

/// Retrieves the parameter vector length from an optimization result handle.
///
/// # Arguments
///
/// * `handle` - Pointer to an `FfiOptimizationResult` instance
///
/// # Returns
///
/// The length of the optimal parameter vector, or 0 if the handle is null.
#[no_mangle]

pub const extern "C" fn numerical_optimize_get_result_param_len_handle(
    handle : *const FfiOptimizationResult
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

/// Copies the optimal parameter vector from an optimization result handle to a buffer.
///
/// # Arguments
///
/// * `handle` - Pointer to an `FfiOptimizationResult` instance
/// * `buffer` - Pointer to a pre-allocated buffer to receive the parameter data
///
/// # Returns
///
/// `true` if the copy succeeded, `false` if either pointer is null.
///
/// # Safety
///
/// The caller must ensure `buffer` points to an allocation with sufficient capacity
/// to hold the parameter vector (obtainable via `numerical_optimize_get_result_param_len_handle`).
#[no_mangle]

pub const extern "C" fn numerical_optimize_get_result_param_handle(
    handle : *const FfiOptimizationResult,
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

/// Frees an optimization result handle previously allocated by an optimization function.
///
/// # Arguments
///
/// * `handle` - Pointer to an `FfiOptimizationResult` instance to deallocate
///
/// # Safety
///
/// The caller must ensure the handle was previously returned by one of the optimization
/// functions in this module and has not already been freed. Passing a null pointer is safe
/// but has no effect.
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
