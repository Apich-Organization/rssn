//! Handle-based FFI API for numerical ODE solvers.

use std::ptr;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::matrix::Matrix;
use crate::numerical::ode::OdeSolverMethod;
use crate::numerical::ode::{
    self,
};
use crate::symbolic::core::Expr;

/// Solves a system of ODEs and returns the results as a Matrix handle.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_ode_solve(
    funcs: *const *const Expr,
    n_funcs: usize,
    y0: *const f64,
    n_y0: usize,
    x_start: f64,
    x_end: f64,
    num_steps: usize,
    method: i32,
) -> *mut Matrix<f64> {

    if funcs.is_null() || y0.is_null() {

        update_last_error(
            "Null pointer passed to \
             rssn_num_ode_solve"
                .to_string(),
        );

        return ptr::null_mut();
    }

    let method_enum = match method {
        | 0 => OdeSolverMethod::Euler,
        | 1 => OdeSolverMethod::Heun,
        | 2 => {
            OdeSolverMethod::RungeKutta4
        },
        | _ => {

            update_last_error(format!(
                "Invalid ODE solver \
                 method code: {method}"
            ));

            return ptr::null_mut();
        },
    };

    let mut funcs_vec =
        Vec::with_capacity(n_funcs);

    for i in 0 .. n_funcs {

        let f_ptr = *funcs.add(i);

        if f_ptr.is_null() {

            update_last_error(format!(
                "Null function \
                 pointer at index {i}"
            ));

            return ptr::null_mut();
        }

        funcs_vec
            .push((*f_ptr).clone());
    }

    let y0_slice =
        std::slice::from_raw_parts(
            y0, n_y0,
        );

    match ode::solve_ode_system(
        &funcs_vec,
        y0_slice,
        (x_start, x_end),
        num_steps,
        method_enum,
    ) {
        | Ok(results) => {

            let rows = results.len();

            let cols = if rows > 0 {

                results[0].len()
            } else {

                0
            };

            let mut flattened =
                Vec::with_capacity(
                    rows * cols,
                );

            for row in results {

                flattened.extend(row);
            }

            Box::into_raw(Box::new(
                Matrix::new(
                    rows,
                    cols,
                    flattened,
                ),
            ))
        },
        | Err(e) => {

            update_last_error(e);

            ptr::null_mut()
        },
    }
}
