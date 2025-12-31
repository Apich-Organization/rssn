use std::collections::HashMap;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::os::raw::c_int;

use crate::symbolic::core::Expr;
use crate::symbolic::proof;

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

unsafe fn parse_c_str_array(
    arr: *const *const c_char,
    len: usize,
) -> Option<Vec<String>> {

    unsafe {

        if arr.is_null() && len > 0 {

            return None;
        }

        let mut vars =
            Vec::with_capacity(len);

        for i in 0 .. len {

            let ptr = *arr.add(i);

            if ptr.is_null() {

                return None;
            }

            let c_str =
                CStr::from_ptr(ptr);

            match c_str.to_str() {
                | Ok(s) => {

                    vars.push(
                        s.to_string(),
                    );
                },
                | Err(_) => {

                    return None;
                },
            }
        }

        Some(vars)
    }
}

/// Verifies an equation solution (Handle)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_verify_equation_solution_handle(
    equations_ptr: *const *const Expr,
    equations_len: c_int,
    sol_vars_ptr: *const *const c_char,
    sol_exprs_ptr: *const *const Expr,
    sol_len: c_int,
    free_vars_ptr: *const *const c_char,
    free_vars_len: c_int,
) -> bool {

    unsafe {

        if equations_ptr.is_null()
            || sol_vars_ptr.is_null()
            || sol_exprs_ptr.is_null()
        {

            return false;
        }

        let mut equations =
            Vec::with_capacity(
                i64::from(
                    equations_len,
                )
                .try_into()
                .unwrap_or(0),
            );

        for i in 0 .. equations_len {

            let ptr = *equations_ptr
                .add(
                    i64::from(i)
                        .try_into()
                        .unwrap_or(0),
                );

            if ptr.is_null() {

                return false;
            }

            equations
                .push((*ptr).clone());
        }

        let sol_vars =
            match parse_c_str_array(
                sol_vars_ptr,
                i64::from(sol_len)
                    .try_into()
                    .unwrap_or(0),
            ) {
                | Some(v) => v,
                | None => return false,
            };

        let mut solution =
            HashMap::new();

        for i in 0 .. sol_len {

            let expr_ptr =
                *sol_exprs_ptr.add(
                    i64::from(i)
                        .try_into()
                        .unwrap_or(0),
                );

            if expr_ptr.is_null() {

                return false;
            }

            solution.insert(
                sol_vars[i64::from(i)
                    .try_into()
                    .unwrap_or(0)]
                .clone(),
                (*expr_ptr).clone(),
            );
        }

        let free_vars_strings =
            match parse_c_str_array(
                free_vars_ptr,
                i64::from(
                    free_vars_len,
                )
                .try_into()
                .unwrap_or(0),
            ) {
                | Some(v) => v,
                | None => return false,
            };

        let free_vars: Vec<&str> =
        free_vars_strings
            .iter()
            .map(std::string::String::as_str)
            .collect();

        proof::verify_equation_solution(
            &equations,
            &solution,
            &free_vars,
        )
    }
}

/// Verifies an indefinite integral (Handle)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_verify_indefinite_integral_handle(
    integrand_ptr: *const Expr,
    integral_result_ptr: *const Expr,
    var_ptr: *const c_char,
) -> bool {

    unsafe {

        if integrand_ptr.is_null()
            || integral_result_ptr
                .is_null()
            || var_ptr.is_null()
        {

            return false;
        }

        let var = match CStr::from_ptr(
            var_ptr,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(_) => return false,
        };

        proof::verify_indefinite_integral(
        &*integrand_ptr,
        &*integral_result_ptr,
        var,
    )
    }
}

/// Verifies a definite integral (Handle)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_verify_definite_integral_handle(
    integrand_ptr: *const Expr,
    var_ptr: *const c_char,
    lower: f64,
    upper: f64,
    symbolic_result_ptr: *const Expr,
) -> bool {

    unsafe {

        if integrand_ptr.is_null()
            || var_ptr.is_null()
            || symbolic_result_ptr
                .is_null()
        {

            return false;
        }

        let var = match CStr::from_ptr(
            var_ptr,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(_) => return false,
        };

        proof::verify_definite_integral(
            &*integrand_ptr,
            var,
            (lower, upper),
            &*symbolic_result_ptr,
        )
    }
}

/// Verifies an ODE solution (Handle)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_verify_ode_solution_handle(
    ode_ptr: *const Expr,
    solution_ptr: *const Expr,
    func_name_ptr: *const c_char,
    var_ptr: *const c_char,
) -> bool {

    unsafe {

        if ode_ptr.is_null()
            || solution_ptr.is_null()
            || func_name_ptr.is_null()
            || var_ptr.is_null()
        {

            return false;
        }

        let func_name =
            match CStr::from_ptr(
                func_name_ptr,
            )
            .to_str()
            {
                | Ok(s) => s,
                | Err(_) => {

                    return false;
                },
            };

        let var = match CStr::from_ptr(
            var_ptr,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(_) => return false,
        };

        proof::verify_ode_solution(
            &*ode_ptr,
            &*solution_ptr,
            func_name,
            var,
        )
    }
}

/// Verifies a matrix inverse (Handle)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_verify_matrix_inverse_handle(
    original_ptr: *const Expr,
    inverse_ptr: *const Expr,
) -> bool {

    unsafe {

        if original_ptr.is_null()
            || inverse_ptr.is_null()
        {

            return false;
        }

        proof::verify_matrix_inverse(
            &*original_ptr,
            &*inverse_ptr,
        )
    }
}

/// Verifies a derivative (Handle)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_verify_derivative_handle(
    original_func_ptr: *const Expr,
    derivative_func_ptr: *const Expr,
    var_ptr: *const c_char,
) -> bool {

    unsafe {

        if original_func_ptr.is_null()
            || derivative_func_ptr
                .is_null()
            || var_ptr.is_null()
        {

            return false;
        }

        let var = match CStr::from_ptr(
            var_ptr,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(_) => return false,
        };

        proof::verify_derivative(
            &*original_func_ptr,
            &*derivative_func_ptr,
            var,
        )
    }
}

/// Verifies a limit (Handle)
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_verify_limit_handle(
    f_ptr: *const Expr,
    var_ptr: *const c_char,
    target_ptr: *const Expr,
    limit_val_ptr: *const Expr,
) -> bool {

    unsafe {

        if f_ptr.is_null()
            || var_ptr.is_null()
            || target_ptr.is_null()
            || limit_val_ptr.is_null()
        {

            return false;
        }

        let var = match CStr::from_ptr(
            var_ptr,
        )
        .to_str()
        {
            | Ok(s) => s,
            | Err(_) => return false,
        };

        proof::verify_limit(
            &*f_ptr,
            var,
            &*target_ptr,
            &*limit_val_ptr,
        )
    }
}
