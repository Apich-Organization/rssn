use std::collections::HashMap;
use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::proof;

/// Verifies an equation solution using JSON.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_verify_equation_solution(
    equations_json: *const c_char,
    solution_json: *const c_char,
    free_vars_json: *const c_char,
) -> bool {

    let equations: Option<Vec<Expr>> =
        from_json_string(
            equations_json,
        );

    let solution: Option<
        HashMap<String, Expr>,
    > = from_json_string(solution_json);

    let free_vars: Option<Vec<String>> =
        from_json_string(free_vars_json);

    if let (
        Some(eqs),
        Some(sol),
        Some(free),
    ) = (
        equations,
        solution,
        free_vars,
    ) {

        let free_refs: Vec<&str> = free
            .iter()
            .map(std::string::String::as_str)
            .collect();

        proof::verify_equation_solution(
            &eqs,
            &sol,
            &free_refs,
        )
    } else {

        false
    }
}

/// Verifies an indefinite integral using JSON.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_verify_indefinite_integral(
    integrand_json: *const c_char,
    integral_result_json: *const c_char,
    var_json: *const c_char,
) -> bool {

    let integrand: Option<Expr> =
        from_json_string(
            integrand_json,
        );

    let integral_result: Option<Expr> =
        from_json_string(
            integral_result_json,
        );

    let var: Option<String> =
        from_json_string(var_json);

    if let (
        Some(f),
        Some(int),
        Some(v),
    ) = (
        integrand,
        integral_result,
        var,
    ) {

        proof::verify_indefinite_integral(&f, &int, &v)
    } else {

        false
    }
}

/// Verifies a matrix inverse using JSON.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_verify_matrix_inverse(
    original_json: *const c_char,
    inverse_json: *const c_char,
) -> bool {

    let original: Option<Expr> =
        from_json_string(original_json);

    let inverse: Option<Expr> =
        from_json_string(inverse_json);

    if let (Some(orig), Some(inv)) =
        (original, inverse)
    {

        proof::verify_matrix_inverse(
            &orig, &inv,
        )
    } else {

        false
    }
}
