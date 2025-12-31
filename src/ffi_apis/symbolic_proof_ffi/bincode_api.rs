use std::collections::HashMap;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::proof;

/// Verifies an equation solution using Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_verify_equation_solution(
    equations_buf: BincodeBuffer,
    solution_buf: BincodeBuffer,
    free_vars_buf: BincodeBuffer,
) -> bool {

    let equations: Option<Vec<Expr>> =
        from_bincode_buffer(
            &equations_buf,
        );

    let solution: Option<
        HashMap<String, Expr>,
    > = from_bincode_buffer(
        &solution_buf,
    );

    let free_vars: Option<Vec<String>> =
        from_bincode_buffer(
            &free_vars_buf,
        );

    match (
        equations,
        solution,
        free_vars,
    ) {
        | (
            Some(eqs),
            Some(sol),
            Some(free),
        ) => {

            let free_refs: Vec<&str> = free
            .iter()
            .map(std::string::String::as_str)
            .collect();

            proof::verify_equation_solution(
            &eqs,
            &sol,
            &free_refs,
        )
        },
        | _ => false,
    }
}

/// Verifies an indefinite integral using Bincode.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_verify_indefinite_integral(
    integrand_buf: BincodeBuffer,
    integral_result_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
) -> bool {

    let integrand: Option<Expr> =
        from_bincode_buffer(
            &integrand_buf,
        );

    let integral_result: Option<Expr> =
        from_bincode_buffer(
            &integral_result_buf,
        );

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    match (
        integrand,
        integral_result,
        var,
    ) { (
        Some(f),
        Some(int),
        Some(v),
    ) => {

        proof::verify_indefinite_integral(&f, &int, &v)
    } _ => {

        false
    }}
}
