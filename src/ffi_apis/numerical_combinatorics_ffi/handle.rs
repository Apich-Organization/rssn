//! Handle-based FFI API for numerical combinatorics.

use std::slice;

use crate::ffi_apis::ffi_api::update_last_error;
use crate::numerical::combinatorics;

/// Computes the factorial of n.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_factorial(
    n: u64,
    result: *mut f64,
) -> i32 {

    if result.is_null() {

        update_last_error(
            "Null pointer passed to \
             rssn_num_comb_factorial"
                .to_string(),
        );

        return -1;
    }

    *result =
        combinatorics::factorial(n);

    0
}

/// Computes the number of permutations P(n, k).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_permutations(
    n: u64,
    k: u64,
    result: *mut f64,
) -> i32 {

    if result.is_null() {

        update_last_error("Null pointer passed to rssn_num_comb_permutations".to_string());

        return -1;
    }

    *result =
        combinatorics::permutations(
            n, k,
        );

    0
}

/// Computes the number of combinations C(n, k).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_combinations(
    n: u64,
    k: u64,
    result: *mut f64,
) -> i32 {

    if result.is_null() {

        update_last_error("Null pointer passed to rssn_num_comb_combinations".to_string());

        return -1;
    }

    *result =
        combinatorics::combinations(
            n, k,
        );

    0
}

/// Solves a linear recurrence relation numerically.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_solve_recurrence(
    coeffs: *const f64,
    coeffs_len: usize,
    initial_conditions: *const f64,
    initial_len: usize,
    target_n: usize,
    result: *mut f64,
) -> i32 {

    if coeffs.is_null()
        || initial_conditions.is_null()
        || result.is_null()
    {

        update_last_error("Null pointer passed to rssn_num_comb_solve_recurrence".to_string());

        return -1;
    }

    let coeffs_slice =
        slice::from_raw_parts(
            coeffs,
            coeffs_len,
        );

    let initial_slice =
        slice::from_raw_parts(
            initial_conditions,
            initial_len,
        );

    match combinatorics::solve_recurrence_numerical(
        coeffs_slice,
        initial_slice,
        target_n,
    ) {
        | Ok(val) => {

            *result = val;

            0
        },
        | Err(e) => {

            update_last_error(e);

            -1
        },
    }
}

/// Computes the Stirling numbers of the second kind S(n, k).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_stirling_second(
    n: u64,
    k: u64,
    result: *mut f64,
) -> i32 {

    if result.is_null() {

        update_last_error("Null pointer passed to rssn_num_comb_stirling_second".to_string());

        return -1;
    }

    *result =
        combinatorics::stirling_second(
            n, k,
        );

    0
}

/// Computes the Bell number B(n).
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_bell(
    n: u64,
    result: *mut f64,
) -> i32 {

    if result.is_null() {

        update_last_error(
            "Null pointer passed to \
             rssn_num_comb_bell"
                .to_string(),
        );

        return -1;
    }

    *result = combinatorics::bell(n);

    0
}

/// Computes the Catalan number `C_n`.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_catalan(
    n: u64,
    result: *mut f64,
) -> i32 {

    if result.is_null() {

        update_last_error(
            "Null pointer passed to \
             rssn_num_comb_catalan"
                .to_string(),
        );

        return -1;
    }

    *result = combinatorics::catalan(n);

    0
}

/// Computes the rising factorial.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_rising_factorial(
    x: f64,
    n: u64,
    result: *mut f64,
) -> i32 {

    if result.is_null() {

        update_last_error("Null pointer passed to rssn_num_comb_rising_factorial".to_string());

        return -1;
    }

    *result =
        combinatorics::rising_factorial(
            x, n,
        );

    0
}

/// Computes the falling factorial.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_num_comb_falling_factorial(
    x: f64,
    n: u64,
    result: *mut f64,
) -> i32 {

    if result.is_null() {

        update_last_error("Null pointer passed to rssn_num_comb_falling_factorial".to_string());

        return -1;
    }

    *result = combinatorics::falling_factorial(x, n);

    0
}
