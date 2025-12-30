use std::os::raw::c_char;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::stats_inference::{
    self,
};

// For JSON, we can return the full HypothesisTest struct serialized.

/// Performs a one-sample t-test.

///

/// Takes JSON strings representing `Vec<Expr>` (data) and `Expr` (target mean).

/// Returns a JSON string representing the `HypothesisTest` result.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_one_sample_t_test(
    data_json: *const c_char,
    target_mean_json: *const c_char,
) -> *mut c_char {

    let data: Option<Vec<Expr>> =
        from_json_string(data_json);

    let target: Option<Expr> =
        from_json_string(
            target_mean_json,
        );

    match (data, target)
    { (Some(data), Some(target)) => {

        let result = stats_inference::one_sample_t_test_symbolic(&data, &target);

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Performs a two-sample t-test.

///

/// Takes JSON strings representing `Vec<Expr>` (two data sets) and `Expr` (hypothesized difference in means).

/// Returns a JSON string representing the `HypothesisTest` result.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_two_sample_t_test(
    data1_json: *const c_char,
    data2_json: *const c_char,
    mu_diff_json: *const c_char,
) -> *mut c_char {

    let data1: Option<Vec<Expr>> =
        from_json_string(data1_json);

    let data2: Option<Vec<Expr>> =
        from_json_string(data2_json);

    let diff: Option<Expr> =
        from_json_string(mu_diff_json);

    match (data1, data2, diff)
    { (
        Some(d1),
        Some(d2),
        Some(diff),
    ) => {

        let result = stats_inference::two_sample_t_test_symbolic(&d1, &d2, &diff);

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Performs a z-test.

///

/// Takes JSON strings representing `Vec<Expr>` (data), `Expr` (target mean), and `Expr` (population standard deviation).

/// Returns a JSON string representing the `HypothesisTest` result.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_z_test(
    data_json: *const c_char,
    target_mean_json: *const c_char,
    pop_std_dev_json: *const c_char,
) -> *mut c_char {

    let data: Option<Vec<Expr>> =
        from_json_string(data_json);

    let target: Option<Expr> =
        from_json_string(
            target_mean_json,
        );

    let sigma: Option<Expr> =
        from_json_string(
            pop_std_dev_json,
        );

    match (data, target, sigma)
    { (
        Some(data),
        Some(target),
        Some(sigma),
    ) => {

        let result = stats_inference::z_test_symbolic(
            &data,
            &target,
            &sigma,
        );

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}
