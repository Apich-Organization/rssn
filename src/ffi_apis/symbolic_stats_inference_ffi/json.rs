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

#[no_mangle]

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

    if let (Some(data), Some(target)) =
        (data, target)
    {

        let result = stats_inference::one_sample_t_test_symbolic(&data, &target);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Performs a two-sample t-test.

///

/// Takes JSON strings representing `Vec<Expr>` (two data sets) and `Expr` (hypothesized difference in means).

/// Returns a JSON string representing the `HypothesisTest` result.

#[no_mangle]

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

    if let (
        Some(d1),
        Some(d2),
        Some(diff),
    ) = (data1, data2, diff)
    {

        let result = stats_inference::two_sample_t_test_symbolic(&d1, &d2, &diff);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Performs a z-test.

///

/// Takes JSON strings representing `Vec<Expr>` (data), `Expr` (target mean), and `Expr` (population standard deviation).

/// Returns a JSON string representing the `HypothesisTest` result.

#[no_mangle]

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

    if let (
        Some(data),
        Some(target),
        Some(sigma),
    ) = (data, target, sigma)
    {

        let result = stats_inference::z_test_symbolic(
            &data,
            &target,
            &sigma,
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
