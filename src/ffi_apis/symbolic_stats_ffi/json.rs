//! JSON-based FFI API for symbolic statistics functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::stats;

/// Computes the symbolic mean of a set of expressions using JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_mean(
    data_json: *const c_char
) -> *mut c_char {

    let data: Option<Vec<Expr>> =
        from_json_string(data_json);

    if let Some(d) = data {

        let result = stats::mean(&d);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic variance of a set of expressions using JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_variance(
    data_json: *const c_char
) -> *mut c_char {

    let data: Option<Vec<Expr>> =
        from_json_string(data_json);

    if let Some(d) = data {

        let result =
            stats::variance(&d);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic standard deviation of a set of expressions using JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_std_dev(
    data_json: *const c_char
) -> *mut c_char {

    let data: Option<Vec<Expr>> =
        from_json_string(data_json);

    if let Some(d) = data {

        let result = stats::std_dev(&d);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the symbolic covariance of two sets of expressions using JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_covariance(
    data1_json: *const c_char,
    data2_json: *const c_char,
) -> *mut c_char {

    let data1: Option<Vec<Expr>> =
        from_json_string(data1_json);

    let data2: Option<Vec<Expr>> =
        from_json_string(data2_json);

    match (data1, data2)
    { (Some(d1), Some(d2)) => {

        let result =
            stats::covariance(&d1, &d2);

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the symbolic Pearson correlation coefficient using JSON.
#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_correlation(
    data1_json: *const c_char,
    data2_json: *const c_char,
) -> *mut c_char {

    let data1: Option<Vec<Expr>> =
        from_json_string(data1_json);

    let data2: Option<Vec<Expr>> =
        from_json_string(data2_json);

    match (data1, data2)
    { (Some(d1), Some(d2)) => {

        let result = stats::correlation(
            &d1, &d2,
        );

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}
