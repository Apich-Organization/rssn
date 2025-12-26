//! Handle-based FFI API for numerical statistics.

use crate::numerical::stats;
use std::slice;

/// Computes the mean of an array.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_mean(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::mean(slice)
}

/// Computes the variance of an array.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_variance(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::variance(slice)
}

/// Computes the standard deviation of an array.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_std_dev(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::std_dev(slice)
}

/// Computes the geometric mean of an array.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_geometric_mean(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::geometric_mean(slice)
}

/// Computes the harmonic mean of an array.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_harmonic_mean(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::harmonic_mean(slice)
}

/// Computes the range of an array.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_range(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::range(slice)
}

/// Computes the coefficient of variation of an array.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_cv(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::coefficient_of_variation(
        slice,
    )
}

/// Computes the standard error of an array.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_standard_error(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::standard_error(slice)
}

/// Computes the Shannon entropy of a probability distribution.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_shannon_entropy(
    data: *const f64,
    len: usize,
) -> f64 {

    if data.is_null() || len == 0 {

        return f64::NAN;
    }

    let slice = slice::from_raw_parts(
        data, len,
    );

    stats::shannon_entropy(slice)
}

/// Computes the covariance of two arrays.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_covariance(
    data1: *const f64,
    len1: usize,
    data2: *const f64,
    len2: usize,
) -> f64 {

    if data1.is_null()
        || data2.is_null()
        || len1 == 0
        || len1 != len2
    {

        return f64::NAN;
    }

    let slice1 = slice::from_raw_parts(
        data1, len1,
    );

    let slice2 = slice::from_raw_parts(
        data2, len2,
    );

    stats::covariance(slice1, slice2)
}

/// Computes the Pearson correlation coefficient of two arrays.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_correlation(
    data1: *const f64,
    len1: usize,
    data2: *const f64,
    len2: usize,
) -> f64 {

    if data1.is_null()
        || data2.is_null()
        || len1 == 0
        || len1 != len2
    {

        return f64::NAN;
    }

    let slice1 = slice::from_raw_parts(
        data1, len1,
    );

    let slice2 = slice::from_raw_parts(
        data2, len2,
    );

    stats::correlation(slice1, slice2)
}

/// Performs a two-sample t-test.
/// Returns t-statistic via out_t and p-value via out_p.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_two_sample_t_test(
    sample1: *const f64,
    len1: usize,
    sample2: *const f64,
    len2: usize,
    out_t: *mut f64,
    out_p: *mut f64,
) -> i32 {

    if sample1.is_null()
        || sample2.is_null()
        || out_t.is_null()
        || out_p.is_null()
    {

        return -1;
    }

    let s1 = slice::from_raw_parts(
        sample1, len1,
    );

    let s2 = slice::from_raw_parts(
        sample2, len2,
    );

    let (t, p) =
        stats::two_sample_t_test(
            s1, s2,
        );

    *out_t = t;

    *out_p = p;

    0
}

/// Performs Welch's t-test.
/// Returns t-statistic via out_t and p-value via out_p.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_welch_t_test(
    sample1: *const f64,
    len1: usize,
    sample2: *const f64,
    len2: usize,
    out_t: *mut f64,
    out_p: *mut f64,
) -> i32 {

    if sample1.is_null()
        || sample2.is_null()
        || out_t.is_null()
        || out_p.is_null()
    {

        return -1;
    }

    let s1 = slice::from_raw_parts(
        sample1, len1,
    );

    let s2 = slice::from_raw_parts(
        sample2, len2,
    );

    let (t, p) =
        stats::welch_t_test(s1, s2);

    *out_t = t;

    *out_p = p;

    0
}

/// Performs a chi-squared test.
/// Returns chi-squared statistic via out_chi and p-value via out_p.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_chi_squared_test(
    observed: *const f64,
    expected: *const f64,
    len: usize,
    out_chi: *mut f64,
    out_p: *mut f64,
) -> i32 {

    if observed.is_null()
        || expected.is_null()
        || out_chi.is_null()
        || out_p.is_null()
    {

        return -1;
    }

    let obs = slice::from_raw_parts(
        observed, len,
    );

    let exp = slice::from_raw_parts(
        expected, len,
    );

    let (chi, p) =
        stats::chi_squared_test(
            obs, exp,
        );

    *out_chi = chi;

    *out_p = p;

    0
}

/// Performs simple linear regression.
/// Returns slope via out_slope and intercept via out_intercept.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_linear_regression(
    x: *const f64,
    y: *const f64,
    len: usize,
    out_slope: *mut f64,
    out_intercept: *mut f64,
) -> i32 {

    if x.is_null()
        || y.is_null()
        || out_slope.is_null()
        || out_intercept.is_null()
    {

        return -1;
    }

    let xs =
        slice::from_raw_parts(x, len);

    let ys =
        slice::from_raw_parts(y, len);

    let data: Vec<(f64, f64)> = xs
        .iter()
        .zip(ys.iter())
        .map(|(&a, &b)| (a, b))
        .collect();

    let (slope, intercept) =
        stats::simple_linear_regression(
            &data,
        );

    *out_slope = slope;

    *out_intercept = intercept;

    0
}
