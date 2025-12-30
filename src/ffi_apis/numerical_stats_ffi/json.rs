//! JSON-based FFI API for numerical statistics.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::stats;

#[derive(Deserialize)]

struct DataInput {
    data: Vec<f64>,
}

#[derive(Deserialize)]

struct TwoDataInput {
    data1: Vec<f64>,
    data2: Vec<f64>,
}

#[derive(Deserialize)]

struct RegressionInput {
    x: Vec<f64>,
    y: Vec<f64>,
}

#[derive(Serialize)]

struct RegressionOutput {
    slope: f64,
    intercept: f64,
}

#[derive(Serialize)]

struct TestOutput {
    statistic: f64,
    p_value: f64,
}

/// Computes the arithmetic mean (average) of a dataset via JSON serialization.
///
/// The arithmetic mean is defined as μ = (1/n) Σxᵢ.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data`: Array of numerical values
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the mean value μ.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_mean_json(
    input: *const c_char
) -> *mut c_char {

    let input : DataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        stats::mean(&input.data);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the sample variance of a dataset via JSON serialization.
///
/// The sample variance is defined as s² = (1/(n-1)) Σ(xᵢ - μ)².
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data`: Array of numerical values
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the variance s².
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_variance_json(
    input: *const c_char
) -> *mut c_char {

    let input : DataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        stats::variance(&input.data);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the sample standard deviation of a dataset via JSON serialization.
///
/// The standard deviation is defined as s = √(s²) where s² is the sample variance.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data`: Array of numerical values
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the standard deviation s.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_std_dev_json(
    input: *const c_char
) -> *mut c_char {

    let input : DataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        stats::std_dev(&input.data);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the geometric mean of a dataset via JSON serialization.
///
/// The geometric mean is defined as (∏xᵢ)^(1/n), useful for quantities with
/// multiplicative relationships (e.g., growth rates).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data`: Array of positive numerical values
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the geometric mean.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_geometric_mean_json(
    input: *const c_char
) -> *mut c_char {

    let input : DataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = stats::geometric_mean(
        &input.data,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the harmonic mean of a dataset via JSON serialization.
///
/// The harmonic mean is defined as n / Σ(1/xᵢ), useful for averaging rates
/// and ratios.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data`: Array of positive numerical values
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the harmonic mean.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_harmonic_mean_json(
    input: *const c_char
) -> *mut c_char {

    let input : DataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = stats::harmonic_mean(
        &input.data,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the covariance between two datasets via JSON serialization.
///
/// The sample covariance is defined as Cov(X,Y) = (1/(n-1)) Σ(xᵢ - μₓ)(yᵢ - μᵧ).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data1`: First dataset
///   - `data2`: Second dataset (must have same length as data1)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the covariance Cov(X,Y).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_covariance_json(
    input: *const c_char
) -> *mut c_char {

    let input : TwoDataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = stats::covariance(
        &input.data1,
        &input.data2,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes the Pearson correlation coefficient between two datasets via JSON serialization.
///
/// The Pearson correlation is defined as ρ = Cov(X,Y) / (σₓσᵧ), ranging from -1 to 1.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data1`: First dataset
///   - `data2`: Second dataset (must have same length as data1)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// the correlation coefficient ρ ∈ [-1, 1].
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_correlation_json(
    input: *const c_char
) -> *mut c_char {

    let input : TwoDataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = stats::correlation(
        &input.data1,
        &input.data2,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Performs a two-sample t-test for equal means via JSON serialization.
///
/// Tests the null hypothesis that two independent samples have equal means,
/// assuming equal variances (pooled variance).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data1`: First sample
///   - `data2`: Second sample
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<TestOutput, String>` with
/// `statistic` (t-statistic) and `p_value`.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_two_sample_t_test_json(
    input: *const c_char
) -> *mut c_char {

    let input : TwoDataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<TestOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let (t, p) =
        stats::two_sample_t_test(
            &input.data1,
            &input.data2,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(TestOutput {
                    statistic: t,
                    p_value: p,
                }),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Performs Welch's t-test for unequal variances via JSON serialization.
///
/// Tests the null hypothesis that two independent samples have equal means,
/// without assuming equal variances (Welch-Satterthwaite correction).
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data1`: First sample
///   - `data2`: Second sample
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<TestOutput, String>` with
/// `statistic` (t-statistic) and `p_value`.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_welch_t_test_json(
    input: *const c_char
) -> *mut c_char {

    let input : TwoDataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<TestOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let (t, p) = stats::welch_t_test(
        &input.data1,
        &input.data2,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(TestOutput {
                    statistic: t,
                    p_value: p,
                }),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Performs a chi-squared goodness-of-fit test via JSON serialization.
///
/// Tests whether observed frequencies match expected frequencies according to
/// the test statistic χ² = Σ(Oᵢ - Eᵢ)² / Eᵢ.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data1`: Observed frequencies
///   - `data2`: Expected frequencies
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<TestOutput, String>` with
/// `statistic` (χ²-statistic) and `p_value`.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_chi_squared_test_json(
    input: *const c_char
) -> *mut c_char {

    let input : TwoDataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<TestOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let (chi, p) =
        stats::chi_squared_test(
            &input.data1,
            &input.data2,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(TestOutput {
                    statistic: chi,
                    p_value: p,
                }),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes simple linear regression using least squares method via JSON serialization.
///
/// Fits a line y = mx + b to the data by minimizing Σ(yᵢ - (mxᵢ + b))².
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `x`: Independent variable values
///   - `y`: Dependent variable values (must have same length as x)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<RegressionOutput, String>` with
/// `slope` (m) and `intercept` (b).
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_linear_regression_json(
    input: *const c_char
) -> *mut c_char {

    let input : RegressionInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<RegressionOutput, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let data: Vec<(f64, f64)> = input
        .x
        .iter()
        .zip(input.y.iter())
        .map(|(&a, &b)| (a, b))
        .collect();

    let (slope, intercept) =
        stats::simple_linear_regression(
            &data,
        );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(
                    RegressionOutput {
                        slope,
                        intercept,
                    },
                ),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes standardized z-scores for a dataset via JSON serialization.
///
/// The z-score is defined as z = (x - μ) / σ, representing the number of
/// standard deviations each value is from the mean.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data`: Array of numerical values
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<Vec<f64>, String>` with
/// the vector of z-scores.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_z_scores_json(
    input: *const c_char
) -> *mut c_char {

    let input : DataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<Vec<f64>, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result =
        stats::z_scores(&input.data);

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}

/// Computes Shannon entropy of a probability distribution via JSON serialization.
///
/// Shannon entropy is defined as H(X) = -Σ p(xᵢ) log₂(p(xᵢ)), measuring
/// the average information content or uncertainty in the distribution.
///
/// # Arguments
///
/// * `input` - A JSON string pointer containing:
///   - `data`: Probability distribution (values should sum to 1)
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `FfiResult<f64, String>` with
/// Shannon entropy H(X) in bits.
///
/// # Safety
///
/// This function is unsafe because it receives a raw C string pointer that must be
/// valid, null-terminated UTF-8. The caller must free the returned pointer.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_num_stats_shannon_entropy_json(
    input: *const c_char
) -> *mut c_char {

    let input : DataInput = match from_json_string(input) {
        | Some(i) => i,
        | None => {
            return to_c_string(
                serde_json::to_string(
                    &FfiResult::<f64, String> {
                        ok : None,
                        err : Some("Invalid JSON".to_string()),
                    },
                )
                .unwrap(),
            )
        },
    };

    let result = stats::shannon_entropy(
        &input.data,
    );

    to_c_string(
        serde_json::to_string(
            &FfiResult {
                ok: Some(result),
                err: None::<String>,
            },
        )
        .unwrap(),
    )
}
