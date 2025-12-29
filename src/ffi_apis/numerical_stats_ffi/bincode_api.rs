//! Bincode-based FFI API for numerical statistics.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
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

/// Computes the arithmetic mean (average) of a dataset using bincode serialization.
///
/// The arithmetic mean is defined as μ = (1/n) Σxᵢ.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `DataInput` with:
///   - `data`: Vector of numerical values
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The mean value μ
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_mean_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : DataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result =
        stats::mean(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the sample variance of a dataset using bincode serialization.
///
/// The sample variance is defined as s² = (1/(n-1)) Σ(xᵢ - μ)².
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `DataInput` with:
///   - `data`: Vector of numerical values
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The variance s²
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_variance_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : DataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result =
        stats::variance(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the sample standard deviation of a dataset using bincode serialization.
///
/// The standard deviation is defined as s = √(s²) where s² is the sample variance.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `DataInput` with:
///   - `data`: Vector of numerical values
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The standard deviation s
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_std_dev_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : DataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result =
        stats::std_dev(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the covariance between two datasets using bincode serialization.
///
/// The sample covariance is defined as Cov(X,Y) = (1/(n-1)) Σ(xᵢ - μₓ)(yᵢ - μᵧ).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TwoDataInput` with:
///   - `data1`: First dataset
///   - `data2`: Second dataset (must have same length as data1)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The covariance Cov(X,Y)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_covariance_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : TwoDataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result = stats::covariance(
        &input.data1,
        &input.data2,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes the Pearson correlation coefficient between two datasets using bincode serialization.
///
/// The Pearson correlation is defined as ρ = Cov(X,Y) / (σₓσᵧ), ranging from -1 to 1.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TwoDataInput` with:
///   - `data1`: First dataset
///   - `data2`: Second dataset (must have same length as data1)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: The correlation coefficient ρ ∈ [-1, 1]
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_correlation_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : TwoDataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result = stats::correlation(
        &input.data1,
        &input.data2,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Performs a two-sample t-test for equal means using bincode serialization.
///
/// Tests the null hypothesis that two independent samples have equal means,
/// assuming equal variances (pooled variance).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TwoDataInput` with:
///   - `data1`: First sample
///   - `data2`: Second sample
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<TestOutput, String>` with either:
/// - `ok`: Object containing `statistic` (t-statistic) and `p_value`
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_two_sample_t_test_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : TwoDataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<TestOutput, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let (t, p) =
        stats::two_sample_t_test(
            &input.data1,
            &input.data2,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(TestOutput {
            statistic: t,
            p_value: p,
        }),
        err: None::<String>,
    })
}

/// Performs Welch's t-test for unequal variances using bincode serialization.
///
/// Tests the null hypothesis that two independent samples have equal means,
/// without assuming equal variances (Welch-Satterthwaite correction).
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TwoDataInput` with:
///   - `data1`: First sample
///   - `data2`: Second sample
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<TestOutput, String>` with either:
/// - `ok`: Object containing `statistic` (t-statistic) and `p_value`
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_welch_t_test_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : TwoDataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<TestOutput, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let (t, p) = stats::welch_t_test(
        &input.data1,
        &input.data2,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(TestOutput {
            statistic: t,
            p_value: p,
        }),
        err: None::<String>,
    })
}

/// Performs a chi-squared goodness-of-fit test using bincode serialization.
///
/// Tests whether observed frequencies match expected frequencies according to
/// the test statistic χ² = Σ(Oᵢ - Eᵢ)² / Eᵢ.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `TwoDataInput` with:
///   - `data1`: Observed frequencies
///   - `data2`: Expected frequencies
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<TestOutput, String>` with either:
/// - `ok`: Object containing `statistic` (χ²-statistic) and `p_value`
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_chi_squared_test_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : TwoDataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<TestOutput, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let (chi, p) =
        stats::chi_squared_test(
            &input.data1,
            &input.data2,
        );

    to_bincode_buffer(&FfiResult {
        ok: Some(TestOutput {
            statistic: chi,
            p_value: p,
        }),
        err: None::<String>,
    })
}

/// Computes simple linear regression using least squares method via bincode serialization.
///
/// Fits a line y = mx + b to the data by minimizing Σ(yᵢ - (mxᵢ + b))².
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `RegressionInput` with:
///   - `x`: Independent variable values
///   - `y`: Dependent variable values (must have same length as x)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<RegressionOutput, String>` with either:
/// - `ok`: Object containing `slope` (m) and `intercept` (b)
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_linear_regression_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : RegressionInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<RegressionOutput, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
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

    to_bincode_buffer(&FfiResult {
        ok: Some(RegressionOutput {
            slope,
            intercept,
        }),
        err: None::<String>,
    })
}

/// Computes standardized z-scores for a dataset using bincode serialization.
///
/// The z-score is defined as z = (x - μ) / σ, representing the number of
/// standard deviations each value is from the mean.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `DataInput` with:
///   - `data`: Vector of numerical values
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<Vec<f64>, String>` with either:
/// - `ok`: Vector of z-scores
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_z_scores_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : DataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<f64>, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result =
        stats::z_scores(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

/// Computes Shannon entropy of a probability distribution using bincode serialization.
///
/// Shannon entropy is defined as H(X) = -Σ p(xᵢ) log₂(p(xᵢ)), measuring
/// the average information content or uncertainty in the distribution.
///
/// # Arguments
///
/// * `buffer` - A bincode-encoded buffer containing `DataInput` with:
///   - `data`: Probability distribution (values should sum to 1)
///
/// # Returns
///
/// A bincode-encoded buffer containing `FfiResult<f64, String>` with either:
/// - `ok`: Shannon entropy H(X) in bits
/// - `err`: Error message if input is invalid
///
/// # Safety
///
/// This function is unsafe because it receives a raw bincode buffer that must be
/// valid and properly encoded.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_shannon_entropy_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input : DataInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode".to_string()),
                },
            )
        },
    };

    let result = stats::shannon_entropy(
        &input.data,
    );

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}
