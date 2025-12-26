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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_mean_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: DataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        stats::mean(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_variance_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: DataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        stats::variance(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_std_dev_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: DataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        stats::std_dev(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_covariance_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoDataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_correlation_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoDataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_two_sample_t_test_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoDataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    TestOutput,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_welch_t_test_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoDataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    TestOutput,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_chi_squared_test_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoDataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    TestOutput,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_linear_regression_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: RegressionInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    RegressionOutput,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
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

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_z_scores_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: DataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    Vec<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
            },
        };

    let result =
        stats::z_scores(&input.data);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_shannon_entropy_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: DataInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    f64,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode".to_string(),
                    ),
                })
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
