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
    data : Vec<f64>,
}

#[derive(Deserialize)]

struct TwoDataInput {
    data1 : Vec<f64>,
    data2 : Vec<f64>,
}

#[derive(Deserialize)]

struct RegressionInput {
    x : Vec<f64>,
    y : Vec<f64>,
}

#[derive(Serialize)]

struct RegressionOutput {
    slope : f64,
    intercept : f64,
}

#[derive(Serialize)]

struct TestOutput {
    statistic : f64,
    p_value : f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_mean_json(input : *const c_char) -> *mut c_char {

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

    let result = stats::mean(&input.data);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_variance_json(input : *const c_char) -> *mut c_char {

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

    let result = stats::variance(&input.data);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_std_dev_json(input : *const c_char) -> *mut c_char {

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

    let result = stats::std_dev(&input.data);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_geometric_mean_json(input : *const c_char) -> *mut c_char {

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

    let result = stats::geometric_mean(&input.data);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_harmonic_mean_json(input : *const c_char) -> *mut c_char {

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

    let result = stats::harmonic_mean(&input.data);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_covariance_json(input : *const c_char) -> *mut c_char {

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
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_correlation_json(input : *const c_char) -> *mut c_char {

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
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_two_sample_t_test_json(
    input : *const c_char
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

    let (t, p) = stats::two_sample_t_test(
        &input.data1,
        &input.data2,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(TestOutput {
                statistic : t,
                p_value : p,
            }),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_welch_t_test_json(input : *const c_char) -> *mut c_char {

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
        serde_json::to_string(&FfiResult {
            ok : Some(TestOutput {
                statistic : t,
                p_value : p,
            }),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_chi_squared_test_json(
    input : *const c_char
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

    let (chi, p) = stats::chi_squared_test(
        &input.data1,
        &input.data2,
    );

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(TestOutput {
                statistic : chi,
                p_value : p,
            }),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_linear_regression_json(
    input : *const c_char
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

    let data : Vec<(f64, f64)> = input
        .x
        .iter()
        .zip(input.y.iter())
        .map(|(&a, &b)| (a, b))
        .collect();

    let (slope, intercept) = stats::simple_linear_regression(&data);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(RegressionOutput {
                slope,
                intercept,
            }),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_z_scores_json(input : *const c_char) -> *mut c_char {

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

    let result = stats::z_scores(&input.data);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_stats_shannon_entropy_json(input : *const c_char) -> *mut c_char {

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

    let result = stats::shannon_entropy(&input.data);

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(result),
            err : None::<String>,
        })
        .unwrap(),
    )
}
