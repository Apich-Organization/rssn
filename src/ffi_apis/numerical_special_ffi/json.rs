//! JSON-based FFI API for numerical special functions.

use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::special;

#[derive(Deserialize)]

struct SingleInput {
    x : f64,
}

#[derive(Deserialize)]

struct TwoInput {
    a : f64,
    b : f64,
}

#[derive(Deserialize)]

struct ThreeInput {
    x : f64,
    a : f64,
    b : f64,
}

#[derive(Deserialize)]

struct PolyInput {
    n : u32,
    x : f64,
}

#[derive(Deserialize)]

struct IntInput {
    n : u64,
}

#[derive(Deserialize)]

struct BinomialInput {
    n : u64,
    k : u64,
}

// Gamma functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_gamma_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::gamma_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_ln_gamma_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::ln_gamma_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_digamma_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::digamma_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// Beta functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_beta_json(input : *const c_char) -> *mut c_char {

    let input : TwoInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::beta_numerical(input.a, input.b)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_regularized_beta_json(
    input : *const c_char
) -> *mut c_char {

    let input : ThreeInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(
                special::regularized_beta(
                    input.x,
                    input.a,
                    input.b,
                ),
            ),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// Error functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_erf_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::erf_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_erfc_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::erfc_numerical(input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// Bessel functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_bessel_j0_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::bessel_j0(
                input.x,
            )),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_bessel_j1_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::bessel_j1(
                input.x,
            )),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// Orthogonal polynomials
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_legendre_p_json(input : *const c_char) -> *mut c_char {

    let input : PolyInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::legendre_p(
                input.n,
                input.x,
            )),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_chebyshev_t_json(input : *const c_char) -> *mut c_char {

    let input : PolyInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::chebyshev_t(input.n, input.x)),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_hermite_h_json(input : *const c_char) -> *mut c_char {

    let input : PolyInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::hermite_h(
                input.n,
                input.x,
            )),
            err : None::<String>,
        })
        .unwrap(),
    )
}

// Other special functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_factorial_json(input : *const c_char) -> *mut c_char {

    let input : IntInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::factorial(
                input.n,
            )),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_binomial_json(input : *const c_char) -> *mut c_char {

    let input : BinomialInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::binomial(
                input.n,
                input.k,
            )),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_sigmoid_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::sigmoid(
                input.x,
            )),
            err : None::<String>,
        })
        .unwrap(),
    )
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_sinc_json(input : *const c_char) -> *mut c_char {

    let input : SingleInput = match from_json_string(input) {
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

    to_c_string(
        serde_json::to_string(&FfiResult {
            ok : Some(special::sinc(
                input.x,
            )),
            err : None::<String>,
        })
        .unwrap(),
    )
}
