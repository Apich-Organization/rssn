//! Bincode-based FFI API for numerical special functions.

use crate::ffi_apis::common::{
    from_bincode_buffer,
    to_bincode_buffer,
    BincodeBuffer,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::special;
use serde::Deserialize;

#[derive(Deserialize)]

struct SingleInput {
    x: f64,
}

#[derive(Deserialize)]

struct TwoInput {
    a: f64,
    b: f64,
}

#[derive(Deserialize)]

struct ThreeInput {
    x: f64,
    a: f64,
    b: f64,
}

#[derive(Deserialize)]

struct PolyInput {
    n: u32,
    x: f64,
}

#[derive(Deserialize)]

struct IntInput {
    n: u64,
}

#[derive(Deserialize)]

struct BinomialInput {
    n: u64,
    k: u64,
}

// Gamma functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_gamma_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::gamma_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_ln_gamma_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::ln_gamma_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_digamma_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::digamma_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

// Beta functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_beta_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: TwoInput = match from_bincode_buffer(&buffer)
    {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                f64,
                String,
            > {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::beta_numerical(
                input.a, input.b,
            ),
        ),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_regularized_beta_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: ThreeInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::regularized_beta(
                input.x, input.a,
                input.b,
            ),
        ),
        err: None::<String>,
    })
}

// Error functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_erf_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::erf_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_erfc_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::erfc_numerical(
                input.x,
            ),
        ),
        err: None::<String>,
    })
}

// Bessel functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_bessel_j0_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(special::bessel_j0(
            input.x,
        )),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_bessel_j1_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(special::bessel_j1(
            input.x,
        )),
        err: None::<String>,
    })
}

// Orthogonal polynomials
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_legendre_p_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: PolyInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(special::legendre_p(
            input.n, input.x,
        )),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_chebyshev_t_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: PolyInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(
            special::chebyshev_t(
                input.n, input.x,
            ),
        ),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_hermite_h_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: PolyInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(special::hermite_h(
            input.n, input.x,
        )),
        err: None::<String>,
    })
}

// Other special functions
#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_factorial_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: IntInput = match from_bincode_buffer(&buffer)
    {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                f64,
                String,
            > {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        },
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(special::factorial(
            input.n,
        )),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_binomial_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: BinomialInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(special::binomial(
            input.n, input.k,
        )),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_sigmoid_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(special::sigmoid(
            input.x,
        )),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_special_sinc_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: SingleInput =
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

    to_bincode_buffer(&FfiResult {
        ok: Some(special::sinc(
            input.x,
        )),
        err: None::<String>,
    })
}
