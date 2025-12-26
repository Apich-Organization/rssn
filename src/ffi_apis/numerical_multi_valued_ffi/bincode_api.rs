//! Bincode-based FFI API for numerical multi-valued functions.

use crate::ffi_apis::common::{
    from_bincode_buffer,
    to_bincode_buffer,
    BincodeBuffer,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::multi_valued;
use crate::symbolic::core::Expr;
use num_complex::Complex;
use serde::{
    Deserialize,
    Serialize,
};

#[derive(Deserialize)]

struct NewtonInput {
    f: Expr,
    f_prime: Expr,
    start_re: f64,
    start_im: f64,
    tolerance: f64,
    max_iter: usize,
}

#[derive(Serialize)]

struct ComplexResult {
    re: f64,
    im: f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_newton_method_complex_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: NewtonInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    ComplexResult,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode input".to_string(),
                    ),
                })
            },
        };

    let start_point = Complex::new(
        input.start_re,
        input.start_im,
    );

    match multi_valued::newton_method_complex(
        &input.f,
        &input.f_prime,
        start_point,
        input.tolerance,
        input.max_iter,
    ) {
        | Some(root) => {

            let res = ComplexResult {
                re: root.re,
                im: root.im,
            };

            to_bincode_buffer(&FfiResult {
                ok: Some(res),
                err: None::<String>,
            })
        },
        | None => {
            to_bincode_buffer(&FfiResult::<
                ComplexResult,
                String,
            > {
                ok: None,
                err: Some(
                    "Newton's method failed to converge"
                        .to_string(),
                ),
            })
        },
    }
}

#[derive(Deserialize)]

struct LogSqrtInput {
    re: f64,
    im: f64,
    k: i32,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_log_k_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: LogSqrtInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    ComplexResult,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode input".to_string(),
                    ),
                })
            },
        };

    let z = Complex::new(
        input.re, input.im,
    );

    let res =
        multi_valued::complex_log_k(
            z, input.k,
        );

    let out = ComplexResult {
        re: res.re,
        im: res.im,
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(out),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_sqrt_k_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: LogSqrtInput =
        match from_bincode_buffer(&buffer) {
            | Some(i) => i,
            | None => {
                return to_bincode_buffer(&FfiResult::<
                    ComplexResult,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Invalid Bincode input".to_string(),
                    ),
                })
            },
        };

    let z = Complex::new(
        input.re, input.im,
    );

    let res =
        multi_valued::complex_sqrt_k(
            z, input.k,
        );

    let out = ComplexResult {
        re: res.re,
        im: res.im,
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(out),
        err: None::<String>,
    })
}

#[derive(Deserialize)]

struct PowInput {
    z_re: f64,
    z_im: f64,
    w_re: f64,
    w_im: f64,
    k: i32,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_mv_complex_pow_k_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: PowInput = match from_bincode_buffer(&buffer)
    {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(&FfiResult::<
                ComplexResult,
                String,
            > {
                ok: None,
                err: Some(
                    "Invalid Bincode input".to_string(),
                ),
            })
        },
    };

    let z = Complex::new(
        input.z_re, input.z_im,
    );

    let w = Complex::new(
        input.w_re, input.w_im,
    );

    let res =
        multi_valued::complex_pow_k(
            z, w, input.k,
        );

    let out = ComplexResult {
        re: res.re,
        im: res.im,
    };

    to_bincode_buffer(&FfiResult {
        ok: Some(out),
        err: None::<String>,
    })
}
