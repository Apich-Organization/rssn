//! Bincode-based FFI API for numerical functional analysis.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::functional_analysis;

#[derive(Deserialize)]

struct PointsInput {
    points: Vec<(f64, f64)>,
}

#[derive(Deserialize)]

struct InnerProductInput {
    f: Vec<(f64, f64)>,
    g: Vec<(f64, f64)>,
}

#[derive(Deserialize)]

struct GramSchmidtInput {
    basis: Vec<Vec<(f64, f64)>>,
}

/// Computes the L2 norm of a function (represented by a series of points) using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_l2_norm_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : PointsInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    let res =
        functional_analysis::l2_norm(
            &input.points,
        );

    let ffi_res = FfiResult {
        ok: Some(res),
        err: None::<String>,
    };

    to_bincode_buffer(&ffi_res)
}

/// Computes the inner product of two functions (represented by series of points) using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_inner_product_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : InnerProductInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    match functional_analysis::inner_product(&input.f, &input.g) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult {
                ok : Some(res),
                err : None::<String>,
            })
        },
        | Err(e) => {
            to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok : None,
                    err : Some(e),
                },
            )
        },
    }
}

/// Applies the Gram-Schmidt orthonormalization process to a set of basis functions using bincode for serialization.

#[no_mangle]

pub unsafe extern "C" fn rssn_num_fa_gram_schmidt_bincode(

    buffer: BincodeBuffer

) -> BincodeBuffer {

    let input : GramSchmidtInput = match from_bincode_buffer(&buffer) {
        | Some(i) => i,
        | None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<(f64, f64)>>, String> {
                    ok : None,
                    err : Some("Invalid Bincode input".to_string()),
                },
            )
        },
    };

    match functional_analysis::gram_schmidt(&input.basis) {
        | Ok(res) => {
            to_bincode_buffer(&FfiResult {
                ok : Some(res),
                err : None::<String>,
            })
        },
        | Err(e) => {
            to_bincode_buffer(
                &FfiResult::<Vec<Vec<(f64, f64)>>, String> {
                    ok : None,
                    err : Some(e),
                },
            )
        },
    }
}
