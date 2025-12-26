//! Bincode-based FFI API for numerical differential geometry.

use crate::ffi_apis::common::{
    from_bincode_buffer,
    to_bincode_buffer,
    BincodeBuffer,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::differential_geometry;
use crate::symbolic::coordinates::CoordinateSystem;
use serde::{
    Deserialize,
    Serialize,
};

#[derive(Deserialize)]

struct DgPointInput {
    system: CoordinateSystem,
    point: Vec<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_metric_tensor_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: DgPointInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<f64>>, String> {
                    ok: None,
                    err: Some("Invalid Bincode input".to_string()),
                },
            )
        }
    };

    match differential_geometry::metric_tensor_at_point(
        input.system,
        &input.point,
    ) {
        Ok(res) => {

            let ffi_res = FfiResult {
                ok: Some(res),
                err: None::<String>,
            };

            to_bincode_buffer(&ffi_res)
        }
        Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<Vec<Vec<f64>>>,
                err: Some(e),
            };

            to_bincode_buffer(&ffi_res)
        }
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_christoffel_symbols_bincode(
    buffer: BincodeBuffer
) -> BincodeBuffer {

    let input: DgPointInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<Vec<f64>>>, String> {
                    ok: None,
                    err: Some("Invalid Bincode input".to_string()),
                },
            )
        }
    };

    match differential_geometry::christoffel_symbols(
        input.system,
        &input.point,
    ) {
        Ok(res) => {

            let ffi_res = FfiResult {
                ok: Some(res),
                err: None::<String>,
            };

            to_bincode_buffer(&ffi_res)
        }
        Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<Vec<Vec<Vec<f64>>>>,
                err: Some(e),
            };

            to_bincode_buffer(&ffi_res)
        }
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_ricci_tensor_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: DgPointInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<Vec<Vec<f64>>, String> {
                    ok: None,
                    err: Some("Invalid Bincode input".to_string()),
                },
            )
        }
    };

    match differential_geometry::ricci_tensor(
        input.system,
        &input.point,
    ) {
        Ok(res) => {

            let ffi_res = FfiResult {
                ok: Some(res),
                err: None::<String>,
            };

            to_bincode_buffer(&ffi_res)
        }
        Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<Vec<Vec<f64>>>,
                err: Some(e),
            };

            to_bincode_buffer(&ffi_res)
        }
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_ricci_scalar_bincode(buffer: BincodeBuffer) -> BincodeBuffer {

    let input: DgPointInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(
                &FfiResult::<f64, String> {
                    ok: None,
                    err: Some("Invalid Bincode input".to_string()),
                },
            )
        }
    };

    match differential_geometry::ricci_scalar(
        input.system,
        &input.point,
    ) {
        Ok(res) => {

            let ffi_res = FfiResult {
                ok: Some(res),
                err: None::<String>,
            };

            to_bincode_buffer(&ffi_res)
        }
        Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<f64>,
                err: Some(e),
            };

            to_bincode_buffer(&ffi_res)
        }
    }
}
