//! JSON-based FFI API for numerical differential geometry.

use std::os::raw::c_char;

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_c_string;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::differential_geometry;
use crate::symbolic::coordinates::CoordinateSystem;

#[derive(Deserialize)]

struct DgPointInput {
    system : CoordinateSystem,
    point : Vec<f64>,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_metric_tensor_json(
    input_json : *const c_char
) -> *mut c_char {

    let input: DgPointInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<Vec<f64>>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON input"
                                .to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    match differential_geometry::metric_tensor_at_point(
        input.system,
        &input.point,
    ) {
        | Ok(res) => {

            let ffi_res = FfiResult {
                ok: Some(res),
                err: None::<String>,
            };

            to_c_string(
                serde_json::to_string(&ffi_res).unwrap(),
            )
        },
        | Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<Vec<Vec<f64>>>,
                err: Some(e),
            };

            to_c_string(
                serde_json::to_string(&ffi_res).unwrap(),
            )
        },
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_christoffel_symbols_json(
    input_json : *const c_char
) -> *mut c_char {

    let input: DgPointInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<Vec<Vec<f64>>>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON input"
                                .to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    match differential_geometry::christoffel_symbols(
        input.system,
        &input.point,
    ) {
        | Ok(res) => {

            let ffi_res = FfiResult {
                ok: Some(res),
                err: None::<String>,
            };

            to_c_string(
                serde_json::to_string(&ffi_res).unwrap(),
            )
        },
        | Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<Vec<Vec<Vec<f64>>>>,
                err: Some(e),
            };

            to_c_string(
                serde_json::to_string(&ffi_res).unwrap(),
            )
        },
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_ricci_tensor_json(
    input_json : *const c_char
) -> *mut c_char {

    let input: DgPointInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        Vec<Vec<f64>>,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON input"
                                .to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    match differential_geometry::ricci_tensor(
        input.system,
        &input.point,
    ) {
        | Ok(res) => {

            let ffi_res = FfiResult {
                ok: Some(res),
                err: None::<String>,
            };

            to_c_string(
                serde_json::to_string(&ffi_res).unwrap(),
            )
        },
        | Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<Vec<Vec<f64>>>,
                err: Some(e),
            };

            to_c_string(
                serde_json::to_string(&ffi_res).unwrap(),
            )
        },
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_dg_ricci_scalar_json(
    input_json : *const c_char
) -> *mut c_char {

    let input: DgPointInput =
        match from_json_string(input_json) {
            | Some(i) => i,
            | None => {
                return to_c_string(
                    serde_json::to_string(&FfiResult::<
                        f64,
                        String,
                    > {
                        ok: None,
                        err: Some(
                            "Invalid JSON input"
                                .to_string(),
                        ),
                    })
                    .unwrap(),
                )
            },
        };

    match differential_geometry::ricci_scalar(
        input.system,
        &input.point,
    ) {
        | Ok(res) => {

            let ffi_res = FfiResult {
                ok: Some(res),
                err: None::<String>,
            };

            to_c_string(
                serde_json::to_string(&ffi_res).unwrap(),
            )
        },
        | Err(e) => {

            let ffi_res = FfiResult {
                ok: None::<f64>,
                err: Some(e),
            };

            to_c_string(
                serde_json::to_string(&ffi_res).unwrap(),
            )
        },
    }
}
