//! Bincode-based FFI API for numerical coordinate transformations.

use serde::Deserialize;
use serde::Serialize;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::coordinates as nc;
use crate::symbolic::coordinates::CoordinateSystem;

#[derive(Deserialize)]

struct CoordinateTransformRequest {
    point: Vec<f64>,
    from: CoordinateSystem,
    to: CoordinateSystem,
}

fn decode<
    T: for<'de> Deserialize<'de>,
>(
    data: *const u8,
    len: usize,
) -> Option<T> {

    if data.is_null() {

        return None;
    }

    let slice = unsafe {

        std::slice::from_raw_parts(
            data, len,
        )
    };

    bincode_next::serde::decode_from_slice(
        slice,
        bincode_next::config::standard(),
    )
    .ok()
    .map(|(v, _)| v)
}

fn encode<T: Serialize>(
    val: &T
) -> BincodeBuffer {

    match bincode_next::serde::encode_to_vec(
        val,
        bincode_next::config::standard(),
    ) {
        | Ok(bytes) => BincodeBuffer::from_vec(bytes),
        | Err(_) => BincodeBuffer::empty(),
    }
}

/// Transforms a point via Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_coord_transform_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: CoordinateTransformRequest =
        match decode(data, len) {
            | Some(r) => r,
            | None => {
                return encode(&FfiResult::<
                    Vec<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Bincode decode error".to_string(),
                    ),
                })
            },
        };

    match nc::transform_point(
        &req.point,
        req.from,
        req.to,
    ) {
        | Ok(res) => {
            encode(&FfiResult::<
                Vec<f64>,
                String,
            > {
                ok: Some(res),
                err: None,
            })
        },
        | Err(e) => {
            encode(&FfiResult::<
                Vec<f64>,
                String,
            > {
                ok: None,
                err: Some(e),
            })
        },
    }
}

/// Transforms a point (pure numerical) via Bincode.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_coord_transform_pure_bincode(
    data: *const u8,
    len: usize,
) -> BincodeBuffer {

    let req: CoordinateTransformRequest =
        match decode(data, len) {
            | Some(r) => r,
            | None => {
                return encode(&FfiResult::<
                    Vec<f64>,
                    String,
                > {
                    ok: None,
                    err: Some(
                        "Bincode decode error".to_string(),
                    ),
                })
            },
        };

    match nc::transform_point_pure(
        &req.point,
        req.from,
        req.to,
    ) {
        | Ok(res) => {
            encode(&FfiResult::<
                Vec<f64>,
                String,
            > {
                ok: Some(res),
                err: None,
            })
        },
        | Err(e) => {
            encode(&FfiResult::<
                Vec<f64>,
                String,
            > {
                ok: None,
                err: Some(e),
            })
        },
    }
}
