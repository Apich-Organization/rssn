//! Bincode-based FFI API for numerical computer graphics functions.

use crate::ffi_apis::common::{from_bincode_buffer, to_bincode_buffer, BincodeBuffer};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::computer_graphics;
use serde::{Deserialize, Serialize};

#[derive(Deserialize)]

struct Vector3DInput {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Deserialize)]

struct TwoVectors3DInput {
    v1: Vector3DInput,
    v2: Vector3DInput,
}

#[derive(Serialize)]

struct Vector3DOutput {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Deserialize)]

struct AngleInput {
    angle: f64,
}

#[derive(Deserialize)]

struct TransformInput {
    dx: f64,
    dy: f64,
    dz: f64,
}

#[derive(Deserialize)]

struct QuaternionInput {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Serialize)]

struct QuaternionOutput {
    w: f64,
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Deserialize)]

struct TwoQuaternionsInput {
    q1: QuaternionInput,
    q2: QuaternionInput,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_dot_product_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: TwoVectors3DInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<f64, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let v1 = computer_graphics::Vector3D::new(input.v1.x, input.v1.y, input.v1.z);

    let v2 = computer_graphics::Vector3D::new(input.v2.x, input.v2.y, input.v2.z);

    let result = computer_graphics::dot_product(&v1, &v2);

    to_bincode_buffer(&FfiResult {
        ok: Some(result),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_cross_product_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: TwoVectors3DInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vector3DOutput, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let v1 = computer_graphics::Vector3D::new(input.v1.x, input.v1.y, input.v1.z);

    let v2 = computer_graphics::Vector3D::new(input.v2.x, input.v2.y, input.v2.z);

    let result = computer_graphics::cross_product(&v1, &v2);

    to_bincode_buffer(&FfiResult {
        ok: Some(Vector3DOutput {
            x: result.x,
            y: result.y,
            z: result.z,
        }),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_normalize_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: Vector3DInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vector3DOutput, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let v = computer_graphics::Vector3D::new(input.x, input.y, input.z);

    let result = v.normalize();

    to_bincode_buffer(&FfiResult {
        ok: Some(Vector3DOutput {
            x: result.x,
            y: result.y,
            z: result.z,
        }),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_rotation_matrix_x_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: AngleInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<f64>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let matrix = computer_graphics::rotation_matrix_x(input.angle);

    to_bincode_buffer(&FfiResult {
        ok: Some(matrix.data()),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_translation_matrix_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: TransformInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<Vec<f64>, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let matrix = computer_graphics::translation_matrix(input.dx, input.dy, input.dz);

    to_bincode_buffer(&FfiResult {
        ok: Some(matrix.data()),
        err: None::<String>,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_graphics_quaternion_multiply_bincode(
    buffer: BincodeBuffer,
) -> BincodeBuffer {

    let input: TwoQuaternionsInput = match from_bincode_buffer(&buffer) {
        Some(i) => i,
        None => {
            return to_bincode_buffer(&FfiResult::<QuaternionOutput, String> {
                ok: None,
                err: Some("Invalid Bincode".to_string()),
            })
        }
    };

    let q1 = computer_graphics::Quaternion::new(input.q1.w, input.q1.x, input.q1.y, input.q1.z);

    let q2 = computer_graphics::Quaternion::new(input.q2.w, input.q2.x, input.q2.y, input.q2.z);

    let result = q1.multiply(&q2);

    to_bincode_buffer(&FfiResult {
        ok: Some(QuaternionOutput {
            w: result.w,
            x: result.x,
            y: result.y,
            z: result.z,
        }),
        err: None::<String>,
    })
}
