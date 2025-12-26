//! Bincode-based FFI API for symbolic vector calculus functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;
use crate::symbolic::vector_calculus::*;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]

struct LineIntegralScalarInput {
    scalar_field: Expr,
    curve: ParametricCurve,
}

#[derive(Serialize, Deserialize)]

struct LineIntegralVectorInput {
    vector_field: Vector,
    curve: ParametricCurve,
}

#[derive(Serialize, Deserialize)]

struct SurfaceIntegralInput {
    vector_field: Vector,
    surface: ParametricSurface,
}

#[derive(Serialize, Deserialize)]

struct VolumeIntegralInput {
    scalar_field: Expr,
    volume: Volume,
}

/// Computes the line integral of a scalar field (Bincode).
#[no_mangle]

pub extern "C" fn rssn_line_integral_scalar_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr as *mut u8,
        len: input_len,
    };

    let input: Option<LineIntegralScalarInput> = from_bincode_buffer(&input_buffer);

    let input = match input {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = line_integral_scalar(&input.scalar_field, &input.curve);

    to_bincode_buffer(&result)
}

/// Computes the line integral of a vector field (Bincode).
#[no_mangle]

pub extern "C" fn rssn_line_integral_vector_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr as *mut u8,
        len: input_len,
    };

    let input: Option<LineIntegralVectorInput> = from_bincode_buffer(&input_buffer);

    let input = match input {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = line_integral_vector(&input.vector_field, &input.curve);

    to_bincode_buffer(&result)
}

/// Computes the surface integral (flux) of a vector field (Bincode).
#[no_mangle]

pub extern "C" fn rssn_surface_integral_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr as *mut u8,
        len: input_len,
    };

    let input: Option<SurfaceIntegralInput> = from_bincode_buffer(&input_buffer);

    let input = match input {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = surface_integral(&input.vector_field, &input.surface);

    to_bincode_buffer(&result)
}

/// Computes the volume integral of a scalar field (Bincode).
#[no_mangle]

pub extern "C" fn rssn_volume_integral_bincode(
    input_ptr: *const u8,
    input_len: usize,
) -> BincodeBuffer {

    let input_buffer = BincodeBuffer {
        data: input_ptr as *mut u8,
        len: input_len,
    };

    let input: Option<VolumeIntegralInput> = from_bincode_buffer(&input_buffer);

    let input = match input {
        Some(i) => i,
        None => return BincodeBuffer::empty(),
    };

    let result = volume_integral(&input.scalar_field, &input.volume);

    to_bincode_buffer(&result)
}
