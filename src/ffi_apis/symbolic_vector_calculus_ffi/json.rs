//! JSON-based FFI API for symbolic vector calculus functions.

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;
use crate::symbolic::vector_calculus::*;
use serde::{Deserialize, Serialize};
use std::os::raw::c_char;

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

/// Computes the line integral of a scalar field (JSON).
#[no_mangle]

pub extern "C" fn rssn_line_integral_scalar_json(input_json: *const c_char) -> *mut c_char {

    let input: Option<LineIntegralScalarInput> = from_json_string(input_json);

    let input = match input {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = line_integral_scalar(&input.scalar_field, &input.curve);

    to_json_string(&result)
}

/// Computes the line integral of a vector field (JSON).
#[no_mangle]

pub extern "C" fn rssn_line_integral_vector_json(input_json: *const c_char) -> *mut c_char {

    let input: Option<LineIntegralVectorInput> = from_json_string(input_json);

    let input = match input {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = line_integral_vector(&input.vector_field, &input.curve);

    to_json_string(&result)
}

/// Computes the surface integral (flux) of a vector field (JSON).
#[no_mangle]

pub extern "C" fn rssn_surface_integral_json(input_json: *const c_char) -> *mut c_char {

    let input: Option<SurfaceIntegralInput> = from_json_string(input_json);

    let input = match input {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = surface_integral(&input.vector_field, &input.surface);

    to_json_string(&result)
}

/// Computes the volume integral of a scalar field (JSON).
#[no_mangle]

pub extern "C" fn rssn_volume_integral_json(input_json: *const c_char) -> *mut c_char {

    let input: Option<VolumeIntegralInput> = from_json_string(input_json);

    let input = match input {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = volume_integral(&input.scalar_field, &input.volume);

    to_json_string(&result)
}
