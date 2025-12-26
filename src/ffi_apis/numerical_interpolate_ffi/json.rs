//! JSON-based FFI API for numerical interpolation.

use crate::ffi_apis::common::{
    from_json_string,
    to_c_string,
};
use crate::ffi_apis::ffi_api::FfiResult;
use crate::numerical::interpolate;
use crate::numerical::polynomial::Polynomial;
use serde::{
    Deserialize,
    Serialize,
};
use std::os::raw::c_char;

#[derive(Deserialize)]

struct LagrangeInput {
    points: Vec<(f64, f64)>,
}

#[derive(Deserialize)]

struct CubicSplineInput {
    points: Vec<(f64, f64)>,
    x_eval: f64,
}

#[derive(Deserialize)]

struct BezierInput {
    control_points: Vec<Vec<f64>>,
    t: f64,
}

#[derive(Deserialize)]

struct BSplineInput {
    control_points: Vec<Vec<f64>>,
    degree: usize,
    knots: Vec<f64>,
    t: f64,
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_lagrange_interpolation_json(
    input_ptr: *const c_char
) -> *mut c_char {

    let input: LagrangeInput = match from_json_string(input_ptr) {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = interpolate::lagrange_interpolation(&input.points);

    let ffi_result = match result {
        Ok(poly) => {
            FfiResult {
                ok: Some(poly),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    to_c_string(serde_json::to_string(&ffi_result).unwrap())
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_cubic_spline_interpolation_json(
    input_ptr: *const c_char
) -> *mut c_char {

    let input: CubicSplineInput = match from_json_string(input_ptr) {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = interpolate::cubic_spline_interpolation(&input.points);

    let ffi_result = match result {
        Ok(spline) => {

            let val = spline(input.x_eval);

            FfiResult {
                ok: Some(val),
                err: None,
            }
        }
        Err(e) => {
            FfiResult {
                ok: None,
                err: Some(e),
            }
        }
    };

    to_c_string(serde_json::to_string(&ffi_result).unwrap())
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_bezier_curve_json(input_ptr: *const c_char) -> *mut c_char {

    let input: BezierInput = match from_json_string(input_ptr) {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = interpolate::bezier_curve(
        &input.control_points,
        input.t,
    );

    let ffi_result = FfiResult::<Vec<f64>, String> {
        ok: Some(result),
        err: None,
    };

    to_c_string(serde_json::to_string(&ffi_result).unwrap())
}

#[no_mangle]

pub unsafe extern "C" fn rssn_num_b_spline_json(input_ptr: *const c_char) -> *mut c_char {

    let input: BSplineInput = match from_json_string(input_ptr) {
        Some(i) => i,
        None => return std::ptr::null_mut(),
    };

    let result = interpolate::b_spline(
        &input.control_points,
        input.degree,
        &input.knots,
        input.t,
    );

    let ffi_result = match result {
        Some(p) => {
            FfiResult {
                ok: Some(p),
                err: None,
            }
        }
        None => {
            FfiResult {
                ok: None,
                err: Some("Invalid B-spline parameters".to_string()),
            }
        }
    };

    to_c_string(serde_json::to_string(&ffi_result).unwrap())
}
