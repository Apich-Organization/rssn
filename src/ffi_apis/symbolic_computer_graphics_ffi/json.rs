//! JSON-based FFI API for computer graphics operations.
//!
//! This module provides JSON string-based FFI functions for 2D/3D transformations
//! and projections, enabling language-agnostic integration.

use crate::symbolic::core::Expr;
use crate::symbolic::computer_graphics::{
    translation_2d, translation_3d, rotation_2d, rotation_3d_x, rotation_3d_y, rotation_3d_z,
    scaling_2d, scaling_3d,
};
use crate::ffi_apis::common::*;
use std::os::raw::c_char;

/// Generates a 3x3 2D translation matrix via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_translation_2d(
    tx_json: *const c_char,
    ty_json: *const c_char
) -> *mut c_char {
    let tx: Option<Expr> = from_json_string(tx_json);
    let ty: Option<Expr> = from_json_string(ty_json);
    if let (Some(tx), Some(ty)) = (tx, ty) {
        to_json_string(&translation_2d(tx, ty))
    } else {
        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D translation matrix via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_translation_3d(
    tx_json: *const c_char,
    ty_json: *const c_char,
    tz_json: *const c_char
) -> *mut c_char {
    let tx: Option<Expr> = from_json_string(tx_json);
    let ty: Option<Expr> = from_json_string(ty_json);
    let tz: Option<Expr> = from_json_string(tz_json);
    if let (Some(tx), Some(ty), Some(tz)) = (tx, ty, tz) {
        to_json_string(&translation_3d(tx, ty, tz))
    } else {
        std::ptr::null_mut()
    }
}

/// Generates a 3x3 2D rotation matrix via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_rotation_2d(angle_json: *const c_char) -> *mut c_char {
    let angle: Option<Expr> = from_json_string(angle_json);
    if let Some(a) = angle {
        to_json_string(&rotation_2d(a))
    } else {
        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D rotation matrix around X-axis via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_rotation_3d_x(angle_json: *const c_char) -> *mut c_char {
    let angle: Option<Expr> = from_json_string(angle_json);
    if let Some(a) = angle {
        to_json_string(&rotation_3d_x(a))
    } else {
        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D rotation matrix around Y-axis via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_rotation_3d_y(angle_json: *const c_char) -> *mut c_char {
    let angle: Option<Expr> = from_json_string(angle_json);
    if let Some(a) = angle {
        to_json_string(&rotation_3d_y(a))
    } else {
        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D rotation matrix around Z-axis via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_rotation_3d_z(angle_json: *const c_char) -> *mut c_char {
    let angle: Option<Expr> = from_json_string(angle_json);
    if let Some(a) = angle {
        to_json_string(&rotation_3d_z(a))
    } else {
        std::ptr::null_mut()
    }
}

/// Generates a 3x3 2D scaling matrix via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_scaling_2d(
    sx_json: *const c_char,
    sy_json: *const c_char
) -> *mut c_char {
    let sx: Option<Expr> = from_json_string(sx_json);
    let sy: Option<Expr> = from_json_string(sy_json);
    if let (Some(sx), Some(sy)) = (sx, sy) {
        to_json_string(&scaling_2d(sx, sy))
    } else {
        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D scaling matrix via JSON interface.
#[no_mangle]
pub unsafe extern "C" fn rssn_json_scaling_3d(
    sx_json: *const c_char,
    sy_json: *const c_char,
    sz_json: *const c_char
) -> *mut c_char {
    let sx: Option<Expr> = from_json_string(sx_json);
    let sy: Option<Expr> = from_json_string(sy_json);
    let sz: Option<Expr> = from_json_string(sz_json);
    if let (Some(sx), Some(sy), Some(sz)) = (sx, sy, sz) {
        to_json_string(&scaling_3d(sx, sy, sz))
    } else {
        std::ptr::null_mut()
    }
}
