//! Handle-based FFI API for symbolic vector calculus functions.

use crate::ffi_apis::common::*;
use crate::input::parser::parse_expr;
use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;
use crate::symbolic::vector_calculus::*;
use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::sync::Arc;

// Helper function to parse expression from string
fn parse_expr_from_str(s: &str) -> Option<Expr> {
    match parse_expr(s) {
        Ok(("", expr)) => Some(expr),
        _ => None,
    }
}

// ===== ParametricCurve Handle Functions =====

/// Creates a new ParametricCurve.
#[no_mangle]
pub extern "C" fn rssn_parametric_curve_new(
    r_x: *const c_char,
    r_y: *const c_char,
    r_z: *const c_char,
    t_var: *const c_char,
    t_lower: *const c_char,
    t_upper: *const c_char,
) -> *mut ParametricCurve {
    if r_x.is_null() || r_y.is_null() || r_z.is_null() || t_var.is_null() || t_lower.is_null() || t_upper.is_null() {
        return std::ptr::null_mut();
    }

    unsafe {
        let r_x_str = match CStr::from_ptr(r_x).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let r_y_str = match CStr::from_ptr(r_y).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let r_z_str = match CStr::from_ptr(r_z).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let t_var_str = match CStr::from_ptr(t_var).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let t_lower_str = match CStr::from_ptr(t_lower).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let t_upper_str = match CStr::from_ptr(t_upper).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };

        let r_x_expr: Expr = match parse_expr_from_str(r_x_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let r_y_expr: Expr = match parse_expr_from_str(r_y_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let r_z_expr: Expr = match parse_expr_from_str(r_z_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let t_lower_expr: Expr = match parse_expr_from_str(t_lower_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let t_upper_expr: Expr = match parse_expr_from_str(t_upper_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };

        let curve = ParametricCurve {
            r: Vector::new(r_x_expr, r_y_expr, r_z_expr),
            t_var: t_var_str.to_string(),
            t_bounds: (t_lower_expr, t_upper_expr),
        };

        Box::into_raw(Box::new(curve))
    }
}

/// Frees a ParametricCurve handle.
#[no_mangle]
pub extern "C" fn rssn_parametric_curve_free(curve: *mut ParametricCurve) {
    if !curve.is_null() {
        unsafe {
            let _ = Box::from_raw(curve);
        }
    }
}

// ===== ParametricSurface Handle Functions =====

/// Creates a new ParametricSurface.
#[no_mangle]
pub extern "C" fn rssn_parametric_surface_new(
    r_x: *const c_char,
    r_y: *const c_char,
    r_z: *const c_char,
    u_var: *const c_char,
    u_lower: *const c_char,
    u_upper: *const c_char,
    v_var: *const c_char,
    v_lower: *const c_char,
    v_upper: *const c_char,
) -> *mut ParametricSurface {
    if r_x.is_null() || r_y.is_null() || r_z.is_null() || u_var.is_null() || u_lower.is_null() 
        || u_upper.is_null() || v_var.is_null() || v_lower.is_null() || v_upper.is_null() {
        return std::ptr::null_mut();
    }

    unsafe {
        let r_x_str = match CStr::from_ptr(r_x).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let r_y_str = match CStr::from_ptr(r_y).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let r_z_str = match CStr::from_ptr(r_z).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let u_var_str = match CStr::from_ptr(u_var).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let u_lower_str = match CStr::from_ptr(u_lower).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let u_upper_str = match CStr::from_ptr(u_upper).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let v_var_str = match CStr::from_ptr(v_var).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let v_lower_str = match CStr::from_ptr(v_lower).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let v_upper_str = match CStr::from_ptr(v_upper).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };

        let r_x_expr: Expr = match parse_expr_from_str(r_x_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let r_y_expr: Expr = match parse_expr_from_str(r_y_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let r_z_expr: Expr = match parse_expr_from_str(r_z_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let u_lower_expr: Expr = match parse_expr_from_str(u_lower_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let u_upper_expr: Expr = match parse_expr_from_str(u_upper_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let v_lower_expr: Expr = match parse_expr_from_str(v_lower_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let v_upper_expr: Expr = match parse_expr_from_str(v_upper_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };

        let surface = ParametricSurface {
            r: Vector::new(r_x_expr, r_y_expr, r_z_expr),
            u_var: u_var_str.to_string(),
            u_bounds: (u_lower_expr, u_upper_expr),
            v_var: v_var_str.to_string(),
            v_bounds: (v_lower_expr, v_upper_expr),
        };

        Box::into_raw(Box::new(surface))
    }
}

/// Frees a ParametricSurface handle.
#[no_mangle]
pub extern "C" fn rssn_parametric_surface_free(surface: *mut ParametricSurface) {
    if !surface.is_null() {
        unsafe {
            let _ = Box::from_raw(surface);
        }
    }
}

// ===== Volume Handle Functions =====

/// Creates a new Volume.
#[no_mangle]
pub extern "C" fn rssn_volume_new(
    z_lower: *const c_char,
    z_upper: *const c_char,
    y_lower: *const c_char,
    y_upper: *const c_char,
    x_lower: *const c_char,
    x_upper: *const c_char,
    x_var: *const c_char,
    y_var: *const c_char,
    z_var: *const c_char,
) -> *mut Volume {
    if z_lower.is_null() || z_upper.is_null() || y_lower.is_null() || y_upper.is_null() 
        || x_lower.is_null() || x_upper.is_null() || x_var.is_null() || y_var.is_null() || z_var.is_null() {
        return std::ptr::null_mut();
    }

    unsafe {
        let z_lower_str = match CStr::from_ptr(z_lower).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let z_upper_str = match CStr::from_ptr(z_upper).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let y_lower_str = match CStr::from_ptr(y_lower).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let y_upper_str = match CStr::from_ptr(y_upper).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let x_lower_str = match CStr::from_ptr(x_lower).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let x_upper_str = match CStr::from_ptr(x_upper).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let x_var_str = match CStr::from_ptr(x_var).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let y_var_str = match CStr::from_ptr(y_var).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let z_var_str = match CStr::from_ptr(z_var).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };

        let z_lower_expr: Expr = match parse_expr_from_str(z_lower_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let z_upper_expr: Expr = match parse_expr_from_str(z_upper_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let y_lower_expr: Expr = match parse_expr_from_str(y_lower_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let y_upper_expr: Expr = match parse_expr_from_str(y_upper_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let x_lower_expr: Expr = match parse_expr_from_str(x_lower_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let x_upper_expr: Expr = match parse_expr_from_str(x_upper_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };

        let volume = Volume {
            z_bounds: (z_lower_expr, z_upper_expr),
            y_bounds: (y_lower_expr, y_upper_expr),
            x_bounds: (x_lower_expr, x_upper_expr),
            vars: (x_var_str.to_string(), y_var_str.to_string(), z_var_str.to_string()),
        };

        Box::into_raw(Box::new(volume))
    }
}

/// Frees a Volume handle.
#[no_mangle]
pub extern "C" fn rssn_volume_free(volume: *mut Volume) {
    if !volume.is_null() {
        unsafe {
            let _ = Box::from_raw(volume);
        }
    }
}

// ===== Vector Calculus Operations =====

/// Computes the line integral of a scalar field along a curve.
#[no_mangle]
pub extern "C" fn rssn_line_integral_scalar(
    scalar_field: *const c_char,
    curve: *const ParametricCurve,
) -> *mut c_char {
    if scalar_field.is_null() || curve.is_null() {
        return std::ptr::null_mut();
    }

    unsafe {
        let field_str = match CStr::from_ptr(scalar_field).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };

        let field_expr: Expr = match parse_expr_from_str(field_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };

        let curve_ref = &*curve;
        let result = line_integral_scalar(&field_expr, curve_ref);
        let result_str = format!("{}", result);

        match CString::new(result_str) {
            Ok(c_str) => c_str.into_raw(),
            Err(_) => std::ptr::null_mut(),
        }
    }
}

/// Computes the line integral of a vector field along a curve.
#[no_mangle]
pub extern "C" fn rssn_line_integral_vector(
    field_x: *const c_char,
    field_y: *const c_char,
    field_z: *const c_char,
    curve: *const ParametricCurve,
) -> *mut c_char {
    if field_x.is_null() || field_y.is_null() || field_z.is_null() || curve.is_null() {
        return std::ptr::null_mut();
    }

    unsafe {
        let field_x_str = match CStr::from_ptr(field_x).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let field_y_str = match CStr::from_ptr(field_y).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let field_z_str = match CStr::from_ptr(field_z).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };

        let field_x_expr: Expr = match parse_expr_from_str(field_x_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let field_y_expr: Expr = match parse_expr_from_str(field_y_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let field_z_expr: Expr = match parse_expr_from_str(field_z_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };

        let field = Vector::new(field_x_expr, field_y_expr, field_z_expr);
        let curve_ref = &*curve;
        let result = line_integral_vector(&field, curve_ref);
        let result_str = format!("{}", result);

        match CString::new(result_str) {
            Ok(c_str) => c_str.into_raw(),
            Err(_) => std::ptr::null_mut(),
        }
    }
}

/// Computes the surface integral (flux) of a vector field.
#[no_mangle]
pub extern "C" fn rssn_surface_integral(
    field_x: *const c_char,
    field_y: *const c_char,
    field_z: *const c_char,
    surface: *const ParametricSurface,
) -> *mut c_char {
    if field_x.is_null() || field_y.is_null() || field_z.is_null() || surface.is_null() {
        return std::ptr::null_mut();
    }

    unsafe {
        let field_x_str = match CStr::from_ptr(field_x).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let field_y_str = match CStr::from_ptr(field_y).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };
        let field_z_str = match CStr::from_ptr(field_z).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };

        let field_x_expr: Expr = match parse_expr_from_str(field_x_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let field_y_expr: Expr = match parse_expr_from_str(field_y_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };
        let field_z_expr: Expr = match parse_expr_from_str(field_z_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };

        let field = Vector::new(field_x_expr, field_y_expr, field_z_expr);
        let surface_ref = &*surface;
        let result = surface_integral(&field, surface_ref);
        let result_str = format!("{}", result);

        match CString::new(result_str) {
            Ok(c_str) => c_str.into_raw(),
            Err(_) => std::ptr::null_mut(),
        }
    }
}

/// Computes the volume integral of a scalar field.
#[no_mangle]
pub extern "C" fn rssn_volume_integral(
    scalar_field: *const c_char,
    volume: *const Volume,
) -> *mut c_char {
    if scalar_field.is_null() || volume.is_null() {
        return std::ptr::null_mut();
    }

    unsafe {
        let field_str = match CStr::from_ptr(scalar_field).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        };
            Ok(s) => s,
            None => return std::ptr::null_mut(),
        };

        let field_expr: Expr = match parse_expr_from_str(field_str) {
            Some(e) => e,
            None => return std::ptr::null_mut(),
        };

        let volume_ref = &*volume;
        let result = volume_integral(&field_expr, volume_ref);
        let result_str = format!("{}", result);

        match CString::new(result_str) {
            Ok(c_str) => c_str.into_raw(),
            Err(_) => std::ptr::null_mut(),
        }
    }
}
