use std::ffi::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::coordinates::*;
use crate::symbolic::core::Expr;

/// Transforms a point between coordinate systems using JSON-encoded coordinates.
///
/// The point is represented as a JSON-encoded `Vec<Expr>` (e.g., \([x,y,z]\)), and the
/// transformation applies the appropriate coordinate mapping.
///
/// # Arguments
///
/// * `point_json` - C string pointer to JSON encoding a `Vec<Expr>` for the point in the `from` system.
/// * `from` - Source [`CoordinateSystem`].
/// * `to` - Target [`CoordinateSystem`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Vec<Expr>` for the point in the `to` system,
/// or null on deserialization or transformation failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and returns
/// ownership of a heap-allocated C string.
#[no_mangle]

pub extern "C" fn rssn_json_transform_point(
    point_json: *const c_char,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut c_char {

    let point: Option<Vec<Expr>> =
        from_json_string(point_json);

    if let Some(p) = point {

        match transform_point(
            &p, from, to,
        ) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Transforms a scalar expression between coordinate systems using JSON serialization.
///
/// # Arguments
///
/// * `expr_json` - C string pointer to JSON encoding an `Expr` in the `from` system.
/// * `from` - Source [`CoordinateSystem`].
/// * `to` - Target [`CoordinateSystem`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the transformed expression in
/// the `to` system, or null on deserialization or transformation failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and returns
/// ownership of a heap-allocated C string.
#[no_mangle]

pub extern "C" fn rssn_json_transform_expression(
    expr_json: *const c_char,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut c_char {

    let expr: Option<Expr> =
        from_json_string(expr_json);

    if let Some(e) = expr {

        match transform_expression(
            &e, from, to,
        ) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Returns the metric tensor for a given coordinate system as JSON-encoded `Expr`.
///
/// # Arguments
///
/// * `system` - [`CoordinateSystem`] for which to compute the metric tensor.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the metric tensor, or null on error.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point and returns
/// ownership of a heap-allocated C string.
#[no_mangle]

pub extern "C" fn rssn_json_coordinates_get_metric_tensor(
    system: CoordinateSystem
) -> *mut c_char {

    match get_metric_tensor(system) {
        | Ok(result) => {
            to_json_string(&result)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Transforms contravariant vector components between coordinate systems using JSON serialization.
///
/// # Arguments
///
/// * `comps_json` - C string pointer to JSON encoding a `Vec<Expr>` of contravariant components
///   in the `from` system.
/// * `from` - Source [`CoordinateSystem`].
/// * `to` - Target [`CoordinateSystem`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Vec<Expr>` with components in the `to` system,
/// or null on deserialization or transformation failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and returns
/// ownership of a heap-allocated C string.
#[no_mangle]

pub extern "C" fn rssn_json_transform_contravariant_vector(
    comps_json: *const c_char,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut c_char {

    let comps: Option<Vec<Expr>> =
        from_json_string(comps_json);

    if let Some(c) = comps {

        match transform_contravariant_vector(&c, from, to) {
            | Ok(result) => to_json_string(&result),
            | Err(_) => std::ptr::null_mut(),
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Transforms covariant vector components between coordinate systems using JSON serialization.
///
/// # Arguments
///
/// * `comps_json` - C string pointer to JSON encoding a `Vec<Expr>` of covariant components
///   in the `from` system.
/// * `from` - Source [`CoordinateSystem`].
/// * `to` - Target [`CoordinateSystem`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Vec<Expr>` with components in the `to` system,
/// or null on deserialization or transformation failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and returns
/// ownership of a heap-allocated C string.
#[no_mangle]

pub extern "C" fn rssn_json_transform_covariant_vector(
    comps_json: *const c_char,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut c_char {

    let comps: Option<Vec<Expr>> =
        from_json_string(comps_json);

    if let Some(c) = comps {

        match transform_covariant_vector(
            &c, from, to,
        ) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the divergence of a vector field in a given coordinate system using JSON serialization.
///
/// # Arguments
///
/// * `comps_json` - C string pointer to JSON encoding a `Vec<Expr>` of vector components.
/// * `from` - [`CoordinateSystem`] with respect to which the divergence is taken.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the divergence, or null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and returns
/// ownership of a heap-allocated C string.
#[no_mangle]

pub extern "C" fn rssn_json_transform_divergence(
    comps_json: *const c_char,
    from: CoordinateSystem,
) -> *mut c_char {

    let comps: Option<Vec<Expr>> =
        from_json_string(comps_json);

    if let Some(c) = comps {

        match transform_divergence(
            &c, from,
        ) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

/// Computes the curl of a vector field in a given coordinate system using JSON serialization.
///
/// # Arguments
///
/// * `comps_json` - C string pointer to JSON encoding a `Vec<Expr>` of vector components.
/// * `from` - [`CoordinateSystem`] with respect to which the curl is taken.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Vec<Expr>` for the curl, or null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and returns
/// ownership of a heap-allocated C string.
#[no_mangle]

pub extern "C" fn rssn_json_transform_curl(
    comps_json: *const c_char,
    from: CoordinateSystem,
) -> *mut c_char {

    let comps: Option<Vec<Expr>> =
        from_json_string(comps_json);

    if let Some(c) = comps {

        match transform_curl(&c, from) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

/// Computes the gradient of a scalar field and transforms it between coordinate systems
/// using JSON serialization.
///
/// # Arguments
///
/// * `scalar_json` - C string pointer to JSON encoding an `Expr` for the scalar field.
/// * `vars_json` - C string pointer to JSON encoding a `Vec<String>` of coordinate variables.
/// * `from` - Source [`CoordinateSystem`].
/// * `to` - Target [`CoordinateSystem`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Vec<Expr>` for the gradient components in
/// the `to` system, or null on deserialization or transformation failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and returns
/// ownership of a heap-allocated C string.
#[no_mangle]

pub extern "C" fn rssn_json_transform_gradient(
    scalar_json: *const c_char,
    vars_json: *const c_char,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut c_char {

    let scalar: Option<Expr> =
        from_json_string(scalar_json);

    let vars: Option<Vec<String>> =
        from_json_string(vars_json);

    if let (Some(s), Some(v)) =
        (scalar, vars)
    {

        match transform_gradient(
            &s, &v, from, to,
        ) {
            | Ok(result) => {
                to_json_string(&result)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}
