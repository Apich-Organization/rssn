use std::ffi::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::vector::curl;
use crate::symbolic::vector::divergence;
use crate::symbolic::vector::gradient;
use crate::symbolic::vector::Vector;

/// Computes the magnitude of a vector.

///

/// Takes a JSON string representing a `Vector`.

/// Returns a JSON string representing the `Expr` of its magnitude.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_vector_magnitude(
    v_json: *const c_char
) -> *mut c_char {

    let v: Option<Vector> =
        from_json_string(v_json);

    if let Some(vector) = v {

        let result = vector.magnitude();

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the dot product of two vectors.

///

/// Takes two JSON strings representing `Vector` objects.

/// Returns a JSON string representing the `Expr` of their dot product.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_vector_dot(
    v1_json: *const c_char,
    v2_json: *const c_char,
) -> *mut c_char {

    let v1: Option<Vector> =
        from_json_string(v1_json);

    let v2: Option<Vector> =
        from_json_string(v2_json);

    match (v1, v2)
    { (Some(vec1), Some(vec2)) => {

        let result = vec1.dot(&vec2);

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the cross product of two vectors.

///

/// Takes two JSON strings representing `Vector` objects.

/// Returns a JSON string representing the `Vector` of their cross product.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_vector_cross(
    v1_json: *const c_char,
    v2_json: *const c_char,
) -> *mut c_char {

    let v1: Option<Vector> =
        from_json_string(v1_json);

    let v2: Option<Vector> =
        from_json_string(v2_json);

    match (v1, v2)
    { (Some(vec1), Some(vec2)) => {

        let result = vec1.cross(&vec2);

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Normalizes a vector.

///

/// Takes a JSON string representing a `Vector`.

/// Returns a JSON string representing the normalized `Vector`.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_vector_normalize(
    v_json: *const c_char
) -> *mut c_char {

    let v: Option<Vector> =
        from_json_string(v_json);

    if let Some(vector) = v {

        let result = vector.normalize();

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the gradient of a scalar field.

///

/// Takes a JSON string representing an `Expr` (scalar field) and three C-style strings for the variable names (x, y, z).

/// Returns a JSON string representing the `Vector` of the gradient.

/// Computes the gradient of a scalar field.

///

/// Takes a JSON string representing an `Expr` (scalar field) and three C-style strings for the variable names (x, y, z).

/// Returns a JSON string representing the `Vector` of the gradient.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_vector_gradient(
    scalar_field_json: *const c_char,
    x_var: *const c_char,
    y_var: *const c_char,
    z_var: *const c_char,
) -> *mut c_char {

    let scalar_field: Option<Expr> =
        from_json_string(
            scalar_field_json,
        );

    let x_str = unsafe {

        if x_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                x_var,
            )
            .to_str()
            .ok()
        }
    };

    let y_str = unsafe {

        if y_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                y_var,
            )
            .to_str()
            .ok()
        }
    };

    let z_str = unsafe {

        if z_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                z_var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        scalar_field,
        x_str,
        y_str,
        z_str,
    ) { (
        Some(field),
        Some(x),
        Some(y),
        Some(z),
    ) => {

        let result =
            gradient(&field, (x, y, z));

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the divergence of a vector field.

///

/// Takes a JSON string representing a `Vector` and three C-style strings for the variable names (x, y, z).

/// Returns a JSON string representing the `Expr` of the divergence.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_vector_divergence(
    v_json: *const c_char,
    x_var: *const c_char,
    y_var: *const c_char,
    z_var: *const c_char,
) -> *mut c_char {

    let v: Option<Vector> =
        from_json_string(v_json);

    let x_str = unsafe {

        if x_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                x_var,
            )
            .to_str()
            .ok()
        }
    };

    let y_str = unsafe {

        if y_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                y_var,
            )
            .to_str()
            .ok()
        }
    };

    let z_str = unsafe {

        if z_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                z_var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        v, x_str, y_str, z_str,
    ) { (
        Some(vector),
        Some(x),
        Some(y),
        Some(z),
    ) => {

        let result = divergence(
            &vector,
            (x, y, z),
        );

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}

/// Computes the curl of a 3D vector field.

///

/// Takes a JSON string representing a `Vector` and three C-style strings for the variable names (x, y, z).

/// Returns a JSON string representing the `Vector` of the curl.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_json_vector_curl(
    v_json: *const c_char,
    x_var: *const c_char,
    y_var: *const c_char,
    z_var: *const c_char,
) -> *mut c_char {

    let v: Option<Vector> =
        from_json_string(v_json);

    let x_str = unsafe {

        if x_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                x_var,
            )
            .to_str()
            .ok()
        }
    };

    let y_str = unsafe {

        if y_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                y_var,
            )
            .to_str()
            .ok()
        }
    };

    let z_str = unsafe {

        if z_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(
                z_var,
            )
            .to_str()
            .ok()
        }
    };

    match (
        v, x_str, y_str, z_str,
    ) { (
        Some(vector),
        Some(x),
        Some(y),
        Some(z),
    ) => {

        let result =
            curl(&vector, (x, y, z));

        to_json_string(&result)
    } _ => {

        std::ptr::null_mut()
    }}
}
