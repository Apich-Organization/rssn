use std::ffi::c_char;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::vector::{Vector, gradient, divergence, curl};

/// Computes the magnitude of a vector.

///

/// Takes a JSON string representing a `Vector`.

/// Returns a JSON string representing the `Expr` of its magnitude.

#[no_mangle]

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

#[no_mangle]

pub extern "C" fn rssn_json_vector_dot(
    v1_json: *const c_char,
    v2_json: *const c_char,
) -> *mut c_char {

    let v1: Option<Vector> =
        from_json_string(v1_json);

    let v2: Option<Vector> =
        from_json_string(v2_json);

    if let (Some(vec1), Some(vec2)) =
        (v1, v2)
    {

        let result = vec1.dot(&vec2);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the cross product of two vectors.

///

/// Takes two JSON strings representing `Vector` objects.

/// Returns a JSON string representing the `Vector` of their cross product.

#[no_mangle]

pub extern "C" fn rssn_json_vector_cross(
    v1_json: *const c_char,
    v2_json: *const c_char,
) -> *mut c_char {

    let v1: Option<Vector> =
        from_json_string(v1_json);

    let v2: Option<Vector> =
        from_json_string(v2_json);

    if let (Some(vec1), Some(vec2)) =
        (v1, v2)
    {

        let result = vec1.cross(&vec2);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Normalizes a vector.

///

/// Takes a JSON string representing a `Vector`.

/// Returns a JSON string representing the normalized `Vector`.

#[no_mangle]

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

#[no_mangle]

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

    if let (
        Some(field),
        Some(x),
        Some(y),
        Some(z),
    ) = (
        scalar_field,
        x_str,
        y_str,
        z_str,
    ) {

        let result =
            gradient(&field, (x, y, z));

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the divergence of a vector field.

///

/// Takes a JSON string representing a `Vector` and three C-style strings for the variable names (x, y, z).

/// Returns a JSON string representing the `Expr` of the divergence.

#[no_mangle]

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

    if let (
        Some(vector),
        Some(x),
        Some(y),
        Some(z),
    ) = (
        v, x_str, y_str, z_str,
    ) {

        let result = divergence(
            &vector,
            (x, y, z),
        );

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the curl of a 3D vector field.

///

/// Takes a JSON string representing a `Vector` and three C-style strings for the variable names (x, y, z).

/// Returns a JSON string representing the `Vector` of the curl.

#[no_mangle]

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

    if let (
        Some(vector),
        Some(x),
        Some(y),
        Some(z),
    ) = (
        v, x_str, y_str, z_str,
    ) {

        let result =
            curl(&vector, (x, y, z));

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
