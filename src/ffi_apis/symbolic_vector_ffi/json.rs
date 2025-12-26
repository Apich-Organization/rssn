use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::vector::*;
use std::ffi::c_char;

#[no_mangle]

pub extern "C" fn rssn_json_vector_magnitude(v_json: *const c_char) -> *mut c_char {

    let v: Option<Vector> = from_json_string(v_json);

    if let Some(vector) = v {

        let result = vector.magnitude();

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub extern "C" fn rssn_json_vector_dot(
    v1_json: *const c_char,
    v2_json: *const c_char,
) -> *mut c_char {

    let v1: Option<Vector> = from_json_string(v1_json);

    let v2: Option<Vector> = from_json_string(v2_json);

    if let (Some(vec1), Some(vec2)) = (v1, v2) {

        let result = vec1.dot(&vec2);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub extern "C" fn rssn_json_vector_cross(
    v1_json: *const c_char,
    v2_json: *const c_char,
) -> *mut c_char {

    let v1: Option<Vector> = from_json_string(v1_json);

    let v2: Option<Vector> = from_json_string(v2_json);

    if let (Some(vec1), Some(vec2)) = (v1, v2) {

        let result = vec1.cross(&vec2);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub extern "C" fn rssn_json_vector_normalize(v_json: *const c_char) -> *mut c_char {

    let v: Option<Vector> = from_json_string(v_json);

    if let Some(vector) = v {

        let result = vector.normalize();

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub extern "C" fn rssn_json_vector_gradient(
    scalar_field_json: *const c_char,
    x_var: *const c_char,
    y_var: *const c_char,
    z_var: *const c_char,
) -> *mut c_char {

    let scalar_field: Option<Expr> = from_json_string(scalar_field_json);

    let x_str = unsafe {

        if x_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(x_var)
                .to_str()
                .ok()
        }
    };

    let y_str = unsafe {

        if y_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(y_var)
                .to_str()
                .ok()
        }
    };

    let z_str = unsafe {

        if z_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(z_var)
                .to_str()
                .ok()
        }
    };

    if let (Some(field), Some(x), Some(y), Some(z)) = (
        scalar_field,
        x_str,
        y_str,
        z_str,
    ) {

        let result = gradient(&field, (x, y, z));

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub extern "C" fn rssn_json_vector_divergence(
    v_json: *const c_char,
    x_var: *const c_char,
    y_var: *const c_char,
    z_var: *const c_char,
) -> *mut c_char {

    let v: Option<Vector> = from_json_string(v_json);

    let x_str = unsafe {

        if x_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(x_var)
                .to_str()
                .ok()
        }
    };

    let y_str = unsafe {

        if y_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(y_var)
                .to_str()
                .ok()
        }
    };

    let z_str = unsafe {

        if z_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(z_var)
                .to_str()
                .ok()
        }
    };

    if let (Some(vector), Some(x), Some(y), Some(z)) = (
        v, x_str, y_str, z_str,
    ) {

        let result = divergence(&vector, (x, y, z));

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub extern "C" fn rssn_json_vector_curl(
    v_json: *const c_char,
    x_var: *const c_char,
    y_var: *const c_char,
    z_var: *const c_char,
) -> *mut c_char {

    let v: Option<Vector> = from_json_string(v_json);

    let x_str = unsafe {

        if x_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(x_var)
                .to_str()
                .ok()
        }
    };

    let y_str = unsafe {

        if y_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(y_var)
                .to_str()
                .ok()
        }
    };

    let z_str = unsafe {

        if z_var.is_null() {

            None
        } else {

            std::ffi::CStr::from_ptr(z_var)
                .to_str()
                .ok()
        }
    };

    if let (Some(vector), Some(x), Some(y), Some(z)) = (
        v, x_str, y_str, z_str,
    ) {

        let result = curl(&vector, (x, y, z));

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
