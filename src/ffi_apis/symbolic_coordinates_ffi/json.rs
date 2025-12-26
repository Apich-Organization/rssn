use std::ffi::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::coordinates::*;
use crate::symbolic::core::Expr;

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
