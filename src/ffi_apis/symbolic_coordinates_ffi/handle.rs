use crate::symbolic::coordinates::*;
use crate::symbolic::core::Expr;

#[no_mangle]

pub extern "C" fn rssn_transform_point_handle(
    point: *const Vec<Expr>,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Vec<Expr> {

    let point_ref = unsafe {

        &*point
    };

    match transform_point(point_ref, from, to) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn rssn_transform_expression_handle(
    expr: *const Expr,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    match transform_expression(expr_ref, from, to) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn rssn_coordinates_get_metric_tensor_handle(system: CoordinateSystem) -> *mut Expr {

    match get_metric_tensor(system) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn rssn_transform_contravariant_vector_handle(
    comps: *const Vec<Expr>,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Vec<Expr> {

    let comps_ref = unsafe {

        &*comps
    };

    match transform_contravariant_vector(comps_ref, from, to) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn rssn_transform_covariant_vector_handle(
    comps: *const Vec<Expr>,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Vec<Expr> {

    let comps_ref = unsafe {

        &*comps
    };

    match transform_covariant_vector(comps_ref, from, to) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn rssn_transform_divergence_handle(
    comps: *const Vec<Expr>,
    from: CoordinateSystem,
) -> *mut Expr {

    let comps_ref = unsafe {

        &*comps
    };

    match transform_divergence(comps_ref, from) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn rssn_transform_curl_handle(
    comps: *const Vec<Expr>,
    from: CoordinateSystem,
) -> *mut Vec<Expr> {

    let comps_ref = unsafe {

        &*comps
    };

    match transform_curl(comps_ref, from) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub extern "C" fn rssn_transform_gradient_handle(
    scalar: *const Expr,
    vars: *const Vec<String>,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Vec<Expr> {

    let scalar_ref = unsafe {

        &*scalar
    };

    let vars_ref = unsafe {

        &*vars
    };

    match transform_gradient(
        scalar_ref, vars_ref, from, to,
    ) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}
