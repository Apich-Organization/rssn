use crate::ffi_apis::common::*;
use crate::symbolic::coordinates::*;
use crate::symbolic::core::Expr;

#[no_mangle]
pub extern "C" fn rssn_bincode_transform_point(
    point_buf: BincodeBuffer,
    from_buf: BincodeBuffer,
    to_buf: BincodeBuffer,
) -> BincodeBuffer {
    let point: Option<Vec<Expr>> = from_bincode_buffer(&point_buf);
    let from: Option<CoordinateSystem> = from_bincode_buffer(&from_buf);
    let to: Option<CoordinateSystem> = from_bincode_buffer(&to_buf);

    if let (Some(p), Some(f), Some(t)) = (point, from, to) {
        match transform_point(&p, f, t) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_transform_expression(
    expr_buf: BincodeBuffer,
    from_buf: BincodeBuffer,
    to_buf: BincodeBuffer,
) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    let from: Option<CoordinateSystem> = from_bincode_buffer(&from_buf);
    let to: Option<CoordinateSystem> = from_bincode_buffer(&to_buf);

    if let (Some(e), Some(f), Some(t)) = (expr, from, to) {
        match transform_expression(&e, f, t) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_coordinates_get_metric_tensor(
    system_buf: BincodeBuffer,
) -> BincodeBuffer {
    let system: Option<CoordinateSystem> = from_bincode_buffer(&system_buf);

    if let Some(s) = system {
        match get_metric_tensor(s) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_transform_contravariant_vector(
    comps_buf: BincodeBuffer,
    from_buf: BincodeBuffer,
    to_buf: BincodeBuffer,
) -> BincodeBuffer {
    let comps: Option<Vec<Expr>> = from_bincode_buffer(&comps_buf);
    let from: Option<CoordinateSystem> = from_bincode_buffer(&from_buf);
    let to: Option<CoordinateSystem> = from_bincode_buffer(&to_buf);

    if let (Some(c), Some(f), Some(t)) = (comps, from, to) {
        match transform_contravariant_vector(&c, f, t) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_transform_covariant_vector(
    comps_buf: BincodeBuffer,
    from_buf: BincodeBuffer,
    to_buf: BincodeBuffer,
) -> BincodeBuffer {
    let comps: Option<Vec<Expr>> = from_bincode_buffer(&comps_buf);
    let from: Option<CoordinateSystem> = from_bincode_buffer(&from_buf);
    let to: Option<CoordinateSystem> = from_bincode_buffer(&to_buf);

    if let (Some(c), Some(f), Some(t)) = (comps, from, to) {
        match transform_covariant_vector(&c, f, t) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_transform_divergence(
    comps_buf: BincodeBuffer,
    from_buf: BincodeBuffer,
) -> BincodeBuffer {
    let comps: Option<Vec<Expr>> = from_bincode_buffer(&comps_buf);
    let from: Option<CoordinateSystem> = from_bincode_buffer(&from_buf);

    if let (Some(c), Some(f)) = (comps, from) {
        match transform_divergence(&c, f) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_transform_curl(
    comps_buf: BincodeBuffer,
    from_buf: BincodeBuffer,
) -> BincodeBuffer {
    let comps: Option<Vec<Expr>> = from_bincode_buffer(&comps_buf);
    let from: Option<CoordinateSystem> = from_bincode_buffer(&from_buf);

    if let (Some(c), Some(f)) = (comps, from) {
        match transform_curl(&c, f) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_transform_gradient(
    scalar_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
    from_buf: BincodeBuffer,
    to_buf: BincodeBuffer,
) -> BincodeBuffer {
    let scalar: Option<Expr> = from_bincode_buffer(&scalar_buf);
    let vars: Option<Vec<String>> = from_bincode_buffer(&vars_buf);
    let from: Option<CoordinateSystem> = from_bincode_buffer(&from_buf);
    let to: Option<CoordinateSystem> = from_bincode_buffer(&to_buf);

    if let (Some(s), Some(v), Some(f), Some(t)) = (scalar, vars, from, to) {
        match transform_gradient(&s, &v, f, t) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}
