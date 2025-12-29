use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::coordinates::{CoordinateSystem, transform_point, transform_expression, get_metric_tensor, transform_contravariant_vector, transform_covariant_vector, transform_divergence, transform_curl, transform_gradient};
use crate::symbolic::core::Expr;

/// Transforms a point from one coordinate system to another using bincode-serialized inputs.

///

/// Takes `BincodeBuffer`s for the point (`Vec<Expr>`), the source `CoordinateSystem`,

/// and the target `CoordinateSystem`.

/// Returns a `BincodeBuffer` containing the transformed point (`Vec<Expr>`).

#[no_mangle]

pub extern "C" fn rssn_bincode_transform_point(
    point_buf: BincodeBuffer,

    from_buf: BincodeBuffer,

    to_buf: BincodeBuffer,
) -> BincodeBuffer {

    let point: Option<Vec<Expr>> =
        from_bincode_buffer(&point_buf);

    let from: Option<CoordinateSystem> =
        from_bincode_buffer(&from_buf);

    let to: Option<CoordinateSystem> =
        from_bincode_buffer(&to_buf);

    if let (Some(p), Some(f), Some(t)) =
        (point, from, to)
    {

        match transform_point(&p, f, t)
        {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Transforms an expression from one coordinate system to another using bincode-serialized inputs.

///

/// Takes `BincodeBuffer`s for the expression (`Expr`), the source `CoordinateSystem`,

/// and the target `CoordinateSystem`.

/// Returns a `BincodeBuffer` containing the transformed expression (`Expr`).

#[no_mangle]

pub extern "C" fn rssn_bincode_transform_expression(
    expr_buf: BincodeBuffer,

    from_buf: BincodeBuffer,

    to_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let from: Option<CoordinateSystem> =
        from_bincode_buffer(&from_buf);

    let to: Option<CoordinateSystem> =
        from_bincode_buffer(&to_buf);

    if let (Some(e), Some(f), Some(t)) =
        (expr, from, to)
    {

        match transform_expression(
            &e, f, t,
        ) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Gets the metric tensor for a given coordinate system using bincode-serialized input.

///

/// Takes a `BincodeBuffer` containing the `CoordinateSystem`.

/// Returns a `BincodeBuffer` containing the metric tensor (`Vec<Vec<Expr>>`).

#[no_mangle]

pub extern "C" fn rssn_bincode_coordinates_get_metric_tensor(
    system_buf: BincodeBuffer
) -> BincodeBuffer {

    let system: Option<
        CoordinateSystem,
    > = from_bincode_buffer(
        &system_buf,
    );

    if let Some(s) = system {

        match get_metric_tensor(s) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Transforms a contravariant vector from one coordinate system to another using bincode-serialized inputs.

///

/// Takes `BincodeBuffer`s for the vector components (`Vec<Expr>`), the source `CoordinateSystem`,

/// and the target `CoordinateSystem`.

/// Returns a `BincodeBuffer` containing the transformed contravariant vector components (`Vec<Expr>`).

#[no_mangle]

pub extern "C" fn rssn_bincode_transform_contravariant_vector(
    comps_buf: BincodeBuffer,

    from_buf: BincodeBuffer,

    to_buf: BincodeBuffer,
) -> BincodeBuffer {

    let comps: Option<Vec<Expr>> =
        from_bincode_buffer(&comps_buf);

    let from: Option<CoordinateSystem> =
        from_bincode_buffer(&from_buf);

    let to: Option<CoordinateSystem> =
        from_bincode_buffer(&to_buf);

    if let (Some(c), Some(f), Some(t)) =
        (comps, from, to)
    {

        match transform_contravariant_vector(&c, f, t) {
            | Ok(result) => to_bincode_buffer(&result),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Transforms a covariant vector from one coordinate system to another using bincode-serialized inputs.

///

/// Takes `BincodeBuffer`s for the vector components (`Vec<Expr>`), the source `CoordinateSystem`,

/// and the target `CoordinateSystem`.

/// Returns a `BincodeBuffer` containing the transformed covariant vector components (`Vec<Expr>`).

#[no_mangle]

pub extern "C" fn rssn_bincode_transform_covariant_vector(
    comps_buf: BincodeBuffer,

    from_buf: BincodeBuffer,

    to_buf: BincodeBuffer,
) -> BincodeBuffer {

    let comps: Option<Vec<Expr>> =
        from_bincode_buffer(&comps_buf);

    let from: Option<CoordinateSystem> =
        from_bincode_buffer(&from_buf);

    let to: Option<CoordinateSystem> =
        from_bincode_buffer(&to_buf);

    if let (Some(c), Some(f), Some(t)) =
        (comps, from, to)
    {

        match transform_covariant_vector(
            &c, f, t,
        ) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Transforms the divergence of a vector field from one coordinate system to another using bincode-serialized inputs.

///

/// Takes `BincodeBuffer`s for the vector components (`Vec<Expr>`) and the source `CoordinateSystem`.

/// Returns a `BincodeBuffer` containing the transformed divergence (`Expr`).

#[no_mangle]

pub extern "C" fn rssn_bincode_transform_divergence(
    comps_buf: BincodeBuffer,

    from_buf: BincodeBuffer,
) -> BincodeBuffer {

    let comps: Option<Vec<Expr>> =
        from_bincode_buffer(&comps_buf);

    let from: Option<CoordinateSystem> =
        from_bincode_buffer(&from_buf);

    if let (Some(c), Some(f)) =
        (comps, from)
    {

        match transform_divergence(
            &c, f,
        ) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Transforms the curl of a vector field from one coordinate system to another using bincode-serialized inputs.

///

/// Takes `BincodeBuffer`s for the vector components (`Vec<Expr>`) and the source `CoordinateSystem`.

/// Returns a `BincodeBuffer` containing the transformed curl (`Vec<Expr>`).

#[no_mangle]

pub extern "C" fn rssn_bincode_transform_curl(
    comps_buf: BincodeBuffer,

    from_buf: BincodeBuffer,
) -> BincodeBuffer {

    let comps: Option<Vec<Expr>> =
        from_bincode_buffer(&comps_buf);

    let from: Option<CoordinateSystem> =
        from_bincode_buffer(&from_buf);

    if let (Some(c), Some(f)) =
        (comps, from)
    {

        match transform_curl(&c, f) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Transforms the gradient of a scalar function from one coordinate system to another using bincode-serialized inputs.

///

/// Takes `BincodeBuffer`s for the scalar function (`Expr`), the variables (`Vec<String>`),

/// the source `CoordinateSystem`, and the target `CoordinateSystem`.

/// Returns a `BincodeBuffer` containing the transformed gradient (`Vec<Expr>`).

#[no_mangle]

pub extern "C" fn rssn_bincode_transform_gradient(
    scalar_buf: BincodeBuffer,

    vars_buf: BincodeBuffer,

    from_buf: BincodeBuffer,

    to_buf: BincodeBuffer,
) -> BincodeBuffer {

    let scalar: Option<Expr> =
        from_bincode_buffer(
            &scalar_buf,
        );

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    let from: Option<CoordinateSystem> =
        from_bincode_buffer(&from_buf);

    let to: Option<CoordinateSystem> =
        from_bincode_buffer(&to_buf);

    if let (
        Some(s),
        Some(v),
        Some(f),
        Some(t),
    ) = (
        scalar,
        vars,
        from,
        to,
    ) {

        match transform_gradient(
            &s, &v, f, t,
        ) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}
