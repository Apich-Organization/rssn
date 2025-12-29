use crate::symbolic::coordinates::{CoordinateSystem, transform_point, transform_expression, get_metric_tensor, transform_contravariant_vector, transform_covariant_vector, transform_divergence, transform_curl, transform_gradient};
use crate::symbolic::core::Expr;

/// Transforms a point between coordinate systems and returns its components in the target system.
///
/// The point is represented as a vector of symbolic expressions (e.g., \(x,y,z\)), and the
/// transformation applies the appropriate coordinate mapping (Cartesian, polar, spherical, etc.).
///
/// # Arguments
///
/// * `point` - Pointer to a `Vec<Expr>` containing the point coordinates in the `from` system.
/// * `from` - Source [`CoordinateSystem`] in which `point` is expressed.
/// * `to` - Target [`CoordinateSystem`] to which the point is transformed.
///
/// # Returns
///
/// A newly allocated `Vec<Expr>` pointer with the point coordinates in the `to` system, or
/// null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer and returns ownership
/// of a heap-allocated vector to the caller.
#[no_mangle]

pub extern "C" fn rssn_transform_point_handle(
    point: *const Vec<Expr>,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Vec<Expr> {

    let point_ref = unsafe {

        &*point
    };

    match transform_point(
        point_ref,
        from,
        to,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Transforms a scalar expression between coordinate systems.
///
/// This replaces the variables of `expr` according to the mapping between the `from` and
/// `to` coordinate systems, yielding an equivalent symbolic expression in the target system.
///
/// # Arguments
///
/// * `expr` - Pointer to an `Expr` representing the scalar field in the `from` system.
/// * `from` - Source [`CoordinateSystem`] of the variables in `expr`.
/// * `to` - Target [`CoordinateSystem`] to which the expression is transformed.
///
/// # Returns
///
/// A newly allocated `Expr` pointer representing the transformed expression, or null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw `Expr` pointer and returns
/// ownership of a heap-allocated `Expr` to the caller.
#[no_mangle]

pub extern "C" fn rssn_transform_expression_handle(
    expr: *const Expr,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Expr {

    let expr_ref = unsafe {

        &*expr
    };

    match transform_expression(
        expr_ref,
        from,
        to,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Returns the metric tensor for a given coordinate system as a symbolic expression.
///
/// The metric tensor encodes the inner product and volume element for the coordinate
/// system, and is typically represented as a matrix-valued `Expr`.
///
/// # Arguments
///
/// * `system` - [`CoordinateSystem`] for which to compute the metric tensor.
///
/// # Returns
///
/// A newly allocated `Expr` pointer representing the metric tensor, or null on failure.
///
/// # Safety
///
/// This function is unsafe at the FFI boundary; the returned pointer must be eventually
/// freed by the caller using the appropriate deallocation routine.
#[no_mangle]

pub extern "C" fn rssn_coordinates_get_metric_tensor_handle(
    system: CoordinateSystem
) -> *mut Expr {

    match get_metric_tensor(system) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Transforms contravariant vector components between coordinate systems.
///
/// Contravariant components transform with the Jacobian of the coordinate change, corresponding
/// to vector components with upper indices.
///
/// # Arguments
///
/// * `comps` - Pointer to a `Vec<Expr>` containing contravariant components in the `from` system.
/// * `from` - Source [`CoordinateSystem`].
/// * `to` - Target [`CoordinateSystem`].
///
/// # Returns
///
/// A newly allocated `Vec<Expr>` pointer with contravariant components in the `to` system, or
/// null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer to a vector and returns
/// ownership of a heap-allocated vector to the caller.
#[no_mangle]

pub extern "C" fn rssn_transform_contravariant_vector_handle(
    comps: *const Vec<Expr>,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Vec<Expr> {

    let comps_ref = unsafe {

        &*comps
    };

    match transform_contravariant_vector(
        comps_ref,
        from,
        to,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Transforms covariant vector components between coordinate systems.
///
/// Covariant components transform with the inverse Jacobian and correspond to components
/// with lower indices (1-forms).
///
/// # Arguments
///
/// * `comps` - Pointer to a `Vec<Expr>` containing covariant components in the `from` system.
/// * `from` - Source [`CoordinateSystem`].
/// * `to` - Target [`CoordinateSystem`].
///
/// # Returns
///
/// A newly allocated `Vec<Expr>` pointer with covariant components in the `to` system, or
/// null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer to a vector and returns
/// ownership of a heap-allocated vector to the caller.
#[no_mangle]

pub extern "C" fn rssn_transform_covariant_vector_handle(
    comps: *const Vec<Expr>,
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> *mut Vec<Expr> {

    let comps_ref = unsafe {

        &*comps
    };

    match transform_covariant_vector(
        comps_ref,
        from,
        to,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Computes the divergence of a vector field in a given coordinate system.
///
/// The divergence is computed using the metric and Christoffel symbols associated with
/// `from`, yielding a scalar `Expr` representing \(\nabla \cdot \vec{v}\).
///
/// # Arguments
///
/// * `comps` - Pointer to a `Vec<Expr>` containing the vector components in the `from` system.
/// * `from` - [`CoordinateSystem`] with respect to which the divergence is taken.
///
/// # Returns
///
/// A newly allocated `Expr` pointer representing the divergence, or null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw vector pointer and returns
/// ownership of a heap-allocated `Expr` to the caller.
#[no_mangle]

pub extern "C" fn rssn_transform_divergence_handle(
    comps: *const Vec<Expr>,
    from: CoordinateSystem,
) -> *mut Expr {

    let comps_ref = unsafe {

        &*comps
    };

    match transform_divergence(
        comps_ref,
        from,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Computes the curl of a vector field in a given coordinate system.
///
/// The curl is computed using the metric and Levi-Civita tensor of the `from` system,
/// yielding a vector-valued `Expr` representing \(\nabla \times \vec{v}\).
///
/// # Arguments
///
/// * `comps` - Pointer to a `Vec<Expr>` containing the vector components in the `from` system.
/// * `from` - [`CoordinateSystem`] with respect to which the curl is taken.
///
/// # Returns
///
/// A newly allocated `Vec<Expr>` pointer representing the curl, or null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw vector pointer and returns
/// ownership of a heap-allocated vector to the caller.
#[no_mangle]

pub extern "C" fn rssn_transform_curl_handle(
    comps: *const Vec<Expr>,
    from: CoordinateSystem,
) -> *mut Vec<Expr> {

    let comps_ref = unsafe {

        &*comps
    };

    match transform_curl(
        comps_ref,
        from,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Computes the gradient of a scalar field and transforms it between coordinate systems.
///
/// The gradient is first formed with respect to the variables `vars` in the `from` system,
/// then mapped into the `to` system as a vector of symbolic components.
///
/// # Arguments
///
/// * `scalar` - Pointer to an `Expr` representing the scalar field.
/// * `vars` - Pointer to a `Vec<String>` listing the coordinate variables.
/// * `from` - Source [`CoordinateSystem`].
/// * `to` - Target [`CoordinateSystem`].
///
/// # Returns
///
/// A newly allocated `Vec<Expr>` pointer representing the gradient components in the `to`
/// system, or null on failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers and returns ownership
/// of a heap-allocated vector to the caller.
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
        scalar_ref,
        vars_ref,
        from,
        to,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}
