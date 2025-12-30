//! Handle-based FFI API for physics BEM functions.

use std::ptr;

use crate::physics::physics_bem::BoundaryCondition;
use crate::physics::physics_bem::{
    self,
};

/// Solves a 2D Laplace problem using BEM and returns the results as a flat array.
/// The `bcs_type` array should be 0 for Potential and 1 for Flux.
/// # Safety
/// This function is unsafe because it dereferences pointers.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_physics_bem_solve_laplace_2d(
    points_x: *const f64,
    points_y: *const f64,
    bcs_type: *const i32,
    bcs_value: *const f64,
    n: usize,
    out_u: *mut f64,
    out_q: *mut f64,
) -> i32 { unsafe {

    if points_x.is_null()
        || points_y.is_null()
        || bcs_type.is_null()
        || bcs_value.is_null()
        || out_u.is_null()
        || out_q.is_null()
    {

        return -1;
    }

    let points_x_slice =
        std::slice::from_raw_parts(
            points_x,
            n,
        );

    let points_y_slice =
        std::slice::from_raw_parts(
            points_y,
            n,
        );

    let bcs_type_slice =
        std::slice::from_raw_parts(
            bcs_type,
            n,
        );

    let bcs_value_slice =
        std::slice::from_raw_parts(
            bcs_value,
            n,
        );

    let points: Vec<(f64, f64)> =
        points_x_slice
            .iter()
            .zip(points_y_slice.iter())
            .map(|(&x, &y)| (x, y))
            .collect();

    let bcs : Vec<BoundaryCondition<f64>> = bcs_type_slice
        .iter()
        .zip(bcs_value_slice.iter())
        .map(|(&t, &v)| {
            if t == 0 {

                BoundaryCondition::Potential(v)
            } else {

                BoundaryCondition::Flux(v)
            }
        })
        .collect();

    match physics_bem::solve_laplace_bem_2d(&points, &bcs) {
        | Ok((u, q)) => {

            ptr::copy_nonoverlapping(u.as_ptr(), out_u, n);

            ptr::copy_nonoverlapping(q.as_ptr(), out_q, n);

            0
        },
        | Err(_) => -2,
    }
}}
