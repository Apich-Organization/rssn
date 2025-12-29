//! JSON-based FFI API for computer graphics operations.
//!
//! This module provides JSON string-based FFI functions for 2D/3D transformations
//! and projections, enabling language-agnostic integration.

use std::os::raw::c_char;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::computer_graphics::reflection_2d;
use crate::symbolic::computer_graphics::reflection_3d;
use crate::symbolic::computer_graphics::rotation_2d;
use crate::symbolic::computer_graphics::rotation_3d_x;
use crate::symbolic::computer_graphics::rotation_3d_y;
use crate::symbolic::computer_graphics::rotation_3d_z;
use crate::symbolic::computer_graphics::rotation_axis_angle;
use crate::symbolic::computer_graphics::scaling_2d;
use crate::symbolic::computer_graphics::scaling_3d;
use crate::symbolic::computer_graphics::shear_2d;
use crate::symbolic::computer_graphics::translation_2d;
use crate::symbolic::computer_graphics::translation_3d;
use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;

/// Generates a 3x3 2D translation matrix via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_translation_2d(
    tx_json: *const c_char,
    ty_json: *const c_char,
) -> *mut c_char {

    let tx: Option<Expr> =
        from_json_string(tx_json);

    let ty: Option<Expr> =
        from_json_string(ty_json);

    if let (Some(tx), Some(ty)) =
        (tx, ty)
    {

        to_json_string(&translation_2d(
            tx, ty,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D translation matrix via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_translation_3d(
    tx_json: *const c_char,
    ty_json: *const c_char,
    tz_json: *const c_char,
) -> *mut c_char {

    let tx: Option<Expr> =
        from_json_string(tx_json);

    let ty: Option<Expr> =
        from_json_string(ty_json);

    let tz: Option<Expr> =
        from_json_string(tz_json);

    if let (
        Some(tx),
        Some(ty),
        Some(tz),
    ) = (tx, ty, tz)
    {

        to_json_string(&translation_3d(
            tx, ty, tz,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 3x3 2D rotation matrix via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rotation_2d(
    angle_json: *const c_char
) -> *mut c_char {

    let angle: Option<Expr> =
        from_json_string(angle_json);

    if let Some(a) = angle {

        to_json_string(&rotation_2d(a))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D rotation matrix around X-axis via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rotation_3d_x(
    angle_json: *const c_char
) -> *mut c_char {

    let angle: Option<Expr> =
        from_json_string(angle_json);

    if let Some(a) = angle {

        to_json_string(&rotation_3d_x(
            a,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D rotation matrix around Y-axis via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rotation_3d_y(
    angle_json: *const c_char
) -> *mut c_char {

    let angle: Option<Expr> =
        from_json_string(angle_json);

    if let Some(a) = angle {

        to_json_string(&rotation_3d_y(
            a,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D rotation matrix around Z-axis via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rotation_3d_z(
    angle_json: *const c_char
) -> *mut c_char {

    let angle: Option<Expr> =
        from_json_string(angle_json);

    if let Some(a) = angle {

        to_json_string(&rotation_3d_z(
            a,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 3x3 2D scaling matrix via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_scaling_2d(
    sx_json: *const c_char,
    sy_json: *const c_char,
) -> *mut c_char {

    let sx: Option<Expr> =
        from_json_string(sx_json);

    let sy: Option<Expr> =
        from_json_string(sy_json);

    if let (Some(sx), Some(sy)) =
        (sx, sy)
    {

        to_json_string(&scaling_2d(
            sx, sy,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D scaling matrix via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_scaling_3d(
    sx_json: *const c_char,
    sy_json: *const c_char,
    sz_json: *const c_char,
) -> *mut c_char {

    let sx: Option<Expr> =
        from_json_string(sx_json);

    let sy: Option<Expr> =
        from_json_string(sy_json);

    let sz: Option<Expr> =
        from_json_string(sz_json);

    if let (
        Some(sx),
        Some(sy),
        Some(sz),
    ) = (sx, sy, sz)
    {

        to_json_string(&scaling_3d(
            sx, sy, sz,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 3x3 2D shear matrix via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_shear_2d(
    shx_json: *const c_char,
    shy_json: *const c_char,
) -> *mut c_char {

    let shx: Option<Expr> =
        from_json_string(shx_json);

    let shy: Option<Expr> =
        from_json_string(shy_json);

    if let (Some(shx), Some(shy)) =
        (shx, shy)
    {

        to_json_string(&shear_2d(
            shx, shy,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 3x3 2D reflection matrix via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_reflection_2d(
    angle_json: *const c_char
) -> *mut c_char {

    let angle: Option<Expr> =
        from_json_string(angle_json);

    if let Some(a) = angle {

        to_json_string(&reflection_2d(
            a,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D reflection matrix via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_reflection_3d(
    nx_json: *const c_char,
    ny_json: *const c_char,
    nz_json: *const c_char,
) -> *mut c_char {

    let nx: Option<Expr> =
        from_json_string(nx_json);

    let ny: Option<Expr> =
        from_json_string(ny_json);

    let nz: Option<Expr> =
        from_json_string(nz_json);

    if let (
        Some(nx),
        Some(ny),
        Some(nz),
    ) = (nx, ny, nz)
    {

        to_json_string(&reflection_3d(
            nx, ny, nz,
        ))
    } else {

        std::ptr::null_mut()
    }
}

/// Generates a 4x4 3D rotation around arbitrary axis via JSON interface.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_rotation_axis_angle(
    axis_x_json: *const c_char,
    axis_y_json: *const c_char,
    axis_z_json: *const c_char,
    angle_json: *const c_char,
) -> *mut c_char {

    let ax: Option<Expr> =
        from_json_string(axis_x_json);

    let ay: Option<Expr> =
        from_json_string(axis_y_json);

    let az: Option<Expr> =
        from_json_string(axis_z_json);

    let angle: Option<Expr> =
        from_json_string(angle_json);

    if let (
        Some(ax),
        Some(ay),
        Some(az),
        Some(a),
    ) = (ax, ay, az, angle)
    {

        let axis =
            Vector::new(ax, ay, az);

        to_json_string(
            &rotation_axis_angle(
                &axis, a,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}
