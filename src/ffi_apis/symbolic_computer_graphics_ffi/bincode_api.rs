//! Bincode-based FFI API for computer graphics operations.
//!
//! This module provides binary serialization-based FFI functions for 2D/3D transformations
//! and projections, offering efficient binary data interchange.

use crate::ffi_apis::common::*;
use crate::symbolic::computer_graphics::{
    reflection_2d,
    reflection_3d,
    rotation_2d,
    rotation_3d_x,
    rotation_3d_y,
    rotation_3d_z,
    rotation_axis_angle,
    scaling_2d,
    scaling_3d,
    shear_2d,
    translation_2d,
    translation_3d,
};
use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;

/// Generates a 3x3 2D translation matrix via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_translation_2d(
    tx_buf: BincodeBuffer,
    ty_buf: BincodeBuffer,
) -> BincodeBuffer {

    let tx: Option<Expr> =
        from_bincode_buffer(&tx_buf);

    let ty: Option<Expr> =
        from_bincode_buffer(&ty_buf);

    if let (Some(tx), Some(ty)) =
        (tx, ty)
    {

        to_bincode_buffer(
            &translation_2d(tx, ty),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 4x4 3D translation matrix via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_translation_3d(
    tx_buf: BincodeBuffer,
    ty_buf: BincodeBuffer,
    tz_buf: BincodeBuffer,
) -> BincodeBuffer {

    let tx: Option<Expr> =
        from_bincode_buffer(&tx_buf);

    let ty: Option<Expr> =
        from_bincode_buffer(&ty_buf);

    let tz: Option<Expr> =
        from_bincode_buffer(&tz_buf);

    if let (
        Some(tx),
        Some(ty),
        Some(tz),
    ) = (tx, ty, tz)
    {

        to_bincode_buffer(
            &translation_3d(tx, ty, tz),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 3x3 2D rotation matrix via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_rotation_2d(
    angle_buf: BincodeBuffer
) -> BincodeBuffer {

    let angle: Option<Expr> =
        from_bincode_buffer(&angle_buf);

    if let Some(a) = angle {

        to_bincode_buffer(&rotation_2d(
            a,
        ))
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 4x4 3D rotation matrix around X-axis via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_rotation_3d_x(
    angle_buf: BincodeBuffer
) -> BincodeBuffer {

    let angle: Option<Expr> =
        from_bincode_buffer(&angle_buf);

    if let Some(a) = angle {

        to_bincode_buffer(
            &rotation_3d_x(a),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 4x4 3D rotation matrix around Y-axis via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_rotation_3d_y(
    angle_buf: BincodeBuffer
) -> BincodeBuffer {

    let angle: Option<Expr> =
        from_bincode_buffer(&angle_buf);

    if let Some(a) = angle {

        to_bincode_buffer(
            &rotation_3d_y(a),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 4x4 3D rotation matrix around Z-axis via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_rotation_3d_z(
    angle_buf: BincodeBuffer
) -> BincodeBuffer {

    let angle: Option<Expr> =
        from_bincode_buffer(&angle_buf);

    if let Some(a) = angle {

        to_bincode_buffer(
            &rotation_3d_z(a),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 3x3 2D scaling matrix via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_scaling_2d(
    sx_buf: BincodeBuffer,
    sy_buf: BincodeBuffer,
) -> BincodeBuffer {

    let sx: Option<Expr> =
        from_bincode_buffer(&sx_buf);

    let sy: Option<Expr> =
        from_bincode_buffer(&sy_buf);

    if let (Some(sx), Some(sy)) =
        (sx, sy)
    {

        to_bincode_buffer(&scaling_2d(
            sx, sy,
        ))
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 4x4 3D scaling matrix via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_scaling_3d(
    sx_buf: BincodeBuffer,
    sy_buf: BincodeBuffer,
    sz_buf: BincodeBuffer,
) -> BincodeBuffer {

    let sx: Option<Expr> =
        from_bincode_buffer(&sx_buf);

    let sy: Option<Expr> =
        from_bincode_buffer(&sy_buf);

    let sz: Option<Expr> =
        from_bincode_buffer(&sz_buf);

    if let (
        Some(sx),
        Some(sy),
        Some(sz),
    ) = (sx, sy, sz)
    {

        to_bincode_buffer(&scaling_3d(
            sx, sy, sz,
        ))
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 3x3 2D shear matrix via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_shear_2d(
    shx_buf: BincodeBuffer,
    shy_buf: BincodeBuffer,
) -> BincodeBuffer {

    let shx: Option<Expr> =
        from_bincode_buffer(&shx_buf);

    let shy: Option<Expr> =
        from_bincode_buffer(&shy_buf);

    if let (Some(shx), Some(shy)) =
        (shx, shy)
    {

        to_bincode_buffer(&shear_2d(
            shx, shy,
        ))
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 3x3 2D reflection matrix via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_reflection_2d(
    angle_buf: BincodeBuffer
) -> BincodeBuffer {

    let angle: Option<Expr> =
        from_bincode_buffer(&angle_buf);

    if let Some(a) = angle {

        to_bincode_buffer(
            &reflection_2d(a),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 4x4 3D reflection matrix via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_reflection_3d(
    nx_buf: BincodeBuffer,
    ny_buf: BincodeBuffer,
    nz_buf: BincodeBuffer,
) -> BincodeBuffer {

    let nx: Option<Expr> =
        from_bincode_buffer(&nx_buf);

    let ny: Option<Expr> =
        from_bincode_buffer(&ny_buf);

    let nz: Option<Expr> =
        from_bincode_buffer(&nz_buf);

    if let (
        Some(nx),
        Some(ny),
        Some(nz),
    ) = (nx, ny, nz)
    {

        to_bincode_buffer(
            &reflection_3d(nx, ny, nz),
        )
    } else {

        BincodeBuffer::empty()
    }
}

/// Generates a 4x4 3D rotation around arbitrary axis via Bincode interface.
#[no_mangle]

pub extern "C" fn rssn_bincode_rotation_axis_angle(
    axis_x_buf: BincodeBuffer,
    axis_y_buf: BincodeBuffer,
    axis_z_buf: BincodeBuffer,
    angle_buf: BincodeBuffer,
) -> BincodeBuffer {

    let ax: Option<Expr> =
        from_bincode_buffer(
            &axis_x_buf,
        );

    let ay: Option<Expr> =
        from_bincode_buffer(
            &axis_y_buf,
        );

    let az: Option<Expr> =
        from_bincode_buffer(
            &axis_z_buf,
        );

    let angle: Option<Expr> =
        from_bincode_buffer(&angle_buf);

    if let (
        Some(ax),
        Some(ay),
        Some(az),
        Some(a),
    ) = (ax, ay, az, angle)
    {

        let axis =
            Vector::new(ax, ay, az);

        to_bincode_buffer(
            &rotation_axis_angle(
                &axis, a,
            ),
        )
    } else {

        BincodeBuffer::empty()
    }
}
