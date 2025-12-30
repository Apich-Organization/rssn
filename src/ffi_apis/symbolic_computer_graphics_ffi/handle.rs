//! Handle-based FFI API for computer graphics operations.
//!
//! This module provides C-compatible FFI functions for 2D/3D transformations,
//! projections, Bezier curves, B-splines, and polygon mesh manipulation.

use crate::symbolic::computer_graphics::orthographic_projection;
use crate::symbolic::computer_graphics::perspective_projection;
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
use crate::symbolic::computer_graphics::BezierCurve;
use crate::symbolic::computer_graphics::PolygonMesh;
use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;

/// Generates a 3x3 2D translation matrix.
///
/// # Safety
/// All Expr pointers must be valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_translation_2d(
    tx: *const Expr,
    ty: *const Expr,
) -> *mut Expr { unsafe {

    if tx.is_null() || ty.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        translation_2d(
            (*tx).clone(),
            (*ty).clone(),
        ),
    ))
}}

/// Generates a 4x4 3D translation matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_translation_3d(
    tx: *const Expr,
    ty: *const Expr,
    tz: *const Expr,
) -> *mut Expr { unsafe {

    if tx.is_null()
        || ty.is_null()
        || tz.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        translation_3d(
            (*tx).clone(),
            (*ty).clone(),
            (*tz).clone(),
        ),
    ))
}}

/// Generates a 3x3 2D rotation matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rotation_2d(
    angle: *const Expr
) -> *mut Expr { unsafe {

    if angle.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        rotation_2d((*angle).clone()),
    ))
}}

/// Generates a 4x4 3D rotation matrix around the X-axis.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rotation_3d_x(
    angle: *const Expr
) -> *mut Expr { unsafe {

    if angle.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        rotation_3d_x((*angle).clone()),
    ))
}}

/// Generates a 4x4 3D rotation matrix around the Y-axis.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rotation_3d_y(
    angle: *const Expr
) -> *mut Expr { unsafe {

    if angle.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        rotation_3d_y((*angle).clone()),
    ))
}}

/// Generates a 4x4 3D rotation matrix around the Z-axis.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rotation_3d_z(
    angle: *const Expr
) -> *mut Expr { unsafe {

    if angle.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        rotation_3d_z((*angle).clone()),
    ))
}}

/// Generates a 3x3 2D scaling matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_scaling_2d(
    sx: *const Expr,
    sy: *const Expr,
) -> *mut Expr { unsafe {

    if sx.is_null() || sy.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        scaling_2d(
            (*sx).clone(),
            (*sy).clone(),
        ),
    ))
}}

/// Generates a 4x4 3D scaling matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_scaling_3d(
    sx: *const Expr,
    sy: *const Expr,
    sz: *const Expr,
) -> *mut Expr { unsafe {

    if sx.is_null()
        || sy.is_null()
        || sz.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        scaling_3d(
            (*sx).clone(),
            (*sy).clone(),
            (*sz).clone(),
        ),
    ))
}}

/// Generates a 3x3 2D shear matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_shear_2d(
    shx: *const Expr,
    shy: *const Expr,
) -> *mut Expr { unsafe {

    if shx.is_null() || shy.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(shear_2d(
        (*shx).clone(),
        (*shy).clone(),
    )))
}}

/// Generates a 3x3 2D reflection matrix across a line.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_reflection_2d(
    angle: *const Expr
) -> *mut Expr { unsafe {

    if angle.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        reflection_2d((*angle).clone()),
    ))
}}

/// Generates a 4x4 3D reflection matrix across a plane.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_reflection_3d(
    nx: *const Expr,
    ny: *const Expr,
    nz: *const Expr,
) -> *mut Expr { unsafe {

    if nx.is_null()
        || ny.is_null()
        || nz.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        reflection_3d(
            (*nx).clone(),
            (*ny).clone(),
            (*nz).clone(),
        ),
    ))
}}

/// Generates a 4x4 3D rotation matrix around an arbitrary axis.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_rotation_axis_angle(
    axis: *const Vector,
    angle: *const Expr,
) -> *mut Expr { unsafe {

    if axis.is_null() || angle.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        rotation_axis_angle(
            &*axis,
            (*angle).clone(),
        ),
    ))
}}

/// Generates a 4x4 perspective projection matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_perspective_projection(
    fovy: *const Expr,
    aspect: *const Expr,
    near: *const Expr,
    far: *const Expr,
) -> *mut Expr { unsafe {

    if fovy.is_null()
        || aspect.is_null()
        || near.is_null()
        || far.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        perspective_projection(
            (*fovy).clone(),
            &*aspect,
            (*near).clone(),
            (*far).clone(),
        ),
    ))
}}

/// Generates a 4x4 orthographic projection matrix.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_orthographic_projection(
    left: *const Expr,
    right: *const Expr,
    bottom: *const Expr,
    top: *const Expr,
    near: *const Expr,
    far: *const Expr,
) -> *mut Expr { unsafe {

    if left.is_null()
        || right.is_null()
        || bottom.is_null()
        || top.is_null()
        || near.is_null()
        || far.is_null()
    {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        orthographic_projection(
            (*left).clone(),
            (*right).clone(),
            (*bottom).clone(),
            (*top).clone(),
            (*near).clone(),
            (*far).clone(),
        ),
    ))
}}

/// Creates a new Bezier curve from control points.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bezier_curve_new(
    points: *const Vector,
    count: usize,
) -> *mut BezierCurve { unsafe {

    if points.is_null() || count == 0 {

        return std::ptr::null_mut();
    }

    let control_points: Vec<Vector> =
        std::slice::from_raw_parts(
            points,
            count,
        )
        .to_vec();

    let degree = count - 1;

    Box::into_raw(Box::new(
        BezierCurve {
            control_points,
            degree,
        },
    ))
}}

/// Evaluates a Bezier curve at parameter t.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bezier_curve_evaluate(
    curve: *const BezierCurve,
    t: *const Expr,
) -> *mut Vector { unsafe {

    if curve.is_null() || t.is_null() {

        return std::ptr::null_mut();
    }

    let result = (*curve).evaluate(&*t);

    Box::into_raw(Box::new(result))
}}

/// Computes the derivative (tangent) of a Bezier curve at parameter t.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bezier_curve_derivative(
    curve: *const BezierCurve,
    t: *const Expr,
) -> *mut Vector { unsafe {

    if curve.is_null() || t.is_null() {

        return std::ptr::null_mut();
    }

    let result =
        (*curve).derivative(&*t);

    Box::into_raw(Box::new(result))
}}

/// Splits a Bezier curve at parameter t into two curves.
/// Returns left curve. Use `rssn_bezier_curve_split_right` for the right curve.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bezier_curve_split_left(
    curve: *const BezierCurve,
    t: *const Expr,
) -> *mut BezierCurve { unsafe {

    if curve.is_null() || t.is_null() {

        return std::ptr::null_mut();
    }

    let (left, _right) =
        (*curve).split(&*t);

    Box::into_raw(Box::new(left))
}}

/// Splits a Bezier curve at parameter t into two curves.
/// Returns right curve.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bezier_curve_split_right(
    curve: *const BezierCurve,
    t: *const Expr,
) -> *mut BezierCurve { unsafe {

    if curve.is_null() || t.is_null() {

        return std::ptr::null_mut();
    }

    let (_left, right) =
        (*curve).split(&*t);

    Box::into_raw(Box::new(right))
}}

/// Frees a Bezier curve.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bezier_curve_free(
    curve: *mut BezierCurve
) { unsafe {

    if !curve.is_null() {

        drop(Box::from_raw(curve));
    }
}}

/// Creates a new polygon mesh.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_polygon_mesh_new(
    vertices: *const Vector,
    vertex_count: usize,
) -> *mut PolygonMesh { unsafe {

    if vertices.is_null()
        || vertex_count == 0
    {

        return std::ptr::null_mut();
    }

    let verts: Vec<Vector> =
        std::slice::from_raw_parts(
            vertices,
            vertex_count,
        )
        .to_vec();

    Box::into_raw(Box::new(
        PolygonMesh::new(verts, vec![]),
    ))
}}

/// Triangulates a polygon mesh.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_polygon_mesh_triangulate(
    mesh: *const PolygonMesh
) -> *mut PolygonMesh { unsafe {

    if mesh.is_null() {

        return std::ptr::null_mut();
    }

    let triangulated =
        (*mesh).triangulate();

    Box::into_raw(Box::new(
        triangulated,
    ))
}}

/// Frees a polygon mesh.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_polygon_mesh_free(
    mesh: *mut PolygonMesh
) { unsafe {

    if !mesh.is_null() {

        drop(Box::from_raw(mesh));
    }
}}

/// Frees a Vector.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_vector_free(
    vec: *mut Vector
) { unsafe {

    if !vec.is_null() {

        drop(Box::from_raw(vec));
    }
}}
