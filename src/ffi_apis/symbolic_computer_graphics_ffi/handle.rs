//! Handle-based FFI API for computer graphics operations.
//!
//! This module provides C-compatible FFI functions for 2D/3D transformations,
//! projections, Bezier curves, B-splines, and polygon mesh manipulation.

use crate::symbolic::core::Expr;
use crate::symbolic::vector::Vector;
use crate::symbolic::computer_graphics::{
    translation_2d, translation_3d, rotation_2d, rotation_3d_x, rotation_3d_y, rotation_3d_z,
    scaling_2d, scaling_3d, perspective_projection, orthographic_projection, look_at,
    BezierCurve, BSplineCurve, Polygon, PolygonMesh,
};

/// Generates a 3x3 2D translation matrix.
///
/// # Safety
/// All Expr pointers must be valid.
#[no_mangle]
pub unsafe extern "C" fn rssn_translation_2d(tx: *const Expr, ty: *const Expr) -> *mut Expr {
    if tx.is_null() || ty.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(translation_2d((*tx).clone(), (*ty).clone())))
}

/// Generates a 4x4 3D translation matrix.
#[no_mangle]
pub unsafe extern "C" fn rssn_translation_3d(
    tx: *const Expr, ty: *const Expr, tz: *const Expr
) -> *mut Expr {
    if tx.is_null() || ty.is_null() || tz.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(translation_3d((*tx).clone(), (*ty).clone(), (*tz).clone())))
}

/// Generates a 3x3 2D rotation matrix.
#[no_mangle]
pub unsafe extern "C" fn rssn_rotation_2d(angle: *const Expr) -> *mut Expr {
    if angle.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(rotation_2d((*angle).clone())))
}

/// Generates a 4x4 3D rotation matrix around the X-axis.
#[no_mangle]
pub unsafe extern "C" fn rssn_rotation_3d_x(angle: *const Expr) -> *mut Expr {
    if angle.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(rotation_3d_x((*angle).clone())))
}

/// Generates a 4x4 3D rotation matrix around the Y-axis.
#[no_mangle]
pub unsafe extern "C" fn rssn_rotation_3d_y(angle: *const Expr) -> *mut Expr {
    if angle.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(rotation_3d_y((*angle).clone())))
}

/// Generates a 4x4 3D rotation matrix around the Z-axis.
#[no_mangle]
pub unsafe extern "C" fn rssn_rotation_3d_z(angle: *const Expr) -> *mut Expr {
    if angle.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(rotation_3d_z((*angle).clone())))
}

/// Generates a 3x3 2D scaling matrix.
#[no_mangle]
pub unsafe extern "C" fn rssn_scaling_2d(sx: *const Expr, sy: *const Expr) -> *mut Expr {
    if sx.is_null() || sy.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(scaling_2d((*sx).clone(), (*sy).clone())))
}

/// Generates a 4x4 3D scaling matrix.
#[no_mangle]
pub unsafe extern "C" fn rssn_scaling_3d(
    sx: *const Expr, sy: *const Expr, sz: *const Expr
) -> *mut Expr {
    if sx.is_null() || sy.is_null() || sz.is_null() { return std::ptr::null_mut(); }
    Box::into_raw(Box::new(scaling_3d((*sx).clone(), (*sy).clone(), (*sz).clone())))
}

/// Generates a 4x4 perspective projection matrix.
#[no_mangle]
pub unsafe extern "C" fn rssn_perspective_projection(
    fovy: *const Expr,
    aspect: *const Expr,
    near: *const Expr,
    far: *const Expr
) -> *mut Expr {
    if fovy.is_null() || aspect.is_null() || near.is_null() || far.is_null() {
        return std::ptr::null_mut();
    }
    Box::into_raw(Box::new(perspective_projection(
        (*fovy).clone(), &*aspect, (*near).clone(), (*far).clone()
    )))
}

/// Generates a 4x4 orthographic projection matrix.
#[no_mangle]
pub unsafe extern "C" fn rssn_orthographic_projection(
    left: *const Expr, right: *const Expr,
    bottom: *const Expr, top: *const Expr,
    near: *const Expr, far: *const Expr
) -> *mut Expr {
    if left.is_null() || right.is_null() || bottom.is_null() || 
       top.is_null() || near.is_null() || far.is_null() {
        return std::ptr::null_mut();
    }
    Box::into_raw(Box::new(orthographic_projection(
        (*left).clone(), (*right).clone(), (*bottom).clone(),
        (*top).clone(), (*near).clone(), (*far).clone()
    )))
}

/// Creates a new Bezier curve from control points.
#[no_mangle]
pub unsafe extern "C" fn rssn_bezier_curve_new(
    points: *const Vector,
    count: usize
) -> *mut BezierCurve {
    if points.is_null() || count == 0 { return std::ptr::null_mut(); }
    let control_points: Vec<Vector> = std::slice::from_raw_parts(points, count)
        .iter()
        .cloned()
        .collect();
    let degree = count - 1;
    Box::into_raw(Box::new(BezierCurve { control_points, degree }))
}

/// Evaluates a Bezier curve at parameter t.
#[no_mangle]
pub unsafe extern "C" fn rssn_bezier_curve_evaluate(
    curve: *const BezierCurve,
    t: *const Expr
) -> *mut Vector {
    if curve.is_null() || t.is_null() { return std::ptr::null_mut(); }
    let result = (*curve).evaluate(&*t);
    Box::into_raw(Box::new(result))
}

/// Frees a Bezier curve.
#[no_mangle]
pub unsafe extern "C" fn rssn_bezier_curve_free(curve: *mut BezierCurve) {
    if !curve.is_null() { drop(Box::from_raw(curve)); }
}

/// Creates a new polygon mesh.
#[no_mangle]
pub unsafe extern "C" fn rssn_polygon_mesh_new(
    vertices: *const Vector,
    vertex_count: usize
) -> *mut PolygonMesh {
    if vertices.is_null() || vertex_count == 0 { return std::ptr::null_mut(); }
    let verts: Vec<Vector> = std::slice::from_raw_parts(vertices, vertex_count)
        .iter()
        .cloned()
        .collect();
    Box::into_raw(Box::new(PolygonMesh::new(verts, vec![])))
}

/// Frees a polygon mesh.
#[no_mangle]
pub unsafe extern "C" fn rssn_polygon_mesh_free(mesh: *mut PolygonMesh) {
    if !mesh.is_null() { drop(Box::from_raw(mesh)); }
}

/// Frees a Vector.
#[no_mangle]
pub unsafe extern "C" fn rssn_vector_free(vec: *mut Vector) {
    if !vec.is_null() { drop(Box::from_raw(vec)); }
}
