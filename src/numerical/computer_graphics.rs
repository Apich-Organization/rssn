//! # Numerical 3D Computer Graphics
//!
//! This module provides fundamental structs and functions for 3D computer graphics,
//! including representations for points and vectors, standard vector operations
//! (dot product, cross product), and functions to generate 4x4 transformation
//! matrices for translation, scaling, and rotation.

use crate::numerical::matrix::Matrix;
use std::ops::{Add, Div, Mul, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

// Basic vector operations
impl Vector3D {
    /// Computes the magnitude (length) of the vector.
    ///
    /// # Returns
    /// The magnitude as an `f64`.
    pub fn magnitude(&self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Normalizes the vector to have a magnitude of 1.
    ///
    /// If the vector's magnitude is zero, it returns the original vector unchanged.
    ///
    /// # Returns
    /// A new `Vector3D` representing the normalized vector.
    pub fn normalize(&self) -> Self {
        let mag = self.magnitude();
        if mag == 0.0 {
            *self
        } else {
            *self / mag
        }
    }
}

impl Add for Vector3D {
    type Output = Self;
    /// Performs vector addition component-wise.
    fn add(self, rhs: Self) -> Self {
        Vector3D {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for Vector3D {
    type Output = Self;
    /// Performs vector subtraction component-wise.
    fn sub(self, rhs: Self) -> Self {
        Vector3D {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Mul<f64> for Vector3D {
    type Output = Self;
    /// Performs scalar multiplication of the vector.
    fn mul(self, rhs: f64) -> Self {
        Vector3D {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Div<f64> for Vector3D {
    type Output = Self;
    /// Performs scalar division of the vector.
    fn div(self, rhs: f64) -> Self {
        Vector3D {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

/// Computes the dot product of two `Vector3D`s.
///
/// The dot product (or scalar product) is a scalar value that represents
/// the projection of one vector onto another. It is defined as `v1.x*v2.x + v1.y*v2.y + v1.z*v2.z`.
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// The scalar dot product as an `f64`.
pub fn dot_product(v1: &Vector3D, v2: &Vector3D) -> f64 {
    v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

/// Computes the cross product of two `Vector3D`s.
///
/// The cross product (or vector product) results in a new vector that is
/// perpendicular to both input vectors. Its magnitude is equal to the area
/// of the parallelogram that the two vectors form.
///
/// # Arguments
/// * `v1` - The first vector.
/// * `v2` - The second vector.
///
/// # Returns
/// A new `Vector3D` representing the cross product.
pub fn cross_product(v1: &Vector3D, v2: &Vector3D) -> Vector3D {
    Vector3D {
        x: v1.y * v2.z - v1.z * v2.y,
        y: v1.z * v2.x - v1.x * v2.z,
        z: v1.x * v2.y - v1.y * v2.x,
    }
}

// Transformation matrix generation
/// Generates a 4x4 translation matrix.
///
/// This matrix can be used to move objects in 3D space by `dx`, `dy`, and `dz`.
///
/// # Arguments
/// * `dx` - Translation along the x-axis.
/// * `dy` - Translation along the y-axis.
/// * `dz` - Translation along the z-axis.
///
/// # Returns
/// A `Matrix<f64>` representing the translation transformation.
pub fn translation_matrix(dx: f64, dy: f64, dz: f64) -> Matrix<f64> {
    Matrix::new(
        4,
        4,
        vec![
            1.0, 0.0, 0.0, dx, 0.0, 1.0, 0.0, dy, 0.0, 0.0, 1.0, dz, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 scaling matrix.
///
/// This matrix can be used to scale objects in 3D space by `sx`, `sy`, and `sz`.
///
/// # Arguments
/// * `sx` - Scaling factor along the x-axis.
/// * `sy` - Scaling factor along the y-axis.
/// * `sz` - Scaling factor along the z-axis.
///
/// # Returns
/// A `Matrix<f64>` representing the scaling transformation.
pub fn scaling_matrix(sx: f64, sy: f64, sz: f64) -> Matrix<f64> {
    Matrix::new(
        4,
        4,
        vec![
            sx, 0.0, 0.0, 0.0, 0.0, sy, 0.0, 0.0, 0.0, 0.0, sz, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 rotation matrix around the X-axis.
///
/// # Arguments
/// * `angle_rad` - The rotation angle in radians.
///
/// # Returns
/// A `Matrix<f64>` representing the rotation transformation around the X-axis.
pub fn rotation_matrix_x(angle_rad: f64) -> Matrix<f64> {
    let (s, c) = angle_rad.sin_cos();
    Matrix::new(
        4,
        4,
        vec![
            1.0, 0.0, 0.0, 0.0, 0.0, c, -s, 0.0, 0.0, s, c, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

// ... (rotation_matrix_y, rotation_matrix_z would follow a similar pattern)
