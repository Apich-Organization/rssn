//! # Numerical 3D Computer Graphics
//!
//! This module provides fundamental structs and functions for 2D and 3D computer graphics,
//! including representations for points and vectors, standard vector operations
//! (dot product, cross product), and functions to generate transformation matrices
//! for translation, scaling, rotation, projections, and camera views.
//!
//! ## Features
//!
//! ### Basic Types
//! - `Point2D` - A 2D point
//! - `Point3D` - A 3D point
//! - `Vector2D` - A 2D vector with standard operations
//! - `Vector3D` - A 3D vector with standard operations
//! - `Color` - RGB/RGBA color representation
//!
//! ### Vector Operations
//! - `dot_product`, `cross_product` - Standard vector products
//! - `reflect`, `refract` - Reflection and refraction vectors
//! - `lerp`, `slerp` - Linear and spherical interpolation
//!
//! ### Transformation Matrices
//! - Translation, scaling, rotation (X, Y, Z, arbitrary axis)
//! - Shearing transformations
//! - Perspective and orthographic projections
//! - Look-at view matrix
//!
//! ### Curves
//! - Bezier curves (quadratic, cubic, general)
//! - B-spline curves
//! - Catmull-Rom splines
//!
//! ### Ray Tracing Primitives
//! - Ray-sphere intersection
//! - Ray-plane intersection
//! - Ray-triangle intersection with reusable geometry and parallel ray batches
//!
//! ### Quaternions
//! - Quaternion representation and operations
//! - Rotation using quaternions
//!
//! ## Examples
//!
//! ```
//! use rssn::kernels::computer_graphics::*;
//!
//! let v1 = Vector3D {
//!     x: 1.0,
//!     y: 0.0,
//!     z: 0.0,
//! };
//!
//! let v2 = Vector3D {
//!     x: 0.0,
//!     y: 1.0,
//!     z: 0.0,
//! };
//!
//! let cross = cross_product(&v1, &v2);
//!
//! assert_eq!(cross.z, 1.0);
//! ```

use std::ops::Add;
use std::ops::Div;
use std::ops::Mul;
use std::ops::Neg;
use std::ops::Sub;

use num_rational::BigRational;
use num_traits::One;
use num_traits::Signed;
use num_traits::ToPrimitive;
use num_traits::Zero;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::matrix::Matrix;

// ============================================================================
// Basic Types: Points and Vectors
// ============================================================================

/// A 2D point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Point2D {
    /// The x-coordinate of the point.
    pub x: f64,
    /// The y-coordinate of the point.
    pub y: f64,
}

impl Point2D {
    /// Creates a new 2D point.
    #[must_use]
    pub const fn new(
        x: f64,
        y: f64,
    ) -> Self {
        Self { x, y }
    }

    /// Distance to another point.
    #[must_use]
    pub fn distance_to(
        &self,
        other: &Self,
    ) -> f64 {
        let dx = self.x - other.x;

        let dy = self.y - other.y;

        dx.hypot(dy)
    }
}

/// A 3D point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Point3D {
    /// The x-coordinate of the point.
    pub x: f64,
    /// The y-coordinate of the point.
    pub y: f64,
    /// The z-coordinate of the point.
    pub z: f64,
}

impl Point3D {
    /// Creates a new 3D point.
    #[must_use]
    pub const fn new(
        x: f64,
        y: f64,
        z: f64,
    ) -> Self {
        Self { x, y, z }
    }

    /// Distance to another point.
    #[must_use]
    pub fn distance_to(
        &self,
        other: &Self,
    ) -> f64 {
        let dx = self.x - other.x;

        let dy = self.y - other.y;

        let dz = self.z - other.z;

        dz.mul_add(dz, dx.mul_add(dx, dy * dy)).sqrt()
    }

    /// Returns the point as a vector from the origin.
    #[must_use]
    pub const fn to_vector(&self) -> Vector3D {
        Vector3D {
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }
}

/// A 2D vector.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Vector2D {
    /// The x-component of the vector.
    pub x: f64,
    /// The y-component of the vector.
    pub y: f64,
}

impl Vector2D {
    /// Creates a new 2D vector.
    #[must_use]
    pub const fn new(
        x: f64,
        y: f64,
    ) -> Self {
        Self { x, y }
    }

    /// Computes the magnitude (length) of the vector.
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        self.x.hypot(self.y)
    }

    /// Normalizes the vector to have a magnitude of 1.
    #[must_use]
    pub fn normalize(&self) -> Self {
        let mag = self.magnitude();

        if mag == 0.0 {
            *self
        } else {
            *self / mag
        }
    }

    /// Rotates the vector by an angle (in radians).
    #[must_use]
    pub fn rotate(
        &self,
        angle: f64,
    ) -> Self {
        let (s, c) = angle.sin_cos();

        Self {
            x: self.x.mul_add(c, -(self.y * s)),
            y: self.x.mul_add(s, self.y * c),
        }
    }

    /// Returns the perpendicular vector (90 degrees counter-clockwise).
    #[must_use]
    pub const fn perpendicular(&self) -> Self {
        Self {
            x: -self.y,
            y: self.x,
        }
    }
}

impl Add for Vector2D {
    type Output = Self;

    fn add(
        self,
        rhs: Self,
    ) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub for Vector2D {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Mul<f64> for Vector2D {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Div<f64> for Vector2D {
    type Output = Self;

    fn div(
        self,
        rhs: f64,
    ) -> Self {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}

impl Neg for Vector2D {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
        }
    }
}

/// A 3D vector.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Vector3D {
    /// The x-component of the vector.
    pub x: f64,
    /// The y-component of the vector.
    pub y: f64,
    /// The z-component of the vector.
    pub z: f64,
}

impl Vector3D {
    /// Creates a new 3D vector.
    #[must_use]
    pub const fn new(
        x: f64,
        y: f64,
        z: f64,
    ) -> Self {
        Self { x, y, z }
    }

    /// Computes the magnitude (length) of the vector.
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        self.z
            .mul_add(self.z, self.x.mul_add(self.x, self.y * self.y))
            .sqrt()
    }

    /// Returns the squared magnitude (avoids sqrt).
    #[must_use]
    pub fn magnitude_squared(&self) -> f64 {
        self.z
            .mul_add(self.z, self.x.mul_add(self.x, self.y * self.y))
    }

    /// Normalizes the vector to have a magnitude of 1.
    #[must_use]
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

    fn add(
        self,
        rhs: Self,
    ) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl Sub for Vector3D {
    type Output = Self;

    fn sub(
        self,
        rhs: Self,
    ) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Mul<f64> for Vector3D {
    type Output = Self;

    fn mul(
        self,
        rhs: f64,
    ) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl Div<f64> for Vector3D {
    type Output = Self;

    fn div(
        self,
        rhs: f64,
    ) -> Self {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

impl Neg for Vector3D {
    type Output = Self;

    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

// ============================================================================
// Color
// ============================================================================

/// RGBA color with values in [0, 1].
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Color {
    /// The red component (range [0, 1]).
    pub r: f64,
    /// The green component (range [0, 1]).
    pub g: f64,
    /// The blue component (range [0, 1]).
    pub b: f64,
    /// The alpha component (range [0, 1]).
    pub a: f64,
}

impl Color {
    /// Black color.
    pub const BLACK: Self = Self {
        r: 0.0,
        g: 0.0,
        b: 0.0,
        a: 1.0,
    };
    /// Blue color.
    pub const BLUE: Self = Self {
        r: 0.0,
        g: 0.0,
        b: 1.0,
        a: 1.0,
    };
    /// Green color.
    pub const GREEN: Self = Self {
        r: 0.0,
        g: 1.0,
        b: 0.0,
        a: 1.0,
    };
    /// Red color.
    pub const RED: Self = Self {
        r: 1.0,
        g: 0.0,
        b: 0.0,
        a: 1.0,
    };
    /// White color.
    pub const WHITE: Self = Self {
        r: 1.0,
        g: 1.0,
        b: 1.0,
        a: 1.0,
    };

    /// Creates a new RGBA color.
    #[must_use]
    pub const fn new(
        r: f64,
        g: f64,
        b: f64,
        a: f64,
    ) -> Self {
        Self { r, g, b, a }
    }

    /// Creates an RGB color with alpha = 1.
    #[must_use]
    pub const fn rgb(
        r: f64,
        g: f64,
        b: f64,
    ) -> Self {
        Self { r, g, b, a: 1.0 }
    }

    /// Clamps color values to [0, 1].
    #[must_use]
    pub const fn clamp(&self) -> Self {
        Self {
            r: self.r.clamp(0.0, 1.0),
            g: self.g.clamp(0.0, 1.0),
            b: self.b.clamp(0.0, 1.0),
            a: self.a.clamp(0.0, 1.0),
        }
    }

    /// Linearly interpolates between two colors.
    #[must_use]
    pub fn lerp(
        &self,
        other: &Self,
        t: f64,
    ) -> Self {
        Self {
            r: (other.r - self.r).mul_add(t, self.r),
            g: (other.g - self.g).mul_add(t, self.g),
            b: (other.b - self.b).mul_add(t, self.b),
            a: (other.a - self.a).mul_add(t, self.a),
        }
    }
}

// ============================================================================
// Vector Operations
// ============================================================================

/// Computes the dot product of two 2D vectors.
#[must_use]
pub fn dot_product_2d(
    v1: &Vector2D,
    v2: &Vector2D,
) -> f64 {
    v1.x.mul_add(v2.x, v1.y * v2.y)
}

/// Computes the dot product of two 3D vectors.
#[must_use]
pub fn dot_product(
    v1: &Vector3D,
    v2: &Vector3D,
) -> f64 {
    v1.z.mul_add(v2.z, v1.x.mul_add(v2.x, v1.y * v2.y))
}

/// Computes the dot product of two 3D vectors; the explicitly named 3D sibling of
/// [`dot_product_2d`] (identical to [`dot_product`]).
#[must_use]
pub fn dot_product_3d(
    v1: &Vector3D,
    v2: &Vector3D,
) -> f64 {
    dot_product(v1, v2)
}

/// Computes the dot product of two vectors of any dimension.
///
/// # Errors
/// Returns an error when the slices have different lengths.
pub fn dot_product_nd(
    v1: &[f64],
    v2: &[f64],
) -> Result<f64, String> {
    if v1.len() != v2.len() {
        return Err(format!(
            "dot product of vectors of different dimension ({} and {})",
            v1.len(),
            v2.len()
        ));
    }

    Ok(v1
        .iter()
        .zip(v2)
        .fold(0.0, |acc, (a, b)| a.mul_add(*b, acc)))
}

/// Inner product that works for 2D vectors, 3D vectors and (via [`dot_product_nd`])
/// slices alike, so that code can be written once for every dimension.
pub trait InnerProduct {
    /// The Euclidean inner product `self . other`.
    fn inner(
        &self,
        other: &Self,
    ) -> f64;
}

impl InnerProduct for Vector2D {
    fn inner(
        &self,
        other: &Self,
    ) -> f64 {
        dot_product_2d(self, other)
    }
}

impl InnerProduct for Vector3D {
    fn inner(
        &self,
        other: &Self,
    ) -> f64 {
        dot_product(self, other)
    }
}

/// Dimension-generic dot product: dispatches on the vector type.
#[must_use]
pub fn dot<V: InnerProduct>(
    v1: &V,
    v2: &V,
) -> f64 {
    v1.inner(v2)
}

/// Computes the (scalar) cross product `x1 y2 - y1 x2` of two 2D vectors, the
/// `z` component of the 3D cross product of their embeddings.
#[must_use]
pub fn cross_product_2d(
    v1: &Vector2D,
    v2: &Vector2D,
) -> f64 {
    v1.x.mul_add(v2.y, -(v1.y * v2.x))
}

/// Computes the cross product of two 3D vectors.
#[must_use]
pub fn cross_product(
    v1: &Vector3D,
    v2: &Vector3D,
) -> Vector3D {
    Vector3D {
        x: v1.y.mul_add(v2.z, -(v1.z * v2.y)),
        y: v1.z.mul_add(v2.x, -(v1.x * v2.z)),
        z: v1.x.mul_add(v2.y, -(v1.y * v2.x)),
    }
}

/// Computes the reflection of a vector about a normal.
///
/// R = I - 2(N·I)N, where I is the incident vector and N is the normal.
#[must_use]
pub fn reflect(
    incident: &Vector3D,
    normal: &Vector3D,
) -> Vector3D {
    let dot = dot_product(incident, normal);

    *incident - *normal * (2.0 * dot)
}

/// Computes the refraction of a vector through a surface.
///
/// Uses Snell's law with `eta` = n1/n2 (ratio of refractive indices).
/// Returns `None` for total internal reflection.
#[must_use]
pub fn refract(
    incident: &Vector3D,
    normal: &Vector3D,
    eta: f64,
) -> Option<Vector3D> {
    let i = incident.normalize();

    let n = normal.normalize();

    let cos_i = -dot_product(&i, &n);

    let sin2_t = eta * eta * cos_i.mul_add(-cos_i, 1.0);

    if sin2_t > 1.0 {
        return None; // Total internal reflection
    }

    let cos_t = (1.0 - sin2_t).sqrt();

    Some(i * eta + n * eta.mul_add(cos_i, -cos_t))
}

/// Linear interpolation between two vectors.
#[must_use]
pub fn lerp(
    v1: &Vector3D,
    v2: &Vector3D,
    t: f64,
) -> Vector3D {
    *v1 + (*v2 - *v1) * t
}

/// Spherical linear interpolation (slerp) between two vectors.
///
/// Useful for smooth rotation interpolation.
#[must_use]
pub fn slerp(
    v1: &Vector3D,
    v2: &Vector3D,
    t: f64,
) -> Vector3D {
    let dot = dot_product(v1, v2).clamp(-1.0, 1.0);

    let theta = dot.acos();

    if theta.abs() < 1e-10 {
        return lerp(v1, v2, t);
    }

    let sin_theta = theta.sin();

    let s1 = ((1.0 - t) * theta).sin() / sin_theta;

    let s2 = (t * theta).sin() / sin_theta;

    *v1 * s1 + *v2 * s2
}

/// Computes the angle between two vectors in radians.
#[must_use]
pub fn angle_between(
    v1: &Vector3D,
    v2: &Vector3D,
) -> f64 {
    let dot = dot_product(v1, v2);

    let mags = v1.magnitude() * v2.magnitude();

    if mags == 0.0 {
        0.0
    } else {
        (dot / mags).clamp(-1.0, 1.0).acos()
    }
}

/// Projects vector `a` onto vector `b`.
#[must_use]
pub fn project(
    a: &Vector3D,
    b: &Vector3D,
) -> Vector3D {
    let b_mag_sq = b.magnitude_squared();

    if b_mag_sq == 0.0 {
        Vector3D::new(0.0, 0.0, 0.0)
    } else {
        *b * (dot_product(a, b) / b_mag_sq)
    }
}

// ============================================================================
// Transformation Matrices
// ============================================================================

/// Generates a 4x4 translation matrix.
#[must_use]
pub fn translation_matrix(
    dx: f64,
    dy: f64,
    dz: f64,
) -> Matrix<f64> {
    Matrix::new(
        4,
        4,
        vec![
            1.0, 0.0, 0.0, dx, 0.0, 1.0, 0.0, dy, 0.0, 0.0, 1.0, dz, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 scaling matrix.
#[must_use]
pub fn scaling_matrix(
    sx: f64,
    sy: f64,
    sz: f64,
) -> Matrix<f64> {
    Matrix::new(
        4,
        4,
        vec![
            sx, 0.0, 0.0, 0.0, 0.0, sy, 0.0, 0.0, 0.0, 0.0, sz, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 uniform scaling matrix.
#[must_use]
pub fn uniform_scaling_matrix(scale: f64) -> Matrix<f64> {
    scaling_matrix(scale, scale, scale)
}

/// Generates a 4x4 rotation matrix around the X-axis.
#[must_use]
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

/// Generates a 4x4 rotation matrix around the Y-axis.
#[must_use]
pub fn rotation_matrix_y(angle_rad: f64) -> Matrix<f64> {
    let (s, c) = angle_rad.sin_cos();

    Matrix::new(
        4,
        4,
        vec![
            c, 0.0, s, 0.0, 0.0, 1.0, 0.0, 0.0, -s, 0.0, c, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 rotation matrix around the Z-axis.
#[must_use]
pub fn rotation_matrix_z(angle_rad: f64) -> Matrix<f64> {
    let (s, c) = angle_rad.sin_cos();

    Matrix::new(
        4,
        4,
        vec![
            c, -s, 0.0, 0.0, s, c, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 rotation matrix around an arbitrary axis.
///
/// Uses the Rodrigues' rotation formula.
#[must_use]
pub fn rotation_matrix_axis(
    axis: &Vector3D,
    angle_rad: f64,
) -> Matrix<f64> {
    let n = axis.normalize();

    let (s, c) = angle_rad.sin_cos();

    let t = 1.0 - c;

    Matrix::new(
        4,
        4,
        vec![
            (t * n.x).mul_add(n.x, c),
            (t * n.x).mul_add(n.y, -(s * n.z)),
            (t * n.x).mul_add(n.z, s * n.y),
            0.0,
            (t * n.x).mul_add(n.y, s * n.z),
            (t * n.y).mul_add(n.y, c),
            (t * n.y).mul_add(n.z, -(s * n.x)),
            0.0,
            (t * n.x).mul_add(n.z, -(s * n.y)),
            (t * n.y).mul_add(n.z, s * n.x),
            (t * n.z).mul_add(n.z, c),
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        ],
    )
}

/// Generates a 4x4 shearing matrix.
///
/// # Arguments
/// * `xy` - Shear X by Y
/// * `xz` - Shear X by Z
/// * `yx` - Shear Y by X
/// * `yz` - Shear Y by Z
/// * `zx` - Shear Z by X
/// * `zy` - Shear Z by Y
#[must_use]
pub fn shearing_matrix(
    xy: f64,
    xz: f64,
    yx: f64,
    yz: f64,
    zx: f64,
    zy: f64,
) -> Matrix<f64> {
    Matrix::new(
        4,
        4,
        vec![
            1.0, xy, xz, 0.0, yx, 1.0, yz, 0.0, zx, zy, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

// ============================================================================
// Projection Matrices
// ============================================================================

/// Generates a 4x4 perspective projection matrix.
///
/// # Arguments
/// * `fov_y` - Field of view in the y-direction (radians).
/// * `aspect` - Aspect ratio (width / height).
/// * `near` - Distance to the near clipping plane.
/// * `far` - Distance to the far clipping plane.
#[must_use]
pub fn perspective_matrix(
    fov_y: f64,
    aspect: f64,
    near: f64,
    far: f64,
) -> Matrix<f64> {
    let f = 1.0 / (fov_y / 2.0).tan();

    let nf = 1.0 / (near - far);

    Matrix::new(
        4,
        4,
        vec![
            f / aspect,
            0.0,
            0.0,
            0.0,
            0.0,
            f,
            0.0,
            0.0,
            0.0,
            0.0,
            (far + near) * nf,
            2.0 * far * near * nf,
            0.0,
            0.0,
            -1.0,
            0.0,
        ],
    )
}

/// Generates a 4x4 orthographic projection matrix.
#[must_use]
pub fn orthographic_matrix(
    left: f64,
    right: f64,
    bottom: f64,
    top: f64,
    near: f64,
    far: f64,
) -> Matrix<f64> {
    let rl = 1.0 / (right - left);

    let tb = 1.0 / (top - bottom);

    let fn_ = 1.0 / (far - near);

    Matrix::new(
        4,
        4,
        vec![
            2.0 * rl,
            0.0,
            0.0,
            -(right + left) * rl,
            0.0,
            2.0 * tb,
            0.0,
            -(top + bottom) * tb,
            0.0,
            0.0,
            -2.0 * fn_,
            -(far + near) * fn_,
            0.0,
            0.0,
            0.0,
            1.0,
        ],
    )
}

/// Generates a 4x4 look-at view matrix.
///
/// # Arguments
/// * `eye` - Camera position.
/// * `center` - Point the camera is looking at.
/// * `up` - Up direction vector.
#[must_use]
pub fn look_at_matrix(
    eye: &Vector3D,
    center: &Vector3D,
    up: &Vector3D,
) -> Matrix<f64> {
    let f = (*center - *eye).normalize();

    let s = cross_product(&f, up).normalize();

    let u = cross_product(&s, &f);

    Matrix::new(
        4,
        4,
        vec![
            s.x,
            s.y,
            s.z,
            -dot_product(&s, eye),
            u.x,
            u.y,
            u.z,
            -dot_product(&u, eye),
            -f.x,
            -f.y,
            -f.z,
            dot_product(&f, eye),
            0.0,
            0.0,
            0.0,
            1.0,
        ],
    )
}

/// Identity 4x4 matrix.
#[must_use]
pub fn identity_matrix() -> Matrix<f64> {
    Matrix::new(
        4,
        4,
        vec![
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

// ============================================================================
// Quaternions
// ============================================================================

/// A quaternion for 3D rotations.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Quaternion {
    /// The real component.
    pub w: f64,
    /// The i component.
    pub x: f64,
    /// The j component.
    pub y: f64,
    /// The k component.
    pub z: f64,
}

impl Quaternion {
    /// Creates a new quaternion.
    #[must_use]
    pub const fn new(
        w: f64,
        x: f64,
        y: f64,
        z: f64,
    ) -> Self {
        Self { w, x, y, z }
    }

    /// Creates an identity quaternion (no rotation).
    #[must_use]
    pub const fn identity() -> Self {
        Self {
            w: 1.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Creates a quaternion from an axis-angle representation.
    #[must_use]
    pub fn from_axis_angle(
        axis: &Vector3D,
        angle: f64,
    ) -> Self {
        let half_angle = angle / 2.0;

        let (s, c) = half_angle.sin_cos();

        let n = axis.normalize();

        Self {
            w: c,
            x: n.x * s,
            y: n.y * s,
            z: n.z * s,
        }
    }

    /// Creates a quaternion from Euler angles (ZYX order).
    #[must_use]
    pub fn from_euler(
        roll: f64,
        pitch: f64,
        yaw: f64,
    ) -> Self {
        let (sr, cr) = (roll / 2.0).sin_cos();

        let (sp, cp) = (pitch / 2.0).sin_cos();

        let (sy, cy) = (yaw / 2.0).sin_cos();

        Self {
            w: (cr * cp).mul_add(cy, sr * sp * sy),
            x: (sr * cp).mul_add(cy, -(cr * sp * sy)),
            y: (cr * sp).mul_add(cy, sr * cp * sy),
            z: (cr * cp).mul_add(sy, -(sr * sp * cy)),
        }
    }

    /// Returns the magnitude of the quaternion.
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        self.z
            .mul_add(
                self.z,
                self.y
                    .mul_add(self.y, self.w.mul_add(self.w, self.x * self.x)),
            )
            .sqrt()
    }

    /// Normalizes the quaternion.
    #[must_use]
    pub fn normalize(&self) -> Self {
        let mag = self.magnitude();

        if mag == 0.0 {
            *self
        } else {
            Self {
                w: self.w / mag,
                x: self.x / mag,
                y: self.y / mag,
                z: self.z / mag,
            }
        }
    }

    /// Returns the conjugate of the quaternion.
    #[must_use]
    pub const fn conjugate(&self) -> Self {
        Self {
            w: self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    /// Returns the inverse of the quaternion.
    #[must_use]
    pub fn inverse(&self) -> Self {
        let mag_sq = self.z.mul_add(
            self.z,
            self.y
                .mul_add(self.y, self.w.mul_add(self.w, self.x * self.x)),
        );

        let conj = self.conjugate();

        Self {
            w: conj.w / mag_sq,
            x: conj.x / mag_sq,
            y: conj.y / mag_sq,
            z: conj.z / mag_sq,
        }
    }

    /// Multiplies two quaternions.
    #[must_use]
    pub fn multiply(
        &self,
        other: &Self,
    ) -> Self {
        Self {
            w: self.z.mul_add(
                -other.z,
                self.y
                    .mul_add(-other.y, self.w.mul_add(other.w, -(self.x * other.x))),
            ),
            x: self.z.mul_add(
                -other.y,
                self.y
                    .mul_add(other.z, self.w.mul_add(other.x, self.x * other.w)),
            ),
            y: self.z.mul_add(
                other.x,
                self.y
                    .mul_add(other.w, self.w.mul_add(other.y, -(self.x * other.z))),
            ),
            z: self.z.mul_add(
                other.w,
                self.y
                    .mul_add(-other.x, self.w.mul_add(other.z, self.x * other.y)),
            ),
        }
    }

    /// Rotates a vector by this quaternion.
    #[must_use]
    pub fn rotate_vector(
        &self,
        v: &Vector3D,
    ) -> Vector3D {
        let q_v = Self::new(0.0, v.x, v.y, v.z);

        let result = self.multiply(&q_v).multiply(&self.conjugate());

        Vector3D::new(result.x, result.y, result.z)
    }

    /// Converts the quaternion to a 4x4 rotation matrix.
    #[must_use]
    pub fn to_matrix(&self) -> Matrix<f64> {
        let q = self.normalize();

        let xx = q.x * q.x;

        let yy = q.y * q.y;

        let zz = q.z * q.z;

        let xy = q.x * q.y;

        let xz = q.x * q.z;

        let yz = q.y * q.z;

        let wx = q.w * q.x;

        let wy = q.w * q.y;

        let wz = q.w * q.z;

        Matrix::new(
            4,
            4,
            vec![
                2.0f64.mul_add(-(yy + zz), 1.0),
                2.0 * (xy - wz),
                2.0 * (xz + wy),
                0.0,
                2.0 * (xy + wz),
                2.0f64.mul_add(-(xx + zz), 1.0),
                2.0 * (yz - wx),
                0.0,
                2.0 * (xz - wy),
                2.0 * (yz + wx),
                2.0f64.mul_add(-(xx + yy), 1.0),
                0.0,
                0.0,
                0.0,
                0.0,
                1.0,
            ],
        )
    }

    /// Spherical linear interpolation between two quaternions.
    #[must_use]
    pub fn slerp(
        &self,
        other: &Self,
        t: f64,
    ) -> Self {
        let dot = self.z.mul_add(
            other.z,
            self.y
                .mul_add(other.y, self.w.mul_add(other.w, self.x * other.x)),
        );

        // If dot is negative, negate one quaternion to take the shorter path
        let (q2, dot) = if dot < 0.0 {
            (Self::new(-other.w, -other.x, -other.y, -other.z), -dot)
        } else {
            (*other, dot)
        };

        // Use linear interpolation for very close quaternions
        if dot > 0.9995 {
            return Self::new(
                t.mul_add(q2.w - self.w, self.w),
                t.mul_add(q2.x - self.x, self.x),
                t.mul_add(q2.y - self.y, self.y),
                t.mul_add(q2.z - self.z, self.z),
            )
            .normalize();
        }

        let theta_0 = dot.acos();

        let theta = theta_0 * t;

        let sin_theta = theta.sin();

        let sin_theta_0 = theta_0.sin();

        // Standard slerp weight sin(theta_0 - theta) / sin(theta_0), since
        // q(t) = sin((1-t) theta_0)/sin(theta_0) q0 + sin(t theta_0)/sin(theta_0) q1.
        let s0 = (theta_0 - theta).sin() / sin_theta_0;

        let s1 = sin_theta / sin_theta_0;

        Self::new(
            s0.mul_add(self.w, s1 * q2.w),
            s0.mul_add(self.x, s1 * q2.x),
            s0.mul_add(self.y, s1 * q2.y),
            s0.mul_add(self.z, s1 * q2.z),
        )
    }
}

// ============================================================================
// Ray Tracing Primitives
// ============================================================================

/// A ray defined by an origin point and direction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Ray {
    /// The origin of the ray.
    pub origin: Point3D,
    /// The direction of the ray (should be normalized).
    pub direction: Vector3D,
}

impl Ray {
    /// Creates a new ray.
    #[must_use]
    pub const fn new(
        origin: Point3D,
        direction: Vector3D,
    ) -> Self {
        Self { origin, direction }
    }

    /// Returns the point at parameter t along the ray.
    #[must_use]
    pub const fn at(
        &self,
        t: f64,
    ) -> Point3D {
        Point3D {
            x: t.mul_add(self.direction.x, self.origin.x),
            y: t.mul_add(self.direction.y, self.origin.y),
            z: t.mul_add(self.direction.z, self.origin.z),
        }
    }
}

/// A sphere defined by center and radius.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Sphere {
    /// The center of the sphere.
    pub center: Point3D,
    /// The radius of the sphere.
    pub radius: f64,
}

impl Sphere {
    /// Creates a new sphere.
    #[must_use]
    pub const fn new(
        center: Point3D,
        radius: f64,
    ) -> Self {
        Self { center, radius }
    }
}

/// Result of a ray intersection.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Intersection {
    /// Parameter t where ray intersects.
    pub t: f64,
    /// Point of intersection.
    pub point: Point3D,
    /// Normal at intersection point.
    pub normal: Vector3D,
}

/// Computes ray-sphere intersection.
///
/// Returns the closest intersection with t > 0, or `None` if no intersection.
#[must_use]
pub fn ray_sphere_intersection(
    ray: &Ray,
    sphere: &Sphere,
) -> Option<Intersection> {
    let oc = Vector3D::new(
        ray.origin.x - sphere.center.x,
        ray.origin.y - sphere.center.y,
        ray.origin.z - sphere.center.z,
    );

    let a = dot_product(&ray.direction, &ray.direction);

    let b = 2.0 * dot_product(&oc, &ray.direction);

    let c = sphere.radius.mul_add(-sphere.radius, dot_product(&oc, &oc));

    let discriminant = b * b - 4.0 * a * c;

    if discriminant < 0.0 {
        return None;
    }

    let sqrt_d = discriminant.sqrt();

    let t1 = (-b - sqrt_d) / (2.0 * a);

    let t2 = (-b + sqrt_d) / (2.0 * a);

    let t = if t1 > 1e-6 {
        t1
    } else if t2 > 1e-6 {
        t2
    } else {
        return None;
    };

    let point = ray.at(t);

    let normal = Vector3D::new(
        (point.x - sphere.center.x) / sphere.radius,
        (point.y - sphere.center.y) / sphere.radius,
        (point.z - sphere.center.z) / sphere.radius,
    );

    Some(Intersection { t, point, normal })
}

/// A plane defined by a point and normal.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Plane {
    /// A point on the plane.
    pub point: Point3D,
    /// The normal vector of the plane (should be normalized).
    pub normal: Vector3D,
}

impl Plane {
    /// Creates a new plane.
    #[must_use]
    pub const fn new(
        point: Point3D,
        normal: Vector3D,
    ) -> Self {
        Self { point, normal }
    }
}

/// Computes ray-plane intersection.
#[must_use]
pub fn ray_plane_intersection(
    ray: &Ray,
    plane: &Plane,
) -> Option<Intersection> {
    let denom = dot_product(&plane.normal, &ray.direction);

    if denom.abs() < 1e-6 {
        return None; // Ray is parallel to plane
    }

    let p0_l0 = Vector3D::new(
        plane.point.x - ray.origin.x,
        plane.point.y - ray.origin.y,
        plane.point.z - ray.origin.z,
    );

    let t = dot_product(&p0_l0, &plane.normal) / denom;

    if t < 1e-6 {
        return None;
    }

    let point = ray.at(t);

    let normal = if denom < 0.0 {
        plane.normal
    } else {
        -plane.normal
    };

    Some(Intersection { t, point, normal })
}

const fn vector_is_finite(v: Vector3D) -> bool {
    v.x.is_finite() && v.y.is_finite() && v.z.is_finite()
}

const fn largest_component(v: Vector3D) -> f64 {
    v.x.abs().max(v.y.abs()).max(v.z.abs())
}

fn difference_of_products(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
) -> f64 {
    let product = c * d;
    a.mul_add(b, -product) + (-c).mul_add(d, product)
}

const fn project_to_triangle_plane(
    v: Vector3D,
    axis: usize,
) -> Vector2D {
    match axis {
        | 0 => Vector2D::new(v.y, v.z),
        | 1 => Vector2D::new(v.z, v.x),
        | _ => Vector2D::new(v.x, v.y),
    }
}

#[derive(Debug, Clone, Copy)]
enum TriangleProjection {
    Inverse([Vector2D; 2]),
    Scaled {
        edges: [Vector2D; 2],
        scales: [f64; 2],
        area: f64,
    },
}

impl TriangleProjection {
    fn new(
        edges: [Vector2D; 2],
        scales: [f64; 2],
        area: f64,
    ) -> Self {
        let inv = scales.map(|scale| (1.0 / area) / scale);
        let inverse = [
            Vector2D::new(edges[1].y, -edges[1].x),
            Vector2D::new(-edges[0].y, edges[0].x),
        ];
        let weights = [inverse[0] * inv[0], inverse[1] * inv[1]];
        let representable = inverse.iter().zip(weights).all(|(unit, weight)| {
            (unit.x == 0.0 || weight.x.is_normal()) && (unit.y == 0.0 || weight.y.is_normal())
        });
        if inv.iter().all(|v| v.is_normal()) && representable {
            Self::Inverse(weights)
        } else {
            Self::Scaled { edges, scales, area }
        }
    }

    fn coordinates(
        &self,
        point: Vector2D,
    ) -> [f64; 2] {
        let coordinate = |i: usize| {
            match self {
                | Self::Inverse(weights) => dot_product_2d(&point, &weights[i]),
                | Self::Scaled { edges, scales, area } => {
                    let value = if i == 0 {
                        difference_of_products(point.x, edges[1].y, point.y, edges[1].x)
                    } else {
                        difference_of_products(edges[0].x, point.y, edges[0].y, point.x)
                    };
                    value / scales[i] / area
                },
            }
        };
        [coordinate(0), coordinate(1)]
    }
}

fn exact_vector(v: Vector3D) -> Option<[BigRational; 3]> {
    Some([
        BigRational::from_float(v.x)?,
        BigRational::from_float(v.y)?,
        BigRational::from_float(v.z)?,
    ])
}

fn exact_triple(
    a: &[BigRational; 3],
    b: &[BigRational; 3],
    c: &[BigRational; 3],
) -> BigRational {
    &a[0] * (&b[1] * &c[2] - &b[2] * &c[1])
        + &a[1] * (&b[2] * &c[0] - &b[0] * &c[2])
        + &a[2] * (&b[0] * &c[1] - &b[1] * &c[0])
}

#[derive(Clone, Copy)]
struct TriangleRayProjection {
    plane_factor: f64,
    absolute_denominator: f64,
    uncertainty: f64,
    direction_scale: Option<f64>,
    projected_direction: Vector2D,
}

const fn valid_ray_range(
    min_t: f64,
    max_t: f64,
) -> bool {
    min_t.is_finite() && min_t >= 0.0 && max_t >= min_t
}

/// A triangle prepared for repeated ray intersections.
#[derive(Debug, Clone, Copy)]
pub struct PreparedTriangle {
    vertices: [Vector3D; 3],
    normal: Vector3D,
    axis: usize,
    area_scale: f64,
    projection: TriangleProjection,
}

impl PreparedTriangle {
    /// Prepares the triangle or returns `None` for invalid or collinear vertices.
    #[must_use]
    pub fn new(
        v0: &Point3D,
        v1: &Point3D,
        v2: &Point3D,
    ) -> Option<Self> {
        let vertices = [v0.to_vector(), v1.to_vector(), v2.to_vector()];
        let edge1 = vertices[1] - vertices[0];
        let edge2 = vertices[2] - vertices[0];
        let scales = [largest_component(edge1), largest_component(edge2)];
        if !vector_is_finite(edge1) || !vector_is_finite(edge2) || scales.contains(&0.0) {
            return None;
        }
        // Normalize each edge separately so widely differing edge lengths do
        // not make the cross product overflow or its smaller edge disappear.
        let edge1 = edge1 / scales[0];
        let edge2 = edge2 / scales[1];
        let cross = Vector3D::new(
            difference_of_products(edge1.y, edge2.z, edge1.z, edge2.y),
            difference_of_products(edge1.z, edge2.x, edge1.x, edge2.z),
            difference_of_products(edge1.x, edge2.y, edge1.y, edge2.x),
        );
        let dominant = largest_component(cross);
        if dominant == 0.0 {
            return None;
        }
        let (axis, projected_area) = if cross.x.abs() >= dominant {
            (0, cross.x)
        } else if cross.y.abs() >= dominant {
            (1, cross.y)
        } else {
            (2, cross.z)
        };
        let edges = [
            project_to_triangle_plane(edge1, axis),
            project_to_triangle_plane(edge2, axis),
        ];
        Some(Self {
            vertices,
            normal: (cross / dominant).normalize(),
            axis,
            area_scale: dominant.min(1.0),
            projection: TriangleProjection::new(edges, scales, projected_area),
        })
    }

    #[cold]
    fn exact_intersection(
        &self,
        ray: &Ray,
        min_t: f64,
        max_t: f64,
    ) -> Option<Intersection> {
        let anchor = exact_vector(self.vertices[0])?;
        let p1 = exact_vector(self.vertices[1])?;
        let p2 = exact_vector(self.vertices[2])?;
        let origin = exact_vector(ray.origin.to_vector())?;
        let direction = exact_vector(ray.direction)?;
        let first = std::array::from_fn(|i| &anchor[i] - &p1[i]);
        let second = std::array::from_fn(|i| &anchor[i] - &p2[i]);
        let rhs = std::array::from_fn(|i| &anchor[i] - &origin[i]);
        let det = exact_triple(&direction, &first, &second);
        if det.is_zero() {
            return None;
        }
        let u = exact_triple(&direction, &rhs, &second) / &det;
        let v = exact_triple(&direction, &first, &rhs) / &det;
        if u.is_negative() || v.is_negative() || u + v > BigRational::one() {
            return None;
        }
        let t = exact_triple(&rhs, &first, &second) / det;
        if t < BigRational::from_float(min_t)?
            || (max_t.is_finite() && t > BigRational::from_float(max_t)?)
        {
            return None;
        }
        self.hit_at(ray, t.to_f64()?)
    }

    fn hit_at(
        &self,
        ray: &Ray,
        t: f64,
    ) -> Option<Intersection> {
        let point = ray.at(t);
        vector_is_finite(point.to_vector()).then_some(Intersection {
            t,
            point,
            normal: self.normal,
        })
    }

    #[inline]
    fn ray_projection(
        &self,
        raw_direction: Vector3D,
    ) -> Option<TriangleRayProjection> {
        let size = largest_component(raw_direction);
        let (direction, direction_scale, angular_scale) = if (1e-150..=1e150).contains(&size) {
            (raw_direction, None, size)
        } else {
            if size == 0.0 || !size.is_finite() {
                return None;
            }
            (raw_direction / size, Some(size), 1.0)
        };
        let denominator = dot_product(&self.normal, &direction);
        let abs_den = denominator.abs();
        if abs_den > 16.0 * f64::EPSILON * angular_scale {
            Some(TriangleRayProjection {
                plane_factor: -1.0 / denominator,
                absolute_denominator: abs_den,
                uncertainty: 64.0 * f64::EPSILON * angular_scale,
                direction_scale,
                projected_direction: project_to_triangle_plane(direction, self.axis),
            })
        } else {
            None
        }
    }

    #[inline]
    fn intersect_projected<const FIXED_PROJECTION: bool>(
        &self,
        ray: &Ray,
        projection: TriangleRayProjection,
        min_t: f64,
        max_t: f64,
    ) -> Option<Intersection> {
        let offset = ray.origin.to_vector() - self.vertices[0];
        let travel = dot_product(&offset, &self.normal) * projection.plane_factor;
        let t = projection
            .direction_scale
            .map_or(travel, |scale| travel / scale);
        if !t.is_finite() {
            return None;
        }
        let near_limit = |bound: f64| {
            bound.is_finite()
                && (t - bound).abs() * projection.absolute_denominator
                    <= projection.uncertainty * t.abs().max(bound)
        };
        if near_limit(min_t) || near_limit(max_t) {
            return self.exact_intersection(ray, min_t, max_t);
        }
        if !(min_t..=max_t).contains(&t) {
            return None;
        }
        let o = project_to_triangle_plane(offset, self.axis);
        let d = projection.projected_direction;
        let point = if FIXED_PROJECTION {
            o
        } else {
            Vector2D::new(travel.mul_add(d.x, o.x), travel.mul_add(d.y, o.y))
        };
        let [u, v] = self.projection.coordinates(point);
        let sum = u + v;
        if !sum.is_finite() {
            return self.exact_intersection(ray, min_t, max_t);
        }
        let lowest = u.min(v).min(1.0 - sum);
        if lowest.abs() * projection.absolute_denominator * self.area_scale
            <= projection.uncertainty
        {
            return self.exact_intersection(ray, min_t, max_t);
        }
        if lowest < 0.0 {
            return None;
        }
        self.hit_at(ray, t)
    }

    /// Intersects a ray over the inclusive parameter interval `[min_t, max_t]`.
    ///
    /// The lower bound is finite and nonnegative. The upper bound may be
    /// positive infinity. The returned normal follows the vertex winding.
    /// Zero directions and rays numerically parallel to the plane have no hit.
    /// Near a triangle edge or range endpoint, the signs are checked over
    /// the exact binary-float inputs.
    #[must_use]
    #[inline]
    pub fn intersect(
        &self,
        ray: &Ray,
        min_t: f64,
        max_t: f64,
    ) -> Option<Intersection> {
        if !valid_ray_range(min_t, max_t) {
            return None;
        }
        self.intersect_projected::<false>(ray, self.ray_projection(ray.direction)?, min_t, max_t)
    }

    /// Intersects parallel rays with a common direction and parameter range.
    ///
    /// Results follow the order of `origins`, using the same bounds as
    /// [`Self::intersect`]. Rays along the projection axis are culled against
    /// the triangle's projected bounding box. An invalid shared setting
    /// gives no hits.
    #[must_use]
    pub fn intersect_parallel(
        &self,
        origins: &[Point3D],
        direction: Vector3D,
        min_t: f64,
        max_t: f64,
    ) -> Vec<Option<Intersection>> {
        let projection = valid_ray_range(min_t, max_t)
            .then(|| self.ray_projection(direction))
            .flatten();
        let Some(projection) = projection else {
            return vec![None; origins.len()];
        };
        let projected_direction = project_to_triangle_plane(direction, self.axis);
        if projected_direction.x == 0.0 && projected_direction.y == 0.0 {
            // Along the projection axis, the ray cannot move into the
            // triangle's projected bounding box from outside it.
            let [a, b, c] = self
                .vertices
                .map(|v| project_to_triangle_plane(v, self.axis));
            let x = a.x.min(b.x).min(c.x)..=a.x.max(b.x).max(c.x);
            let y = a.y.min(b.y).min(c.y)..=a.y.max(b.y).max(c.y);
            return origins
                .iter()
                .map(|&origin| {
                    let point = project_to_triangle_plane(origin.to_vector(), self.axis);
                    if x.contains(&point.x) && y.contains(&point.y) {
                        self.intersect_projected::<true>(
                            &Ray::new(origin, direction),
                            projection,
                            min_t,
                            max_t,
                        )
                    } else {
                        None
                    }
                })
                .collect();
        }
        origins
            .iter()
            .map(|&origin| {
                self.intersect_projected::<false>(
                    &Ray::new(origin, direction),
                    projection,
                    min_t,
                    max_t,
                )
            })
            .collect()
    }
}

enum TriangleCandidate {
    Hit(Intersection),
    Miss,
    Prepare,
}

// Ordinary one-shot queries use the Möller–Trumbore determinants directly.
// Exceptional scales and uncertain signs share the prepared path.
#[inline]
fn direct_triangle_candidate(
    ray: &Ray,
    vertices: [&Point3D; 3],
    min_t: f64,
) -> TriangleCandidate {
    let [v0, v1, v2] = vertices.map(Point3D::to_vector);
    let first = v1 - v0;
    let second = v2 - v0;
    let across = cross_product(&ray.direction, &second);
    let det = dot_product(&first, &across);
    if !(1e-6..=1e150).contains(&det.abs()) {
        return TriangleCandidate::Prepare;
    }
    let inverse = 1.0 / det;
    let offset = ray.origin.to_vector() - v0;
    let u = inverse * dot_product(&offset, &across);
    let uncertainty = 64.0 * f64::EPSILON;
    if u.is_finite() && (u < -uncertainty || u > 1.0 + uncertainty) {
        return TriangleCandidate::Miss;
    }
    let across = cross_product(&offset, &first);
    let v = inverse * dot_product(&ray.direction, &across);
    let sum = u + v;
    if !sum.is_finite() {
        return TriangleCandidate::Prepare;
    }
    let lowest = u.min(v).min(1.0 - sum);
    if lowest.abs() <= uncertainty {
        return TriangleCandidate::Prepare;
    }
    if lowest < 0.0 {
        return TriangleCandidate::Miss;
    }
    let t = inverse * dot_product(&second, &across);
    if !t.is_finite() || (t - min_t).abs() <= uncertainty * t.abs().max(min_t) {
        return TriangleCandidate::Prepare;
    }
    if t < min_t {
        return TriangleCandidate::Miss;
    }
    let normal = cross_product(&first, &second);
    let magnitude = normal.magnitude();
    if !magnitude.is_normal()
        || det.abs() <= 16.0 * f64::EPSILON * magnitude * largest_component(ray.direction)
    {
        return TriangleCandidate::Prepare;
    }
    let point = ray.at(t);
    if !vector_is_finite(point.to_vector()) {
        return TriangleCandidate::Prepare;
    }
    TriangleCandidate::Hit(Intersection {
        t,
        point,
        normal: normal / magnitude,
    })
}

/// Intersects a ray with a triangle using a minimum ray parameter of `1e-6`.
///
/// Use [`PreparedTriangle`] to reuse triangle geometry or specify another
/// parameter range. The returned normal follows the vertex winding.
#[must_use]
pub fn ray_triangle_intersection(
    ray: &Ray,
    v0: &Point3D,
    v1: &Point3D,
    v2: &Point3D,
) -> Option<Intersection> {
    const MIN_T: f64 = 1e-6;
    match direct_triangle_candidate(ray, [v0, v1, v2], MIN_T) {
        | TriangleCandidate::Hit(hit) => Some(hit),
        | TriangleCandidate::Miss => None,
        | TriangleCandidate::Prepare => {
            PreparedTriangle::new(v0, v1, v2)?.intersect(ray, MIN_T, f64::INFINITY)
        },
    }
}

// ============================================================================
// Curves
// ============================================================================

/// Evaluates a quadratic Bezier curve at parameter t.
///
/// # Arguments
/// * `p0`, `p1`, `p2` - Control points
/// * `t` - Parameter in [0, 1]
#[must_use]
pub fn bezier_quadratic(
    p0: &Point3D,
    p1: &Point3D,
    p2: &Point3D,
    t: f64,
) -> Point3D {
    let t2 = t * t;

    let mt = 1.0 - t;

    let mt2 = mt * mt;

    Point3D {
        x: t2.mul_add(p2.x, mt2 * p0.x + 2.0 * mt * t * p1.x),
        y: t2.mul_add(p2.y, mt2 * p0.y + 2.0 * mt * t * p1.y),
        z: t2.mul_add(p2.z, mt2 * p0.z + 2.0 * mt * t * p1.z),
    }
}

/// Evaluates a cubic Bezier curve at parameter t.
///
/// # Arguments
/// * `p0`, `p1`, `p2`, `p3` - Control points
/// * `t` - Parameter in [0, 1]
#[must_use]
pub fn bezier_cubic(
    p0: &Point3D,
    p1: &Point3D,
    p2: &Point3D,
    p3: &Point3D,
    t: f64,
) -> Point3D {
    let t2 = t * t;

    let t3 = t2 * t;

    let mt = 1.0 - t;

    let mt2 = mt * mt;

    let mt3 = mt2 * mt;

    Point3D {
        x: t3.mul_add(
            p3.x,
            (3.0 * mt * t2).mul_add(p2.x, mt3 * p0.x + 3.0 * mt2 * t * p1.x),
        ),
        y: t3.mul_add(
            p3.y,
            (3.0 * mt * t2).mul_add(p2.y, mt3 * p0.y + 3.0 * mt2 * t * p1.y),
        ),
        z: t3.mul_add(
            p3.z,
            (3.0 * mt * t2).mul_add(p2.z, mt3 * p0.z + 3.0 * mt2 * t * p1.z),
        ),
    }
}

/// Evaluates a Catmull-Rom spline at parameter t.
///
/// The spline passes through p1 and p2, with p0 and p3 as guide points.
#[must_use]
pub fn catmull_rom(
    p0: &Point3D,
    p1: &Point3D,
    p2: &Point3D,
    p3: &Point3D,
    t: f64,
) -> Point3D {
    let t2 = t * t;

    let t3 = t2 * t;

    Point3D {
        x: f64::midpoint(
            (4.0f64.mul_add(p2.x, 2.0f64.mul_add(p0.x, -(5.0 * p1.x))) - p3.x)
                .mul_add(t2, 2.0f64.mul_add(p1.x, (-p0.x + p2.x) * t)),
            t3.mul_add(
                3.0f64.mul_add(-p2.x, 3.0f64.mul_add(p1.x, -p0.x)) + p3.x,
                0.0,
            ),
        ),
        y: f64::midpoint(
            (4.0f64.mul_add(p2.y, 2.0f64.mul_add(p0.y, -(5.0 * p1.y))) - p3.y)
                .mul_add(t2, 2.0f64.mul_add(p1.y, (-p0.y + p2.y) * t)),
            t3.mul_add(
                3.0f64.mul_add(-p2.y, 3.0f64.mul_add(p1.y, -p0.y)) + p3.y,
                0.0,
            ),
        ),
        z: f64::midpoint(
            (4.0f64.mul_add(p2.z, 2.0f64.mul_add(p0.z, -(5.0 * p1.z))) - p3.z)
                .mul_add(t2, 2.0f64.mul_add(p1.z, (-p0.z + p2.z) * t)),
            t3.mul_add(
                3.0f64.mul_add(-p2.z, 3.0f64.mul_add(p1.z, -p0.z)) + p3.z,
                0.0,
            ),
        ),
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Converts degrees to radians.
#[must_use]
pub const fn degrees_to_radians(degrees: f64) -> f64 {
    degrees.to_radians()
}

/// Converts radians to degrees.
#[must_use]
pub const fn radians_to_degrees(radians: f64) -> f64 {
    radians.to_degrees()
}

/// Applies a 4x4 transformation matrix to a 3D point.
#[must_use]
pub fn transform_point(
    matrix: &Matrix<f64>,
    point: &Point3D,
) -> Point3D {
    let x = matrix.get(0, 0) * point.x
        + matrix.get(0, 1) * point.y
        + matrix.get(0, 2) * point.z
        + matrix.get(0, 3);

    let y = matrix.get(1, 0) * point.x
        + matrix.get(1, 1) * point.y
        + matrix.get(1, 2) * point.z
        + matrix.get(1, 3);

    let z = matrix.get(2, 0) * point.x
        + matrix.get(2, 1) * point.y
        + matrix.get(2, 2) * point.z
        + matrix.get(2, 3);

    let w = matrix.get(3, 0) * point.x
        + matrix.get(3, 1) * point.y
        + matrix.get(3, 2) * point.z
        + matrix.get(3, 3);

    if w.abs() > 1e-10 && (w - 1.0).abs() > f64::EPSILON {
        Point3D::new(x / w, y / w, z / w)
    } else {
        Point3D::new(x, y, z)
    }
}

/// Applies a 4x4 transformation matrix to a 3D vector (ignores translation).
#[must_use]
pub fn transform_vector(
    matrix: &Matrix<f64>,
    vector: &Vector3D,
) -> Vector3D {
    let x = matrix.get(0, 0) * vector.x + matrix.get(0, 1) * vector.y + matrix.get(0, 2) * vector.z;

    let y = matrix.get(1, 0) * vector.x + matrix.get(1, 1) * vector.y + matrix.get(1, 2) * vector.z;

    let z = matrix.get(2, 0) * vector.x + matrix.get(2, 1) * vector.y + matrix.get(2, 2) * vector.z;

    Vector3D::new(x, y, z)
}

/// Computes the barycentric coordinates of a point in a triangle.
///
/// Returns (u, v, w) where point = u*v0 + v*v1 + w*v2.
#[must_use]
pub fn barycentric_coordinates(
    point: &Point3D,
    v0: &Point3D,
    v1: &Point3D,
    v2: &Point3D,
) -> (f64, f64, f64) {
    let v0v1 = Vector3D::new(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);

    let v0v2 = Vector3D::new(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);

    let v0p = Vector3D::new(point.x - v0.x, point.y - v0.y, point.z - v0.z);

    let d00 = dot_product(&v0v1, &v0v1);

    let d01 = dot_product(&v0v1, &v0v2);

    let d11 = dot_product(&v0v2, &v0v2);

    let d20 = dot_product(&v0p, &v0v1);

    let d21 = dot_product(&v0p, &v0v2);

    let denom = d00.mul_add(d11, -(d01 * d01));

    if denom.abs() < 1e-10 {
        return (1.0, 0.0, 0.0);
    }

    let v = d11.mul_add(d20, -(d01 * d21)) / denom;

    let w = d00.mul_add(d21, -(d01 * d20)) / denom;

    let u = 1.0 - v - w;

    (u, v, w)
}
