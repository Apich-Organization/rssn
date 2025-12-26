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
//! - Ray-triangle intersection
//!
//! ### Quaternions
//! - Quaternion representation and operations
//! - Rotation using quaternions
//!
//! ## Examples
//!
//! ```
//! 
//! use rssn::numerical::computer_graphics::*;
//!
//! let v1 = Vector3D {
//!     x : 1.0,
//!     y : 0.0,
//!     z : 0.0,
//! };
//!
//! let v2 = Vector3D {
//!     x : 0.0,
//!     y : 1.0,
//!     z : 0.0,
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

use serde::Deserialize;
use serde::Serialize;

use crate::numerical::matrix::Matrix;

// ============================================================================
// Basic Types: Points and Vectors
// ============================================================================

/// A 2D point.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Point2D {
    pub x : f64,
    pub y : f64,
}

impl Point2D {
    /// Creates a new 2D point.
    #[must_use]

    pub const fn new(
        x : f64,
        y : f64,
    ) -> Self {

        Self {
            x,
            y,
        }
    }

    /// Distance to another point.
    #[must_use]

    pub fn distance_to(
        &self,
        other : &Self,
    ) -> f64 {

        let dx = self.x - other.x;

        let dy = self.y - other.y;

        dx.hypot(dy)
    }
}

/// A 3D point.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Point3D {
    pub x : f64,
    pub y : f64,
    pub z : f64,
}

impl Point3D {
    /// Creates a new 3D point.
    #[must_use]

    pub const fn new(
        x : f64,
        y : f64,
        z : f64,
    ) -> Self {

        Self {
            x,
            y,
            z,
        }
    }

    /// Distance to another point.
    #[must_use]

    pub fn distance_to(
        &self,
        other : &Self,
    ) -> f64 {

        let dx = self.x - other.x;

        let dy = self.y - other.y;

        let dz = self.z - other.z;

        (dx.mul_add(dx, dy * dy)
            + dz * dz)
            .sqrt()
    }

    /// Returns the point as a vector from the origin.
    #[must_use]

    pub const fn to_vector(
        &self
    ) -> Vector3D {

        Vector3D {
            x : self.x,
            y : self.y,
            z : self.z,
        }
    }
}

/// A 2D vector.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Vector2D {
    pub x : f64,
    pub y : f64,
}

impl Vector2D {
    /// Creates a new 2D vector.
    #[must_use]

    pub const fn new(
        x : f64,
        y : f64,
    ) -> Self {

        Self {
            x,
            y,
        }
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
        angle : f64,
    ) -> Self {

        let (s, c) = angle.sin_cos();

        Self {
            x : self.x.mul_add(
                c,
                -(self.y * s),
            ),
            y : self
                .x
                .mul_add(s, self.y * c),
        }
    }

    /// Returns the perpendicular vector (90 degrees counter-clockwise).
    #[must_use]

    pub const fn perpendicular(
        &self
    ) -> Self {

        Self {
            x : -self.y,
            y : self.x,
        }
    }
}

impl Add for Vector2D {
    type Output = Self;

    fn add(
        self,
        rhs : Self,
    ) -> Self {

        Self {
            x : self.x + rhs.x,
            y : self.y + rhs.y,
        }
    }
}

impl Sub for Vector2D {
    type Output = Self;

    fn sub(
        self,
        rhs : Self,
    ) -> Self {

        Self {
            x : self.x - rhs.x,
            y : self.y - rhs.y,
        }
    }
}

impl Mul<f64> for Vector2D {
    type Output = Self;

    fn mul(
        self,
        rhs : f64,
    ) -> Self {

        Self {
            x : self.x * rhs,
            y : self.y * rhs,
        }
    }
}

impl Div<f64> for Vector2D {
    type Output = Self;

    fn div(
        self,
        rhs : f64,
    ) -> Self {

        Self {
            x : self.x / rhs,
            y : self.y / rhs,
        }
    }
}

impl Neg for Vector2D {
    type Output = Self;

    fn neg(self) -> Self {

        Self {
            x : -self.x,
            y : -self.y,
        }
    }
}

/// A 3D vector.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Vector3D {
    pub x : f64,
    pub y : f64,
    pub z : f64,
}

impl Vector3D {
    /// Creates a new 3D vector.
    #[must_use]

    pub const fn new(
        x : f64,
        y : f64,
        z : f64,
    ) -> Self {

        Self {
            x,
            y,
            z,
        }
    }

    /// Computes the magnitude (length) of the vector.
    #[must_use]

    pub fn magnitude(&self) -> f64 {

        (self.x.mul_add(
            self.x,
            self.y * self.y,
        ) + self.z * self.z)
            .sqrt()
    }

    /// Returns the squared magnitude (avoids sqrt).
    #[must_use]

    pub fn magnitude_squared(
        &self
    ) -> f64 {

        self.x.mul_add(
            self.x,
            self.y * self.y,
        ) + self.z * self.z
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
        rhs : Self,
    ) -> Self {

        Self {
            x : self.x + rhs.x,
            y : self.y + rhs.y,
            z : self.z + rhs.z,
        }
    }
}

impl Sub for Vector3D {
    type Output = Self;

    fn sub(
        self,
        rhs : Self,
    ) -> Self {

        Self {
            x : self.x - rhs.x,
            y : self.y - rhs.y,
            z : self.z - rhs.z,
        }
    }
}

impl Mul<f64> for Vector3D {
    type Output = Self;

    fn mul(
        self,
        rhs : f64,
    ) -> Self {

        Self {
            x : self.x * rhs,
            y : self.y * rhs,
            z : self.z * rhs,
        }
    }
}

impl Div<f64> for Vector3D {
    type Output = Self;

    fn div(
        self,
        rhs : f64,
    ) -> Self {

        Self {
            x : self.x / rhs,
            y : self.y / rhs,
            z : self.z / rhs,
        }
    }
}

impl Neg for Vector3D {
    type Output = Self;

    fn neg(self) -> Self {

        Self {
            x : -self.x,
            y : -self.y,
            z : -self.z,
        }
    }
}

// ============================================================================
// Color
// ============================================================================

/// RGBA color with values in [0, 1].
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Color {
    pub r : f64,
    pub g : f64,
    pub b : f64,
    pub a : f64,
}

impl Color {
    /// Black color.

    pub const BLACK : Self = Self {
        r : 0.0,
        g : 0.0,
        b : 0.0,
        a : 1.0,
    };
    /// Blue color.

    pub const BLUE : Self = Self {
        r : 0.0,
        g : 0.0,
        b : 1.0,
        a : 1.0,
    };
    /// Green color.

    pub const GREEN : Self = Self {
        r : 0.0,
        g : 1.0,
        b : 0.0,
        a : 1.0,
    };
    /// Red color.

    pub const RED : Self = Self {
        r : 1.0,
        g : 0.0,
        b : 0.0,
        a : 1.0,
    };
    /// White color.

    pub const WHITE : Self = Self {
        r : 1.0,
        g : 1.0,
        b : 1.0,
        a : 1.0,
    };

    /// Creates a new RGBA color.
    #[must_use]

    pub const fn new(
        r : f64,
        g : f64,
        b : f64,
        a : f64,
    ) -> Self {

        Self {
            r,
            g,
            b,
            a,
        }
    }

    /// Creates an RGB color with alpha = 1.
    #[must_use]

    pub const fn rgb(
        r : f64,
        g : f64,
        b : f64,
    ) -> Self {

        Self {
            r,
            g,
            b,
            a : 1.0,
        }
    }

    /// Clamps color values to [0, 1].
    #[must_use]

    pub const fn clamp(&self) -> Self {

        Self {
            r : self
                .r
                .clamp(0.0, 1.0),
            g : self
                .g
                .clamp(0.0, 1.0),
            b : self
                .b
                .clamp(0.0, 1.0),
            a : self
                .a
                .clamp(0.0, 1.0),
        }
    }

    /// Linearly interpolates between two colors.
    #[must_use]

    pub fn lerp(
        &self,
        other : &Self,
        t : f64,
    ) -> Self {

        Self {
            r : (other.r - self.r)
                .mul_add(t, self.r),
            g : (other.g - self.g)
                .mul_add(t, self.g),
            b : (other.b - self.b)
                .mul_add(t, self.b),
            a : (other.a - self.a)
                .mul_add(t, self.a),
        }
    }
}

// ============================================================================
// Vector Operations
// ============================================================================

/// Computes the dot product of two 2D vectors.
#[must_use]

pub fn dot_product_2d(
    v1 : &Vector2D,
    v2 : &Vector2D,
) -> f64 {

    v1.x.mul_add(v2.x, v1.y * v2.y)
}

/// Computes the dot product of two 3D vectors.
#[must_use]

pub fn dot_product(
    v1 : &Vector3D,
    v2 : &Vector3D,
) -> f64 {

    v1.x.mul_add(v2.x, v1.y * v2.y)
        + v1.z * v2.z
}

/// Computes the cross product of two 3D vectors.
#[must_use]

pub fn cross_product(
    v1 : &Vector3D,
    v2 : &Vector3D,
) -> Vector3D {

    Vector3D {
        x : v1.y.mul_add(
            v2.z,
            -(v1.z * v2.y),
        ),
        y : v1.z.mul_add(
            v2.x,
            -(v1.x * v2.z),
        ),
        z : v1.x.mul_add(
            v2.y,
            -(v1.y * v2.x),
        ),
    }
}

/// Computes the reflection of a vector about a normal.
///
/// R = I - 2(N·I)N, where I is the incident vector and N is the normal.
#[must_use]

pub fn reflect(
    incident : &Vector3D,
    normal : &Vector3D,
) -> Vector3D {

    let dot =
        dot_product(incident, normal);

    *incident - *normal * (2.0 * dot)
}

/// Computes the refraction of a vector through a surface.
///
/// Uses Snell's law with `eta` = n1/n2 (ratio of refractive indices).
/// Returns `None` for total internal reflection.
#[must_use]

pub fn refract(
    incident : &Vector3D,
    normal : &Vector3D,
    eta : f64,
) -> Option<Vector3D> {

    let i = incident.normalize();

    let n = normal.normalize();

    let cos_i = -dot_product(&i, &n);

    let sin2_t = eta
        * eta
        * cos_i.mul_add(-cos_i, 1.0);

    if sin2_t > 1.0 {

        return None; // Total internal reflection
    }

    let cos_t = (1.0 - sin2_t).sqrt();

    Some(
        i * eta
            + n * eta
                .mul_add(cos_i, -cos_t),
    )
}

/// Linear interpolation between two vectors.
#[must_use]

pub fn lerp(
    v1 : &Vector3D,
    v2 : &Vector3D,
    t : f64,
) -> Vector3D {

    *v1 + (*v2 - *v1) * t
}

/// Spherical linear interpolation (slerp) between two vectors.
///
/// Useful for smooth rotation interpolation.
#[must_use]

pub fn slerp(
    v1 : &Vector3D,
    v2 : &Vector3D,
    t : f64,
) -> Vector3D {

    let dot = dot_product(v1, v2)
        .clamp(-1.0, 1.0);

    let theta = dot.acos();

    if theta.abs() < 1e-10 {

        return lerp(v1, v2, t);
    }

    let sin_theta = theta.sin();

    let s1 = ((1.0 - t) * theta).sin()
        / sin_theta;

    let s2 =
        (t * theta).sin() / sin_theta;

    *v1 * s1 + *v2 * s2
}

/// Computes the angle between two vectors in radians.
#[must_use]

pub fn angle_between(
    v1 : &Vector3D,
    v2 : &Vector3D,
) -> f64 {

    let dot = dot_product(v1, v2);

    let mags =
        v1.magnitude() * v2.magnitude();

    if mags == 0.0 {

        0.0
    } else {

        (dot / mags)
            .clamp(-1.0, 1.0)
            .acos()
    }
}

/// Projects vector `a` onto vector `b`.
#[must_use]

pub fn project(
    a : &Vector3D,
    b : &Vector3D,
) -> Vector3D {

    let b_mag_sq =
        b.magnitude_squared();

    if b_mag_sq == 0.0 {

        Vector3D::new(0.0, 0.0, 0.0)
    } else {

        *b * (dot_product(a, b)
            / b_mag_sq)
    }
}

// ============================================================================
// Transformation Matrices
// ============================================================================

/// Generates a 4x4 translation matrix.
#[must_use]

pub fn translation_matrix(
    dx : f64,
    dy : f64,
    dz : f64,
) -> Matrix<f64> {

    Matrix::new(
        4,
        4,
        vec![
            1.0, 0.0, 0.0, dx, 0.0,
            1.0, 0.0, dy, 0.0, 0.0,
            1.0, dz, 0.0, 0.0, 0.0,
            1.0,
        ],
    )
}

/// Generates a 4x4 scaling matrix.
#[must_use]

pub fn scaling_matrix(
    sx : f64,
    sy : f64,
    sz : f64,
) -> Matrix<f64> {

    Matrix::new(
        4,
        4,
        vec![
            sx, 0.0, 0.0, 0.0, 0.0, sy,
            0.0, 0.0, 0.0, 0.0, sz,
            0.0, 0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 uniform scaling matrix.
#[must_use]

pub fn uniform_scaling_matrix(
    scale : f64
) -> Matrix<f64> {

    scaling_matrix(scale, scale, scale)
}

/// Generates a 4x4 rotation matrix around the X-axis.
#[must_use]

pub fn rotation_matrix_x(
    angle_rad : f64
) -> Matrix<f64> {

    let (s, c) = angle_rad.sin_cos();

    Matrix::new(
        4,
        4,
        vec![
            1.0, 0.0, 0.0, 0.0, 0.0, c,
            -s, 0.0, 0.0, s, c, 0.0,
            0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 rotation matrix around the Y-axis.
#[must_use]

pub fn rotation_matrix_y(
    angle_rad : f64
) -> Matrix<f64> {

    let (s, c) = angle_rad.sin_cos();

    Matrix::new(
        4,
        4,
        vec![
            c, 0.0, s, 0.0, 0.0, 1.0,
            0.0, 0.0, -s, 0.0, c, 0.0,
            0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 rotation matrix around the Z-axis.
#[must_use]

pub fn rotation_matrix_z(
    angle_rad : f64
) -> Matrix<f64> {

    let (s, c) = angle_rad.sin_cos();

    Matrix::new(
        4,
        4,
        vec![
            c, -s, 0.0, 0.0, s, c, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        ],
    )
}

/// Generates a 4x4 rotation matrix around an arbitrary axis.
///
/// Uses the Rodrigues' rotation formula.
#[must_use]

pub fn rotation_matrix_axis(
    axis : &Vector3D,
    angle_rad : f64,
) -> Matrix<f64> {

    let n = axis.normalize();

    let (s, c) = angle_rad.sin_cos();

    let t = 1.0 - c;

    Matrix::new(
        4,
        4,
        vec![
            (t * n.x).mul_add(n.x, c),
            (t * n.x).mul_add(
                n.y,
                -(s * n.z),
            ),
            (t * n.x)
                .mul_add(n.z, s * n.y),
            0.0,
            (t * n.x)
                .mul_add(n.y, s * n.z),
            (t * n.y).mul_add(n.y, c),
            (t * n.y).mul_add(
                n.z,
                -(s * n.x),
            ),
            0.0,
            (t * n.x).mul_add(
                n.z,
                -(s * n.y),
            ),
            (t * n.y)
                .mul_add(n.z, s * n.x),
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
    xy : f64,
    xz : f64,
    yx : f64,
    yz : f64,
    zx : f64,
    zy : f64,
) -> Matrix<f64> {

    Matrix::new(
        4,
        4,
        vec![
            1.0, xy, xz, 0.0, yx, 1.0,
            yz, 0.0, zx, zy, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
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
    fov_y : f64,
    aspect : f64,
    near : f64,
    far : f64,
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
    left : f64,
    right : f64,
    bottom : f64,
    top : f64,
    near : f64,
    far : f64,
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
    eye : &Vector3D,
    center : &Vector3D,
    up : &Vector3D,
) -> Matrix<f64> {

    let f =
        (*center - *eye).normalize();

    let s = cross_product(&f, up)
        .normalize();

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

pub fn identity_matrix() -> Matrix<f64>
{

    Matrix::new(
        4,
        4,
        vec![
            1.0, 0.0, 0.0, 0.0, 0.0,
            1.0, 0.0, 0.0, 0.0, 0.0,
            1.0, 0.0, 0.0, 0.0, 0.0,
            1.0,
        ],
    )
}

// ============================================================================
// Quaternions
// ============================================================================

/// A quaternion for 3D rotations.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Quaternion {
    pub w : f64,
    pub x : f64,
    pub y : f64,
    pub z : f64,
}

impl Quaternion {
    /// Creates a new quaternion.
    #[must_use]

    pub const fn new(
        w : f64,
        x : f64,
        y : f64,
        z : f64,
    ) -> Self {

        Self {
            w,
            x,
            y,
            z,
        }
    }

    /// Creates an identity quaternion (no rotation).
    #[must_use]

    pub const fn identity() -> Self {

        Self {
            w : 1.0,
            x : 0.0,
            y : 0.0,
            z : 0.0,
        }
    }

    /// Creates a quaternion from an axis-angle representation.
    #[must_use]

    pub fn from_axis_angle(
        axis : &Vector3D,
        angle : f64,
    ) -> Self {

        let half_angle = angle / 2.0;

        let (s, c) =
            half_angle.sin_cos();

        let n = axis.normalize();

        Self {
            w : c,
            x : n.x * s,
            y : n.y * s,
            z : n.z * s,
        }
    }

    /// Creates a quaternion from Euler angles (ZYX order).
    #[must_use]

    pub fn from_euler(
        roll : f64,
        pitch : f64,
        yaw : f64,
    ) -> Self {

        let (sr, cr) =
            (roll / 2.0).sin_cos();

        let (sp, cp) =
            (pitch / 2.0).sin_cos();

        let (sy, cy) =
            (yaw / 2.0).sin_cos();

        Self {
            w : (cr * cp).mul_add(
                cy,
                sr * sp * sy,
            ),
            x : (sr * cp).mul_add(
                cy,
                -(cr * sp * sy),
            ),
            y : (cr * sp).mul_add(
                cy,
                sr * cp * sy,
            ),
            z : (cr * cp).mul_add(
                sy,
                -(sr * sp * cy),
            ),
        }
    }

    /// Returns the magnitude of the quaternion.
    #[must_use]

    pub fn magnitude(&self) -> f64 {

        (self.w.mul_add(
            self.w,
            self.x * self.x,
        ) + self.y * self.y
            + self.z * self.z)
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
                w : self.w / mag,
                x : self.x / mag,
                y : self.y / mag,
                z : self.z / mag,
            }
        }
    }

    /// Returns the conjugate of the quaternion.
    #[must_use]

    pub const fn conjugate(
        &self
    ) -> Self {

        Self {
            w : self.w,
            x : -self.x,
            y : -self.y,
            z : -self.z,
        }
    }

    /// Returns the inverse of the quaternion.
    #[must_use]

    pub fn inverse(&self) -> Self {

        let mag_sq = self.w.mul_add(
            self.w,
            self.x * self.x,
        ) + self.y
            * self.y
            + self.z * self.z;

        let conj = self.conjugate();

        Self {
            w : conj.w / mag_sq,
            x : conj.x / mag_sq,
            y : conj.y / mag_sq,
            z : conj.z / mag_sq,
        }
    }

    /// Multiplies two quaternions.
    #[must_use]

    pub fn multiply(
        &self,
        other : &Self,
    ) -> Self {

        Self {
            w : self.w.mul_add(
                other.w,
                -(self.x * other.x),
            ) - self.y * other.y
                - self.z * other.z,
            x : self.w.mul_add(
                other.x,
                self.x * other.w,
            ) + self.y * other.z
                - self.z * other.y,
            y : self.w.mul_add(
                other.y,
                -(self.x * other.z),
            ) + self.y * other.w
                + self.z * other.x,
            z : self.w.mul_add(
                other.z,
                self.x * other.y,
            ) - self.y * other.x
                + self.z * other.w,
        }
    }

    /// Rotates a vector by this quaternion.
    #[must_use]

    pub fn rotate_vector(
        &self,
        v : &Vector3D,
    ) -> Vector3D {

        let q_v = Self::new(
            0.0, v.x, v.y, v.z,
        );

        let result = self
            .multiply(&q_v)
            .multiply(
                &self.conjugate(),
            );

        Vector3D::new(
            result.x,
            result.y,
            result.z,
        )
    }

    /// Converts the quaternion to a 4x4 rotation matrix.
    #[must_use]

    pub fn to_matrix(
        &self
    ) -> Matrix<f64> {

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
                2.0f64.mul_add(
                    -(yy + zz),
                    1.0,
                ),
                2.0 * (xy - wz),
                2.0 * (xz + wy),
                0.0,
                2.0 * (xy + wz),
                2.0f64.mul_add(
                    -(xx + zz),
                    1.0,
                ),
                2.0 * (yz - wx),
                0.0,
                2.0 * (xz - wy),
                2.0 * (yz + wx),
                2.0f64.mul_add(
                    -(xx + yy),
                    1.0,
                ),
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
        other : &Self,
        t : f64,
    ) -> Self {

        let dot = self.w.mul_add(
            other.w,
            self.x * other.x,
        ) + self.y * other.y
            + self.z * other.z;

        // If dot is negative, negate one quaternion to take the shorter path
        let (q2, dot) = if dot < 0.0 {

            (
                Self::new(
                    -other.w,
                    -other.x,
                    -other.y,
                    -other.z,
                ),
                -dot,
            )
        } else {

            (*other, dot)
        };

        // Use linear interpolation for very close quaternions
        if dot > 0.9995 {

            return Self::new(
                t.mul_add(
                    q2.w - self.w,
                    self.w,
                ),
                t.mul_add(
                    q2.x - self.x,
                    self.x,
                ),
                t.mul_add(
                    q2.y - self.y,
                    self.y,
                ),
                t.mul_add(
                    q2.z - self.z,
                    self.z,
                ),
            )
            .normalize();
        }

        let theta_0 = dot.acos();

        let theta = theta_0 * t;

        let sin_theta = theta.sin();

        let sin_theta_0 = theta_0.sin();

        let s0 = (theta_0 - theta)
            .cos()
            - dot * sin_theta
                / sin_theta_0;

        let s1 =
            sin_theta / sin_theta_0;

        Self::new(
            s0.mul_add(
                self.w,
                s1 * q2.w,
            ),
            s0.mul_add(
                self.x,
                s1 * q2.x,
            ),
            s0.mul_add(
                self.y,
                s1 * q2.y,
            ),
            s0.mul_add(
                self.z,
                s1 * q2.z,
            ),
        )
    }
}

// ============================================================================
// Ray Tracing Primitives
// ============================================================================

/// A ray defined by an origin point and direction.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Ray {
    pub origin : Point3D,
    pub direction : Vector3D,
}

impl Ray {
    /// Creates a new ray.
    #[must_use]

    pub const fn new(
        origin : Point3D,
        direction : Vector3D,
    ) -> Self {

        Self {
            origin,
            direction,
        }
    }

    /// Returns the point at parameter t along the ray.
    #[must_use]

    pub fn at(
        &self,
        t : f64,
    ) -> Point3D {

        Point3D {
            x : t.mul_add(
                self.direction.x,
                self.origin.x,
            ),
            y : t.mul_add(
                self.direction.y,
                self.origin.y,
            ),
            z : t.mul_add(
                self.direction.z,
                self.origin.z,
            ),
        }
    }
}

/// A sphere defined by center and radius.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Sphere {
    pub center : Point3D,
    pub radius : f64,
}

impl Sphere {
    /// Creates a new sphere.
    #[must_use]

    pub const fn new(
        center : Point3D,
        radius : f64,
    ) -> Self {

        Self {
            center,
            radius,
        }
    }
}

/// Result of a ray intersection.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Intersection {
    /// Parameter t where ray intersects.
    pub t : f64,
    /// Point of intersection.
    pub point : Point3D,
    /// Normal at intersection point.
    pub normal : Vector3D,
}

/// Computes ray-sphere intersection.
///
/// Returns the closest intersection with t > 0, or `None` if no intersection.
#[must_use]

pub fn ray_sphere_intersection(
    ray : &Ray,
    sphere : &Sphere,
) -> Option<Intersection> {

    let oc = Vector3D::new(
        ray.origin.x - sphere.center.x,
        ray.origin.y - sphere.center.y,
        ray.origin.z - sphere.center.z,
    );

    let a = dot_product(
        &ray.direction,
        &ray.direction,
    );

    let b = 2.0
        * dot_product(
            &oc,
            &ray.direction,
        );

    let c = sphere
        .radius
        .mul_add(
            -sphere.radius,
            dot_product(&oc, &oc),
        );

    let discriminant =
        b * b - 4.0 * a * c;

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
        (point.x - sphere.center.x)
            / sphere.radius,
        (point.y - sphere.center.y)
            / sphere.radius,
        (point.z - sphere.center.z)
            / sphere.radius,
    );

    Some(Intersection {
        t,
        point,
        normal,
    })
}

/// A plane defined by a point and normal.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Serialize,
    Deserialize,
)]

pub struct Plane {
    pub point : Point3D,
    pub normal : Vector3D,
}

impl Plane {
    /// Creates a new plane.
    #[must_use]

    pub const fn new(
        point : Point3D,
        normal : Vector3D,
    ) -> Self {

        Self {
            point,
            normal,
        }
    }
}

/// Computes ray-plane intersection.
#[must_use]

pub fn ray_plane_intersection(
    ray : &Ray,
    plane : &Plane,
) -> Option<Intersection> {

    let denom = dot_product(
        &plane.normal,
        &ray.direction,
    );

    if denom.abs() < 1e-6 {

        return None; // Ray is parallel to plane
    }

    let p0_l0 = Vector3D::new(
        plane.point.x - ray.origin.x,
        plane.point.y - ray.origin.y,
        plane.point.z - ray.origin.z,
    );

    let t = dot_product(
        &p0_l0,
        &plane.normal,
    ) / denom;

    if t < 1e-6 {

        return None;
    }

    let point = ray.at(t);

    let normal = if denom < 0.0 {

        plane.normal
    } else {

        -plane.normal
    };

    Some(Intersection {
        t,
        point,
        normal,
    })
}

/// Computes ray-triangle intersection using Möller–Trumbore algorithm.
#[must_use]

pub fn ray_triangle_intersection(
    ray : &Ray,
    v0 : &Point3D,
    v1 : &Point3D,
    v2 : &Point3D,
) -> Option<Intersection> {

    let edge1 = Vector3D::new(
        v1.x - v0.x,
        v1.y - v0.y,
        v1.z - v0.z,
    );

    let edge2 = Vector3D::new(
        v2.x - v0.x,
        v2.y - v0.y,
        v2.z - v0.z,
    );

    let h = cross_product(
        &ray.direction,
        &edge2,
    );

    let a = dot_product(&edge1, &h);

    if a.abs() < 1e-6 {

        return None;
    }

    let f = 1.0 / a;

    let s = Vector3D::new(
        ray.origin.x - v0.x,
        ray.origin.y - v0.y,
        ray.origin.z - v0.z,
    );

    let u = f * dot_product(&s, &h);

    if !(0.0 ..= 1.0).contains(&u) {

        return None;
    }

    let q = cross_product(&s, &edge1);

    let v = f * dot_product(
        &ray.direction,
        &q,
    );

    if v < 0.0 || u + v > 1.0 {

        return None;
    }

    let t = f * dot_product(&edge2, &q);

    if t < 1e-6 {

        return None;
    }

    let point = ray.at(t);

    let normal =
        cross_product(&edge1, &edge2)
            .normalize();

    Some(Intersection {
        t,
        point,
        normal,
    })
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
    p0 : &Point3D,
    p1 : &Point3D,
    p2 : &Point3D,
    t : f64,
) -> Point3D {

    let t2 = t * t;

    let mt = 1.0 - t;

    let mt2 = mt * mt;

    Point3D {
        x : t2.mul_add(
            p2.x,
            mt2 * p0.x
                + 2.0 * mt * t * p1.x,
        ),
        y : t2.mul_add(
            p2.y,
            mt2 * p0.y
                + 2.0 * mt * t * p1.y,
        ),
        z : t2.mul_add(
            p2.z,
            mt2 * p0.z
                + 2.0 * mt * t * p1.z,
        ),
    }
}

/// Evaluates a cubic Bezier curve at parameter t.
///
/// # Arguments
/// * `p0`, `p1`, `p2`, `p3` - Control points
/// * `t` - Parameter in [0, 1]
#[must_use]

pub fn bezier_cubic(
    p0 : &Point3D,
    p1 : &Point3D,
    p2 : &Point3D,
    p3 : &Point3D,
    t : f64,
) -> Point3D {

    let t2 = t * t;

    let t3 = t2 * t;

    let mt = 1.0 - t;

    let mt2 = mt * mt;

    let mt3 = mt2 * mt;

    Point3D {
        x : (3.0 * mt * t2).mul_add(
            p2.x,
            mt3 * p0.x
                + 3.0 * mt2 * t * p1.x,
        ) + t3 * p3.x,
        y : (3.0 * mt * t2).mul_add(
            p2.y,
            mt3 * p0.y
                + 3.0 * mt2 * t * p1.y,
        ) + t3 * p3.y,
        z : (3.0 * mt * t2).mul_add(
            p2.z,
            mt3 * p0.z
                + 3.0 * mt2 * t * p1.z,
        ) + t3 * p3.z,
    }
}

/// Evaluates a Catmull-Rom spline at parameter t.
///
/// The spline passes through p1 and p2, with p0 and p3 as guide points.
#[must_use]

pub fn catmull_rom(
    p0 : &Point3D,
    p1 : &Point3D,
    p2 : &Point3D,
    p3 : &Point3D,
    t : f64,
) -> Point3D {

    let t2 = t * t;

    let t3 = t2 * t;

    Point3D {
        x : 0.5
            * (2.0f64.mul_add(
                p1.x,
                (-p0.x + p2.x) * t,
            ) + (2.0f64.mul_add(
                p0.x,
                -(5.0 * p1.x),
            ) + 4.0 * p2.x
                - p3.x)
                * t2
                + (3.0f64.mul_add(
                    p1.x, -p0.x,
                ) - 3.0 * p2.x
                    + p3.x)
                    * t3),
        y : 0.5
            * (2.0f64.mul_add(
                p1.y,
                (-p0.y + p2.y) * t,
            ) + (2.0f64.mul_add(
                p0.y,
                -(5.0 * p1.y),
            ) + 4.0 * p2.y
                - p3.y)
                * t2
                + (3.0f64.mul_add(
                    p1.y, -p0.y,
                ) - 3.0 * p2.y
                    + p3.y)
                    * t3),
        z : 0.5
            * (2.0f64.mul_add(
                p1.z,
                (-p0.z + p2.z) * t,
            ) + (2.0f64.mul_add(
                p0.z,
                -(5.0 * p1.z),
            ) + 4.0 * p2.z
                - p3.z)
                * t2
                + (3.0f64.mul_add(
                    p1.z, -p0.z,
                ) - 3.0 * p2.z
                    + p3.z)
                    * t3),
    }
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Converts degrees to radians.
#[must_use]

pub fn degrees_to_radians(
    degrees : f64
) -> f64 {

    degrees.to_radians()
}

/// Converts radians to degrees.
#[must_use]

pub fn radians_to_degrees(
    radians : f64
) -> f64 {

    radians.to_degrees()
}

/// Applies a 4x4 transformation matrix to a 3D point.
#[must_use]

pub fn transform_point(
    matrix : &Matrix<f64>,
    point : &Point3D,
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

    if w.abs() > 1e-10 && w != 1.0 {

        Point3D::new(
            x / w,
            y / w,
            z / w,
        )
    } else {

        Point3D::new(x, y, z)
    }
}

/// Applies a 4x4 transformation matrix to a 3D vector (ignores translation).
#[must_use]

pub fn transform_vector(
    matrix : &Matrix<f64>,
    vector : &Vector3D,
) -> Vector3D {

    let x = matrix.get(0, 0) * vector.x
        + matrix.get(0, 1) * vector.y
        + matrix.get(0, 2) * vector.z;

    let y = matrix.get(1, 0) * vector.x
        + matrix.get(1, 1) * vector.y
        + matrix.get(1, 2) * vector.z;

    let z = matrix.get(2, 0) * vector.x
        + matrix.get(2, 1) * vector.y
        + matrix.get(2, 2) * vector.z;

    Vector3D::new(x, y, z)
}

/// Computes the barycentric coordinates of a point in a triangle.
///
/// Returns (u, v, w) where point = u*v0 + v*v1 + w*v2.
#[must_use]

pub fn barycentric_coordinates(
    point : &Point3D,
    v0 : &Point3D,
    v1 : &Point3D,
    v2 : &Point3D,
) -> (f64, f64, f64) {

    let v0v1 = Vector3D::new(
        v1.x - v0.x,
        v1.y - v0.y,
        v1.z - v0.z,
    );

    let v0v2 = Vector3D::new(
        v2.x - v0.x,
        v2.y - v0.y,
        v2.z - v0.z,
    );

    let v0p = Vector3D::new(
        point.x - v0.x,
        point.y - v0.y,
        point.z - v0.z,
    );

    let d00 = dot_product(&v0v1, &v0v1);

    let d01 = dot_product(&v0v1, &v0v2);

    let d11 = dot_product(&v0v2, &v0v2);

    let d20 = dot_product(&v0p, &v0v1);

    let d21 = dot_product(&v0p, &v0v2);

    let denom =
        d00.mul_add(d11, -(d01 * d01));

    if denom.abs() < 1e-10 {

        return (1.0, 0.0, 0.0);
    }

    let v = d11
        .mul_add(d20, -(d01 * d21))
        / denom;

    let w = d00
        .mul_add(d21, -(d01 * d20))
        / denom;

    let u = 1.0 - v - w;

    (u, v, w)
}
