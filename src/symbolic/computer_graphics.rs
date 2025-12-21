//! # Symbolic Computer Graphics
//!
//! This module provides symbolic tools for 2D and 3D computer graphics,
//! including transformations (translation, rotation, scaling), projections
//! (perspective, orthographic), curve representations (Bezier, B-spline),
//! and mesh manipulation.
use crate::symbolic::core::Expr;
use crate::symbolic::elementary::{cos, sin, tan};
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::vector::Vector;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use std::sync::Arc;
/// Generates a 3x3 2D translation matrix.
///
/// This matrix can be used to move objects in 2D space by `tx` and `ty`.
///
/// # Arguments
/// * `tx` - Translation along the x-axis.
/// * `ty` - Translation along the y-axis.
///
/// # Returns
/// An `Expr::Matrix` representing the 2D translation transformation.
#[must_use]
pub fn translation_2d(tx: Expr, ty: Expr) -> Expr {
    Expr::Matrix(vec![
        vec![
            Expr::BigInt(BigInt::one()),
            Expr::BigInt(BigInt::zero()),
            tx,
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
            ty,
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 4x4 3D translation matrix.
///
/// This matrix can be used to move objects in 3D space by `tx`, `ty`, and `tz`.
///
/// # Arguments
/// * `tx` - Translation along the x-axis.
/// * `ty` - Translation along the y-axis.
/// * `tz` - Translation along the z-axis.
///
/// # Returns
/// An `Expr::Matrix` representing the 3D translation transformation.
#[must_use]
pub fn translation_3d(tx: Expr, ty: Expr, tz: Expr) -> Expr {
    Expr::Matrix(vec![
        vec![
            Expr::BigInt(BigInt::one()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            tx,
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
            Expr::BigInt(BigInt::zero()),
            ty,
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
            tz,
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 3x3 2D rotation matrix.
///
/// This matrix rotates objects around the origin in 2D space by a given `angle`.
///
/// # Arguments
/// * `angle` - The rotation angle.
///
/// # Returns
/// An `Expr::Matrix` representing the 2D rotation transformation.
#[must_use]
pub fn rotation_2d(angle: Expr) -> Expr {
    let c = cos(angle.clone());
    let s = sin(angle);
    Expr::Matrix(vec![
        vec![
            c.clone(),
            Expr::Neg(Arc::new(s.clone())),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![s, c, Expr::BigInt(BigInt::zero())],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 4x4 3D rotation matrix around the X-axis.
///
/// # Arguments
/// * `angle` - The rotation angle.
///
/// # Returns
/// An `Expr::Matrix` representing the 3D rotation transformation around the X-axis.
#[must_use]
pub fn rotation_3d_x(angle: Expr) -> Expr {
    let c = cos(angle.clone());
    let s = sin(angle);
    Expr::Matrix(vec![
        vec![
            Expr::BigInt(BigInt::one()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            c.clone(),
            Expr::Neg(Arc::new(s.clone())),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            s,
            c,
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 4x4 3D rotation matrix around the Y-axis.
///
/// # Arguments
/// * `angle` - The rotation angle.
///
/// # Returns
/// An `Expr::Matrix` representing the 3D rotation transformation around the Y-axis.
#[must_use]
pub fn rotation_3d_y(angle: Expr) -> Expr {
    let c = cos(angle.clone());
    let s = sin(angle);
    Expr::Matrix(vec![
        vec![
            c.clone(),
            Expr::BigInt(BigInt::zero()),
            s.clone(),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::Neg(Arc::new(s)),
            Expr::BigInt(BigInt::zero()),
            c,
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 4x4 3D rotation matrix around the Z-axis.
///
/// # Arguments
/// * `angle` - The rotation angle.
///
/// # Returns
/// An `Expr::Matrix` representing the 3D rotation transformation around the Z-axis.
#[must_use]
pub fn rotation_3d_z(angle: Expr) -> Expr {
    let c = cos(angle.clone());
    let s = sin(angle);
    Expr::Matrix(vec![
        vec![
            c.clone(),
            Expr::Neg(Arc::new(s.clone())),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            s,
            c,
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 3x3 2D scaling matrix.
///
/// This matrix scales objects in 2D space by `sx` and `sy`.
///
/// # Arguments
/// * `sx` - Scaling factor along the x-axis.
/// * `sy` - Scaling factor along the y-axis.
///
/// # Returns
/// An `Expr::Matrix` representing the 2D scaling transformation.
#[must_use]
pub fn scaling_2d(sx: Expr, sy: Expr) -> Expr {
    Expr::Matrix(vec![
        vec![
            sx,
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            sy,
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 4x4 3D scaling matrix.
///
/// This matrix scales objects in 3D space by `sx`, `sy`, and `sz`.
///
/// # Arguments
/// * `sx` - Scaling factor along the x-axis.
/// * `sy` - Scaling factor along the y-axis.
/// * `sz` - Scaling factor along the z-axis.
///
/// # Returns
/// An `Expr::Matrix` representing the 3D scaling transformation.
#[must_use]
pub fn scaling_3d(sx: Expr, sy: Expr, sz: Expr) -> Expr {
    Expr::Matrix(vec![
        vec![
            sx,
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            sy,
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            sz,
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 4x4 perspective projection matrix.
///
/// This matrix transforms 3D points into 2D points, simulating how objects
/// appear smaller when they are farther away.
///
/// # Arguments
/// * `fovy` - Field of view in the y-direction.
/// * `aspect` - Aspect ratio of the viewport (width / height).
/// * `near` - Distance to the near clipping plane.
/// * `far` - Distance to the far clipping plane.
///
/// # Returns
/// An `Expr::Matrix` representing the perspective projection.
#[must_use]
pub fn perspective_projection(fovy: Expr, aspect: &Expr, near: Expr, far: Expr) -> Expr {
    let f = tan(Expr::new_div(fovy, Expr::BigInt(BigInt::from(2))));
    let range_inv = Expr::new_div(
        Expr::BigInt(BigInt::one()),
        Expr::new_sub(near.clone(), far.clone()),
    );
    Expr::Matrix(vec![
        vec![
            Expr::Div(
                Arc::new(Expr::BigInt(BigInt::one())),
                Arc::new(Expr::Mul(Arc::new(f.clone()), Arc::new(aspect.clone()))),
            ),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::Div(Arc::new(Expr::BigInt(BigInt::one())), Arc::new(f)),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::Mul(
                Arc::new(Expr::Add(Arc::new(near.clone()), Arc::new(far.clone()))),
                Arc::new(range_inv),
            ),
            Expr::Mul(
                Arc::new(Expr::Mul(
                    Arc::new(Expr::BigInt(BigInt::from(2))),
                    Arc::new(near),
                )),
                Arc::new(far),
            ),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::from(-1)),
            Expr::BigInt(BigInt::zero()),
        ],
    ])
}
/// Generates a 4x4 orthographic projection matrix.
///
/// This matrix transforms 3D points into 2D points without perspective distortion,
/// meaning parallel lines remain parallel.
///
/// # Arguments
/// * `left` - Coordinate of the left vertical clipping plane.
/// * `right` - Coordinate of the right vertical clipping plane.
/// * `bottom` - Coordinate of the bottom horizontal clipping plane.
/// * `top` - Coordinate of the top horizontal clipping plane.
/// * `near` - Distance to the near clipping plane.
/// * `far` - Distance to the far clipping plane.
///
/// # Returns
/// An `Expr::Matrix` representing the orthographic projection.
#[must_use]
pub fn orthographic_projection(
    left: Expr,
    right: Expr,
    bottom: Expr,
    top: Expr,
    near: Expr,
    far: Expr,
) -> Expr {
    let r_l = Expr::new_div(
        Expr::BigInt(BigInt::one()),
        Expr::new_sub(right.clone(), left.clone()),
    );
    let t_b = Expr::new_div(
        Expr::BigInt(BigInt::one()),
        Expr::new_sub(top.clone(), bottom.clone()),
    );
    let f_n = Expr::new_div(
        Expr::BigInt(BigInt::one()),
        Expr::new_sub(far.clone(), near.clone()),
    );
    Expr::Matrix(vec![
        vec![
            Expr::Mul(
                Arc::new(Expr::BigInt(BigInt::from(2))),
                Arc::new(r_l.clone()),
            ),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::Neg(Arc::new(Expr::Mul(
                Arc::new(Expr::Add(Arc::new(right), Arc::new(left))),
                Arc::new(r_l),
            ))),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::Mul(
                Arc::new(Expr::BigInt(BigInt::from(2))),
                Arc::new(t_b.clone()),
            ),
            Expr::BigInt(BigInt::zero()),
            Expr::Neg(Arc::new(Expr::Mul(
                Arc::new(Expr::Add(Arc::new(top), Arc::new(bottom))),
                Arc::new(t_b),
            ))),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::Neg(Arc::new(Expr::Mul(
                Arc::new(Expr::BigInt(BigInt::from(2))),
                Arc::new(f_n.clone()),
            ))),
            Expr::Neg(Arc::new(Expr::Mul(
                Arc::new(Expr::Add(Arc::new(far), Arc::new(near))),
                Arc::new(f_n),
            ))),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Generates a 4x4 "look at" view matrix.
///
/// This matrix transforms world coordinates into view coordinates, effectively
/// positioning and orienting a camera in 3D space.
///
/// # Arguments
/// * `eye` - The position of the camera in world space.
/// * `center` - The point in world space that the camera is looking at.
/// * `up` - The up direction vector in world space.
///
/// # Returns
/// An `Expr::Matrix` representing the view transformation.
#[must_use]
pub fn look_at(eye: &Vector, center: &Vector, up: &Vector) -> Expr {
    let f = (center.clone() - eye.clone()).normalize();
    let s = f.cross(up).normalize();
    let u = s.cross(&f);
    Expr::Matrix(vec![
        vec![
            s.x.clone(),
            s.y.clone(),
            s.z.clone(),
            Expr::Neg(Arc::new(s.dot(eye))),
        ],
        vec![
            u.x.clone(),
            u.y.clone(),
            u.z.clone(),
            Expr::Neg(Arc::new(u.dot(eye))),
        ],
        vec![
            Expr::Neg(Arc::new(f.x.clone())),
            Expr::Neg(Arc::new(f.y.clone())),
            Expr::Neg(Arc::new(f.z.clone())),
            f.dot(eye),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
/// Represents a Bézier curve defined by a set of control points.
///
/// A Bézier curve of degree `n` is defined by `n + 1` control points. The curve
/// starts at the first control point and ends at the last control point,
/// passing smoothly between them based on the Bernstein polynomial basis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BezierCurve {
    /// The control points defining the curve shape.
    pub control_points: Vec<Vector>,
    /// The degree of the curve (number of control points - 1).
    pub degree: usize,
}
impl BezierCurve {
    /// Evaluates the Bezier curve at a given parameter `t`.
    ///
    /// The Bezier curve is defined by a set of control points and the Bernstein polynomials.
    /// `B(t) = Σ_{i=0 to n} [ B_i,n(t) * P_i ]`, where `B_i,n(t)` are the Bernstein polynomials
    /// and `P_i` are the control points.
    ///
    /// # Arguments
    /// * `t` - The parameter `t` (typically in the range [0, 1]).
    ///
    /// # Returns
    /// A `Vector` representing the point on the curve at parameter `t`.
    #[allow(clippy::cast_possible_wrap)]
    #[must_use]
    pub fn evaluate(&self, t: &Expr) -> Vector {
        let n = self.degree as i64;
        let mut result = Vector::new(
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
        );
        for (i, pt) in self.control_points.iter().enumerate() {
            let i_bigint = BigInt::from(i);
            let n_bigint = BigInt::from(n);
            let bernstein = Expr::new_mul(
                Expr::Binomial(
                    Arc::new(Expr::BigInt(n_bigint.clone())),
                    Arc::new(Expr::BigInt(i_bigint.clone())),
                ),
                Expr::new_pow(t.clone(), Expr::BigInt(i_bigint.clone())),
            );
            let bernstein = Expr::new_mul(
                bernstein,
                Expr::new_pow(
                    Expr::new_sub(Expr::BigInt(BigInt::one()), t.clone()),
                    Expr::BigInt(n_bigint - i_bigint),
                ),
            );
            result = result + pt.scalar_mul(&bernstein);
        }
        result
    }

    /// Computes the derivative (tangent vector) of the Bezier curve at parameter `t`.
    ///
    /// The derivative of a Bézier curve of degree n is a Bézier curve of degree n-1
    /// with control points n * (P_{i+1} - `P_i`).
    ///
    /// # Arguments
    /// * `t` - The parameter `t` (typically in the range [0, 1]).
    ///
    /// # Returns
    /// A `Vector` representing the tangent vector at parameter `t`.
    #[must_use]
    pub fn derivative(&self, t: &Expr) -> Vector {
        if self.degree == 0 || self.control_points.len() < 2 {
            return Vector::new(
                Expr::BigInt(BigInt::zero()),
                Expr::BigInt(BigInt::zero()),
                Expr::BigInt(BigInt::zero()),
            );
        }

        // Derivative control points: n * (P_{i+1} - P_i)
        let n = Expr::BigInt(BigInt::from(self.degree as i64));
        let derivative_points: Vec<Vector> = (0..self.control_points.len() - 1)
            .map(|i| {
                let diff = self.control_points[i + 1].clone() - self.control_points[i].clone();
                diff.scalar_mul(&n)
            })
            .collect();

        let derivative_curve = Self {
            control_points: derivative_points,
            degree: self.degree - 1,
        };
        derivative_curve.evaluate(t)
    }

    /// Splits the Bezier curve at parameter `t` using De Casteljau's algorithm.
    ///
    /// This subdivides the curve into two Bézier curves that together represent
    /// the same curve as the original.
    ///
    /// # Arguments
    /// * `t` - The parameter `t` at which to split (typically in range [0, 1]).
    ///
    /// # Returns
    /// A tuple of two `BezierCurve` representing the left and right portions.
    #[must_use]
    pub fn split(&self, t: &Expr) -> (Self, Self) {
        let n = self.control_points.len();
        let mut pyramid: Vec<Vec<Vector>> = vec![self.control_points.clone()];

        // Build the De Casteljau pyramid
        for level in 1..n {
            let prev = &pyramid[level - 1];
            let mut current = Vec::with_capacity(n - level);
            for i in 0..(n - level) {
                let one_minus_t = simplify(&Expr::new_sub(Expr::BigInt(BigInt::one()), t.clone()));
                let left = prev[i].scalar_mul(&one_minus_t);
                let right = prev[i + 1].scalar_mul(t);
                current.push(left + right);
            }
            pyramid.push(current);
        }

        // Left curve: first element of each level
        let left_points: Vec<Vector> = (0..n).map(|i| pyramid[i][0].clone()).collect();

        // Right curve: last element of each level (in reverse)
        let right_points: Vec<Vector> = (0..n).map(|i| pyramid[n - 1 - i][i].clone()).collect();

        (
            Self {
                control_points: left_points,
                degree: self.degree,
            },
            Self {
                control_points: right_points,
                degree: self.degree,
            },
        )
    }
}
/// Represents a B-spline curve with local control.
///
/// B-splines provide local control over curve shape, meaning that moving a control
/// point only affects a local region of the curve. The curve is defined by control
/// points, a knot vector, and a degree.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BSplineCurve {
    /// The control points defining the curve shape.
    pub control_points: Vec<Vector>,
    /// The knot vector defining the parameterization.
    pub knots: Vec<Expr>,
    /// The degree of the B-spline basis functions.
    pub degree: usize,
}
impl BSplineCurve {
    /// Evaluates the B-spline curve at parameter `t` using De Boor's algorithm.
    ///
    /// De Boor's algorithm is a numerically stable and efficient method for evaluating
    /// B-spline curves. It uses a recursive formula to compute the point on the curve.
    ///
    /// # Arguments
    /// * `t` - The parameter `t`.
    ///
    /// # Returns
    /// A `Vector` representing the point on the curve at parameter `t`.
    #[must_use]
    pub fn evaluate(&self, t: &Expr) -> Vector {
        let p = self.degree;
        let k = self.knots.len() - 1 - p - 1;
        let mut d: Vec<Vector> = self.control_points[..=k + p].to_vec();
        for r in 1..=p {
            for j in (r..=k + p).rev() {
                let t_j = &self.knots[j];
                let t_j_p1 = &self.knots[j + p + 1 - r];
                let alpha = simplify(&Expr::new_div(
                    Expr::new_sub(t.clone(), t_j.clone()),
                    Expr::new_sub(t_j_p1.clone(), t_j.clone()),
                ));
                d[j] = d[j - 1].scalar_mul(&simplify(&Expr::new_sub(
                    Expr::BigInt(BigInt::one()),
                    alpha.clone(),
                ))) + d[j].scalar_mul(&alpha);
            }
        }
        d[k + p].clone()
    }
}
/// Represents a single polygon face in a mesh.
/// The `indices` field contains a list of indices that point to vertices
/// in the `vertices` list of a `PolygonMesh`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Polygon {
    pub indices: Vec<usize>,
}
impl Polygon {
    /// Creates a new polygon from a list of vertex indices.
    #[must_use]
    pub const fn new(indices: Vec<usize>) -> Self {
        Self { indices }
    }
}
/// Represents a 3D object as a polygon mesh.
/// A mesh is composed of a list of vertices (3D points) and a list of polygons (faces)
/// that connect those vertices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolygonMesh {
    pub vertices: Vec<Vector>,
    pub polygons: Vec<Polygon>,
}
impl PolygonMesh {
    /// Creates a new polygon mesh from a list of vertices and polygons.
    ///
    /// # Arguments
    /// * `vertices` - A vector of `Vector` representing the 3D points of the mesh.
    /// * `polygons` - A vector of `Polygon` representing the faces of the mesh.
    #[must_use]
    pub const fn new(vertices: Vec<Vector>, polygons: Vec<Polygon>) -> Self {
        Self { vertices, polygons }
    }
    /// Applies a geometric transformation to the entire mesh.
    ///
    /// This function iterates through all vertices of the mesh and applies the given
    /// transformation matrix to each one. The transformation is performed symbolically.
    ///
    /// # Arguments
    /// * `transformation` - An `Expr::Matrix` representing the 4x4 transformation matrix.
    ///
    /// # Returns
    /// A new `PolygonMesh` with the transformed vertices.
    ///
    /// # Panics
    /// Panics if the provided `transformation` is not a matrix.
    pub fn apply_transformation(&self, transformation: &Expr) -> Result<Self, String> {
        if let Expr::Matrix(matrix) = transformation {
            let transformed_vertices = self
                .vertices
                .iter()
                .map(|vertex| {
                    let homogeneous_vertex = Expr::Vector(vec![
                        vertex.x.clone(),
                        vertex.y.clone(),
                        vertex.z.clone(),
                        Expr::BigInt(BigInt::one()),
                    ]);
                    let transformed_homogeneous = simplify(&Expr::new_matrix_vec_mul(
                        Expr::Matrix(matrix.clone()),
                        homogeneous_vertex,
                    ));
                    if let Expr::Vector(vec) = transformed_homogeneous {
                        let w = vec
                            .get(3)
                            .cloned()
                            .unwrap_or_else(|| Expr::BigInt(BigInt::one()));
                        let x = simplify(&Expr::new_div(vec[0].clone(), w.clone()));
                        let y = simplify(&Expr::new_div(vec[1].clone(), w.clone()));
                        let z = simplify(&Expr::new_div(vec[2].clone(), w));
                        Vector::new(x, y, z)
                    } else {
                        vertex.clone()
                    }
                })
                .collect();
            Ok(Self {
                vertices: transformed_vertices,
                polygons: self.polygons.clone(),
            })
        } else {
            Err("Transformation must be an Expr::Matrix".to_string())
        }
    }

    /// Computes the surface normal vectors for each polygon in the mesh.
    ///
    /// The normal is computed using the cross product of two edges of each polygon.
    /// For a triangle with vertices A, B, C, the normal is (B - A) × (C - A).
    ///
    /// # Returns
    /// A vector of `Vector` representing the normal for each polygon.
    #[must_use]
    pub fn compute_normals(&self) -> Vec<Vector> {
        self.polygons
            .iter()
            .filter_map(|poly| {
                if poly.indices.len() >= 3 {
                    let v0 = &self.vertices[poly.indices[0]];
                    let v1 = &self.vertices[poly.indices[1]];
                    let v2 = &self.vertices[poly.indices[2]];
                    let edge1 = v1.clone() - v0.clone();
                    let edge2 = v2.clone() - v0.clone();
                    Some(edge1.cross(&edge2).normalize())
                } else {
                    None
                }
            })
            .collect()
    }

    /// Triangulates all polygons in the mesh into triangles.
    ///
    /// Uses a simple fan triangulation where each polygon with n vertices
    /// is converted into n-2 triangles sharing the first vertex.
    ///
    /// # Returns
    /// A new `PolygonMesh` containing only triangular polygons.
    #[must_use]
    pub fn triangulate(&self) -> Self {
        let triangles: Vec<Polygon> = self
            .polygons
            .iter()
            .flat_map(|poly| {
                if poly.indices.len() <= 3 {
                    vec![poly.clone()]
                } else {
                    // Fan triangulation
                    (1..poly.indices.len() - 1)
                        .map(|i| {
                            Polygon::new(vec![
                                poly.indices[0],
                                poly.indices[i],
                                poly.indices[i + 1],
                            ])
                        })
                        .collect()
                }
            })
            .collect();
        Self {
            vertices: self.vertices.clone(),
            polygons: triangles,
        }
    }
}

/// Generates a 3x3 2D shear matrix.
///
/// This matrix shears objects in 2D space. Shearing displaces each point
/// in a fixed direction by an amount proportional to its perpendicular distance.
///
/// # Arguments
/// * `shx` - Shear factor along the x-axis (y displacement creates x movement).
/// * `shy` - Shear factor along the y-axis (x displacement creates y movement).
///
/// # Returns
/// An `Expr::Matrix` representing the 2D shear transformation.
#[must_use]
pub fn shear_2d(shx: Expr, shy: Expr) -> Expr {
    Expr::Matrix(vec![
        vec![
            Expr::BigInt(BigInt::one()),
            shx,
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            shy,
            Expr::BigInt(BigInt::one()),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}

/// Generates a 3x3 2D reflection matrix across a line through the origin.
///
/// # Arguments
/// * `angle` - The angle of the reflection line from the x-axis.
///
/// # Returns
/// An `Expr::Matrix` representing the 2D reflection transformation.
#[must_use]
pub fn reflection_2d(angle: Expr) -> Expr {
    let two_angle = Expr::new_mul(Expr::Constant(2.0), angle);
    let c = cos(two_angle.clone());
    let s = sin(two_angle);
    Expr::Matrix(vec![
        vec![c.clone(), s.clone(), Expr::BigInt(BigInt::zero())],
        vec![s, Expr::Neg(Arc::new(c)), Expr::BigInt(BigInt::zero())],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}

/// Generates a 4x4 3D reflection matrix across a plane through the origin.
///
/// The plane is defined by its normal vector (nx, ny, nz) which should be normalized.
///
/// # Arguments
/// * `nx` - X component of the plane normal.
/// * `ny` - Y component of the plane normal.
/// * `nz` - Z component of the plane normal.
///
/// # Returns
/// An `Expr::Matrix` representing the 3D reflection transformation.
#[must_use]
pub fn reflection_3d(nx: Expr, ny: Expr, nz: Expr) -> Expr {
    // Reflection matrix: I - 2 * n * n^T
    let two = Expr::Constant(2.0);
    Expr::Matrix(vec![
        vec![
            simplify(&Expr::new_sub(
                Expr::BigInt(BigInt::one()),
                Expr::new_mul(two.clone(), Expr::new_mul(nx.clone(), nx.clone())),
            )),
            simplify(&Expr::new_neg(Expr::new_mul(
                two.clone(),
                Expr::new_mul(nx.clone(), ny.clone()),
            ))),
            simplify(&Expr::new_neg(Expr::new_mul(
                two.clone(),
                Expr::new_mul(nx.clone(), nz.clone()),
            ))),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            simplify(&Expr::new_neg(Expr::new_mul(
                two.clone(),
                Expr::new_mul(ny.clone(), nx.clone()),
            ))),
            simplify(&Expr::new_sub(
                Expr::BigInt(BigInt::one()),
                Expr::new_mul(two.clone(), Expr::new_mul(ny.clone(), ny.clone())),
            )),
            simplify(&Expr::new_neg(Expr::new_mul(
                two.clone(),
                Expr::new_mul(ny.clone(), nz.clone()),
            ))),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            simplify(&Expr::new_neg(Expr::new_mul(
                two.clone(),
                Expr::new_mul(nz.clone(), nx),
            ))),
            simplify(&Expr::new_neg(Expr::new_mul(
                two.clone(),
                Expr::new_mul(nz.clone(), ny),
            ))),
            simplify(&Expr::new_sub(
                Expr::BigInt(BigInt::one()),
                Expr::new_mul(two, Expr::new_mul(nz.clone(), nz)),
            )),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}

/// Generates a 4x4 3D rotation matrix around an arbitrary axis using Rodrigues' formula.
///
/// # Arguments
/// * `axis` - A `Vector` representing the axis of rotation (should be normalized).
/// * `angle` - The rotation angle.
///
/// # Returns
/// An `Expr::Matrix` representing the 3D rotation transformation around the given axis.
#[must_use]
pub fn rotation_axis_angle(axis: &Vector, angle: Expr) -> Expr {
    let c = cos(angle.clone());
    let s = sin(angle);
    let one_minus_c = simplify(&Expr::new_sub(Expr::BigInt(BigInt::one()), c.clone()));

    let ux = axis.x.clone();
    let uy = axis.y.clone();
    let uz = axis.z.clone();

    // Rodrigues' rotation formula
    Expr::Matrix(vec![
        vec![
            simplify(&Expr::new_add(
                c.clone(),
                Expr::new_mul(Expr::new_mul(ux.clone(), ux.clone()), one_minus_c.clone()),
            )),
            simplify(&Expr::new_sub(
                Expr::new_mul(Expr::new_mul(ux.clone(), uy.clone()), one_minus_c.clone()),
                Expr::new_mul(uz.clone(), s.clone()),
            )),
            simplify(&Expr::new_add(
                Expr::new_mul(Expr::new_mul(ux.clone(), uz.clone()), one_minus_c.clone()),
                Expr::new_mul(uy.clone(), s.clone()),
            )),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            simplify(&Expr::new_add(
                Expr::new_mul(Expr::new_mul(uy.clone(), ux.clone()), one_minus_c.clone()),
                Expr::new_mul(uz.clone(), s.clone()),
            )),
            simplify(&Expr::new_add(
                c.clone(),
                Expr::new_mul(Expr::new_mul(uy.clone(), uy.clone()), one_minus_c.clone()),
            )),
            simplify(&Expr::new_sub(
                Expr::new_mul(Expr::new_mul(uy.clone(), uz.clone()), one_minus_c.clone()),
                Expr::new_mul(ux.clone(), s.clone()),
            )),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            simplify(&Expr::new_sub(
                Expr::new_mul(Expr::new_mul(uz.clone(), ux.clone()), one_minus_c.clone()),
                Expr::new_mul(uy.clone(), s.clone()),
            )),
            simplify(&Expr::new_add(
                Expr::new_mul(Expr::new_mul(uz.clone(), uy), one_minus_c.clone()),
                Expr::new_mul(ux, s),
            )),
            simplify(&Expr::new_add(
                c,
                Expr::new_mul(Expr::new_mul(uz.clone(), uz), one_minus_c),
            )),
            Expr::BigInt(BigInt::zero()),
        ],
        vec![
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::zero()),
            Expr::BigInt(BigInt::one()),
        ],
    ])
}
