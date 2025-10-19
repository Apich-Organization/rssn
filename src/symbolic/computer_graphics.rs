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
            Expr::Neg(Arc::new(s.clone())),
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
                Arc::new(range_inv.clone()),
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
#[derive(Clone, Debug, PartialEq)]
pub struct BezierCurve {
    pub control_points: Vec<Vector>,
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
}
#[derive(Clone, Debug, PartialEq)]
pub struct BSplineCurve {
    pub control_points: Vec<Vector>,
    pub knots: Vec<Expr>,
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
    pub fn new(indices: Vec<usize>) -> Self {
        Self { indices }
    }
}
/// Represents a 3D object as a polygon mesh.
/// A mesh is composed of a list of vertices (3D points) and a list of polygons (faces)
/// that connect those vertices.
#[derive(Clone, Debug, PartialEq)]
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
    pub fn new(vertices: Vec<Vector>, polygons: Vec<Polygon>) -> Self {
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
                        let z = simplify(&Expr::new_div(vec[2].clone(), w.clone()));
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
}
