//! # Symbolic Vector Algebra and Calculus
//!
//! This module defines a 3D symbolic vector and implements a range of operations
//! for vector algebra and vector calculus. It includes basic vector arithmetic,
//! dot and cross products, as well as differential operators like gradient,
//! divergence, and curl.

use std::ops::Add;
use std::ops::Sub;

use num_bigint::BigInt;
use num_traits::One;
use serde::Deserialize;
use serde::Serialize;

use crate::symbolic::calculus::differentiate;
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::is_zero;
use crate::symbolic::simplify_dag::simplify;

/// Represents a symbolic vector in 3D space.
#[derive(
    Clone,
    Debug,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
)]

pub struct Vector {
    /// The x-component of the vector.
    pub x: Expr,
    /// The y-component of the vector.
    pub y: Expr,
    /// The z-component of the vector.
    pub z: Expr,
}

impl Vector {
    /// Creates a new symbolic vector with the given components.
    ///
    /// # Arguments
    /// * `x` - The expression for the x-component.
    /// * `y` - The expression for the y-component.
    /// * `z` - The expression for the z-component.
    #[must_use]

    pub const fn new(
        x: Expr,
        y: Expr,
        z: Expr,
    ) -> Self {

        Self {
            x,
            y,
            z,
        }
    }

    /// Computes the magnitude (Euclidean norm) of the vector.
    ///
    /// The magnitude is defined as `||V|| = sqrt(x^2 + y^2 + z^2)`.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic magnitude.
    #[must_use]

    pub fn magnitude(&self) -> Expr {

        simplify(&Expr::new_sqrt(
            Expr::new_add(
                Expr::new_add(
                    Expr::new_pow(
                        self.x.clone(),
                        Expr::new_bigint(BigInt::from(2)),
                    ),
                    Expr::new_pow(
                        self.y.clone(),
                        Expr::new_bigint(BigInt::from(2)),
                    ),
                ),
                Expr::new_pow(
                    self.z.clone(),
                    Expr::new_bigint(BigInt::from(2)),
                ),
            ),
        ))
    }

    /// Computes the dot product of this vector with another vector.
    ///
    /// The dot product is defined as `V1 . V2 = x1*x2 + y1*y2 + z1*z2`.
    ///
    /// # Arguments
    /// * `other` - The other vector to compute the dot product with.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic dot product.
    #[must_use]

    pub fn dot(
        &self,
        other: &Self,
    ) -> Expr {

        simplify(&Expr::new_add(
            Expr::new_add(
                Expr::new_mul(
                    self.x.clone(),
                    other.x.clone(),
                ),
                Expr::new_mul(
                    self.y.clone(),
                    other.y.clone(),
                ),
            ),
            Expr::new_mul(
                self.z.clone(),
                other.z.clone(),
            ),
        ))
    }

    /// Computes the cross product of this vector with another vector.
    ///
    /// The cross product `V1 x V2` results in a new vector that is perpendicular
    /// to both `V1` and `V2`. The components are calculated as:
    /// - `x = y1*z2 - z1*y2`
    /// - `y = z1*x2 - x1*z2`
    /// - `z = x1*y2 - y1*x2`
    ///
    /// # Arguments
    /// * `other` - The other vector to compute the cross product with.
    ///
    /// # Returns
    /// A new `Vector` representing the symbolic cross product.
    #[must_use]

    pub fn cross(
        &self,
        other: &Self,
    ) -> Self {

        let x_comp =
            simplify(&Expr::new_sub(
                Expr::new_mul(
                    self.y.clone(),
                    other.z.clone(),
                ),
                Expr::new_mul(
                    self.z.clone(),
                    other.y.clone(),
                ),
            ));

        let y_comp =
            simplify(&Expr::new_sub(
                Expr::new_mul(
                    self.z.clone(),
                    other.x.clone(),
                ),
                Expr::new_mul(
                    self.x.clone(),
                    other.z.clone(),
                ),
            ));

        let z_comp =
            simplify(&Expr::new_sub(
                Expr::new_mul(
                    self.x.clone(),
                    other.y.clone(),
                ),
                Expr::new_mul(
                    self.y.clone(),
                    other.x.clone(),
                ),
            ));

        Self::new(
            x_comp,
            y_comp,
            z_comp,
        )
    }

    /// Normalizes the vector to have a magnitude of 1.
    ///
    /// This is achieved by dividing each component of the vector by its magnitude.
    /// If the vector has a magnitude of zero, it is returned unchanged.
    ///
    /// # Returns
    /// A new `Vector` representing the normalized vector.
    #[must_use]

    pub fn normalize(&self) -> Self {

        let mag = self.magnitude();

        if is_zero(&mag) {

            return self.clone();
        }

        self.scalar_mul(&Expr::new_div(
            Expr::new_bigint(BigInt::one()),
            mag,
        ))
    }

    /// Multiplies the vector by a scalar expression.
    ///
    /// Each component of the vector is multiplied by the given scalar.
    ///
    /// # Arguments
    /// * `scalar` - The `Expr` to multiply the vector by.
    ///
    /// # Returns
    /// A new `Vector` representing the result of the scalar multiplication.
    #[must_use]

    pub fn scalar_mul(
        &self,
        scalar: &Expr,
    ) -> Self {

        Self::new(
            simplify(&Expr::new_mul(
                scalar.clone(),
                self.x.clone(),
            )),
            simplify(&Expr::new_mul(
                scalar.clone(),
                self.y.clone(),
            )),
            simplify(&Expr::new_mul(
                scalar.clone(),
                self.z.clone(),
            )),
        )
    }

    /// Computes the angle between this vector and another vector.
    ///
    /// The angle is calculated using the dot product formula: `cos(theta) = (A . B) / (|A| |B|)`.
    /// Thus, `theta = acos((A . B) / (|A| |B|))`.
    ///
    /// # Arguments
    /// * `other` - The other vector.
    ///
    /// # Returns
    /// An `Expr` representing the symbolic angle in radians.
    #[must_use]

    pub fn angle(
        &self,
        other: &Self,
    ) -> Expr {

        let dot_prod = self.dot(other);

        let mag_self = self.magnitude();

        let mag_other =
            other.magnitude();

        let cos_theta =
            simplify(&Expr::new_div(
                dot_prod,
                Expr::new_mul(
                    mag_self,
                    mag_other,
                ),
            ));

        simplify(&Expr::new_arccos(
            cos_theta,
        ))
    }

    /// Computes the projection of this vector onto another vector.
    ///
    /// The projection of A onto B is given by: `proj_B(A) = ((A . B) / |B|^2) * B`.
    ///
    /// # Arguments
    /// * `other` - The vector to project onto (B).
    ///
    /// # Returns
    /// A new `Vector` representing the projection of self onto other.
    #[must_use]

    pub fn project_onto(
        &self,
        other: &Self,
    ) -> Self {

        let dot_prod = self.dot(other);

        let mag_other_sq =
            other.dot(other); // |B|^2 = B . B

        if is_zero(&mag_other_sq) {

            // Projection onto zero vector is zero vector
            return Self::new(
                Expr::new_constant(0.0),
                Expr::new_constant(0.0),
                Expr::new_constant(0.0),
            );
        }

        let scalar =
            simplify(&Expr::new_div(
                dot_prod,
                mag_other_sq,
            ));

        other.scalar_mul(&scalar)
    }

    /// Converts the `Vector` into a `Expr::Vector` variant.
    ///
    /// # Returns
    /// An `Expr` that contains the vector's components.
    #[must_use]

    pub fn to_expr(&self) -> Expr {

        Expr::Vector(vec![
            self.x.clone(),
            self.y.clone(),
            self.z.clone(),
        ])
    }
}

/// Overloads the '+' operator for Vector addition.

impl Add for Vector {
    type Output = Self;

    fn add(
        self,
        other: Self,
    ) -> Self {

        Self::new(
            simplify(&Expr::new_add(
                self.x,
                other.x,
            )),
            simplify(&Expr::new_add(
                self.y,
                other.y,
            )),
            simplify(&Expr::new_add(
                self.z,
                other.z,
            )),
        )
    }
}

/// Overloads the '-' operator for Vector subtraction.

impl Sub for Vector {
    type Output = Self;

    fn sub(
        self,
        other: Self,
    ) -> Self {

        Self::new(
            simplify(&Expr::new_sub(
                self.x,
                other.x,
            )),
            simplify(&Expr::new_sub(
                self.y,
                other.y,
            )),
            simplify(&Expr::new_sub(
                self.z,
                other.z,
            )),
        )
    }
}

/// Computes the gradient of a scalar field `f(x, y, z)`.
///
/// The gradient is a vector field that points in the direction of the greatest rate of
/// increase of the scalar field, and its magnitude is the rate of increase.
/// It is defined as `grad(f) = (df/dx, df/dy, df/dz)`.
///
/// # Arguments
/// * `scalar_field` - An `Expr` representing the scalar function `f`.
/// * `vars` - A tuple of string slices `("x", "y", "z")` representing the variables.
///
/// # Returns
/// A `Vector` representing the symbolic gradient of the scalar field.
#[must_use]

pub fn gradient(
    scalar_field: &Expr,
    vars: (&str, &str, &str),
) -> Vector {

    let df_dx = differentiate(
        scalar_field,
        vars.0,
    );

    let df_dy = differentiate(
        scalar_field,
        vars.1,
    );

    let df_dz = differentiate(
        scalar_field,
        vars.2,
    );

    Vector::new(df_dx, df_dy, df_dz)
}

/// Computes the divergence of a vector field `F = (Fx, Fy, Fz)`.
///
/// The divergence measures the magnitude of a vector field's source or sink at a given point.
/// It is a scalar quantity defined as `div(F) = dFx/dx + dFy/dy + dFz/dz`.
///
/// # Arguments
/// * `vector_field` - The `Vector` representing the vector field `F`.
/// * `vars` - A tuple of string slices `("x", "y", "z")` representing the variables.
///
/// # Returns
/// An `Expr` representing the symbolic divergence of the vector field.
#[must_use]

pub fn divergence(
    vector_field: &Vector,
    vars: (&str, &str, &str),
) -> Expr {

    let d_fx_dx = differentiate(
        &vector_field.x,
        vars.0,
    );

    let d_fy_dy = differentiate(
        &vector_field.y,
        vars.1,
    );

    let d_fz_dz = differentiate(
        &vector_field.z,
        vars.2,
    );

    simplify(&Expr::new_add(
        Expr::new_add(d_fx_dx, d_fy_dy),
        d_fz_dz,
    ))
}

/// Computes the curl of a vector field `F = (Fx, Fy, Fz)`.
///
/// The curl measures the infinitesimal rotation of a 3D vector field.
/// It is a vector field defined as:
/// `curl(F) = (dFz/dy - dFy/dz, dFx/dz - dFz/dx, dFy/dx - dFx/dy)`
///
/// # Arguments
/// * `vector_field` - The `Vector` representing the vector field `F`.
/// * `vars` - A tuple of string slices `("x", "y", "z")` representing the variables.
///
/// # Returns
/// A `Vector` representing the symbolic curl of the vector field.
#[must_use]

pub fn curl(
    vector_field: &Vector,
    vars: (&str, &str, &str),
) -> Vector {

    let d_fz_dy = differentiate(
        &vector_field.z,
        vars.1,
    );

    let d_fy_dz = differentiate(
        &vector_field.y,
        vars.2,
    );

    let d_fx_dz = differentiate(
        &vector_field.x,
        vars.2,
    );

    let d_fz_dx = differentiate(
        &vector_field.z,
        vars.0,
    );

    let d_fy_dx = differentiate(
        &vector_field.y,
        vars.0,
    );

    let d_fx_dy = differentiate(
        &vector_field.x,
        vars.1,
    );

    let x_comp =
        simplify(&Expr::new_sub(
            d_fz_dy,
            d_fy_dz,
        ));

    let y_comp =
        simplify(&Expr::new_sub(
            d_fx_dz,
            d_fz_dx,
        ));

    let z_comp =
        simplify(&Expr::new_sub(
            d_fy_dx,
            d_fx_dy,
        ));

    Vector::new(
        x_comp,
        y_comp,
        z_comp,
    )
}

/// Computes the directional derivative of a scalar field `f` in the direction of a vector `v`.
///
/// The directional derivative represents the rate of change of the function `f`
/// along the direction of `v`. It is calculated as the dot product of the gradient of `f`
/// and the direction vector `v`: `D_v(f) = grad(f) . v`.
///
/// # Arguments
/// * `scalar_field` - The scalar function `f` as an `Expr`.
/// * `direction` - The `Vector` specifying the direction.
/// * `vars` - A tuple of string slices `("x", "y", "z")` representing the variables.
///
/// # Returns
/// An `Expr` representing the symbolic directional derivative.
#[must_use]

pub fn directional_derivative(
    scalar_field: &Expr,
    direction: &Vector,
    vars: (&str, &str, &str),
) -> Expr {

    let grad_f =
        gradient(scalar_field, vars);

    grad_f.dot(direction)
}

/// Computes the partial derivative of a vector field with respect to a single variable.
///
/// This is done by taking the partial derivative of each component of the vector field
/// with respect to the given variable.
///
/// # Arguments
/// * `vector_field` - The `Vector` to differentiate.
/// * `var` - The variable to differentiate with respect to.
///
/// # Returns
/// A new `Vector` where each component is the partial derivative of the original component.
#[must_use]

pub fn partial_derivative_vector(
    vector_field: &Vector,
    var: &str,
) -> Vector {

    Vector::new(
        differentiate(
            &vector_field.x,
            var,
        ),
        differentiate(
            &vector_field.y,
            var,
        ),
        differentiate(
            &vector_field.z,
            var,
        ),
    )
}
