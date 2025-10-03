//! # Symbolic Vector Calculus Operations
//!
//! This module provides functions for performing symbolic vector calculus operations,
//! including line integrals (scalar and vector fields), surface integrals (flux),
//! and volume integrals. It defines structures for `ParametricCurve`, `ParametricSurface`,
//! and `Volume` to represent the domains of integration.

use crate::symbolic::calculus::{definite_integrate, substitute};
use crate::symbolic::core::Expr;
use crate::symbolic::simplify::simplify;
use crate::symbolic::vector::partial_derivative_vector;
use crate::symbolic::vector::Vector;

/// Represents a parameterized curve C given by r(t).
pub struct ParametricCurve {
    /// The vector expression for the curve, e.g., [cos(t), sin(t), t].
    pub r: Vector,
    /// The name of the parameter, e.g., "t".
    pub t_var: String,
    /// The integration bounds for the parameter, e.g., (0, 2*pi).
    pub t_bounds: (Expr, Expr),
}

/// Represents a parameterized surface S given by r(u, v).
pub struct ParametricSurface {
    /// The vector expression for the surface, e.g., [u*cos(v), u*sin(v), v].
    pub r: Vector,
    /// The name of the first parameter, e.g., "u".
    pub u_var: String,
    /// The integration bounds for the first parameter, e.g., (0, 1).
    pub u_bounds: (Expr, Expr),
    /// The name of the second parameter, e.g., "v".
    pub v_var: String,
    /// The integration bounds for the second parameter, e.g., (0, 2*pi).
    pub v_bounds: (Expr, Expr),
}

/// Represents a volume V for triple integration.
/// Defines the integration order as dz dy dx.
pub struct Volume {
    /// The bounds for the innermost integral (dz). Can be expressions in terms of x and y.
    pub z_bounds: (Expr, Expr),
    /// The bounds for the middle integral (dy). Can be expressions in terms of x.
    pub y_bounds: (Expr, Expr),
    /// The bounds for the outermost integral (dx). Must be constants.
    pub x_bounds: (Expr, Expr),
    /// The variable names for (x, y, z).
    pub vars: (String, String, String),
}

/// Computes the line integral of a scalar field `f` along a parameterized curve C.
///
/// The integral is `∫_C f ds = ∫_a^b f(r(t)) ||r'(t)|| dt`.
///
/// # Arguments
/// * `scalar_field` - The scalar field `f` as an `Expr`.
/// * `curve` - The `ParametricCurve` representing the path of integration.
///
/// # Returns
/// An `Expr` representing the symbolic line integral.
pub fn line_integral_scalar(scalar_field: &Expr, curve: &ParametricCurve) -> Expr {
    let r_prime = partial_derivative_vector(&curve.r, &curve.t_var);
    let r_prime_magnitude = r_prime.magnitude();

    let sub = |expr: &Expr| {
        let e1 = substitute(expr, "x", &curve.r.x);
        let e2 = substitute(&e1, "y", &curve.r.y);
        substitute(&e2, "z", &curve.r.z)
    };
    let field_on_curve = sub(scalar_field);

    let integrand = simplify(Expr::Mul(
        Box::new(field_on_curve),
        Box::new(r_prime_magnitude),
    ));
    let integral = definite_integrate(
        &integrand,
        &curve.t_var,
        &curve.t_bounds.0,
        &curve.t_bounds.1,
    );

    simplify(integral)
}

/// Computes the line integral of a vector field `F` along a parameterized curve C (work).
///
/// The integral is `∫_C F · dr = ∫_a^b F(r(t)) · r'(t) dt`.
///
/// # Arguments
/// * `vector_field` - The vector field `F` as a `Vector`.
/// * `curve` - The `ParametricCurve` representing the path of integration.
///
/// # Returns
/// An `Expr` representing the symbolic line integral.
pub fn line_integral_vector(vector_field: &Vector, curve: &ParametricCurve) -> Expr {
    let r_prime = partial_derivative_vector(&curve.r, &curve.t_var);

    let sub = |expr: &Expr| {
        let e1 = substitute(expr, "x", &curve.r.x);
        let e2 = substitute(&e1, "y", &curve.r.y);
        substitute(&e2, "z", &curve.r.z)
    };
    let field_on_curve = Vector::new(
        sub(&vector_field.x),
        sub(&vector_field.y),
        sub(&vector_field.z),
    );

    let integrand = field_on_curve.dot(&r_prime);
    let integral = definite_integrate(
        &integrand,
        &curve.t_var,
        &curve.t_bounds.0,
        &curve.t_bounds.1,
    );

    simplify(integral)
}

/// Computes the surface integral (flux) of a vector field F over a parameterized surface S.
///
/// The integral is `∫∫_S (F · dS) = ∫∫_D F(r(u,v)) · (r_u × r_v) du dv`.
///
/// # Arguments
/// * `field` - The vector field `F` as a `Vector`.
/// * `surface` - The `ParametricSurface` representing the surface of integration.
///
/// # Returns
/// An `Expr` representing the symbolic surface integral.
pub fn surface_integral(field: &Vector, surface: &ParametricSurface) -> Expr {
    let r_u = partial_derivative_vector(&surface.r, &surface.u_var);
    let r_v = partial_derivative_vector(&surface.r, &surface.v_var);
    let normal_vector = r_u.cross(&r_v);

    let sub = |expr: &Expr| {
        let e1 = substitute(expr, "x", &surface.r.x);
        let e2 = substitute(&e1, "y", &surface.r.y);
        substitute(&e2, "z", &surface.r.z)
    };
    let field_on_surface = Vector::new(sub(&field.x), sub(&field.y), sub(&field.z));

    let integrand = field_on_surface.dot(&normal_vector);

    let inner_integral = definite_integrate(
        &integrand,
        &surface.u_var,
        &surface.u_bounds.0,
        &surface.u_bounds.1,
    );
    let outer_integral = definite_integrate(
        &inner_integral,
        &surface.v_var,
        &surface.v_bounds.0,
        &surface.v_bounds.1,
    );

    simplify(outer_integral)
}

/// Computes the volume integral of a scalar field `f` over a defined volume V.
///
/// The integral is `∫∫∫_V f dV = ∫_x1^x2 ∫_y1^y2 ∫_z1^z2 f(x,y,z) dz dy dx`.
///
/// # Arguments
/// * `scalar_field` - The scalar field `f` as an `Expr`.
/// * `volume` - The `Volume` structure defining the integration domain and order.
///
/// # Returns
/// An `Expr` representing the symbolic volume integral.
pub fn volume_integral(scalar_field: &Expr, volume: &Volume) -> Expr {
    let (x_var, y_var, z_var) = (&volume.vars.0, &volume.vars.1, &volume.vars.2);

    let integral_z =
        definite_integrate(scalar_field, z_var, &volume.z_bounds.0, &volume.z_bounds.1);
    let integral_y = definite_integrate(&integral_z, y_var, &volume.y_bounds.0, &volume.y_bounds.1);
    let integral_x = definite_integrate(&integral_y, x_var, &volume.x_bounds.0, &volume.x_bounds.1);

    simplify(integral_x)
}
