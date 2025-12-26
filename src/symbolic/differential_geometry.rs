//! # Differential Geometry
//!
//! This module provides symbolic tools for calculus on manifolds, including the
//! manipulation of differential forms and the application of major integral theorems.
//! It includes implementations for exterior derivatives, wedge products, and symbolic
//! representations of generalized Stokes' theorem, Gauss's theorem, and Green's theorem.

use crate::symbolic::calculus::{definite_integrate, differentiate};
use crate::symbolic::core::Expr;
use crate::symbolic::simplify_dag::simplify;
use crate::symbolic::vector::Vector;
use num_bigint::BigInt;
use num_traits::Zero;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

/// Represents a differential k-form.
///
/// A k-form is a mathematical object that can be integrated over a k-dimensional manifold.
/// It is a sum of terms, where each term is a scalar function (coefficient) multiplied by
/// a wedge product of k basis 1-forms (like dx, dy, etc.).
///
/// For example, a 2-form in R^3 could be `f(x,y,z) dx^dy + g(x,y,z) dx^dz`.
///
/// Here, the basis wedge products (e.g., `dx^dy`) are represented by a bitmask (`blade`).
/// If `vars = ["x", "y", "z"]`, then `dx` is `1<<0`, `dy` is `1<<1`, `dz` is `1<<2`.
/// The wedge product `dx^dy` corresponds to the bitmask `(1<<0) | (1<<1) = 3`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]

pub struct DifferentialForm {
    /// A map from the basis wedge product (represented by a bitmask) to its coefficient expression.
    pub terms: std::collections::BTreeMap<u32, Expr>,
}

/// Computes the exterior derivative of a k-form, resulting in a (k+1)-form.
///
/// The exterior derivative `dω` of a k-form `ω` is a central concept in differential geometry.
/// For a 0-form (a scalar function `f`), `df` is its total differential.
/// For a k-form `ω = f dx_I`, `dω = df ^ dx_I`.
///
/// # Arguments
/// * `form` - The differential form `ω` to differentiate.
/// * `vars` - A slice of variable names (e.g., `["x", "y", "z"]`) defining the basis.
///
/// # Returns
/// A new `DifferentialForm` representing `dω`.
#[must_use]

pub fn exterior_derivative(form: &DifferentialForm, vars: &[&str]) -> DifferentialForm {

    let mut result_terms = std::collections::BTreeMap::new();

    for (blade, coeff) in &form.terms {

        for (i, _item) in vars.iter().enumerate() {

            let new_blade = (1 << i) | blade;

            if new_blade.count_ones() != blade.count_ones() + 1 {

                continue;
            }

            let d_coeff = differentiate(coeff, vars[i]);

            let mut sign = 1i64;

            for j in 0..i {

                if (blade >> j) & 1 == 1 {

                    sign *= -1;
                }
            }

            let signed_coeff = if sign == 1 {

                d_coeff
            } else {

                simplify(&Expr::new_neg(d_coeff))
            };

            let entry = result_terms
                .entry(new_blade)
                .or_insert(Expr::BigInt(BigInt::zero()));

            *entry = simplify(&Expr::new_add(entry.clone(), signed_coeff));
        }
    }

    DifferentialForm {
        terms: result_terms,
    }
}

/// Computes the wedge product (exterior product) of two differential forms.
///
/// The wedge product is an antisymmetric, associative product.
/// For basis 1-forms, `dx_i ^ dx_j = -dx_j ^ dx_i`, and `dx_i ^ dx_i = 0`.
///
/// # Arguments
/// * `form1` - The first differential form.
/// * `form2` - The second differential form.
///
/// # Returns
/// A new `DifferentialForm` representing the wedge product.
#[must_use]

pub fn wedge_product(form1: &DifferentialForm, form2: &DifferentialForm) -> DifferentialForm {

    let mut result_terms = std::collections::BTreeMap::new();

    for (blade1, coeff1) in &form1.terms {

        for (blade2, coeff2) in &form2.terms {

            if (blade1 & blade2) != 0 {

                continue;
            }

            let new_blade = blade1 | blade2;

            let mut sign = 1i64;

            let mut temp_blade2 = *blade2;

            while temp_blade2 > 0 {

                let i = temp_blade2.trailing_zeros();

                let swaps = (blade1 >> (i + 1)).count_ones();

                if swaps % 2 != 0 {

                    sign *= -1;
                }

                temp_blade2 &= !(1 << i);
            }

            let new_coeff = simplify(&Expr::new_mul(coeff1.clone(), coeff2.clone()));

            let signed_coeff = if sign == 1 {

                new_coeff
            } else {

                simplify(&Expr::new_neg(new_coeff))
            };

            let entry = result_terms
                .entry(new_blade)
                .or_insert(Expr::BigInt(BigInt::zero()));

            *entry = simplify(&Expr::new_add(entry.clone(), signed_coeff));
        }
    }

    DifferentialForm {
        terms: result_terms,
    }
}

/// Represents the boundary of a domain, denoted as `∂M` for a manifold `M`.
/// This is a symbolic representation used in the integral theorems.
#[must_use]

pub fn boundary(domain: &Expr) -> Expr {

    Expr::Boundary(Arc::new(domain.clone()))
}

/// Represents the generalized Stokes' Theorem.
///
/// The theorem states that the integral of the exterior derivative of a differential form `ω`
/// over some manifold `M` is equal to the integral of `ω` over the boundary of `M` (`∂M`).
/// Formula: `∫_M dω = ∫_{∂M} ω`
/// This function returns a symbolic equation representing the theorem.
#[must_use]

pub fn generalized_stokes_theorem(
    omega: &DifferentialForm,
    manifold: &Expr,
    vars: &[&str],
) -> Expr {

    let d_omega = exterior_derivative(omega, vars);

    let integral_d_omega = Expr::Integral {
        integrand: Arc::new(Expr::Variable(format!("{d_omega:?}"))),
        var: Arc::new(Expr::Variable(manifold.to_string())),
        lower_bound: Arc::new(Expr::Variable("M".to_string())),
        upper_bound: Arc::new(Expr::BigInt(BigInt::zero())),
    };

    let integral_omega = Expr::Integral {
        integrand: Arc::new(Expr::Variable(format!("{omega:?}"))),
        var: Arc::new(Expr::Variable(manifold.to_string())),
        lower_bound: Arc::new(boundary(manifold)),
        upper_bound: Arc::new(Expr::BigInt(BigInt::zero())),
    };

    Expr::Eq(Arc::new(integral_d_omega), Arc::new(integral_omega))
}

/// Represents Gauss's Theorem (Divergence Theorem) as a special case of Stokes' Theorem.
///
/// It relates the flux of a vector field through a closed surface to the divergence
/// of the field in the volume enclosed.
/// Formula: `∫_V (∇ · F) dV = ∮_{∂V} (F · n) dS`
/// This function returns a symbolic equation representing the theorem.
#[must_use]

pub fn gauss_theorem(vector_field: &Vector, volume: &Expr) -> Expr {

    let div_f = super::vector::divergence(vector_field, ("x", "y", "z"));

    let integral_div = Expr::VolumeIntegral {
        scalar_field: Arc::new(div_f),
        volume: Arc::new(volume.clone()),
    };

    let surface_integral = Expr::SurfaceIntegral {
        vector_field: Arc::new(Expr::Variable("F".to_string())),
        surface: Arc::new(boundary(volume)),
    };

    Expr::Eq(Arc::new(integral_div), Arc::new(surface_integral))
}

/// Represents the classical Stokes' Theorem as a special case of the generalized theorem.
///
/// It relates the integral of the curl of a vector field over a surface `S` to the
/// line integral of the vector field over its boundary `∂S`.
/// Formula: `∫_S (∇ × F) · dS = ∮_{∂S} F · dr`
/// This function returns a symbolic equation representing the theorem.
#[must_use]

pub fn stokes_theorem(vector_field: &Vector, surface: &Expr) -> Expr {

    let curl_f = super::vector::curl(vector_field, ("x", "y", "z"));

    let integral_curl = Expr::SurfaceIntegral {
        vector_field: Arc::new(curl_f.to_expr()),
        surface: Arc::new(surface.clone()),
    };

    let line_integral = Expr::Integral {
        integrand: Arc::new(Expr::Variable("F · dr".to_string())),
        var: Arc::new(Expr::Variable("t".to_string())),
        lower_bound: Arc::new(boundary(surface)),
        upper_bound: Arc::new(Expr::BigInt(BigInt::zero())),
    };

    Expr::Eq(Arc::new(integral_curl), Arc::new(line_integral))
}

/// Represents Green's Theorem as a 2D special case of Stokes' Theorem.
///
/// It relates a line integral around a simple closed curve `C` to a double integral
/// over the plane region `D` bounded by `C`.
/// Formula: `∬_D (∂Q/∂x - ∂P/∂y) dA = ∮_{∂D} P dx + Q dy`
/// This function returns a symbolic equation representing the theorem.
#[must_use]

pub fn greens_theorem(p: &Expr, q: &Expr, domain: &Expr) -> Expr {

    let dq_dx = differentiate(q, "x");

    let dp_dy = differentiate(p, "y");

    let integrand_da = simplify(&Expr::new_sub(dq_dx, dp_dy));

    let integral_da = definite_integrate(
        &integrand_da,
        "A",
        &Expr::Domain(format!("{domain}")),
        &Expr::BigInt(BigInt::zero()),
    );

    let integrand_line = Expr::new_add(
        Expr::new_mul(p.clone(), Expr::Variable("dx".to_string())),
        Expr::new_mul(q.clone(), Expr::Variable("dy".to_string())),
    );

    let line_integral = Expr::Integral {
        integrand: Arc::new(integrand_line),
        var: Arc::new(Expr::Variable("t".to_string())),
        lower_bound: Arc::new(boundary(domain)),
        upper_bound: Arc::new(Expr::BigInt(BigInt::zero())),
    };

    Expr::Eq(Arc::new(integral_da), Arc::new(line_integral))
}
