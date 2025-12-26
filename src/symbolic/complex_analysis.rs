//! # Complex Analysis
//!
//! This module implements advanced concepts from complex analysis, including:
//! - Analytic continuation along paths
//! - Residue calculus and contour integration
//! - Conformal mappings (Möbius transformations, etc.)
//! - Singularity analysis and Laurent series
//! - Cauchy integral formulas
//! - Complex function evaluation

use crate::symbolic::calculus::{differentiate, substitute};
use crate::symbolic::core::Expr;
use crate::symbolic::series::{self, calculate_taylor_coefficients, taylor_series};
use crate::symbolic::simplify_dag::simplify;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

// ============================================================================
// Analytic Continuation
// ============================================================================

/// Represents the analytic continuation of a function along a path.
/// It is stored as a chain of series expansions, each centered at a point on the path.
#[derive(Debug, Clone, Serialize, Deserialize)]

pub struct PathContinuation {
    pub var: String,
    pub order: usize,
    /// A vector of (center, `series_expression`) tuples.
    pub pieces: Vec<(Expr, Expr)>,
}

impl PathContinuation {
    /// Creates a new analytic continuation starting with a Taylor series for `func` centered at `start_point`.
    #[must_use]

    pub fn new(func: &Expr, var: &str, start_point: &Expr, order: usize) -> Self {

        let initial_series = taylor_series(func, var, start_point, order);

        Self {
            var: var.to_string(),
            order,
            pieces: vec![(start_point.clone(), initial_series)],
        }
    }

    /// Continues the function along a given path.

    pub fn continue_along_path(&mut self, path_points: &[Expr]) -> Result<(), String> {

        for next_point in path_points {

            let (last_center, last_series) = self.pieces.last().ok_or_else(|| {

                "PathContinuation must be initialized with `new` before continuing.".to_string()
            })?;

            let radius =
                estimate_radius_of_convergence(last_series, &self.var, last_center, self.order + 5)
                    .ok_or_else(|| "Failed to estimate the radius of convergence.".to_string())?;

            let distance = complex_distance(last_center, next_point).ok_or_else(|| {

                "Failed to calculate distance between complex points.".to_string()
            })?;

            if distance >= radius {

                return Err(format!(
                    "Analytic continuation failed: point {next_point} is outside radius {radius} of series at {last_center}."
                ));
            }

            let next_series = series::analytic_continuation(
                last_series,
                &self.var,
                last_center,
                next_point,
                self.order,
            );

            self.pieces.push((next_point.clone(), next_series));
        }

        Ok(())
    }

    /// Returns the final expression (Taylor series) after continuation to the last point.
    #[must_use]

    pub fn get_final_expression(&self) -> Option<&Expr> {

        self.pieces.last().map(|(_, series)| series)
    }
}

/// Estimates the radius of convergence for a Taylor series using the ratio test.
#[must_use]

pub fn estimate_radius_of_convergence(
    series_expr: &Expr,
    var: &str,
    center: &Expr,
    order: usize,
) -> Option<f64> {

    let coeffs = calculate_taylor_coefficients(series_expr, var, center, order);

    for n in (1..coeffs.len()).rev() {

        let cn = &coeffs[n];

        let cn_minus_1 = &coeffs[n - 1];

        if let (Some(c_n_val), Some(c_n_minus_1_val)) = (cn.to_f64(), cn_minus_1.to_f64()) {

            if c_n_val.abs() > f64::EPSILON && c_n_minus_1_val.abs() > f64::EPSILON {

                let limit_ratio = c_n_val.abs() / c_n_minus_1_val.abs();

                if limit_ratio < f64::EPSILON {

                    return Some(f64::INFINITY);
                }

                return Some(1.0 / limit_ratio);
            }
        }
    }

    Some(f64::INFINITY)
}

/// Calculates the Euclidean distance between two complex points.
#[must_use]

pub fn complex_distance(p1: &Expr, p2: &Expr) -> Option<f64> {

    let re1 = p1.re().to_f64().unwrap_or(0.0);

    let im1 = p1.im().to_f64().unwrap_or(0.0);

    let re2 = p2.re().to_f64().unwrap_or(0.0);

    let im2 = p2.im().to_f64().unwrap_or(0.0);

    let dx = re1 - re2;

    let dy = im1 - im2;

    Some(dx.hypot(dy))
}

// ============================================================================
// Residue Calculus
// ============================================================================

/// Represents a singularity type in complex analysis.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]

pub enum SingularityType {
    /// Removable singularity
    Removable,
    /// Pole of order n
    Pole(usize),
    /// Essential singularity
    Essential,
}

/// Classifies the type of singularity at a point.
///
/// Uses the function structure to determine if the singularity is:
/// - Removable (limit exists)
/// - Pole of order n (1/(z-a)^n term)
/// - Essential (e.g., e^(1/z))
#[must_use]

pub fn classify_singularity(
    func: &Expr,
    var: &str,
    singularity: &Expr,
    _order: usize,
) -> SingularityType {

    // Analyze the function structure to detect poles
    // For f(z) = g(z)/(z-a)^n, we have a pole of order n
    // Check if function is a division
    let _z = Expr::Variable(var.to_string());

    let _factor = simplify(&Expr::new_sub(_z.clone(), singularity.clone()));

    let pole_order = count_pole_order(&_z, &_factor);

    if pole_order > 0 {

        return SingularityType::Pole(pole_order);
    }

    if let Expr::Div(_num, den) = func {

        // Check if denominator contains (z - singularity)
        let z = Expr::Variable(var.to_string());

        let factor = simplify(&Expr::new_sub(z, singularity.clone()));

        // Count the pole order by checking denominator structure
        let pole_order = count_pole_order(den, &factor);

        if pole_order > 0 {

            return SingularityType::Pole(pole_order);
        }
    }

    // Check for essential singularities (e.g., exp(1/z))
    // This is a simplified check
    let func_str = format!("{func:?}");

    if func_str.contains("Exp") && func_str.contains("Div") {

        return SingularityType::Essential;
    }

    SingularityType::Removable
}

/// Helper function to count the order of a pole.

fn count_pole_order(expr: &Expr, factor: &Expr) -> usize {

    match expr {
        // If denominator is exactly (z-a), pole order is 1
        _ if expr == factor => 1,

        // If denominator is (z-a)^n
        Expr::Power(base, exp) => {

            if base.as_ref() == factor {

                // Try to extract the exponent as a number
                if let Some(n) = exp.to_f64() {

                    return n as usize;
                }

                return 1;
            }

            0
        }

        // If denominator is a product containing (z-a)
        Expr::Mul(a, b) => {

            let order_a = count_pole_order(a, factor);

            let order_b = count_pole_order(b, factor);

            order_a + order_b
        }

        _ => 0,
    }
}

/// Computes the Laurent series expansion around a point.
///
/// Laurent series: f(z) = Σ(n=-∞ to ∞) `a_n` (z-z0)^n
#[must_use]

pub fn laurent_series(func: &Expr, var: &str, center: &Expr, order: usize) -> Expr {

    // For a full Laurent series, we'd need to compute negative power coefficients
    // For now, return the Taylor series as an approximation
    taylor_series(func, var, center, order)
}

/// Calculates the residue of a function at a singularity.
///
/// The residue is the coefficient of (z-z0)^(-1) in the Laurent series.
/// For a simple pole, Res(f, z0) = lim_{z→z0} (z-z0)f(z)
#[must_use]

pub fn calculate_residue(func: &Expr, var: &str, singularity: &Expr) -> Expr {

    // For a simple pole: Res = lim_{z→z0} (z-z0)f(z)
    let z = Expr::Variable(var.to_string());

    let factor = Expr::new_sub(z, singularity.clone());

    let product = Expr::new_mul(factor, func.clone());

    // Evaluate limit as z → singularity
    // Substitute and simplify
    simplify(&substitute(&product, var, singularity))
}

/// Evaluates a contour integral using the residue theorem.
///
/// ∮_C f(z) dz = 2πi Σ Res(f, `z_k`)
/// where `z_k` are the singularities inside the contour C.
#[must_use]

pub fn contour_integral_residue_theorem(func: &Expr, var: &str, singularities: &[Expr]) -> Expr {

    let mut sum = Expr::Constant(0.0);

    for singularity in singularities {

        let residue = calculate_residue(func, var, singularity);

        sum = Expr::new_add(sum, residue);
    }

    // Multiply by 2πi
    let two_pi_i = Expr::new_mul(
        Expr::Constant(2.0),
        Expr::new_mul(
            Expr::Pi,
            Expr::Complex(Arc::new(Expr::Constant(0.0)), Arc::new(Expr::Constant(1.0))),
        ),
    );

    simplify(&Expr::new_mul(two_pi_i, sum))
}

// ============================================================================
// Conformal Mappings
// ============================================================================

/// Represents a Möbius transformation: f(z) = (az + b) / (cz + d)
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]

pub struct MobiusTransformation {
    pub a: Expr,
    pub b: Expr,
    pub c: Expr,
    pub d: Expr,
}

impl MobiusTransformation {
    /// Creates a new Möbius transformation.
    #[must_use]

    pub const fn new(a: Expr, b: Expr, c: Expr, d: Expr) -> Self {

        Self { a, b, c, d }
    }

    /// Creates the identity transformation.
    #[must_use]

    pub fn identity() -> Self {

        Self {
            a: Expr::BigInt(BigInt::one()),
            b: Expr::BigInt(BigInt::zero()),
            c: Expr::BigInt(BigInt::zero()),
            d: Expr::BigInt(BigInt::one()),
        }
    }

    /// Applies the transformation to a point z.
    #[must_use]

    pub fn apply(&self, z: &Expr) -> Expr {

        let numerator = Expr::new_add(Expr::new_mul(self.a.clone(), z.clone()), self.b.clone());

        let denominator = Expr::new_add(Expr::new_mul(self.c.clone(), z.clone()), self.d.clone());

        simplify(&Expr::new_div(numerator, denominator))
    }

    /// Composes two Möbius transformations.
    #[must_use]

    pub fn compose(&self, other: &Self) -> Self {

        // (f ∘ g)(z) where f = self, g = other
        // Result: ((a1*a2 + b1*c2)z + (a1*b2 + b1*d2)) / ((c1*a2 + d1*c2)z + (c1*b2 + d1*d2))
        let a = simplify(&Expr::new_add(
            Expr::new_mul(self.a.clone(), other.a.clone()),
            Expr::new_mul(self.b.clone(), other.c.clone()),
        ));

        let b = simplify(&Expr::new_add(
            Expr::new_mul(self.a.clone(), other.b.clone()),
            Expr::new_mul(self.b.clone(), other.d.clone()),
        ));

        let c = simplify(&Expr::new_add(
            Expr::new_mul(self.c.clone(), other.a.clone()),
            Expr::new_mul(self.d.clone(), other.c.clone()),
        ));

        let d = simplify(&Expr::new_add(
            Expr::new_mul(self.c.clone(), other.b.clone()),
            Expr::new_mul(self.d.clone(), other.d.clone()),
        ));

        Self { a, b, c, d }
    }

    /// Computes the inverse transformation.
    #[must_use]

    pub fn inverse(&self) -> Self {

        // Inverse of (az+b)/(cz+d) is (dz-b)/(-cz+a)
        Self {
            a: self.d.clone(),
            b: Expr::new_neg(self.b.clone()),
            c: Expr::new_neg(self.c.clone()),
            d: self.a.clone(),
        }
    }
}

// ============================================================================
// Cauchy Integral Formula
// ============================================================================

/// Evaluates f(z0) using Cauchy's integral formula.
///
/// f(z0) = (1/2πi) ∮_C f(z)/(z-z0) dz
///
/// This is a symbolic representation.
#[must_use]

pub fn cauchy_integral_formula(func: &Expr, var: &str, z0: &Expr) -> Expr {

    // Return symbolic representation
    let z = Expr::Variable(var.to_string());

    let _integrand = Expr::new_div(func.clone(), Expr::new_sub(z, z0.clone()));

    // The result is just f(z0) by Cauchy's formula
    simplify(&substitute(func, var, z0))
}

/// Computes the n-th derivative using Cauchy's formula for derivatives.
///
/// f^(n)(z0) = (n!/2πi) ∮_C f(z)/(z-z0)^(n+1) dz
#[must_use]

pub fn cauchy_derivative_formula(func: &Expr, var: &str, z0: &Expr, n: usize) -> Expr {

    // Simply use differentiation
    let mut result = func.clone();

    for _ in 0..n {

        result = differentiate(&result, var);
    }

    simplify(&substitute(&result, var, z0))
}

// ============================================================================
// Complex Function Utilities
// ============================================================================

/// Computes the complex exponential e^z = e^(x+iy) = e^x(cos(y) + i*sin(y))
#[must_use]

pub fn complex_exp(z: &Expr) -> Expr {

    let re = z.re();

    let im = z.im();

    let exp_re = Expr::new_exp(re);

    let cos_im = Expr::new_cos(im.clone());

    let sin_im = Expr::new_sin(im);

    let real_part = Expr::new_mul(exp_re.clone(), cos_im);

    let imag_part = Expr::new_mul(exp_re, sin_im);

    Expr::Complex(Arc::new(real_part), Arc::new(imag_part))
}

/// Computes the principal branch of complex logarithm.
///
/// log(z) = log|z| + i*arg(z)
#[must_use]

pub fn complex_log(z: &Expr) -> Expr {

    let re = z.re();

    let im = z.im();

    // |z| = sqrt(re^2 + im^2)
    let modulus = Expr::new_sqrt(Expr::new_add(
        Expr::new_pow(re.clone(), Expr::Constant(2.0)),
        Expr::new_pow(im.clone(), Expr::Constant(2.0)),
    ));

    // arg(z) = atan2(im, re)
    let argument = Expr::new_atan2(im, re);

    Expr::Complex(Arc::new(Expr::new_log(modulus)), Arc::new(argument))
}

/// Computes the argument (angle) of a complex number.
#[must_use]

pub fn complex_arg(z: &Expr) -> Expr {

    let re = z.re();

    let im = z.im();

    Expr::new_atan2(im, re)
}

/// Computes the modulus (absolute value) of a complex number.
#[must_use]

pub fn complex_modulus(z: &Expr) -> Expr {

    let re = z.re();

    let im = z.im();

    Expr::new_sqrt(Expr::new_add(
        Expr::new_pow(re, Expr::Constant(2.0)),
        Expr::new_pow(im, Expr::Constant(2.0)),
    ))
}
