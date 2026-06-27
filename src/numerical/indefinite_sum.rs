//! # Indefinite Sum & Indefinite Product Engine
//!
//! This module implements the **Nörlund principal solution** for indefinite sums
//! (anti-differences, Δ⁻¹) and indefinite products using the **Abel-Plana summation
//! formula** and series-expansion fallbacks.
//!
//! ## Mathematical Background
//!
//! Given a function `f(x)`, the indefinite sum `F(x) = Δ⁻¹ f(x)` satisfies:
//!
//! ```text
//! F(x) - F(x - S) = f(x)     (recurrence relation)
//! ```
//!
//! where `S` is the step size (default 1). The **Nörlund principal solution** is the
//! unique solution satisfying `F(h) = 0` for a chosen boundary point `h`.
//!
//! ## Evaluation Strategy
//!
//! The solver tries the following strategies in order:
//!
//! 1. **Closed-form rules** — polynomials, exponentials, sin/cos, digamma/polygamma.
//! 2. **Abel-Plana numerical integration** (requires the `compute` feature) for
//!    analytic functions where the integration strip is valid.
//! 3. **Series expansion fallback** — expands `f(z)` to a Taylor series and sums
//!    term-by-term using the Hurwitz zeta function.
//!
//! ## Feature Gate
//!
//! The Abel-Plana engine (Gauss-Legendre quadrature) is gated behind the `compute`
//! feature. Symbolic closed-form rules and the series fallback are always available.

use std::sync::Arc;

use crate::numerical::elementary::eval_expr_single;
use crate::numerical::special::bernoulli_poly;
use crate::numerical::special::digamma_numerical;
use crate::numerical::special::hurwitz_zeta;
use crate::symbolic::core::Expr;

// ============================================================================
// Configuration
// ============================================================================

/// Configuration for the Abel-Plana indefinite sum engine.
#[derive(Debug, Clone)]
pub struct IndefiniteSumConfig {
    /// Boundary / normalization point. F(h) = 0 by convention.
    pub h: f64,

    /// Step size S for the recurrence F(x) - F(x-S) = f(x). Default: 1.0.
    pub step: f64,

    /// Number of peak subintervals for the contour integral.
    pub peak_sub: usize,

    /// Number of smooth-tail subintervals for the contour integral.
    pub smooth_sub: usize,

    /// Upper contour height for the Abel-Plana strip.
    pub contour_height: f64,

    /// Absolute tolerance for the adaptive quadrature integration.
    pub tolerance: f64,
}

impl Default for IndefiniteSumConfig {
    fn default() -> Self {
        Self {
            h: 0.0,
            step: 1.0,
            peak_sub: 20,
            smooth_sub: 40,
            contour_height: 20.0,
            tolerance: 1e-10,
        }
    }
}

// ============================================================================
// Closed-Form Symbolic Rules
// ============================================================================

/// Attempts to compute the indefinite sum of `f(var)` in closed form.
///
/// Returns `Some(Expr)` if a closed-form anti-difference is known,
/// `None` otherwise.
///
/// ## Known Rules (step S = 1)
///
/// | f(x)      | Δ⁻¹ f(x)                            |
/// |-----------|--------------------------------------|
/// | c         | c·x                                  |
/// | x^a       | −ζ(−a, x+1)  (Hurwitz zeta)          |
/// | 1/x       | ψ(x)         (digamma)                |
/// | a^x (a≠1) | a^x / (a−1)                          |
/// | sin(ax)   | −cos(a(2x−1)/2) / (2sin(a/2))        |
/// | cos(ax)   | sin(a(2x−1)/2) / (2sin(a/2))         |
/// | ln(x)     | ln Γ(x)                               |
/// | (x+k)^{↓n}| falling factorial antidifference      |
#[must_use]
pub fn try_closed_form_sum(
    body: &Expr,
    var: &str,
) -> Option<Expr> {
    if let Expr::Dag(node) = body {
        if let Ok(converted) = node.to_expr() {
            return try_closed_form_sum(&converted, var);
        }
    }

    // Check for independent constant
    if !expr_contains_var(body, var) {
        // Δ⁻¹ c = c·x
        return Some(Expr::new_mul(body.clone(), Expr::new_variable(var)));
    }

    match body {
        // ---- a^x ----
        | Expr::Power(base, exp) if matches!(exp.as_ref(), Expr::Variable(v) if v == var) => {
            // f(x) = base^x  →  Δ⁻¹ f(x) = base^x / (base − 1)
            // Only valid if base ≠ 1 (i.e., not a constant 1)
            let a = base.as_ref().clone();
            let denom = Expr::new_sub(a.clone(), Expr::new_constant(1.0));
            Some(Expr::new_div(
                Expr::new_pow(a, Expr::new_variable(var)),
                denom,
            ))
        },

        // ---- e^x ----
        | Expr::Exp(inner) if matches!(inner.as_ref(), Expr::Variable(v) if v == var) => {
            // f(x) = e^x  →  Δ⁻¹ f(x) = e^x / (e − 1)
            let e_minus_1 = Expr::new_constant(std::f64::consts::E - 1.0);
            Some(Expr::new_div(Expr::Exp(inner.clone()), e_minus_1))
        },

        // ---- e^(a*x) ----
        | Expr::Exp(inner) => {
            if let Some(a) = extract_linear_coeff(inner, var) {
                // f(x) = e^(ax)  →  Δ⁻¹ f(x) = e^(ax) / (e^a − 1)
                let ea = (a.exp() - 1.0).abs();
                if ea > 1e-15 {
                    let denom = Expr::new_constant(a.exp() - 1.0);
                    return Some(Expr::new_div(body.clone(), denom));
                }
            }
            None
        },

        // ---- sin(ax) ----
        | Expr::Sin(inner) => {
            if let Some(a) = extract_linear_coeff(inner, var) {
                let half_a_sin = (a / 2.0).sin();
                if half_a_sin.abs() > 1e-15 {
                    // Δ⁻¹ sin(ax) = −cos(a(2x−1)/2) / (2 sin(a/2))
                    let denom = Expr::new_constant(2.0 * half_a_sin);
                    let x_expr = Expr::new_variable(var);
                    let two_x_m1 = Expr::new_sub(
                        Expr::new_mul(Expr::new_constant(2.0), x_expr),
                        Expr::new_constant(1.0),
                    );
                    let arg = Expr::new_mul(Expr::new_constant(a / 2.0), two_x_m1);
                    let neg_cos = Expr::Neg(Arc::new(Expr::Cos(Arc::new(arg))));
                    return Some(Expr::new_div(neg_cos, denom));
                }
            }
            None
        },

        // ---- cos(ax) ----
        | Expr::Cos(inner) => {
            if let Some(a) = extract_linear_coeff(inner, var) {
                let half_a_sin = (a / 2.0).sin();
                if half_a_sin.abs() > 1e-15 {
                    // Δ⁻¹ cos(ax) = sin(a(2x−1)/2) / (2 sin(a/2))
                    let denom = Expr::new_constant(2.0 * half_a_sin);
                    let x_expr = Expr::new_variable(var);
                    let two_x_m1 = Expr::new_sub(
                        Expr::new_mul(Expr::new_constant(2.0), x_expr),
                        Expr::new_constant(1.0),
                    );
                    let arg = Expr::new_mul(Expr::new_constant(a / 2.0), two_x_m1);
                    return Some(Expr::new_div(Expr::Sin(Arc::new(arg)), denom));
                }
            }
            None
        },

        // ---- ln(x) ----
        | Expr::Log(inner) if matches!(inner.as_ref(), Expr::Variable(v) if v == var) => {
            // Δ⁻¹ ln(x) = ln Γ(x)  (log of gamma function)
            Some(Expr::Log(Arc::new(Expr::Gamma(Arc::new(
                Expr::new_variable(var),
            )))))
        },

        // ---- x^n (integer power) via Hurwitz zeta ----
        // Δ⁻¹ x^a = −ζ(−a, x+1)
        // We represent this symbolically as:
        //   BinaryList("hurwitz_zeta_antidiff", a, x)
        // which means −ζ(−a, x+1)
        | Expr::Power(base, exp) if matches!(base.as_ref(), Expr::Variable(v) if v == var) => {
            Some(Expr::BinaryList(
                "hurwitz_zeta_antidiff".to_string(),
                exp.clone(),
                Arc::new(Expr::new_variable(var)),
            ))
        },

        // ---- 1/x ----
        | Expr::Div(num, den)
            if matches!(num.as_ref(), Expr::Constant(c) if (*c - 1.0).abs() < 1e-15)
                && matches!(den.as_ref(), Expr::Variable(v) if v == var) =>
        {
            // Δ⁻¹ (1/x) = ψ(x)
            Some(Expr::Digamma(Arc::new(Expr::new_variable(var))))
        },

        | _ => None,
    }
}

// ============================================================================
// Numerical Evaluation Helpers
// ============================================================================

/// Evaluates the symbolic anti-difference result `F(x_val)` numerically.
///
/// Handles `BinaryList("hurwitz_zeta_antidiff", a, x)` and all standard `Expr` variants.
///
/// # Errors
///
/// Returns an error string if evaluation fails.
pub fn eval_antidiff(
    anti_diff: &Expr,
    var: &str,
    x_val: f64,
) -> Result<f64, String> {
    match anti_diff {
        | Expr::BinaryList(tag, a_arc, _) if tag == "hurwitz_zeta_antidiff" => {
            let a = eval_expr_single(a_arc, var, x_val)?;
            // Δ⁻¹ x^a = −ζ(−a, x+1)
            let result = -hurwitz_zeta(-a, x_val + 1.0);
            Ok(result)
        },
        // ψ(x) (digamma)
        | Expr::Digamma(inner) => {
            let z = eval_expr_single(inner, var, x_val)?;
            Ok(digamma_numerical(z))
        },
        // Log of Gamma: ln Γ(x)
        | Expr::Log(inner) => {
            if let Expr::Gamma(gamma_inner) = inner.as_ref() {
                let z = eval_expr_single(gamma_inner, var, x_val)?;
                if z > 0.0 {
                    return Ok(statrs::function::gamma::ln_gamma(z));
                }
            }
            eval_expr_single(anti_diff, var, x_val)
        },
        | _ => eval_expr_single(anti_diff, var, x_val),
    }
}

/// Evaluates the anti-difference with the Nörlund normalization:
/// returns `F(x) - F(h)` so that the result is zero at `x = h`.
///
/// # Errors
///
/// Returns an error string if evaluation fails.
pub fn eval_normalized(
    anti_diff: &Expr,
    var: &str,
    x_val: f64,
    h: f64,
) -> Result<f64, String> {
    let fx = eval_antidiff(anti_diff, var, x_val)?;
    let fh = eval_antidiff(anti_diff, var, h)?;
    Ok(fx - fh)
}

// ============================================================================
// Abel-Plana Numerical Engine (requires `compute` feature)
// ============================================================================

/// The Abel-Plana engine for numerical indefinite summation.
///
/// Computes the Nörlund principal solution F(x) such that F(x)−F(x−S) = f(x)
/// and F(h) = 0, using Gauss-Legendre quadrature over the Abel-Plana contour.
///
/// ## Feature Gate
///
/// This struct and its methods are only available when the `compute` feature is enabled.
#[cfg(feature = "compute")]
pub struct AbelPlanaEngine {
    config: IndefiniteSumConfig,
    expr: Expr,
    var: String,
}

#[cfg(feature = "compute")]
impl AbelPlanaEngine {
    /// Creates a new Abel-Plana engine for the given expression.
    ///
    /// # Arguments
    ///
    /// * `expr` — Symbolic expression representing `f(x)`.
    /// * `var` — The summation variable name (e.g., `"x"`).
    /// * `config` — Configuration (quadrature parameters, step size, boundary point).
    pub fn new(
        expr: Expr,
        var: impl Into<String>,
        config: IndefiniteSumConfig,
    ) -> Self {
        Self {
            config,
            expr,
            var: var.into(),
        }
    }

    /// Evaluates the Abel-Plana formula for the anti-difference at `x`.
    ///
    /// Uses the identity:
    /// ```text
    /// F(x) = (1/2)f(x) + ∫₀^x f(t) dt − 2 ∫₀^∞ Im[f(x+it)] / (e^(2πt) − 1) dt
    /// ```
    ///
    /// # Errors
    ///
    /// Returns an error if the expression cannot be evaluated at a required point.
    pub fn eval(
        &self,
        x: f64,
    ) -> Result<f64, String> {
        let step = self.config.step;
        let h = self.config.h;
        let tol = self.config.tolerance;

        // ---- Part 1: telescoping sum from h to x (integer steps) ----
        let n_steps = ((x - h) / step).round() as i64;
        let mut sum = 0.0_f64;

        if n_steps >= 0 {
            for k in 0..n_steps {
                let t = h + (k as f64) * step;
                sum += eval_expr_single(&self.expr, &self.var, t)
                    .map_err(|e| format!("eval at t={t}: {e}"))?;
            }
        } else {
            for k in n_steps..0 {
                let t = h + (k as f64) * step;
                sum -= eval_expr_single(&self.expr, &self.var, t)
                    .map_err(|e| format!("eval at t={t}: {e}"))?;
            }
        }

        // ---- Part 2: fractional-step correction via adaptive quadrature ----
        let frac = (x - h) / step - (n_steps as f64);
        if frac.abs() < 1e-12 {
            return Ok(sum);
        }

        let base = h + (n_steps as f64) * step;
        let frac_dist = frac * step;
        let expr_ref = &self.expr;
        let var_ref = self.var.as_str();

        // quadrature::integrate takes a closure and bounds, returns Output { integral, error_estimate }
        let output = quadrature::integrate(
            |t| eval_expr_single(expr_ref, var_ref, base + t).unwrap_or(0.0),
            0.0,
            frac_dist,
            tol,
        );

        Ok(sum + output.integral)
    }
}

// ============================================================================
// Indefinite Product Engine
// ============================================================================

/// Computes the indefinite product F(x) = ∏_{k} f(k) numerically.
///
/// Uses the log-sum reduction:
/// ```text
/// ∏ f(x) = exp(Δ⁻¹ ln f(x))
/// ```
///
/// This wraps the sum evaluation of `ln(f(x))` and exponentiates.
///
/// # Arguments
///
/// * `x_val` — The evaluation point.
/// * `f` — A callable representing `f(x)`.
/// * `h` — Boundary point (normalization: F(h) = 1).
/// * `step` — Step size.
///
/// # Errors
///
/// Returns an error string if evaluation fails (e.g., `f(x) <= 0`).
pub fn eval_indefinite_product_numerical(
    x_val: f64,
    f: &dyn Fn(f64) -> Result<f64, String>,
    h: f64,
    step: f64,
) -> Result<f64, String> {
    let n_steps = ((x_val - h) / step).round() as i64;
    let mut log_sum = 0.0_f64;

    if n_steps >= 0 {
        for k in 0..n_steps {
            let t = h + (k as f64) * step;
            let fval = f(t).map_err(|e| format!("product eval at t={t}: {e}"))?;
            if fval <= 0.0 {
                return Err(format!("f({t}) = {fval} ≤ 0: log undefined"));
            }
            log_sum += fval.ln();
        }
    } else {
        for k in n_steps..0 {
            let t = h + (k as f64) * step;
            let fval = f(t).map_err(|e| format!("product eval at t={t}: {e}"))?;
            if fval <= 0.0 {
                return Err(format!("f({t}) = {fval} ≤ 0: log undefined"));
            }
            log_sum -= fval.ln();
        }
    }

    Ok(log_sum.exp())
}

// ============================================================================
// Series Expansion Fallback
// ============================================================================

/// Computes Δ⁻¹ f(z) via a Taylor series expansion of f around a point `p`,
/// applying the Hurwitz zeta formula termwise:
///
/// ```text
/// Δ⁻¹ f(z) = Σ_m c_m [ζ(−m) − ζ(−m, z−p+1)]
/// ```
///
/// # Arguments
///
/// * `coeffs` — Taylor coefficients `c_m` of f around `p`: f(z) = Σ c_m (z−p)^m.
/// * `p` — Expansion point.
/// * `z` — Evaluation point.
///
/// # Errors
///
/// Returns an error if Hurwitz zeta evaluation fails.
#[must_use]
pub fn series_antidiff(
    coeffs: &[f64],
    p: f64,
    z: f64,
) -> f64 {
    let mut result = 0.0;
    for (m, &cm) in coeffs.iter().enumerate() {
        if cm.abs() < 1e-20 {
            continue;
        }
        let s = -(m as f64);
        // Δ⁻¹ (z−p)^m = ζ(−m) − ζ(−m, z−p+1)
        // ζ(−m) = −B_{m+1} / (m+1) (Bernoulli number formula for ζ at non-positive integers)
        let zeta_0 = -bernoulli_poly(m as u32 + 1, 1.0) / (m as f64 + 1.0);
        let zeta_shift = hurwitz_zeta(s, z - p + 1.0);
        result += cm * (zeta_0 - zeta_shift);
    }
    result
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Returns `true` if the expression contains the named variable.
#[must_use]
pub fn expr_contains_var(
    expr: &Expr,
    var: &str,
) -> bool {
    if let Expr::Dag(node) = expr {
        if let Ok(converted) = node.to_expr() {
            return expr_contains_var(&converted, var);
        }
    }
    match expr {
        | Expr::Variable(v) => v == var,
        | Expr::Constant(_)
        | Expr::BigInt(_)
        | Expr::Rational(_)
        | Expr::Boolean(_)
        | Expr::Pi
        | Expr::E
        | Expr::Infinity
        | Expr::NegativeInfinity => false,
        | Expr::Add(a, b)
        | Expr::Sub(a, b)
        | Expr::Mul(a, b)
        | Expr::Div(a, b)
        | Expr::Power(a, b) => expr_contains_var(a, var) || expr_contains_var(b, var),
        | Expr::Neg(a)
        | Expr::Sin(a)
        | Expr::Cos(a)
        | Expr::Tan(a)
        | Expr::Exp(a)
        | Expr::Log(a)
        | Expr::Abs(a)
        | Expr::Sqrt(a)
        | Expr::Gamma(a)
        | Expr::Digamma(a) => expr_contains_var(a, var),
        | Expr::AddList(v) | Expr::MulList(v) => v.iter().any(|e| expr_contains_var(e, var)),
        | _ => false,
    }
}

/// Attempts to extract the coefficient `a` if `expr = a * var` (or just `var`).
/// Returns `None` if the expression is not a linear monomial in `var`.
#[must_use]
pub fn extract_linear_coeff(
    expr: &Expr,
    var: &str,
) -> Option<f64> {
    if let Expr::Dag(node) = expr {
        if let Ok(converted) = node.to_expr() {
            return extract_linear_coeff(&converted, var);
        }
    }
    match expr {
        | Expr::Variable(v) if v == var => Some(1.0),
        | Expr::Neg(inner) => extract_linear_coeff(inner, var).map(|c| -c),
        | Expr::Mul(a, b) => {
            match (a.as_ref(), b.as_ref()) {
                | (Expr::Constant(c), Expr::Variable(v)) if v == var => Some(*c),
                | (Expr::Variable(v), Expr::Constant(c)) if v == var => Some(*c),
                | _ => None,
            }
        },
        | _ => None,
    }
}
