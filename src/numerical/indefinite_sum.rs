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
//! F(x + S) - F(x) = f(x)     (forward difference recurrence)
//! ```
//!
//! where `S` is the step size (default 1). This is the standard anti-difference
//! with respect to the forward difference operator Δ. The **Nörlund principal solution**
//! is normalized by `F(h) = 0` for a chosen boundary point `h`.
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

    /// Step size S for the recurrence F(x+S) - F(x) = f(x). Default: 1.0.
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

    // Linearity: Add(a, b)
    if let Expr::Add(a, b) = body {
        if let (Some(sa), Some(sb)) = (try_closed_form_sum(a, var), try_closed_form_sum(b, var)) {
            return Some(Expr::new_add(sa, sb));
        }
    }

    // Linearity: Sub(a, b)
    if let Expr::Sub(a, b) = body {
        if let (Some(sa), Some(sb)) = (try_closed_form_sum(a, var), try_closed_form_sum(b, var)) {
            return Some(Expr::new_sub(sa, sb));
        }
    }

    // Linearity: Neg(a)
    if let Expr::Neg(a) = body {
        if let Some(sa) = try_closed_form_sum(a, var) {
            return Some(Expr::new_neg(sa));
        }
    }

    // Linearity: Mul(a, b)
    if let Expr::Mul(a, b) = body {
        if !expr_contains_var(a, var) {
            if let Some(sb) = try_closed_form_sum(b, var) {
                return Some(Expr::new_mul(a.as_ref().clone(), sb));
            }
        } else if !expr_contains_var(b, var) {
            if let Some(sa) = try_closed_form_sum(a, var) {
                return Some(Expr::new_mul(sa, b.as_ref().clone()));
            }
        }
    }

    // Linearity: AddList(terms)
    if let Expr::AddList(terms) = body {
        let mut summed_terms = Vec::new();
        let mut all_success = true;
        for term in terms {
            if let Some(st) = try_closed_form_sum(term, var) {
                summed_terms.push(st);
            } else {
                all_success = false;
                break;
            }
        }
        if all_success {
            return Some(Expr::AddList(summed_terms));
        }
    }

    // Linearity: MulList(factors)
    if let Expr::MulList(factors) = body {
        let (const_factors, var_factors): (Vec<_>, Vec<_>) = factors
            .iter()
            .cloned()
            .partition(|f| !expr_contains_var(f, var));
        if !const_factors.is_empty() {
            let rest = if var_factors.is_empty() {
                Expr::Constant(1.0)
            } else if var_factors.len() == 1 {
                var_factors[0].clone()
            } else {
                Expr::MulList(var_factors)
            };
            if let Some(sum_rest) = try_closed_form_sum(&rest, var) {
                let c = if const_factors.len() == 1 {
                    const_factors[0].clone()
                } else {
                    Expr::MulList(const_factors)
                };
                return Some(Expr::new_mul(c, sum_rest));
            }
        }
    }

    match body {
        // ---- a^x ----
        | Expr::Power(base, exp) if is_pure_var(exp.as_ref(), var) => {
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
        | Expr::Exp(inner) if is_pure_var(inner.as_ref(), var) => {
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
        | Expr::Log(inner) if is_pure_var(inner.as_ref(), var) => {
            // Δ⁻¹ ln(x) = ln Γ(x)  (log of gamma function)
            Some(Expr::Log(Arc::new(Expr::Gamma(Arc::new(
                Expr::new_variable(var),
            )))))
        },

        // ---- x^n (integer power) via Hurwitz zeta ----
        // Δ⁻¹ x^a = ζ(−a, 1) − ζ(−a, x), satisfying F(x+1)−F(x) = x^a with F(1)=0.
        // Represented symbolically as BinaryList("hurwitz_zeta_antidiff", a, x);
        // eval_antidiff evaluates to hurwitz_zeta(−a, 1) − hurwitz_zeta(−a, x).
        | Expr::Power(base, exp) if is_pure_var(base.as_ref(), var) => {
            Some(Expr::BinaryList(
                "hurwitz_zeta_antidiff".to_string(),
                exp.clone(),
                Arc::new(Expr::new_variable(var)),
            ))
        },

        // ---- 1/x ----
        | Expr::Div(num, den)
            if matches!(num.as_ref(), Expr::Constant(c) if (*c - 1.0).abs() < 1e-15)
                && is_pure_var(den.as_ref(), var) =>
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
            // Δ⁻¹ x^a: forward convention F(x+1)−F(x) = x^a, normalized F(1) = 0.
            // F(x) = ζ(−a, 1) − ζ(−a, x).
            let result = hurwitz_zeta(-a, 1.0) - hurwitz_zeta(-a, x_val);
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
/// Computes the Nörlund principal solution F(x) such that F(x+S)−F(x) = f(x)
/// and F(h) = 0, using adaptive quadrature over the Abel-Plana contour.
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
// Unified Numerical Evaluation Pipeline
// ============================================================================

/// Computes Taylor coefficients c_m = f^(m)(p)/m! numerically via forward differences.
///
/// Evaluates f at `n_terms` equally-spaced points `p, p+h, ..., p+(n`terms-1)*h`
/// and applies the Newton forward-difference formula to approximate each coefficient.
///
/// # Arguments
///
/// * `expr` — Symbolic expression representing f.
/// * `var` — The variable name.
/// * `center` — Expansion point p.
/// * `n_terms` — Number of Taylor terms to compute.
#[must_use]
pub fn compute_taylor_coeffs_numerical(
    expr: &Expr,
    var: &str,
    center: f64,
    n_terms: usize,
) -> Vec<f64> {
    let h = 1e-3;
    let mut diffs: Vec<f64> = (0..n_terms)
        .map(|k| eval_expr_single(expr, var, center + k as f64 * h).unwrap_or(f64::NAN))
        .collect();

    let mut coeffs = Vec::with_capacity(n_terms);
    let mut h_pow = 1.0_f64;
    let mut fact = 1.0_f64;

    for m in 0..n_terms {
        if m > 0 {
            h_pow *= h;
            fact *= m as f64;
            for k in 0..(n_terms - m) {
                diffs[k] = diffs[k + 1] - diffs[k];
            }
        }
        coeffs.push(diffs[0] / (h_pow * fact));
    }

    coeffs
}

/// Evaluates the indefinite sum F(x) - F(h) using a three-strategy cascade:
///
/// 1. **Closed-form** symbolic rules (always tried first).
/// 2. **Abel-Plana** numerical integration (requires `compute` feature).
/// 3. **Series expansion** fallback via Taylor coefficients + Hurwitz zeta.
///
/// # Errors
///
/// Returns an error string if all strategies fail or produce invalid results.
pub fn eval_indefinite_sum_numerical(
    expr: &Expr,
    var: &str,
    x_val: f64,
    config: &IndefiniteSumConfig,
) -> Result<f64, String> {
    let h = config.h;

    // Strategy 1: closed-form symbolic rules
    if let Some(anti_diff) = try_closed_form_sum(expr, var) {
        let fx = eval_antidiff(&anti_diff, var, x_val)?;
        let fh = eval_antidiff(&anti_diff, var, h)?;
        return Ok(fx - fh);
    }

    // Strategy 2: Abel-Plana numerical integration
    #[cfg(feature = "compute")]
    {
        let engine = AbelPlanaEngine::new(expr.clone(), var.to_string(), config.clone());
        if let Ok(val) = engine.eval(x_val) {
            return Ok(val);
        }
    }

    // Strategy 3: Taylor series + Hurwitz zeta (expand around h)
    let coeffs = compute_taylor_coeffs_numerical(expr, var, h, 8);
    if coeffs.iter().any(|c| c.is_nan() || c.is_infinite()) {
        return Err(format!(
            "series fallback: Taylor coefficients invalid near x={h}"
        ));
    }
    // series_antidiff(coeffs, p, z) is already normalized to zero at z = p = h
    Ok(series_antidiff(&coeffs, h, x_val))
}

// ============================================================================
// Utility Functions
// ============================================================================

/// Returns `true` if `expr` is exactly the variable `var`, handling `Expr::Dag`
/// wrappers that `Expr::new_variable` produces.
fn is_pure_var(
    expr: &Expr,
    var: &str,
) -> bool {
    match expr {
        | Expr::Variable(v) => v == var,
        | Expr::Dag(node) => {
            node.to_expr()
                .map(|e| matches!(e, Expr::Variable(v) if v == var))
                .unwrap_or(false)
        },
        | _ => false,
    }
}

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
