use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Computes the numerical Taylor series coefficients of a function around a point.
///
/// # Example
/// ```rust
/// use rssn::symbolic::core::Expr;
/// use rssn::numerical::series::taylor_coefficients;
///
/// let x = Expr::new_variable("x");
/// let f = Expr::new_pow(x, Expr::new_constant(2.0)); // f(x) = x^2
/// let coeffs = taylor_coefficients(&f, "x", 0.0, 2).unwrap();
/// // coeffs = [f(0), f'(0), f''(0)/2] = [0, 0, 1]
/// assert!((coeffs[2] - 1.0).abs() < 1e-5);
/// ```

pub fn taylor_coefficients(
    f: &Expr,
    var: &str,
    at_point: f64,
    order: usize,
) -> Result<Vec<f64>, String> {

    let mut coeffs = Vec::with_capacity(order + 1);

    let mut current_f = f.clone();

    let mut factorial = 1.0;

    let mut vars_map = HashMap::new();

    vars_map.insert(var.to_string(), at_point);

    coeffs.push(eval_expr(&current_f, &vars_map)?);

    for i in 1..=order {

        current_f = crate::symbolic::calculus::differentiate(&current_f, var);

        factorial *= i as f64;

        let val = eval_expr(&current_f, &vars_map)? / factorial;

        coeffs.push(val);
    }

    Ok(coeffs)
}

/// Evaluates a power series at a point given its coefficients and center.
///
/// # Example
/// ```rust
/// use rssn::numerical::series::evaluate_power_series;
///
/// let coeffs = vec![1.0, 1.0, 0.5]; // 1 + x + x^2/2 (approx e^x)
/// let val = evaluate_power_series(&coeffs, 0.0, 1.0); // e^1 approx 2.5
/// assert!((val - 2.5).abs() < 1e-5);
/// ```

pub fn evaluate_power_series(coeffs: &[f64], at_point: f64, x: f64) -> f64 {

    let dx = x - at_point;

    let mut sum = 0.0;

    let mut p = 1.0;

    for &c in coeffs {

        sum += c * p;

        p *= dx;
    }

    sum
}

/// Computes the sum of a sequence of expressions for a range of indices.
///
/// # Example
/// ```rust
/// use rssn::symbolic::core::Expr;
/// use rssn::numerical::series::sum_series;
///
/// let n = Expr::new_variable("n");
/// let f = n; // sum n from 1 to 10
/// let sum = sum_series(&f, "n", 1, 10).unwrap();
/// assert!((sum - 55.0).abs() < 1e-5);
/// ```

pub fn sum_series(f: &Expr, var: &str, start: i64, end: i64) -> Result<f64, String> {

    let mut sum = 0.0;

    let mut vars_map = HashMap::new();

    for i in start..=end {

        vars_map.insert(var.to_string(), i as f64);

        sum += eval_expr(f, &vars_map)?;
    }

    Ok(sum)
}
