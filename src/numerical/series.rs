use crate::numerical::elementary::eval_expr;
use crate::symbolic::core::Expr;
use std::collections::HashMap;

/// Computes the numerical Taylor series expansion of a function around a point.
///
/// This function calculates the coefficients of the Taylor series by numerically
/// evaluating the function and its derivatives at the `at_point`. It then returns
/// a closure that can evaluate the resulting Taylor polynomial at any `x`.
///
/// # Arguments
/// * `f` - The symbolic expression for the function `f(x)`.
/// * `var` - The variable `x`.
/// * `at_point` - The point `a` around which to expand the series.
/// * `order` - The maximum order `N` of the series to compute.
///
/// # Returns
/// A `Result` containing a closure `Box<dyn Fn(f64) -> f64>` that evaluates the Taylor polynomial,
/// or an error string if evaluation fails.
pub fn taylor_series_numerical(
    f: &Expr,
    var: &str,
    at_point: f64,
    order: usize,
) -> Result<Box<dyn Fn(f64) -> f64>, String> {
    let mut coeffs = Vec::with_capacity(order + 1);
    let mut current_f = f.clone();
    let mut factorial = 1.0;

    let mut vars_map = HashMap::new();
    vars_map.insert(var.to_string(), at_point);

    // c_0 = f(a)
    coeffs.push(eval_expr(&current_f, &vars_map)?);

    for i in 1..=order {
        // c_i = f^(i)(a) / i!
        current_f = crate::symbolic::calculus::differentiate(&current_f, var);
        factorial *= i as f64;
        let coeff_val = eval_expr(&current_f, &vars_map)? / factorial;
        coeffs.push(coeff_val);
    }

    let a = at_point;
    let taylor_poly = move |x: f64| -> f64 {
        let mut sum = 0.0;
        for (i, coeff) in coeffs.iter().enumerate() {
            sum += coeff * (x - a).powi(i as i32);
        }
        sum
    };

    Ok(Box::new(taylor_poly))
}
