//! Prediction from a fitted generalised linear model.

use crate::kernels::glm::fit::GlmFit;
use crate::kernels::glm::link::Link;
use crate::kernels::glm::link::link_inverse;

/// Predicts the mean response for new covariate rows.
///
/// `x` holds one row of covariates per observation, without the intercept
/// column, exactly as accepted by [`fit`](crate::kernels::glm::fit::fit). The
/// returned values are on the response scale, `g⁻¹(x'β)`.
///
/// # Errors
/// Returns a lowercase message when the model has no coefficients or a row
/// does not have one entry per fitted slope. A non-finite linear predictor
/// or response is also an error.
pub fn predict(
    fit: &GlmFit,
    x: &[Vec<f64>],
    kind: Link,
) -> Result<Vec<f64>, String> {
    predict_coefficients(&fit.beta, x, kind)
}

/// Predicts response means from an intercept followed by fitted coefficients.
///
/// # Errors
/// Fails on an empty coefficient vector, a mismatched row, or a non-finite mean.
pub fn predict_coefficients(
    coefficients: &[f64],
    x: &[Vec<f64>],
    kind: Link,
) -> Result<Vec<f64>, String> {
    let Some((&intercept, slopes)) = coefficients.split_first() else {
        return Err("fit has no coefficients".to_string());
    };
    if coefficients
        .iter()
        .any(|coefficient| !coefficient.is_finite())
    {
        return Err("fit coefficients must be finite".into());
    }
    let mut out = Vec::with_capacity(x.len());
    for (i, row) in x.iter().enumerate() {
        if row.len() != slopes.len() {
            return Err(format!(
                "row {i} has {} entries but the model expects {}",
                row.len(),
                slopes.len(),
            ));
        }
        let mut linear = intercept;
        for (value, coefficient) in row.iter().zip(slopes) {
            linear += value * coefficient;
        }
        if !linear.is_finite() {
            return Err(format!("row {i} produced a non-finite linear predictor"));
        }
        let mean = link_inverse(kind, linear);
        if !mean.is_finite() {
            return Err(format!("row {i} produced a non-finite mean"));
        }
        out.push(mean);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernels::glm::family::Family;
    use crate::kernels::glm::fit::fit;

    #[test]
    fn gaussian_prediction_follows_the_fitted_line() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0]];
        let y = vec![1.0, 3.0, 5.0, 7.0];
        let model = fit(&x, &y, Family::Gaussian, Link::Identity, None, 50, 1e-12).expect("fit");

        let predicted =
            predict(&model, &[vec![4.0], vec![5.0]], Link::Identity).expect("prediction");
        assert!((predicted[0] - 9.0).abs() < 1e-6);
        assert!((predicted[1] - 11.0).abs() < 1e-6);
    }

    #[test]
    fn poisson_prediction_is_on_the_response_scale() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0]];
        let y = vec![2.0, 3.0, 5.0, 8.0];
        let model = fit(&x, &y, Family::Poisson, Link::Log, None, 200, 1e-12).expect("fit");

        let predicted = predict(&model, &[vec![4.0]], Link::Log).expect("prediction");
        let eta = model.beta[0] + 4.0 * model.beta[1];
        assert!((predicted[0] - eta.exp()).abs() < 1e-10);
    }

    #[test]
    fn wrong_row_width_is_rejected() {
        let x = vec![vec![0.0], vec![1.0]];
        let y = vec![1.0, 2.0];
        let model = fit(&x, &y, Family::Gaussian, Link::Identity, None, 50, 1e-12).expect("fit");
        assert!(predict(&model, &[vec![1.0, 2.0]], Link::Identity).is_err());
    }
}

#[cfg(test)]
mod finite_tests {
    use super::*;

    #[test]
    fn saturated_inverse_links_do_not_hide_invalid_linear_arguments() {
        let link = Link::Logit;
        for value in [f64::NAN, f64::INFINITY, f64::MAX] {
            let rows = [vec![value]];
            assert!(predict_coefficients(&[0.0, 2.0], &rows, link).is_err());
        }
        assert!(predict_coefficients(&[f64::NAN], &[], link).is_err());
    }
}
