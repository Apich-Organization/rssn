//! System identification: ARX regression and eigensystem realisation.
//!
//! ARX fits the difference equation
//! `y[k] = -a1 y[k-1] - ... - ana y[k-na] + b1 u[k-delay-1] + ...` by least
//! squares. ERA builds a block Hankel matrix from the impulse response, takes
//! an SVD truncation to the requested order and reads `A`, `B`, `C` from the
//! shift-invariant subspace.

use crate::kernels::dense;

use super::linalg;

/// An identified ARX model.
#[derive(Debug, Clone, PartialEq)]
pub struct ArxModel {
    /// Autoregressive coefficients `a1..ana` (positive sign in the regression
    /// `y[k] + a1 y[k-1] + ... = b1 u[k-delay-1] + ...`).
    pub a: Vec<f64>,
    /// Input coefficients `b1..bnb`.
    pub b: Vec<f64>,
    /// Coefficient of determination `1 - SSE/SST` of the fit.
    pub fit: f64,
    /// One-step prediction residuals on the regression rows.
    pub residuals: Vec<f64>,
}

/// A realised state-space model.
#[derive(Debug, Clone, PartialEq)]
pub struct EraModel {
    /// Realised state matrix.
    pub a: Vec<Vec<f64>>,
    /// Realised input matrix.
    pub b: Vec<Vec<f64>>,
    /// Realised output matrix.
    pub c: Vec<Vec<f64>>,
    /// Realisation order.
    pub order: usize,
}

/// Fits an ARX model to input/output data.
///
/// # Errors
/// Returns an error when the signals differ in length, when both `na` and
/// `nb` are zero, or when there are not enough samples to form a full-rank
/// regression.
pub fn arx(
    u: &[f64],
    y: &[f64],
    na: usize,
    nb: usize,
    delay: usize,
) -> Result<ArxModel, String> {
    if u.len() != y.len() {
        return Err("input and output must have the same length".to_string());
    }
    if na == 0 && nb == 0 {
        return Err("at least one of na and nb must be positive".to_string());
    }
    let start = na.max(delay.saturating_add(nb));
    if start >= y.len() {
        return Err("not enough samples for the requested model order".to_string());
    }
    let rows = y.len() - start;
    let cols = na + nb;
    if rows < cols {
        return Err("not enough samples for the requested model order".to_string());
    }
    let mut regressors = vec![vec![0.0; cols]; rows];
    let mut target = vec![0.0; rows];
    for (row, k) in (start..y.len()).enumerate() {
        for i in 0..na {
            regressors[row][i] = -y[k - 1 - i];
        }
        for j in 0..nb {
            regressors[row][na + j] = u[k - delay - 1 - j];
        }
        target[row] = y[k];
    }
    let theta = linalg::lstsq(&regressors, &target)?;
    let a = theta[..na].to_vec();
    let b = theta[na..].to_vec();
    let mut residuals = Vec::with_capacity(rows);
    let mut sse = 0.0_f64;
    for row in 0..rows {
        let mut prediction = 0.0;
        for (col, &value) in regressors[row].iter().enumerate() {
            prediction += value * theta[col];
        }
        let residual = target[row] - prediction;
        sse += residual * residual;
        residuals.push(residual);
    }
    let mean = target.iter().sum::<f64>() / rows as f64;
    let sst: f64 = target.iter().map(|v| (v - mean) * (v - mean)).sum();
    let fit = if sst > 0.0 {
        1.0 - sse / sst
    } else if sse < 1e-24 {
        1.0
    } else {
        0.0
    };
    Ok(ArxModel { a, b, fit, residuals })
}

/// Estimates the impulse response from an input/output pair by recursive
/// deconvolution, assuming the system starts at rest.
fn deconvolve(
    y: &[Vec<f64>],
    u: &[f64],
) -> Result<Vec<Vec<f64>>, String> {
    if y.is_empty() {
        return Err("output data must not be empty".to_string());
    }
    let outputs = y[0].len();
    if outputs == 0 {
        return Err("each output sample must have at least one entry".to_string());
    }
    if y.iter().any(|row| row.len() != outputs) {
        return Err("output rows must have equal length".to_string());
    }
    if u.len() != y.len() {
        return Err("input and output must have the same length".to_string());
    }
    if u[0] == 0.0 {
        return Err("the first input sample must be nonzero".to_string());
    }
    let samples = y.len();
    let mut h = vec![vec![0.0; outputs]; samples];
    for k in 0..samples {
        let mut value = y[k].clone();
        for j in 0..k {
            let factor = u[k - j];
            for c in 0..outputs {
                value[c] -= h[j][c] * factor;
            }
        }
        for c in 0..outputs {
            h[k][c] = value[c] / u[0];
        }
    }
    Ok(h)
}

/// Realises a state-space model from input/output data by the eigensystem
/// realisation algorithm.
///
/// The impulse response is recovered by deconvolution and arranged in a block
/// Hankel matrix; an SVD truncated to `order` gives the shift-invariant
/// subspace from which `A`, `B` and `C` are read.
///
/// # Errors
/// Returns an error for inconsistent data, a zero first input sample, too few
/// samples for the requested order, or a rank-deficient Hankel matrix.
pub fn era(
    y: &[Vec<f64>],
    u: &[f64],
    order: usize,
) -> Result<EraModel, String> {
    if order == 0 {
        return Err("order must be positive".to_string());
    }
    let h = deconvolve(y, u)?;
    let samples = h.len();
    let outputs = h[0].len();
    let rows_blocks = order;
    if samples < 2 * order + 1 {
        return Err("not enough samples for the requested order".to_string());
    }
    let cols = samples - rows_blocks - 1;
    if cols < order {
        return Err("not enough samples for the requested order".to_string());
    }
    let h0_rows = outputs * rows_blocks;
    let mut h0 = vec![vec![0.0; cols]; h0_rows];
    let mut h1 = vec![vec![0.0; cols]; h0_rows];
    for i in 0..rows_blocks {
        for j in 0..cols {
            for c in 0..outputs {
                h0[i * outputs + c][j] = h[i + j + 1][c];
                h1[i * outputs + c][j] = h[i + j + 2][c];
            }
        }
    }
    let matrix0 = linalg::to_mat(&h0)?;
    let decomposition = dense::svd(&matrix0).map_err(|e| e.to_string())?;
    if decomposition.s.len() < order {
        return Err("the hankel matrix has rank below the requested order".to_string());
    }
    let singular: Vec<f64> = decomposition.s[..order].to_vec();
    let largest = singular.first().copied().unwrap_or(0.0);
    if singular.iter().any(|&v| v <= 1e-12 * largest.max(1.0)) {
        return Err("the hankel matrix has rank below the requested order".to_string());
    }
    let sqrt_s: Vec<f64> = singular.iter().map(|v| v.sqrt()).collect();
    let inv_sqrt: Vec<f64> = sqrt_s.iter().map(|v| 1.0 / v).collect();
    let u_n = &decomposition.u;
    let v_n = &decomposition.v;
    let mut tmp = vec![vec![0.0; cols]; order];
    for a in 0..order {
        for j in 0..cols {
            let mut acc = 0.0;
            for i in 0..h0_rows {
                acc += u_n.at(i, a) * h1[i][j];
            }
            tmp[a][j] = acc;
        }
    }
    let mut projected = vec![vec![0.0; order]; order];
    for a in 0..order {
        for b in 0..order {
            let mut acc = 0.0;
            for j in 0..cols {
                acc += tmp[a][j] * v_n.at(j, b);
            }
            projected[a][b] = acc;
        }
    }
    let mut a = vec![vec![0.0; order]; order];
    for i in 0..order {
        for j in 0..order {
            a[i][j] = inv_sqrt[i] * projected[i][j] * inv_sqrt[j];
        }
    }
    let mut b = vec![vec![0.0; 1]; order];
    for i in 0..order {
        b[i][0] = sqrt_s[i] * v_n.at(0, i);
    }
    let mut c = vec![vec![0.0; order]; outputs];
    for r in 0..outputs {
        for i in 0..order {
            c[r][i] = u_n.at(r, i) * sqrt_s[i];
        }
    }
    Ok(EraModel { a, b, c, order })
}

/// Impulse response of a realised model, `h[0] = 0` and
/// `h[k] = C A^{k-1} B`.
///
/// # Errors
/// Returns an error when the model has no inputs.
pub fn impulse_response(
    model: &EraModel,
    steps: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let outputs = model.c.len();
    let inputs = model.b.first().map_or(0, Vec::len);
    if inputs == 0 {
        return Err("the model must have at least one input".to_string());
    }
    let mut out = vec![vec![0.0; outputs]; steps];
    if steps == 0 {
        return Ok(out);
    }
    let mut excitation = vec![0.0; inputs];
    excitation[0] = 1.0;
    let mut state = linalg::matvec(&model.b, &excitation);
    for k in 1..steps {
        out[k] = linalg::matvec(&model.c, &state);
        state = linalg::matvec(&model.a, &state);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::linalg;
    use super::*;

    #[test]
    fn arx_recovers_known_coefficients() {
        let a_true = [0.6, -0.2];
        let b_true = [0.5];
        let samples = 400;
        let mut u = vec![0.0; samples];
        for k in 0..samples {
            u[k] = (0.3 * k as f64).sin() + 0.5 * (1.1 * k as f64).cos();
        }
        let mut y = vec![0.0; samples];
        for k in 2..samples {
            y[k] = -a_true[0] * y[k - 1] - a_true[1] * y[k - 2] + b_true[0] * u[k - 1];
        }
        let model = arx(&u, &y, 2, 1, 0).unwrap();
        assert!((model.a[0] - a_true[0]).abs() < 1e-6);
        assert!((model.a[1] - a_true[1]).abs() < 1e-6);
        assert!((model.b[0] - b_true[0]).abs() < 1e-6);
        assert!(model.fit > 0.999_999);
    }

    #[test]
    fn era_recovers_poles_of_a_two_state_system() {
        let a = vec![vec![0.9, 0.1], vec![-0.1, 0.8]];
        let b = vec![vec![1.0], vec![0.5]];
        let c = vec![vec![1.0, 0.0]];
        let samples = 40;
        let mut u = vec![0.0; samples];
        u[0] = 1.0;
        let mut x = vec![0.0, 0.0];
        let mut y = vec![vec![0.0]; samples];
        for k in 0..samples {
            let cx = linalg::matvec(&c, &x);
            y[k][0] = cx[0];
            let ax = linalg::matvec(&a, &x);
            let bu = linalg::matvec(&b, &[u[k]]);
            for i in 0..2 {
                x[i] = ax[i] + bu[i];
            }
        }
        let model = era(&y, &u, 2).unwrap();
        let true_poles = linalg::eigenvalues(&a).unwrap();
        let estimated = linalg::eigenvalues(&model.a).unwrap();
        assert_eq!(true_poles.len(), estimated.len());
        for pole in &true_poles {
            let best = estimated
                .iter()
                .map(|e| (e - pole).norm())
                .fold(f64::INFINITY, f64::min);
            assert!(best < 1e-6, "pole {pole:?} not recovered");
        }
    }

    #[test]
    fn impulse_response_has_zero_direct_term() {
        let model = EraModel {
            a: vec![vec![0.5]],
            b: vec![vec![1.0]],
            c: vec![vec![2.0]],
            order: 1,
        };
        let h = impulse_response(&model, 4).unwrap();
        assert!(h[0][0].abs() < 1e-12);
        assert!((h[1][0] - 2.0).abs() < 1e-12);
        assert!((h[2][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn arx_rejects_short_data() {
        let u = vec![1.0, 2.0];
        let y = vec![0.0, 1.0];
        assert!(arx(&u, &y, 3, 3, 0).is_err());
    }
}
