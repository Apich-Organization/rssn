//! Linear Gaussian state-space filtering and Rauch–Tung–Striebel smoothing.

use crate::kernels::dense::Mat;
use crate::kernels::dense::lu_factor;

/// The output of a Kalman filter pass.
#[derive(Debug, Clone, PartialEq)]
pub struct KalmanResult {
    /// Filtered state mean at each observation, `E[x_t | y_1..y_t]`.
    pub filtered: Vec<Vec<f64>>,
    /// Predicted state mean at each observation, `E[x_t | y_1..y_{t-1}]`.
    pub predicted: Vec<Vec<f64>>,
    /// Filtered state covariance at each observation.
    pub covariance: Vec<Vec<Vec<f64>>>,
    /// Gaussian log-likelihood of the observations under the model.
    pub loglik: f64,
}

/// A filter pass with the covariances kept as matrices for the smoother.
struct FilterRun {
    filtered: Vec<Vec<f64>>,
    predicted: Vec<Vec<f64>>,
    filtered_cov: Vec<Mat>,
    predicted_cov: Vec<Mat>,
    loglik: f64,
}

/// Builds a matrix from rows, mapping the shape error to a message.
fn from_rows(
    rows: &[Vec<f64>],
    name: &str,
) -> Result<Mat, String> {
    Mat::from_rows(rows).map_err(|_| format!("{name} rows have unequal lengths"))
}

/// Converts a matrix to row vectors.
fn to_rows(matrix: &Mat) -> Vec<Vec<f64>> {
    let mut rows = Vec::with_capacity(matrix.rows);
    for i in 0..matrix.rows {
        let mut row = Vec::with_capacity(matrix.cols);
        for j in 0..matrix.cols {
            row.push(matrix.at(i, j));
        }
        rows.push(row);
    }
    rows
}

/// Entry-wise sum of two equally shaped matrices.
fn add(
    a: &Mat,
    b: &Mat,
) -> Result<Mat, String> {
    if a.rows != b.rows || a.cols != b.cols {
        return Err("matrix shapes differ in an addition".to_string());
    }
    let mut out = Mat::zeros(a.rows, a.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            out.set(i, j, a.at(i, j) + b.at(i, j));
        }
    }
    Ok(out)
}

/// Entry-wise difference of two equally shaped matrices.
fn sub(
    a: &Mat,
    b: &Mat,
) -> Result<Mat, String> {
    if a.rows != b.rows || a.cols != b.cols {
        return Err("matrix shapes differ in a subtraction".to_string());
    }
    let mut out = Mat::zeros(a.rows, a.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            out.set(i, j, a.at(i, j) - b.at(i, j));
        }
    }
    Ok(out)
}

/// Symmetrises a square matrix, damping the drift of the update formulas.
fn symmetrise(a: &Mat) -> Mat {
    let n = a.rows;
    let mut out = a.clone();
    for i in 0..n {
        for j in (i + 1)..n {
            let mean = 0.5 * a.at(i, j) + 0.5 * a.at(j, i);
            out.set(i, j, mean);
            out.set(j, i, mean);
        }
    }
    out
}

/// Runs the prediction and update steps, keeping full covariance matrices.
fn run_filter(
    y: &[Vec<f64>],
    transition: &Mat,
    observation: &Mat,
    process_cov: &Mat,
    noise_cov: &Mat,
    x0: &[f64],
    p0: &Mat,
) -> Result<FilterRun, String> {
    let state = x0.len();
    if state == 0 {
        return Err("initial state is empty".to_string());
    }
    if transition.rows != state || transition.cols != state {
        return Err("transition matrix must be square with the state dimension".to_string());
    }
    if process_cov.rows != state || process_cov.cols != state {
        return Err("process covariance must match the state dimension".to_string());
    }
    if p0.rows != state || p0.cols != state {
        return Err("initial covariance must match the state dimension".to_string());
    }
    if observation.cols != state {
        return Err("observation matrix must have one column per state".to_string());
    }
    let measured = observation.rows;
    if measured == 0 {
        return Err("observation matrix has no rows".to_string());
    }
    if noise_cov.rows != measured || noise_cov.cols != measured {
        return Err("noise covariance must be square with the observation dimension".to_string());
    }
    if y.is_empty() {
        return Err("observation sequence is empty".to_string());
    }
    for (t, row) in y.iter().enumerate() {
        if row.len() != measured {
            return Err(format!(
                "observation {t} does not match the observation dimension"
            ));
        }
        if row.iter().any(|v| !v.is_finite()) {
            return Err(format!("observation {t} contains a non-finite value"));
        }
    }

    let identity = Mat::identity(state);
    let h_transpose = observation.transpose();
    let f_transpose = transition.transpose();

    let mut predicted_state = transition.matvec(x0);
    let mut predicted_cov = add(
        &transition
            .matmul(p0)
            .and_then(|m| m.matmul(&f_transpose))
            .map_err(|err| format!("state prediction failed: {err}"))?,
        process_cov,
    )?;

    let mut filtered = Vec::with_capacity(y.len());
    let mut predicted = Vec::with_capacity(y.len());
    let mut filtered_cov = Vec::with_capacity(y.len());
    let mut predicted_covariances = Vec::with_capacity(y.len());
    let mut loglik = 0.0;

    for row in y {
        predicted.push(predicted_state.clone());
        predicted_covariances.push(predicted_cov.clone());

        let projected = observation.matvec(&predicted_state);
        let innovation: Vec<f64> = row
            .iter()
            .zip(projected.iter())
            .map(|(observed, expected)| observed - expected)
            .collect();

        let innovation_cov = add(
            &observation
                .matmul(&predicted_cov)
                .and_then(|m| m.matmul(&h_transpose))
                .map_err(|err| format!("innovation covariance failed: {err}"))?,
            noise_cov,
        )?;
        let factor = lu_factor(&innovation_cov)
            .map_err(|err| format!("innovation covariance is singular: {err}"))?;
        let determinant = factor.det();
        if determinant <= 0.0 || !determinant.is_finite() {
            return Err("innovation covariance is not positive definite".to_string());
        }
        let innovation_inverse = factor.inverse();

        let gain = predicted_cov
            .matmul(&h_transpose)
            .and_then(|m| m.matmul(&innovation_inverse))
            .map_err(|err| format!("kalman gain failed: {err}"))?;

        let correction = gain.matvec(&innovation);
        let filtered_state: Vec<f64> = predicted_state
            .iter()
            .zip(correction.iter())
            .map(|(prior, delta)| prior + delta)
            .collect();

        let gain_observation = gain
            .matmul(observation)
            .map_err(|err| format!("covariance update failed: {err}"))?;
        let complement = sub(&identity, &gain_observation)?;
        let updated_cov = symmetrise(
            &complement
                .matmul(&predicted_cov)
                .map_err(|err| format!("covariance update failed: {err}"))?,
        );

        let mut quadratic = 0.0;
        for i in 0..measured {
            for j in 0..measured {
                quadratic += innovation[i] * innovation_inverse.at(i, j) * innovation[j];
            }
        }
        loglik -= 0.5
            * (measured as f64 * (2.0 * std::f64::consts::PI).ln() + determinant.ln() + quadratic);

        filtered.push(filtered_state.clone());
        filtered_cov.push(updated_cov.clone());

        predicted_state = transition.matvec(&filtered_state);
        predicted_cov = add(
            &transition
                .matmul(&updated_cov)
                .and_then(|m| m.matmul(&f_transpose))
                .map_err(|err| format!("state prediction failed: {err}"))?,
            process_cov,
        )?;
        predicted_cov = symmetrise(&predicted_cov);
    }

    Ok(FilterRun {
        filtered,
        predicted,
        filtered_cov,
        predicted_cov: predicted_covariances,
        loglik,
    })
}

/// Runs the Kalman filter for a linear Gaussian state-space model
///
/// `x_t = F x_{t-1} + w_t` with `w_t ~ N(0, Q)` and
/// `y_t = H x_t + v_t` with `v_t ~ N(0, R)`, started from
/// `x_0 ~ N(x0, p0)`.
///
/// # Errors
/// Returns a lowercase message when a matrix has the wrong shape, the
/// observation sequence is empty or has rows of the wrong length, a value is
/// not finite, or the innovation covariance is not positive definite.
pub fn kalman_filter(
    y: &[Vec<f64>],
    transition: &[Vec<f64>],
    observation: &[Vec<f64>],
    process_cov: &[Vec<f64>],
    noise_cov: &[Vec<f64>],
    x0: &[f64],
    p0: &[Vec<f64>],
) -> Result<KalmanResult, String> {
    let transition = from_rows(transition, "transition")?;
    let observation = from_rows(observation, "observation")?;
    let process_cov = from_rows(process_cov, "process covariance")?;
    let noise_cov = from_rows(noise_cov, "noise covariance")?;
    let p0 = from_rows(p0, "initial covariance")?;

    let run = run_filter(
        y,
        &transition,
        &observation,
        &process_cov,
        &noise_cov,
        x0,
        &p0,
    )?;

    Ok(KalmanResult {
        filtered: run.filtered,
        predicted: run.predicted,
        covariance: run.filtered_cov.iter().map(to_rows).collect(),
        loglik: run.loglik,
    })
}

/// Rauch–Tung–Striebel smoothed state means.
///
/// Runs the filter and then the backward recursion
/// `x_t^s = x_t + C_t (x_{t+1}^s - x_{t+1})` with
/// `C_t = P_t F' P_{t+1}^{-1}`, returning one smoothed state per observation.
///
/// # Errors
/// Returns the errors of [`kalman_filter`], plus a message when a predicted
/// covariance is singular and the smoother cannot be formed.
pub fn kalman_smooth(
    y: &[Vec<f64>],
    transition: &[Vec<f64>],
    observation: &[Vec<f64>],
    process_cov: &[Vec<f64>],
    noise_cov: &[Vec<f64>],
    x0: &[f64],
    p0: &[Vec<f64>],
) -> Result<Vec<Vec<f64>>, String> {
    let transition = from_rows(transition, "transition")?;
    let observation = from_rows(observation, "observation")?;
    let process_cov = from_rows(process_cov, "process covariance")?;
    let noise_cov = from_rows(noise_cov, "noise covariance")?;
    let p0 = from_rows(p0, "initial covariance")?;

    let run = run_filter(
        y,
        &transition,
        &observation,
        &process_cov,
        &noise_cov,
        x0,
        &p0,
    )?;

    let f_transpose = transition.transpose();
    let mut smoothed = run.filtered.clone();
    let mut covariances = run.filtered_cov.clone();

    for t in (0..y.len().saturating_sub(1)).rev() {
        let next_predicted_cov = run
            .predicted_cov
            .get(t + 1)
            .ok_or("predicted covariance index out of range")?;
        let factor = lu_factor(next_predicted_cov)
            .map_err(|err| format!("predicted covariance is singular: {err}"))?;
        let inverse = factor.inverse();

        let gain = covariances
            .get(t)
            .ok_or("filtered covariance index out of range")?
            .matmul(&f_transpose)
            .and_then(|m| m.matmul(&inverse))
            .map_err(|err| format!("smoother gain failed: {err}"))?;

        let next_smoothed = smoothed.get(t + 1).ok_or("smoothed index out of range")?;
        let next_predicted = run
            .predicted
            .get(t + 1)
            .ok_or("predicted state index out of range")?;
        let difference: Vec<f64> = next_smoothed
            .iter()
            .zip(next_predicted.iter())
            .map(|(smoothed_value, predicted_value)| smoothed_value - predicted_value)
            .collect();
        let correction = gain.matvec(&difference);

        let current = smoothed.get(t).ok_or("smoothed index out of range")?;
        smoothed[t] = current
            .iter()
            .zip(correction.iter())
            .map(|(value, delta)| value + delta)
            .collect();

        let covariance_gap = sub(
            covariances
                .get(t + 1)
                .ok_or("covariance index out of range")?,
            next_predicted_cov,
        )?;
        let updated = add(
            covariances.get(t).ok_or("covariance index out of range")?,
            &gain
                .matmul(&covariance_gap)
                .and_then(|m| m.matmul(&gain.transpose()))
                .map_err(|err| format!("covariance smoothing failed: {err}"))?,
        )?;
        covariances[t] = symmetrise(&updated);
    }

    Ok(smoothed)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Scalar random walk observed without measurement noise.
    fn random_walk_setup() -> (
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
        Vec<Vec<f64>>,
        Vec<f64>,
        Vec<Vec<f64>>,
    ) {
        let y = vec![vec![1.0], vec![2.0], vec![3.0], vec![4.0], vec![5.0]];
        let transition = vec![vec![1.0]];
        let observation = vec![vec![1.0]];
        let process_cov = vec![vec![0.1]];
        let noise_cov = vec![vec![1e-12]];
        let x0 = vec![0.0];
        let p0 = vec![vec![1.0]];
        (y, transition, observation, process_cov, noise_cov, x0, p0)
    }

    #[test]
    fn exact_observations_are_reproduced() {
        let (y, transition, observation, process_cov, noise_cov, x0, p0) = random_walk_setup();
        let result = kalman_filter(
            &y,
            &transition,
            &observation,
            &process_cov,
            &noise_cov,
            &x0,
            &p0,
        )
        .expect("filter");

        assert_eq!(result.filtered.len(), y.len());
        assert_eq!(result.covariance.len(), y.len());
        for (t, state) in result.filtered.iter().enumerate() {
            assert!((state[0] - y[t][0]).abs() < 1e-6, "t {t}: {}", state[0]);
        }
        assert!(result.loglik.is_finite());
    }

    #[test]
    fn predicted_state_leads_the_filter() {
        let (y, transition, observation, process_cov, noise_cov, x0, p0) = random_walk_setup();
        let result = kalman_filter(
            &y,
            &transition,
            &observation,
            &process_cov,
            &noise_cov,
            &x0,
            &p0,
        )
        .expect("filter");
        // The first prediction is the prior mean, the initial state.
        assert!((result.predicted[0][0] - 0.0).abs() < 1e-12);
        // Later predictions repeat the previous filtered state for a walk.
        for t in 1..y.len() {
            assert!((result.predicted[t][0] - result.filtered[t - 1][0]).abs() < 1e-9);
        }
    }

    #[test]
    fn smoother_matches_the_filter_when_observations_are_exact() {
        let (y, transition, observation, process_cov, noise_cov, x0, p0) = random_walk_setup();
        let smoothed = kalman_smooth(
            &y,
            &transition,
            &observation,
            &process_cov,
            &noise_cov,
            &x0,
            &p0,
        )
        .expect("smoother");
        assert_eq!(smoothed.len(), y.len());
        for (t, state) in smoothed.iter().enumerate() {
            assert!((state[0] - y[t][0]).abs() < 1e-6, "t {t}: {}", state[0]);
        }
    }

    #[test]
    fn two_dimensional_system_stays_finite() {
        // Constant velocity model: position and velocity, position observed.
        let transition = vec![vec![1.0, 1.0], vec![0.0, 1.0]];
        let observation = vec![vec![1.0, 0.0]];
        let process_cov = vec![vec![0.01, 0.0], vec![0.0, 0.01]];
        let noise_cov = vec![vec![0.5]];
        let x0 = vec![0.0, 1.0];
        let p0 = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let y: Vec<Vec<f64>> = (1..=6).map(|t| vec![t as f64]).collect();

        let result = kalman_filter(
            &y,
            &transition,
            &observation,
            &process_cov,
            &noise_cov,
            &x0,
            &p0,
        )
        .expect("filter");
        for state in &result.filtered {
            assert_eq!(state.len(), 2);
            assert!(state.iter().all(|v| v.is_finite()));
        }
        // The velocity estimate settles near one, matching the data.
        let last = result.filtered.last().expect("non-empty");
        assert!((last[1] - 1.0).abs() < 0.2, "velocity {}", last[1]);

        let smoothed = kalman_smooth(
            &y,
            &transition,
            &observation,
            &process_cov,
            &noise_cov,
            &x0,
            &p0,
        )
        .expect("smoother");
        assert_eq!(smoothed.len(), y.len());
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        let y = vec![vec![1.0]];
        let good = vec![vec![1.0]];
        let x0 = vec![0.0];
        assert!(kalman_filter(&[], &good, &good, &good, &good, &x0, &good).is_err());
        assert!(
            kalman_filter(&y, &good, &good, &good, &good, &[], &good).is_err(),
            "empty initial state must fail"
        );
        let wrong = vec![vec![1.0, 0.0]];
        assert!(
            kalman_filter(&y, &good, &wrong, &good, &good, &x0, &good).is_err(),
            "observation width mismatch must fail"
        );
        let bad_y = vec![vec![1.0, 2.0]];
        assert!(kalman_filter(&bad_y, &good, &good, &good, &good, &x0, &good).is_err());
        let ragged = vec![vec![1.0], vec![2.0, 3.0]];
        assert!(kalman_filter(&y, &ragged, &good, &good, &good, &x0, &good).is_err());
    }
}
