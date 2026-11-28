//! Linear Gaussian state-space filtering and Rauch–Tung–Striebel smoothing.

use crate::kernels::dense::DenseError;
use crate::kernels::dense::Mat;
use crate::kernels::dense::cholesky;
use crate::kernels::dense::solve_cholesky_factor_in_place;
use crate::kernels::dense::solve_cholesky_lower_in_place;

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

struct KalmanModel {
    transition: Mat,
    observation: Mat,
    process_cov: Mat,
    noise_cov: Mat,
    p0: Mat,
}

impl KalmanModel {
    fn from_rows(
        transition: &[Vec<f64>],
        observation: &[Vec<f64>],
        process_cov: &[Vec<f64>],
        noise_cov: &[Vec<f64>],
        p0: &[Vec<f64>],
    ) -> Result<Self, String> {
        Ok(Self {
            transition: from_rows(transition, "transition")?,
            observation: from_rows(observation, "observation")?,
            process_cov: covariance_from_rows(process_cov, "process covariance")?,
            noise_cov: covariance_from_rows(noise_cov, "noise covariance")?,
            p0: covariance_from_rows(p0, "initial covariance")?,
        })
    }
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

/// Adds two equally shaped matrices, reusing the left allocation.
fn add(
    mut a: Mat,
    b: &Mat,
) -> Result<Mat, String> {
    if a.rows != b.rows || a.cols != b.cols {
        return Err("matrix shapes differ in an addition".to_string());
    }
    for (left, &right) in a.data.iter_mut().zip(&b.data) {
        *left += right;
    }
    Ok(a)
}

/// Reuses the dense product and attaches the failing step to its error.
fn multiply(
    a: &Mat,
    b: &Mat,
    stage: &str,
) -> Result<Mat, String> {
    a.matmul(b).map_err(|err| format!("{stage} failed: {err}"))
}

fn symmetrise(a: &mut Mat) {
    for i in 0..a.rows {
        for j in (i + 1)..a.rows {
            let mean = a.at(i, j).midpoint(a.at(j, i));
            a.set(i, j, mean);
            a.set(j, i, mean);
        }
    }
}

fn covariance_from_rows(
    rows: &[Vec<f64>],
    name: &str,
) -> Result<Mat, String> {
    let mut cov = from_rows(rows, name)?;
    if cov.rows != cov.cols || !cov.is_finite() {
        return Err(format!("{name} must be square and finite"));
    }
    let mut deviations = crate::kernels::reuse::f64_scratch(cov.rows);
    for i in 0..cov.rows {
        let variance = cov.at(i, i);
        if variance < 0.0 {
            return Err(format!("{name} has a negative variance"));
        }
        deviations.push(variance.sqrt());
    }

    // Check the correlation matrix so a large variance cannot mask an
    // invalid covariance among small components. A diagonal perturbation
    // only on this validation copy admits singular PSD matrices and lets
    // the same checked Cholesky factorization reject negative directions.
    let tolerance = 16.0 * f64::EPSILON * cov.rows as f64;
    let mut correlations = None;
    for i in 0..cov.rows {
        for j in (i + 1)..cov.rows {
            let a = cov.at(i, j);
            let b = cov.at(j, i);
            if (a - b).abs() > 16.0 * f64::EPSILON * a.abs().max(b.abs()) {
                return Err(format!("{name} must be symmetric"));
            }
            let mean = a.midpoint(b);
            cov.set(i, j, mean);
            cov.set(j, i, mean);
            if mean == 0.0 {
                continue;
            }
            let large = deviations[i].max(deviations[j]);
            let small = deviations[i].min(deviations[j]);
            if small == 0.0 {
                return Err(format!("{name} must be positive semidefinite"));
            }
            let correlation = (mean / large) / small;
            if correlation.abs() > 1.0 + tolerance {
                return Err(format!("{name} must be positive semidefinite"));
            }
            let matrix = correlations.get_or_insert_with(|| {
                let mut matrix = Mat::identity(cov.rows);
                for k in 0..cov.rows {
                    matrix.set(k, k, 1.0 + tolerance);
                }
                matrix
            });
            matrix.set(i, j, correlation);
            matrix.set(j, i, correlation);
        }
    }
    if let Some(correlations) = correlations {
        cholesky(&correlations).map_err(|error| {
            match error {
                | DenseError::Singular => format!("{name} must be positive semidefinite"),
                | other => format!("could not validate {name}: {other}"),
            }
        })?;
    }
    Ok(cov)
}

/// Applies the same factored right solve to each independent row.
fn solve_gain(
    mut cross: Mat,
    factor: &Mat,
) -> Mat {
    let cols = cross.cols;
    super::super::par::for_chunks_mut(&mut cross.data, cols, 32, |_, row| {
        solve_cholesky_factor_in_place(factor, row);
    });
    cross
}

/// Joseph covariance update, using the gain from the factored innovation.
fn posterior_covariance(
    gain: &Mat,
    observation: &Mat,
    prior: &Mat,
    noise: &Mat,
) -> Result<Mat, String> {
    const UPDATE: &str = "covariance update";
    let single_observation = gain.cols == 1;
    let mut complement = if single_observation {
        let mut complement = Mat::zeros(gain.rows, observation.cols);
        for (row, &weight) in complement
            .data
            .chunks_exact_mut(observation.cols)
            .zip(&gain.data)
        {
            for (entry, &coefficient) in row.iter_mut().zip(&observation.data) {
                *entry = -weight * coefficient;
            }
        }
        complement
    } else {
        let mut complement = multiply(gain, observation, UPDATE)?;
        for value in &mut complement.data {
            *value = -*value;
        }
        complement
    };
    for i in 0..complement.rows {
        complement.set(i, i, 1.0 + complement.at(i, i));
    }
    let mut covariance = multiply(
        &multiply(&complement, prior, UPDATE)?,
        &complement.transpose(),
        UPDATE,
    )?;
    if single_observation {
        // K R K' is a scaled outer product when the measurement is scalar.
        let variance = noise.at(0, 0);
        for (row, &weight) in covariance
            .data
            .chunks_exact_mut(covariance.cols)
            .zip(&gain.data)
        {
            let scaled_weight = weight * variance;
            if scaled_weight != 0.0 {
                for (entry, &other_weight) in row.iter_mut().zip(&gain.data) {
                    *entry += scaled_weight * other_weight;
                }
            }
        }
    } else {
        let from_noise = multiply(&multiply(gain, noise, UPDATE)?, &gain.transpose(), UPDATE)?;
        covariance = add(covariance, &from_noise)?;
    }
    symmetrise(&mut covariance);
    Ok(covariance)
}

/// Runs the prediction and update steps. Only the smoother needs the
/// predicted covariance history.
fn run_filter(
    y: &[Vec<f64>],
    model: &KalmanModel,
    x0: &[f64],
    for_smoothing: bool,
) -> Result<FilterRun, String> {
    let KalmanModel {
        transition,
        observation,
        process_cov,
        noise_cov,
        p0,
    } = model;
    let state = x0.len();
    if state == 0 || x0.iter().any(|value| !value.is_finite()) {
        return Err("initial state must be nonempty and finite".to_string());
    }
    if transition.rows != state || transition.cols != state || !transition.is_finite() {
        return Err("transition must be finite and square with the state dimension".to_string());
    }
    if process_cov.rows != state || p0.rows != state {
        return Err("process and initial covariances must match the state dimension".to_string());
    }
    let measured = observation.rows;
    if observation.cols != state || measured == 0 || !observation.is_finite() {
        return Err("observation matrix must be finite and have one column per state".to_string());
    }
    if noise_cov.rows != measured {
        return Err("noise covariance must match the observation dimension".to_string());
    }
    if y.is_empty() {
        return Err("observation sequence is empty".to_string());
    }
    for (t, row) in y.iter().enumerate() {
        if row.len() != measured || row.iter().any(|v| !v.is_finite()) {
            return Err(format!(
                "observation {t} has the wrong size or a non-finite value"
            ));
        }
    }

    let h_transpose = observation.transpose();
    let f_transpose = transition.transpose();
    let mut filtered: Vec<Vec<f64>> = Vec::with_capacity(y.len());
    let mut predicted = Vec::with_capacity(y.len());
    let mut filtered_cov: Vec<Mat> = Vec::with_capacity(y.len());
    let mut predicted_cov = if for_smoothing {
        Vec::with_capacity(y.len())
    } else {
        Vec::new()
    };
    let mut innovation = crate::kernels::reuse::f64_scratch(measured);
    innovation.resize(measured, 0.0);
    let mut loglik = 0.0;
    let gaussian_constant = measured as f64 * (2.0 * std::f64::consts::PI).ln();

    for row in y {
        const PREDICT: &str = "state prediction";
        let prior_state = filtered.last().map_or(x0, Vec::as_slice);
        let prior_cov = filtered_cov.last().unwrap_or(p0);
        let predicted_state = transition.matvec(prior_state);
        let mut predicted_matrix = add(
            multiply(
                &multiply(transition, prior_cov, PREDICT)?,
                &f_transpose,
                PREDICT,
            )?,
            process_cov,
        )?;
        symmetrise(&mut predicted_matrix);

        observation.matvec_into(&predicted_state, &mut innovation);
        for (error, &observed) in innovation.iter_mut().zip(row) {
            *error = observed - *error;
        }

        let cross = multiply(
            &predicted_matrix,
            &h_transpose,
            "state/observation covariance",
        )?;
        let mut innovation_cov = add(
            multiply(observation, &cross, "innovation covariance")?,
            noise_cov,
        )?;
        symmetrise(&mut innovation_cov);
        let factor =
            cholesky(&innovation_cov).map_err(|err| format!("innovation covariance: {err}"))?;

        let (mut filtered_state, updated_cov) = if state == 1 && measured == 1 {
            let innovation_std = factor.at(0, 0);
            // Scale both operands by the innovation deviation so a gain
            // outside binary64 does not overflow an in-range correction.
            let correction = (cross.at(0, 0) / innovation_std) * (innovation[0] / innovation_std);

            // In the scalar case P_updated = P_prior R / S. Taking the
            // ratio of deviations avoids both subtracting 1 - K H and
            // underflowing R / S before multiplication by P_prior.
            let ratio = noise_cov.at(0, 0).sqrt() / innovation_std;
            let mut covariance = Mat::zeros(1, 1);
            covariance.set(0, 0, ratio * (ratio * predicted_matrix.at(0, 0)));
            (vec![correction], covariance)
        } else {
            let gain = solve_gain(cross, &factor);
            let covariance =
                posterior_covariance(&gain, observation, &predicted_matrix, noise_cov)?;
            (gain.matvec(&innovation), covariance)
        };
        for (estimate, &prior) in filtered_state.iter_mut().zip(&predicted_state) {
            *estimate += prior;
        }

        if !for_smoothing {
            // L^-1 innovation supplies the quadratic form. The logs of
            // the diagonal give the log determinant of the covariance.
            solve_cholesky_lower_in_place(&factor, &mut innovation);
            let log_sqrt_det = (0..measured).map(|i| factor.at(i, i).ln()).sum::<f64>();
            // Apply the Gaussian half-factor before squaring so a finite
            // likelihood survives when the unscaled squared norm overflows.
            crate::kernels::simd::scale_f64(std::f64::consts::FRAC_1_SQRT_2, &mut innovation);
            let half_quadratic = crate::kernels::simd::norm2_sq_f64(&innovation);
            loglik -= 0.5 * gaussian_constant + log_sqrt_det + half_quadratic;
        }
        if !loglik.is_finite()
            || !updated_cov.is_finite()
            || filtered_state.iter().any(|value| !value.is_finite())
        {
            return Err("filter update produced a non-finite result".to_string());
        }

        filtered.push(filtered_state);
        filtered_cov.push(updated_cov);
        predicted.push(predicted_state);
        if for_smoothing {
            predicted_cov.push(predicted_matrix);
        }
    }

    Ok(FilterRun {
        filtered,
        predicted,
        filtered_cov,
        predicted_cov,
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
/// Returns a diagnostic for invalid shapes, non-finite data, nonsymmetric
/// or indefinite model covariances and an innovation covariance that is not
/// positive definite. Model covariances may be positive semidefinite.
pub fn kalman_filter(
    y: &[Vec<f64>],
    transition: &[Vec<f64>],
    observation: &[Vec<f64>],
    process_cov: &[Vec<f64>],
    noise_cov: &[Vec<f64>],
    x0: &[f64],
    p0: &[Vec<f64>],
) -> Result<KalmanResult, String> {
    let model = KalmanModel::from_rows(transition, observation, process_cov, noise_cov, p0)?;
    let run = run_filter(y, &model, x0, false)?;

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
    let model = KalmanModel::from_rows(transition, observation, process_cov, noise_cov, p0)?;
    let run = run_filter(y, &model, x0, true)?;

    let f_transpose = model.transition.transpose();
    let mut smoothed = run.filtered;
    let mut difference = crate::kernels::reuse::f64_scratch(x0.len());
    difference.resize(x0.len(), 0.0);
    let mut correction = crate::kernels::reuse::f64_scratch(x0.len());
    correction.resize(x0.len(), 0.0);

    // The mean recursion uses each original filtered covariance. There
    // is no need to construct the smoothed covariance sequence.
    for t in (0..y.len().saturating_sub(1)).rev() {
        let factor = cholesky(&run.predicted_cov[t + 1])
            .map_err(|err| format!("predicted covariance: {err}"))?;
        let cross = multiply(&run.filtered_cov[t], &f_transpose, "smoother gain")?;
        let gain = solve_gain(cross, &factor);
        for ((delta, &posterior), &prediction) in difference
            .iter_mut()
            .zip(&smoothed[t + 1])
            .zip(&run.predicted[t + 1])
        {
            *delta = posterior - prediction;
        }
        gain.matvec_into(&difference, &mut correction);
        for (mean, &delta) in smoothed[t].iter_mut().zip(correction.iter()) {
            *mean += delta;
        }
        if smoothed[t].iter().any(|value| !value.is_finite()) {
            return Err("smoothing update produced a non-finite result".to_string());
        }
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
        let y: Vec<Vec<f64>> = (1..=6).map(|t| vec![f64::from(t)]).collect();

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
