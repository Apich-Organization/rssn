//! Prepared Kerzman–Stein quadrature and Cauchy evaluation.
use num_complex::Complex64;
use std::f64::consts::PI;
use std::f64::consts::TAU;

use super::divide_projective;
use super::integrate_pure_segment;
use crate::kernels::exact_eq;
use crate::kernels::geometry::point_in_polygon;
use crate::kernels::geometry::segment_intersection_2d;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::fft_scratch_for;

/// Reusable quadrature for the Szegő kernel and Riemann map of a smooth
/// counterclockwise Jordan boundary.
pub struct SzegoKerzmanSteinSolver {
    pub(crate) points: Vec<Complex64>,
    pub(crate) weights: Vec<Complex64>,
    pub(crate) center: Complex64,
    at_center: Complex64,
    polygon: Vec<[f64; 2]>,
    scale: f64,
}

fn shape(
    points: &[Complex64],
    center: Complex64,
) -> Result<(Vec<[f64; 2]>, f64), String> {
    if points.len() < 16 || !center.is_finite() || points.iter().any(|value| !value.is_finite()) {
        return Err(
            "Szegő solve needs at least 16 finite boundary samples and a finite center".into(),
        );
    }
    let scale = points
        .iter()
        .fold(center.re.abs().max(center.im.abs()), |largest, value| {
            largest.max(value.re.abs()).max(value.im.abs())
        });
    if scale <= 0.0 {
        return Err("Szegő boundary has zero extent".into());
    }
    let origin = center / scale;
    let polygon: Vec<_> = points
        .iter()
        .map(|&value| {
            let coordinate = value / scale - origin;
            [coordinate.re, coordinate.im]
        })
        .collect();
    let area: f64 = polygon
        .iter()
        .zip(polygon.iter().cycle().skip(1))
        .map(|(a, b)| a[0] * b[1] - a[1] * b[0])
        .sum();
    if area <= 0.0 || !point_in_polygon([0.0, 0.0], &polygon) {
        return Err(
            "Szegő boundary needs counterclockwise orientation and an interior center".into(),
        );
    }
    Ok((polygon, scale))
}

fn dense_solve(
    mut matrix: Scratch<'static, Complex64>,
    mut rhs: Vec<Complex64>,
) -> Result<Vec<Complex64>, String> {
    let n = rhs.len();
    for col in 0..n {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let pivot = (col..n)
            .max_by(|&a, &b| {
                matrix[a * n + col]
                    .norm()
                    .total_cmp(&matrix[b * n + col].norm())
            })
            .ok_or("Missing Szegő collocation pivot")?;
        let diagonal = matrix[pivot * n + col];
        if !diagonal.is_finite() || diagonal.norm() <= 1e-14 {
            return Err("Szegő collocation matrix is singular".into());
        }
        if col != pivot {
            rhs.swap(col, pivot);
            for j in col..n {
                matrix.swap(col * n + j, pivot * n + j);
            }
        }
        let (above, below) = matrix.split_at_mut((col + 1) * n);
        let pivot_row = &above[col * n..];
        let (known, rhs_rows) = rhs.split_at_mut(col + 1);
        let value = known[col];
        for (row, b) in below.chunks_exact_mut(n).zip(rhs_rows) {
            let factor = divide_projective(row[col], diagonal);
            for (cell, &entry) in row[col + 1..].iter_mut().zip(&pivot_row[col + 1..]) {
                *cell -= factor * entry;
            }
            *b -= factor * value;
        }
    }
    for i in (0..n).rev() {
        let (above, known) = rhs.split_at_mut(i + 1);
        let dot: Complex64 = matrix[i * n + i + 1..(i + 1) * n]
            .iter()
            .zip(known)
            .map(|(&a, &mut b)| a * b)
            .sum();
        above[i] = divide_projective(above[i] - dot, matrix[i * n + i]);
        if !above[i].is_finite() {
            return Err("Szegő collocation solution exceeds Float range".into());
        }
    }
    Ok(rhs)
}

impl SzegoKerzmanSteinSolver {
    /// Retains the Szegő quadrature for an interior point and smooth Jordan
    /// boundary sampled counterclockwise at equally spaced parameters.
    ///
    /// # Errors
    /// Rejects invalid sampling geometry, nonfinite samples, singular systems,
    /// and cancellation.
    pub fn solve<F, DF>(
        gamma: F,
        dgamma: DF,
        n_quad: usize,
        a: Complex64,
    ) -> Result<Self, String>
    where
        F: Fn(f64) -> Complex64,
        DF: Fn(f64) -> Complex64,
    {
        if n_quad < 16 || n_quad.checked_mul(n_quad).is_none() || !a.is_finite() {
            return Err(
                "Szegő solve needs at least 16 boundary samples and a finite center".into(),
            );
        }
        let dt = TAU / n_quad as f64;
        let mut points = Vec::with_capacity(n_quad);
        let mut derivatives = Vec::with_capacity(n_quad);
        for i in 0..n_quad {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
            let t = i as f64 * dt;
            let point = gamma(t);
            let derivative = dgamma(t);
            if !point.is_finite() || !derivative.is_finite() || derivative == Complex64::ZERO {
                return Err("Szegő boundary position and nonzero derivative must be finite".into());
            }
            points.push(point);
            derivatives.push(derivative);
        }
        Self::prepare(points, &derivatives, a)
    }

    pub(crate) fn prepare(
        points: Vec<Complex64>,
        derivatives: &[Complex64],
        center: Complex64,
    ) -> Result<Self, String> {
        let n = points.len();
        let (polygon, scale) = shape(&points, center)?;
        if derivatives.len() != n || derivatives.iter().any(|v| !v.is_finite()) {
            return Err("Szegő derivative tensor must match the finite boundary".into());
        }
        for (i, &point) in polygon.iter().enumerate() {
            let next = polygon[(i + 1) % n];
            let collapsed = exact_eq(point[0], next[0]) && exact_eq(point[1], next[1]);
            if collapsed || segment_intersection_2d([0.0, 0.0], [0.0, 0.0], point, next).is_some() {
                return Err(
                    "Szegő center is on the sampled boundary or a boundary edge collapses".into(),
                );
            }
            for j in (i + 2)..n {
                if i == 0 && j + 1 == n {
                    continue;
                }
                if segment_intersection_2d(point, next, polygon[j], polygon[(j + 1) % n]).is_some()
                {
                    return Err("Szegő sampled boundary edges intersect".into());
                }
            }
        }
        let dt = TAU / n as f64;
        let inverse_two_pi_i = Complex64::new(0.0, -1.0 / TAU);
        let chart: Vec<_> = polygon.iter().map(|&[x, y]| Complex64::new(x, y)).collect();
        let mut speeds = Vec::with_capacity(n);
        let mut tangents = Vec::with_capacity(n);
        for &dp in derivatives {
            let dp = dp / scale;
            let speed = dp.norm();
            if !speed.is_finite() || speed <= 0.0 {
                return Err("Szegő boundary tangent collapses in the normalized chart".into());
            }
            speeds.push(speed);
            tangents.push(dp / speed);
        }
        let elements = n
            .checked_mul(n)
            .ok_or("Szegő collocation matrix exceeds addressable storage")?;
        let mut matrix = fft_scratch_for(elements);
        matrix.resize(elements, Complex64::ZERO);
        let mut rhs = Vec::with_capacity(n);
        for (i, (&zi, &ti)) in chart.iter().zip(&tangents).enumerate() {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
            rhs.push((divide_projective(ti, zi) * inverse_two_pi_i).conj());
            matrix[i * n + i] = Complex64::ONE;
            // The kernel is skew-Hermitian under arclength weighting.
            // Form every off-diagonal pair once and scatter both weights.
            for (j, (&zj, &tj)) in chart.iter().zip(&tangents).enumerate().skip(i + 1) {
                let inverse = divide_projective(inverse_two_pi_i, zj - zi);
                let kernel = (-ti * inverse).conj() - tj * inverse;
                matrix[i * n + j] = kernel * (speeds[j] * dt);
                matrix[j * n + i] = -kernel.conj() * (speeds[i] * dt);
            }
        }
        let solution = dense_solve(matrix, rhs)?;
        let weights = solution
            .into_iter()
            .zip(tangents)
            .zip(speeds)
            .map(|((value, tangent), speed)| value * tangent * inverse_two_pi_i * (speed * dt))
            .collect();
        Self::finish(points, weights, center, polygon, scale)
    }

    pub(crate) fn from_weights(
        points: Vec<Complex64>,
        weights: Vec<Complex64>,
        center: Complex64,
    ) -> Result<Self, String> {
        let (polygon, scale) = shape(&points, center)?;
        Self::finish(points, weights, center, polygon, scale)
    }

    fn finish(
        points: Vec<Complex64>,
        weights: Vec<Complex64>,
        center: Complex64,
        polygon: Vec<[f64; 2]>,
        scale: f64,
    ) -> Result<Self, String> {
        if weights.len() != points.len() || weights.iter().any(|v| !v.is_finite()) {
            return Err("Szegő quadrature needs finite weights matching its boundary".into());
        }
        let mut solver = Self {
            points,
            weights,
            center,
            at_center: Complex64::ZERO,
            polygon,
            scale,
        };
        solver.at_center = solver.eval_szego(center);
        if !solver.at_center.is_finite() || solver.at_center.re <= 0.0 {
            return Err("Szegő normalization at the center is nonpositive or nonfinite".into());
        }
        Ok(solver)
    }

    fn polygon_point(
        &self,
        point: Complex64,
    ) -> [f64; 2] {
        let relative = point / self.scale - self.center / self.scale;
        [relative.re, relative.im]
    }

    pub(crate) fn checked_point(
        &self,
        point: Complex64,
    ) -> Result<(), String> {
        let position = self.polygon_point(point);
        if !point.is_finite() || !point_in_polygon(position, &self.polygon) {
            return Err("Szegő evaluation point must lie inside the sampled boundary".into());
        }
        if self
            .polygon
            .iter()
            .zip(self.polygon.iter().cycle().skip(1))
            .any(|(&a, &b)| segment_intersection_2d(position, position, a, b).is_some())
        {
            return Err("Szegő evaluation point lies on the sampled boundary".into());
        }
        Ok(())
    }

    /// Evaluates the Szegő kernel at an interior point `w`.
    #[must_use]
    pub fn eval_szego(
        &self,
        w: Complex64,
    ) -> Complex64 {
        self.points
            .iter()
            .zip(&self.weights)
            .map(|(&point, &weight)| {
                let delta = point - w;
                let extent = delta.re.abs().max(delta.im.abs());
                if (1e-145..=1e145).contains(&extent)
                    && weight.re.abs().max(weight.im.abs()) <= 1e145
                {
                    weight / delta
                } else if delta.is_finite() {
                    divide_projective(weight, delta)
                } else {
                    divide_projective(weight, point / self.scale - w / self.scale) / self.scale
                }
            })
            .sum()
    }

    pub(crate) fn eval_map_checked(
        &self,
        point: Complex64,
    ) -> Result<(Complex64, Complex64), String> {
        self.checked_point(point)?;
        let position = self.polygon_point(point);
        if self
            .polygon
            .iter()
            .zip(self.polygon.iter().cycle().skip(1))
            .any(|(&a, &b)| segment_intersection_2d([0.0, 0.0], position, a, b).is_some())
        {
            return Err("Szegő straight evaluation path leaves the sampled boundary".into());
        }
        let coefficient = divide_projective(Complex64::from(TAU), self.at_center);
        let result = integrate_pure_segment(
            |z| {
                let v = self.eval_szego(z);
                coefficient * v * v
            },
            self.center,
            point,
            1e-9,
            128,
        )
        .map_err(|error| error.to_string())?;
        let at_point = self.eval_szego(point);
        if !result.converged || !at_point.is_finite() {
            return Err("Szegő map quadrature failed to meet its tolerance".into());
        }
        Ok((result.value, at_point))
    }

    /// Evaluates the Riemann map and Szegő kernel along a straight interior path.
    /// An invalid path or unresolved quadrature returns nonfinite results.
    #[must_use]
    pub fn eval_map(
        &self,
        point: Complex64,
    ) -> (Complex64, Complex64) {
        let invalid = Complex64::new(f64::NAN, f64::NAN);
        self.eval_map_checked(point).unwrap_or((invalid, invalid))
    }

    /// Evaluates the Bergman kernel `4 * pi * S(z, a)^2`.
    #[must_use]
    pub fn eval_bergman(
        &self,
        point: Complex64,
    ) -> Complex64 {
        let s = self.eval_szego(point);
        4.0 * PI * s * s
    }
}

/// Computes a Riemann map and the Szegő kernel for a smooth counterclockwise
/// Jordan boundary. The segment from the center to the query must stay inside.
///
/// # Errors
/// Rejects invalid boundary data, failed collocation, and unresolved quadrature.
pub fn szego_kerzman_stein_map<F, DF>(
    gamma: F,
    dgamma: DF,
    n_quad: usize,
    center: Complex64,
    point: Complex64,
) -> Result<(Complex64, Complex64), String>
where
    F: Fn(f64) -> Complex64,
    DF: Fn(f64) -> Complex64,
{
    SzegoKerzmanSteinSolver::solve(gamma, dgamma, n_quad, center)?.eval_map_checked(point)
}
