//! # Numerical differential geometry
//!
//! The metric callback returns the covariant matrix `g_ij` in one coordinate
//! chart. The kernel uses fourth-order central differences with relative step
//! `1e-3`. Metric samples must be finite and symmetric. Indefinite metrics are
//! supported when the local chart is nonsingular.
//!
//! Index conventions give `Gamma[k][i][j]` = `Γ^k_{ij}` and
//! `Riemann[r][s][m][n]` = `R^r_{smn}`. The curvature convention is
//! `∂_m Γ^r_{sn} - ∂_n Γ^r_{sm} + Γ^r_{lm} Γ^l_{sn} - Γ^r_{ln} Γ^l_{sm}`.
//! Ricci contracts the first and third indices. A sphere has a positive
//! Ricci scalar with this convention.
//!
//! Curvature computations reuse samples at identical coordinates during one
//! evaluation. The callback defines a single metric at each point and runs
//! on the calling thread. Dense Riemann contractions can run in parallel
//! once the stencil has been evaluated. Chart boundaries, poor conditioning
//! and rapidly varying metrics limit the accuracy of the fixed difference step.

#![allow(
    clippy::indexing_slicing,
    reason = "metric samples are shape-checked before the tensor contractions use chart indices"
)]

use std::collections::HashMap;
use std::rc::Rc;

use crate::kernels::dense::DenseError;
use crate::kernels::dense::Mat;
use crate::kernels::dense::lu_factor;
use crate::kernels::par::try_map_indexed;

type Metric = Vec<Vec<f64>>;

const STEP: f64 = 1e-3;
const SYMMETRY_TOLERANCE: f64 = 1e-10;
const RIEMANN_PARALLEL_DIM: usize = 12;
// Squaring this bound is strictly below f64::MAX.
const FINITE_CONNECTION_PRODUCT_BOUND: f64 = 1e154;

/// Row-major working cube. Dense row `k * dim + i` stores `tensor[k][i]`.
struct Tensor3(Mat);

impl Tensor3 {
    fn zeros(dim: usize) -> Result<Self, String> {
        let rows = dim
            .checked_mul(dim)
            .filter(|rows| rows.checked_mul(dim).is_some())
            .ok_or_else(|| "metric dimension is too large".to_string())?;
        Ok(Self(Mat::zeros(rows, dim)))
    }

    #[allow(
        clippy::arithmetic_side_effects,
        reason = "validated chart indices fit the allocated matrix"
    )]
    const fn row(
        &self,
        k: usize,
        i: usize,
    ) -> usize {
        k * self.0.cols + i
    }

    fn at(
        &self,
        k: usize,
        i: usize,
        j: usize,
    ) -> f64 {
        self.0.at(self.row(k, i), j)
    }

    fn set(
        &mut self,
        k: usize,
        i: usize,
        j: usize,
        value: f64,
    ) {
        let row = self.row(k, i);
        self.0.set(row, j, value);
    }

    #[allow(
        clippy::single_call_fn,
        reason = "only the public Christoffel output needs nested rows"
    )]
    fn into_nested(self) -> Vec<Metric> {
        let dim = self.0.cols;
        let mut values = self.0.data.into_iter();
        (0..dim)
            .map(|_| {
                (0..dim)
                    .map(|_| values.by_ref().take(dim).collect())
                    .collect()
            })
            .collect()
    }
}

fn stencil<T>(
    point: &[f64],
    axis: usize,
    mut sample: impl FnMut(&[f64]) -> Result<T, String>,
) -> Result<([T; 4], f64), String> {
    let coordinate = point[axis];
    let h = STEP * coordinate.abs().max(1.0);
    let mut shifted = point.to_vec();
    let mut at = |offset: f64| {
        let shifted_coordinate = coordinate + offset * h;
        if !shifted_coordinate.is_finite() {
            return Err("metric stencil coordinates must be finite".to_string());
        }
        shifted[axis] = shifted_coordinate;
        sample(&shifted)
    };
    Ok(([at(2.0)?, at(1.0)?, at(-1.0)?, at(-2.0)?], h))
}

const fn derivative(
    [p2, p1, m1, m2]: [f64; 4],
    h: f64,
) -> f64 {
    (8.0 * (p1 - m1) - (p2 - m2)) / (12.0 * h)
}

fn checked_metric(
    g: &impl Fn(&[f64]) -> Metric,
    point: &[f64],
) -> Result<(Metric, f64), String> {
    if point.iter().any(|x| !x.is_finite()) {
        return Err("metric point coordinates must be finite".to_string());
    }
    let metric = g(point);
    if metric.len() != point.len() || metric.iter().any(|row| row.len() != point.len()) {
        return Err(format!(
            "metric must be a {0}x{0} matrix at a point with {0} coordinates",
            point.len()
        ));
    }
    let mut scale = 0.0_f64;
    for value in metric.iter().flatten() {
        if !value.is_finite() {
            return Err("metric tensor entries must be finite".to_string());
        }
        scale = scale.max(value.abs());
    }
    for (i, row) in metric.iter().enumerate() {
        for (j, &value) in row.iter().take(i).enumerate() {
            if (value - metric[j][i]).abs() > SYMMETRY_TOLERANCE * scale {
                return Err("metric tensor must be symmetric".to_string());
            }
        }
    }
    Ok((metric, scale))
}

fn normalize_metric(
    mut metric: Metric,
    scale: f64,
) -> Result<Metric, String> {
    let scaled = |value: f64| {
        let value = value / scale;
        if value.is_finite() {
            Ok(value)
        } else {
            Err("metric stencil is outside the representable range".to_string())
        }
    };
    // Remove small discrepancies admitted by the numerical symmetry check.
    for i in 0..metric.len() {
        metric[i][i] = scaled(metric[i][i])?;
        for j in 0..i {
            let value = scaled(metric[i][j].midpoint(metric[j][i]))?;
            metric[i][j] = value;
            metric[j][i] = value;
        }
    }
    Ok(metric)
}

struct MetricSamples<'a, G> {
    evaluate: &'a G,
    scale: f64,
    cache: Option<HashMap<Vec<u64>, Rc<Metric>>>,
    key: Vec<u64>,
}

impl<'a, G: Fn(&[f64]) -> Metric> MetricSamples<'a, G> {
    fn start(
        evaluate: &'a G,
        point: &[f64],
        reuse: bool,
    ) -> Result<(Self, Rc<Metric>), String> {
        let (metric, scale) = checked_metric(evaluate, point)?;
        let scale = if point.is_empty() {
            1.0
        } else if scale == 0.0 {
            return Err("metric tensor is singular".to_string());
        } else {
            scale
        };
        // A constant rescaling leaves the connection and Ricci tensor
        // unchanged. Apply the opposite factor only to scalar curvature.
        let center = Rc::new(normalize_metric(metric, scale)?);
        let cache = reuse.then(|| {
            let key = point.iter().map(|value| value.to_bits()).collect();
            HashMap::from([(key, Rc::clone(&center))])
        });
        Ok((
            Self {
                evaluate,
                scale,
                cache,
                key: Vec::new(),
            },
            center,
        ))
    }

    fn at(
        &mut self,
        point: &[f64],
    ) -> Result<Rc<Metric>, String> {
        if let Some(cache) = &self.cache {
            self.key.clear();
            self.key.extend(point.iter().map(|value| value.to_bits()));
            if let Some(metric) = cache.get(&self.key) {
                return Ok(Rc::clone(metric));
            }
        }
        let (metric, _) = checked_metric(self.evaluate, point)?;
        let metric = Rc::new(normalize_metric(metric, self.scale)?);
        if let Some(cache) = &mut self.cache {
            cache.insert(self.key.clone(), Rc::clone(&metric));
        }
        Ok(metric)
    }
}

/// Evaluates the metric at `point`, retaining the callback's entries.
///
/// # Errors
///
/// Returns an error for non-finite coordinates or a non-finite, asymmetric
/// or incorrectly shaped metric. Symmetry is checked relative to the largest
/// absolute metric entry.
pub fn metric_tensor_at_point(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<Vec<Vec<f64>>, String> {
    checked_metric(&g, point).map(|(metric, _)| metric)
}

fn connection(
    samples: &mut MetricSamples<'_, impl Fn(&[f64]) -> Metric>,
    point: &[f64],
    metric: &[Vec<f64>],
) -> Result<(Tensor3, Mat), String> {
    let dim = point.len();
    let matrix = Mat::from_rows(metric).map_err(|error| error.to_string())?;
    let inverse = lu_factor(&matrix)
        .map_err(|error| {
            match error {
                | DenseError::Singular => "metric tensor is singular".to_string(),
                | _ => error.to_string(),
            }
        })?
        .inverse();
    if !inverse.is_finite() {
        return Err("metric inverse is outside the representable range".to_string());
    }

    // dg[k][i][j] = ∂_k g_ij
    let mut dg = Tensor3::zeros(dim)?;
    for axis in 0..dim {
        let ([p2, p1, m1, m2], h) = stencil(point, axis, |shifted| samples.at(shifted))?;
        for i in 0..dim {
            for j in i..dim {
                let value = derivative([p2[i][j], p1[i][j], m1[i][j], m2[i][j]], h);
                dg.set(axis, i, j, value);
                dg.set(axis, j, i, value);
            }
        }
    }
    let mut gamma = Tensor3::zeros(dim)?;
    let mut small_terms = [0.0; 6];
    let mut large_terms;
    let terms = if let Some(terms) = small_terms.get_mut(..dim) {
        terms
    } else {
        large_terms = vec![0.0; dim];
        &mut large_terms
    };
    for i in 0..dim {
        for j in i..dim {
            for (m, term) in terms.iter_mut().enumerate() {
                *term = dg.at(j, m, i) + dg.at(i, m, j) - dg.at(m, i, j);
            }
            for (k, inverse_row) in inverse.data.chunks_exact(dim).enumerate() {
                let sum: f64 = inverse_row
                    .iter()
                    .zip(terms.iter())
                    .map(|(entry, term)| entry * term)
                    .sum();
                let value = 0.5 * sum;
                if !value.is_finite() {
                    return Err(
                        "Christoffel symbols are outside the representable range".to_string()
                    );
                }
                gamma.set(k, i, j, value);
                gamma.set(k, j, i, value);
            }
        }
    }
    Ok((gamma, inverse))
}

/// Christoffel symbols of the second kind `Γ^k_{ij}` at `point`, computed as
/// `½ g^{km} (∂_j g_{mi} + ∂_i g_{mj} - ∂_m g_{ij})`.
///
/// # Errors
///
/// Returns an error for an invalid metric sample, a singular metric at
/// `point`, or an unrepresentable stencil, inverse or result.
pub fn christoffel_symbols(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<Vec<Vec<Vec<f64>>>, String> {
    let (mut samples, center) = MetricSamples::start(&g, point, false)?;
    connection(&mut samples, point, &center).map(|(gamma, _)| gamma.into_nested())
}

fn curvature_stencil(
    g: &impl Fn(&[f64]) -> Metric,
    point: &[f64],
    mut accumulate: impl FnMut(usize, [Tensor3; 4], f64),
) -> Result<(Tensor3, Mat, f64), String> {
    let (mut samples, center) = MetricSamples::start(g, point, true)?;
    let (gamma, inverse) = connection(&mut samples, point, &center)?;
    for axis in 0..point.len() {
        let (neighbors, h) = stencil(point, axis, |shifted| {
            let metric = samples.at(shifted)?;
            connection(&mut samples, shifted, &metric).map(|(gamma, _)| gamma)
        })?;
        accumulate(axis, neighbors, h);
    }
    Ok((gamma, inverse, samples.scale))
}

fn connection_products(
    gamma: &Tensor3,
    rho: usize,
    sigma: usize,
    mu: usize,
    nu: usize,
) -> f64 {
    (0..gamma.0.cols)
        .map(|l| {
            gamma.at(rho, l, mu) * gamma.at(l, sigma, nu)
                - gamma.at(rho, l, nu) * gamma.at(l, sigma, mu)
        })
        .sum()
}

/// Riemann curvature tensor `R^ρ_{σμν}` at `point`.
///
/// # Errors
///
/// Returns an error for an invalid or singular metric at a connection stencil
/// point, or an unrepresentable stencil or result.
pub fn riemann_tensor(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<Vec<Vec<Vec<Vec<f64>>>>, String> {
    let dim = point.len();
    let mut d_gamma = Vec::with_capacity(dim);
    let (gamma, _, _) = curvature_stencil(&g, point, |_, [p2, p1, m1, m2], h| {
        let mut differential = p2;
        let neighbors = p1.0.data.into_iter().zip(m1.0.data).zip(m2.0.data);
        for (value, ((p1, m1), m2)) in differential.0.data.iter_mut().zip(neighbors) {
            *value = derivative([*value, p1, m1, m2], h);
        }
        d_gamma.push(differential);
    })?;
    let overflow = || "Riemann tensor is outside the representable range".to_string();
    if d_gamma.iter().any(|tensor| !tensor.0.is_finite()) {
        return Err(overflow());
    }
    // Below this conservative bound every connection product is finite.
    // Otherwise explicitly check the diagonal for overflow before mirroring.
    let check_diagonal = gamma
        .0
        .data
        .iter()
        .any(|value| value.abs() > FINITE_CONNECTION_PRODUCT_BOUND);
    try_map_indexed(dim, RIEMANN_PARALLEL_DIM, |rho| {
        let component = |sigma, mu, nu| {
            let products = connection_products(&gamma, rho, sigma, mu, nu);
            d_gamma[mu].at(rho, sigma, nu) - d_gamma[nu].at(rho, sigma, mu) + products
        };
        (0..dim)
            .map(|sigma| {
                let mut plane = vec![vec![0.0; dim]; dim];
                for mu in 0..dim {
                    if check_diagonal && !component(sigma, mu, mu).is_finite() {
                        return Err(overflow());
                    }
                    for nu in (mu + 1)..dim {
                        let value = component(sigma, mu, nu);
                        if !value.is_finite() {
                            return Err(overflow());
                        }
                        plane[mu][nu] = value;
                        // Exact zero may have a different sign from its negation.
                        plane[nu][mu] = if value == 0.0 {
                            component(sigma, nu, mu)
                        } else {
                            -value
                        };
                    }
                }
                Ok(plane)
            })
            .collect()
    })
}

fn ricci_contraction(
    g: &impl Fn(&[f64]) -> Metric,
    point: &[f64],
) -> Result<(Metric, Mat, f64), String> {
    let dim = point.len();
    let mut ricci = vec![vec![0.0; dim]; dim];
    let (gamma, inverse, scale) = curvature_stencil(g, point, |axis, [p2, p1, m1, m2], h| {
        let d_gamma = |r: usize, s: usize, n: usize| {
            derivative(
                [
                    p2.at(r, s, n),
                    p1.at(r, s, n),
                    m1.at(r, s, n),
                    m2.at(r, s, n),
                ],
                h,
            )
        };
        for (sigma, row) in ricci.iter_mut().enumerate() {
            for (nu, value) in row.iter_mut().enumerate() {
                *value += d_gamma(axis, sigma, nu);
            }
            row[axis] -= (0..dim).map(|mu| d_gamma(mu, sigma, mu)).sum::<f64>();
        }
    })?;
    for (sigma, row) in ricci.iter_mut().enumerate() {
        for (nu, value) in row.iter_mut().enumerate() {
            *value += (0..dim)
                .map(|mu| connection_products(&gamma, mu, sigma, mu, nu))
                .sum::<f64>();
            if !value.is_finite() {
                return Err("Ricci tensor is outside the representable range".to_string());
            }
        }
    }
    Ok((ricci, inverse, scale))
}

/// Ricci tensor `R_{σν} = R^μ_{σμν}` at `point`.
///
/// The derivative and product terms are contracted as they are evaluated.
///
/// # Errors
///
/// Returns an error for an invalid or singular metric at a connection stencil
/// point, or an unrepresentable stencil or result.
pub fn ricci_tensor(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<Vec<Vec<f64>>, String> {
    ricci_contraction(&g, point).map(|(tensor, _, _)| tensor)
}

/// Ricci tensor and its scalar contraction from one metric evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct RicciCurvature {
    /// Covariant Ricci tensor `R_{ij}`.
    pub tensor: Vec<Vec<f64>>,
    /// Scalar curvature `g^{ij} R_{ij}`.
    pub scalar: f64,
}

/// Computes the Ricci tensor and scalar with one curvature stencil and one
/// factorization of the metric at `point`.
///
/// ```
/// use rssn::kernels::differential_geometry::ricci_curvature;
/// let flat = |_: &[f64]| vec![vec![1.0, 0.0], vec![0.0, 1.0]];
/// let curvature = ricci_curvature(flat, &[0.0, 0.0])?;
/// assert_eq!(curvature.scalar, 0.0);
/// # Ok::<(), String>(())
/// ```
///
/// # Errors
///
/// Returns an error for an invalid or singular metric at a connection stencil
/// point, or an unrepresentable stencil or result.
pub fn ricci_curvature(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<RicciCurvature, String> {
    let (tensor, inverse, scale) = ricci_contraction(&g, point)?;
    let scalar = inverse
        .data
        .iter()
        .zip(tensor.iter().flatten())
        .map(|(g, r)| g * r)
        .sum::<f64>()
        / scale;
    if !scalar.is_finite() {
        return Err("Ricci scalar is outside the representable range".to_string());
    }
    Ok(RicciCurvature { tensor, scalar })
}

/// Ricci scalar `R = g^{μν} R_{μν}` at `point`.
///
/// # Errors
///
/// Returns an error for an invalid or singular metric at a connection stencil
/// point, or an unrepresentable stencil or result.
pub fn ricci_scalar(
    g: impl Fn(&[f64]) -> Vec<Vec<f64>>,
    point: &[f64],
) -> Result<f64, String> {
    ricci_curvature(g, point).map(|result| result.scalar)
}
