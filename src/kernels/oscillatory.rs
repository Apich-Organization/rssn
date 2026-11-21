//! Highly oscillatory integrals by Levin collocation.
//!
//! `I = integral from a to b of f(x) exp(i w g(x)) dx` with a large frequency
//! `w` defeats every rule that samples the integrand: the number of points
//! must grow with `w`. Levin's observation is that if `p` solves the
//! *non-oscillatory* differential equation `p' + i w g' p = f`, then
//! `(p exp(i w g))' = f exp(i w g)` and
//!
//! ```text
//! I = p(b) exp(i w g(b)) - p(a) exp(i w g(a)).
//! ```
//!
//! [`levin`] finds `p` by collocation at Chebyshev points, where the slowly
//! varying particular solution is well approximated by a polynomial of low
//! degree *whatever* `w` is; the cost is one small dense solve, and the
//! accuracy improves as `w` grows. The oscillatory rule needs `g'` to keep
//! one sign on the interval. [`levin_adaptive`] compares successive rules
//! and may fail to converge near a stationary point.
//!
//! [`LevinPlan`] prepares endpoint weights for a fixed interval and phase.
//! A new amplitude then takes one weighted pass over the nodes. Slowly
//! varying phases use Gauss–Legendre quadrature on the same interval.

#![allow(
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    clippy::cast_precision_loss,
    clippy::many_single_char_names,
    clippy::module_name_repetitions,
    clippy::single_call_fn,
    clippy::missing_docs_in_private_items,
    clippy::needless_range_loop,
    reason = "small dense collocation matrices indexed within their own extents"
)]

use crate::kernels::dense;
use crate::kernels::dense::Mat;
use crate::kernels::numeric_common::scale_product;
use crate::kernels::quadrature::gauss_legendre_rule;
use num_complex::Complex64;

/// A complex integral with an error estimate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Oscillatory {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
    /// An estimate of the absolute error (the difference to the next coarser
    /// collocation); `NAN` from [`levin`], which makes no estimate.
    pub error: f64,
}

fn chebyshev_points(n: usize) -> Vec<f64> {
    let top = (n - 1) as f64;
    (0..n)
        .map(|j| (std::f64::consts::PI * j as f64 / top).cos())
        .collect()
}

// Trefethen's differentiation matrix on [-1, 1], with the negative
// sum of the off-diagonal entries on each diagonal.
fn chebyshev_differentiation(
    x: &[f64],
    mut record: impl FnMut(usize, usize, f64),
) {
    let n = x.len();
    let c = |j: usize| {
        if j == 0 || j == n - 1 {
            2.0
        } else {
            1.0
        }
    };
    for i in 0..n {
        let mut sum = 0.0;
        for j in 0..n {
            if i != j {
                let sign = if (i + j) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                let value = c(i) / c(j) * sign / (x[i] - x[j]);
                record(i, j, value);
                sum += value;
            }
        }
        record(i, i, -sum);
    }
}

#[derive(Clone, Copy, Debug)]
struct Interval {
    start: f64,
    end: f64,
    extent: f64,
    half_factor: f64,
}

impl Interval {
    fn new(
        start: f64,
        end: f64,
    ) -> Self {
        let width = (end - start).abs();
        let (extent, half_factor) = if width.is_finite() {
            (width, 0.5)
        } else {
            ((0.5 * end - 0.5 * start).abs(), 1.0)
        };
        Self {
            start,
            end,
            extent,
            half_factor,
        }
    }

    fn at(
        self,
        node: f64,
    ) -> f64 {
        self.start * f64::midpoint(1.0, -node) + self.end * f64::midpoint(1.0, node)
    }

    fn orientation(self) -> f64 {
        if self.end < self.start {
            -1.0
        } else {
            1.0
        }
    }

    fn slope(
        self,
        derivative: f64,
        frequency: f64,
    ) -> Result<f64, String> {
        if !derivative.is_finite() {
            return Err("the phase derivative must be finite on the interval".into());
        }
        let product = frequency * derivative;
        let candidate = product * (self.extent * self.half_factor);
        let magnitude = if product.is_normal() && candidate.is_normal() {
            candidate
        } else {
            scale_product(
                derivative,
                &[frequency.abs(), self.extent, self.half_factor],
                &[],
            )
            .map_err(|_| "the scaled phase derivative is not representable")?
                * frequency.signum()
        };
        Ok(magnitude * self.orientation())
    }
}

/// Prepared quadrature for repeated amplitudes with one interval and phase.
///
/// Preparation solves the transposed Levin system for its two endpoint
/// functionals. Integration then applies the resulting weights in linear
/// time. Slowly varying phases, including zero frequency, use a Gauss rule
/// of the same size. Reversing the interval reverses the integral.
#[derive(Debug)]
pub struct LevinPlan {
    interval: Interval,
    nodes: Vec<f64>,
    weights: Vec<Complex64>,
    weight_scale: f64,
}

impl LevinPlan {
    /// Prepares `∫ f(x) exp(i w g(x)) dx` on the interval from `a` to `b`.
    ///
    /// `dg` is the phase derivative. Highly oscillatory phases should have
    /// no stationary points on the interval. Empty intervals and zero
    /// frequency do not require phase or derivative evaluations.
    ///
    /// # Errors
    /// Fewer than two points, non-finite parameters or phase samples, an
    /// unrepresentable phase angle or slope, a failed collocation solve,
    /// or cancellation by the enclosing execution scope.
    #[allow(
        clippy::float_cmp,
        reason = "An exactly empty oriented interval has zero integral."
    )]
    pub fn new(
        g: &dyn Fn(f64) -> f64,
        dg: &dyn Fn(f64) -> f64,
        a: f64,
        b: f64,
        w: f64,
        n: usize,
    ) -> Result<Self, String> {
        if n < 2 {
            return Err("Levin collocation needs at least two points".into());
        }
        if !(a.is_finite() && b.is_finite() && w.is_finite()) {
            return Err("the interval and the frequency must be finite".into());
        }
        let interval = Interval::new(a, b);
        if a == b {
            return Ok(Self {
                interval,
                nodes: Vec::new(),
                weights: Vec::new(),
                weight_scale: 1.0,
            });
        }
        let dimension = n.checked_mul(2).ok_or("Levin system size overflows")?;
        let size = dimension
            .checked_mul(dimension)
            .filter(|&size| size <= isize::MAX as usize / size_of::<f64>())
            .ok_or("Levin system size overflows")?;
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let phase = |x| -> Result<Complex64, String> {
            if w == 0.0 {
                return Ok(Complex64::ONE);
            }
            let angle = w * g(x);
            if !angle.is_finite() {
                return Err("the phase angle is not finite on the interval".into());
            }
            Ok(Complex64::from_polar(1.0, angle))
        };
        let mut collocation = None;
        if w != 0.0 {
            let points = chebyshev_points(n);
            let mut nodes: Vec<_> = points.iter().map(|&point| interval.at(point)).collect();
            nodes[0] = b;
            nodes[n - 1] = a;
            let slopes = nodes
                .iter()
                .map(|&x| interval.slope(dg(x), w))
                .collect::<Result<Vec<_>, _>>()?;
            let max_slope = slopes
                .iter()
                .fold(0.0_f64, |max, slope| max.max(slope.abs()));
            // With 2n polynomial moments, the Gauss rule resolves a phase
            // sufficiently small relative to the number of nodes.
            let direct_limit = 0.5 * n as f64 * f64::EPSILON.powf(0.5 / n as f64);
            if max_slope > direct_limit {
                collocation = Some((points, nodes, slopes));
            }
        }
        let (nodes, mut weights) = if let Some((points, nodes, slopes)) = collocation {
            let weights = adjoint_weights(&points, &slopes, phase(a)?, phase(b)?, size)?;
            (nodes, weights)
        } else {
            let (nodes, weights) = gauss_legendre_rule(n);
            let nodes: Vec<_> = nodes.into_iter().map(|x| interval.at(x)).collect();
            let weights = nodes
                .iter()
                .zip(weights)
                .map(|(&x, weight)| Ok(phase(x)? * weight))
                .collect::<Result<Vec<_>, String>>()?;
            (nodes, weights)
        };
        let weight_scale = weights.iter().fold(0.0_f64, |max, weight| {
            max.max(weight.re.abs()).max(weight.im.abs())
        });
        if weight_scale > 0.0 {
            for weight in &mut weights {
                *weight /= weight_scale;
            }
        }
        Ok(Self {
            interval,
            nodes,
            weights,
            weight_scale,
        })
    }

    /// The physical nodes, in the order required by [`Self::integrate_samples`].
    #[must_use]
    pub fn nodes(&self) -> &[f64] {
        &self.nodes
    }

    pub(crate) const fn heap_bytes(&self) -> usize {
        self.nodes.capacity() * size_of::<f64>() + self.weights.capacity() * size_of::<Complex64>()
    }

    /// Integrates an amplitude sampled once at each prepared node.
    ///
    /// # Errors
    /// A non-finite sample, an unrepresentable result, or cancellation.
    pub fn integrate(
        &self,
        f: &dyn Fn(f64) -> f64,
    ) -> Result<Oscillatory, String> {
        self.integrate_values(self.nodes.iter().map(|&x| f(x)))
    }

    /// Integrates amplitude values already sampled at [`Self::nodes`].
    ///
    /// # Errors
    /// The wrong sample count, a non-finite sample, an unrepresentable result,
    /// or cancellation.
    pub fn integrate_samples(
        &self,
        samples: &[f64],
    ) -> Result<Oscillatory, String> {
        if samples.len() != self.nodes.len() {
            return Err("one amplitude sample is required at each Levin node".into());
        }
        self.integrate_values(samples.iter().copied())
    }

    fn integrate_values(
        &self,
        samples: impl Iterator<Item = f64>,
    ) -> Result<Oscillatory, String> {
        if !self.nodes.is_empty() {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        }
        let mut sum = Complex64::ZERO;
        let mut amplitude_scale = 0.0_f64;
        for (sample, &weight) in samples.zip(&self.weights) {
            if !sample.is_finite() {
                return Err("the integrand is not finite at a Levin node".into());
            }
            let magnitude = sample.abs();
            if magnitude > amplitude_scale {
                sum *= amplitude_scale / magnitude;
                amplitude_scale = magnitude;
            }
            if amplitude_scale > 0.0 {
                sum += weight * (sample / amplitude_scale);
            }
        }
        let restore = |value| -> Result<f64, String> {
            if amplitude_scale == 0.0 || self.weight_scale == 0.0 {
                return Ok(0.0);
            }
            scale_product(
                value,
                &[
                    amplitude_scale,
                    self.weight_scale,
                    self.interval.extent,
                    self.interval.half_factor,
                ],
                &[],
            )
            .map(|part| part * self.interval.orientation())
            .map_err(|_| "the oscillatory integral is not representable".into())
        };
        Ok(Oscillatory {
            re: restore(sum.re)?,
            im: restore(sum.im)?,
            error: f64::NAN,
        })
    }
}

// For A p = [f, 0] and endpoint functional l, solve Aᵀ w = l
// once and evaluate lᵀ p as wᵀ [f, 0] for every new amplitude.
fn adjoint_weights(
    points: &[f64],
    slopes: &[f64],
    phase_a: Complex64,
    phase_b: Complex64,
    size: usize,
) -> Result<Vec<Complex64>, String> {
    let n = points.len();
    let dimension = 2 * n;
    let mut data = Vec::new();
    data.try_reserve_exact(size)
        .map_err(|_| "Levin matrix allocation failed")?;
    data.resize(size, 0.0);
    let mut transpose = Mat {
        rows: dimension,
        cols: dimension,
        data,
    };
    chebyshev_differentiation(points, |i, j, d| {
        transpose.set(j, i, d);
        transpose.set(n + j, n + i, d);
    });
    for (i, &slope) in slopes.iter().enumerate() {
        transpose.set(n + i, i, -slope);
        transpose.set(i, n + i, slope);
    }
    let mut endpoint = vec![0.0; dimension];
    endpoint[0] = phase_b.re;
    endpoint[n - 1] = -phase_a.re;
    endpoint[n] = -phase_b.im;
    endpoint[dimension - 1] = phase_a.im;
    let weights = match dense::solve(&transpose, &endpoint) {
        | Ok(weights) if weights.iter().all(|weight| weight.is_finite()) => weights,
        | Err(dense::DenseError::Cancelled) => {
            return Err("collocation solve: computation cancelled".into());
        },
        | _ => {
            dense::lstsq_svd(&transpose, &endpoint)
                .map_err(|error| format!("collocation solve: {error}"))?
        },
    };
    if weights.iter().any(|weight| !weight.is_finite()) {
        return Err("the collocation solution is not finite".into());
    }
    // A quarter-turn maps the real endpoint functional to the imaginary one
    // and commutes with this complex-linear system, so one solve supplies both.
    let (real, rotated) = weights.split_at(n);
    Ok(real
        .iter()
        .zip(rotated)
        .map(|(&re, &im)| Complex64::new(re, -im))
        .collect())
}

/// `integral from a to b of f(x) exp(i w g(x)) dx` by Levin collocation at
/// `n` Chebyshev points (`n >= 2`); `dg` is the derivative of `g`.
///
/// # Errors
/// A message when `n < 2`, the data is not finite, or the collocation system
/// cannot be solved.
pub fn levin(
    f: &dyn Fn(f64) -> f64,
    g: &dyn Fn(f64) -> f64,
    dg: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
    w: f64,
    n: usize,
) -> Result<Oscillatory, String> {
    LevinPlan::new(g, dg, a, b, w, n)?.integrate(f)
}

/// [`levin`] with the number of points doubled until two successive results
/// agree to `tol` (absolute, or relative to the result when that is larger),
/// starting from 8 points and stopping at `max_points`.
///
/// # Errors
/// A message when no pair of successive collocations agrees: typically a
/// stationary point of the phase, or an integrand that is itself
/// oscillatory.
#[allow(clippy::too_many_arguments)] // flat signature mirroring `levin` plus the two stopping controls
pub fn levin_adaptive(
    f: &dyn Fn(f64) -> f64,
    g: &dyn Fn(f64) -> f64,
    dg: &dyn Fn(f64) -> f64,
    a: f64,
    b: f64,
    w: f64,
    tol: f64,
    max_points: usize,
) -> Result<Oscillatory, String> {
    let mut n = 8;
    let mut previous = levin(f, g, dg, a, b, w, n)?;
    let mut best = f64::INFINITY;
    while n < max_points {
        n = (2 * n).min(max_points.max(9));
        let next = levin(f, g, dg, a, b, w, n)?;
        let error = (next.re - previous.re).hypot(next.im - previous.im);
        let scale = next.re.hypot(next.im).max(1.0);
        best = best.min(error);
        if error <= tol * scale {
            return Ok(Oscillatory { error, ..next });
        }
        previous = next;
    }
    Err(format!(
        "Levin collocation did not converge with {max_points} points (best agreement {best:e}); \
         the phase may have a stationary point on the interval"
    ))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

    use super::*;

    #[test]
    fn the_differentiation_matrix_is_exact_on_polynomials() {
        let x = chebyshev_points(6);
        let mut d = [[0.0; 6]; 6];
        chebyshev_differentiation(&x, |i, j, value| d[i][j] = value);
        // d/dx of x^3 at the nodes.
        for i in 0..6 {
            let got: f64 = (0..6).map(|j| d[i][j] * x[j].powi(3)).sum();
            assert!((got - 3.0 * x[i] * x[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn linear_phase_matches_the_closed_forms() {
        for w in [10.0, 1.0e3, 1.0e6] {
            // integral of exp(i w x) over [0, 1] = (exp(i w) - 1) / (i w).
            let r = levin(&|_| 1.0, &|x| x, &|_| 1.0, 0.0, 1.0, w, 8).unwrap();
            let (s, c) = w.sin_cos();
            assert!(
                (r.re - s / w).abs() < 1e-12 && (r.im - (1.0 - c) / w).abs() < 1e-12,
                "w = {w}"
            );
            // integral of x exp(i w x) = exp(i w) (1/(i w) + 1/w^2) - 1/w^2.
            let r = levin(&|x| x, &|x| x, &|_| 1.0, 0.0, 1.0, w, 8).unwrap();
            let re = s / w + (c - 1.0) / (w * w);
            let im = -c / w + s / (w * w);
            assert!(
                (r.re - re).abs() < 1e-12 && (r.im - im).abs() < 1e-12,
                "w = {w}"
            );
        }
    }

    #[test]
    fn a_nonlinear_phase_costs_the_same_at_any_frequency() {
        // integral of 2 x exp(i w x^2) over [1, 2] = (exp(4 i w) - exp(i w)) / (i w).
        for w in [50.0, 5.0e4] {
            let r = levin_adaptive(
                &|x| 2.0 * x,
                &|x| x * x,
                &|x| 2.0 * x,
                1.0,
                2.0,
                w,
                1e-11,
                64,
            )
            .unwrap();
            let re = ((4.0 * w).sin() - w.sin()) / w;
            let im = (w.cos() - (4.0 * w).cos()) / w;
            assert!(
                (r.re - re).abs() < 1e-11 && (r.im - im).abs() < 1e-11,
                "w = {w}: {r:?}"
            );
        }
        // A smooth amplitude: compare with a brute-force Simpson rule.
        let w = 200.0;
        let f = |x: f64| 1.0 / (1.0 + x);
        let r = levin_adaptive(
            &f,
            &|x| x + 0.5 * x * x,
            &|x| 1.0 + x,
            0.0,
            1.0,
            w,
            1e-10,
            128,
        )
        .unwrap();
        let panels = 400_000;
        let h = 1.0 / f64::from(panels);
        let (mut re, mut im) = (0.0, 0.0);
        for k in 0..=panels {
            let x = f64::from(k) * h;
            let weight = if k == 0 || k == panels {
                1.0
            } else if k % 2 == 1 {
                4.0
            } else {
                2.0
            };
            let phase = w * (x + 0.5 * x * x);
            re += weight * f(x) * phase.cos();
            im += weight * f(x) * phase.sin();
        }
        assert!((r.re - re * h / 3.0).abs() < 1e-9 && (r.im - im * h / 3.0).abs() < 1e-9);
    }

    #[test]
    fn failures_are_reported() {
        // A stationary point at 0: collocation does not settle.
        assert!(
            levin_adaptive(
                &|_| 1.0,
                &|x| x * x,
                &|x| 2.0 * x,
                -1.0,
                1.0,
                400.0,
                1e-10,
                64
            )
            .is_err()
        );
        assert!(levin(&|_| 1.0, &|x| x, &|_| 1.0, 0.0, 1.0, 1.0, 1).is_err());
        assert!(levin(&|_| f64::NAN, &|x| x, &|_| 1.0, 0.0, 1.0, 1.0, 8).is_err());
        assert!(levin(&|_| 1.0, &|x| x, &|_| 1.0, 0.0, f64::INFINITY, 1.0, 8).is_err());
        let empty = levin(&|_| 1.0, &|x| x, &|_| 1.0, 2.0, 2.0, 5.0, 8).unwrap();
        assert_eq!((empty.re, empty.im), (0.0, 0.0));
    }
}
