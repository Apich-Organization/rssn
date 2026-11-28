//! Checked complex quadrature for oriented contour paths and symmetric folds.

use std::cell::Cell;
use std::cell::RefCell;
use std::collections::HashMap;
use std::f64::consts::TAU;

use num_complex::Complex64;

use super::ContourSegment;
use super::I;
use crate::kernels::quadrature::GkError;
use crate::kernels::quadrature::GkRule;
use crate::kernels::quadrature::gk_panel_complex;
use crate::kernels::quadrature::integrate_gk_breakpoints;

/// Estimate, termination and actual complex callback count for a contour integral.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContourIntegral {
    /// Estimated contour integral.
    pub value: Complex64,
    /// Combined absolute error estimate from the embedded quadrature rules.
    pub error: f64,
    /// Number of complex integrand calls, including both sides of any fold.
    pub evaluations: usize,
    /// Whether refinement met the requested error estimate in both components.
    pub converged: bool,
}

/// Failure to construct or evaluate a checked contour quadrature.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContourQuadratureError {
    /// A segment, interval, or pole has invalid or unrepresentable geometry.
    InvalidGeometry,
    /// The checked Gauss–Kronrod engine reported a failure.
    Quadrature(GkError),
}

impl From<GkError> for ContourQuadratureError {
    fn from(value: GkError) -> Self {
        Self::Quadrature(value)
    }
}

impl std::fmt::Display for ContourQuadratureError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::InvalidGeometry => {
                f.write_str("contour geometry is invalid or unrepresentable")
            },
            | Self::Quadrature(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for ContourQuadratureError {}

fn complex_partition(
    pieces: usize,
    sample: impl Fn(usize, f64) -> (Complex64, usize),
    tolerance: f64,
    max_panels: usize,
) -> Result<ContourIntegral, GkError> {
    const MAX_PANELS: usize = 4096;
    if pieces > max_panels || max_panels == 0 || max_panels > MAX_PANELS {
        return Err(GkError::PanelBudget);
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(GkError::InvalidTolerance);
    }
    if crate::cancellation::is_cancelled() {
        return Err(GkError::Cancelled);
    }
    let upper = pieces as f64;
    let breaks: Vec<f64> = (1..pieces).map(|i| i as f64).collect();
    let evaluations = Cell::new(0);
    let failed = Cell::new(false);
    let evaluate = |t: f64| {
        if failed.get() || crate::cancellation::is_cancelled() {
            return Complex64::new(f64::NAN, f64::NAN);
        }
        let piece = (t.floor() as usize).min(pieces.saturating_sub(1));
        let (value, calls) = sample(piece, t - piece as f64);
        evaluations.set(evaluations.get() + calls);
        if value.is_finite() {
            value
        } else {
            failed.set(true);
            Complex64::new(f64::NAN, f64::NAN)
        }
    };
    // Common single panels need neither a hash table nor separately scaled
    // adaptive partitions. Keep their samples if either component must refine.
    let samples = RefCell::new(HashMap::new());
    if pieces == 1 {
        let first = RefCell::new([(0_u64, Complex64::ZERO); 15]);
        let size = Cell::new(0_usize);
        let (value, [error_re, error_im], _) = gk_panel_complex(0.0, 1.0, GkRule::G7K15, |t| {
            let value = evaluate(t);
            first.borrow_mut()[size.get()] = (t.to_bits(), value);
            size.set(size.get() + 1);
            value
        });
        if crate::cancellation::is_cancelled() {
            return Err(GkError::Cancelled);
        }
        let error = error_re.hypot(error_im);
        if !value.is_finite() || !error.is_finite() {
            return Err(GkError::NonFiniteEstimate);
        }
        let converged = error_re <= tolerance * 0.5 && error_im <= tolerance * 0.5;
        if converged || max_panels == 1 {
            return Ok(ContourIntegral {
                value,
                error,
                evaluations: evaluations.get(),
                converged,
            });
        }
        samples
            .borrow_mut()
            .extend(first.into_inner().into_iter().take(size.get()));
    }
    // Each callback result can be reused by both component partitions.
    let at = |t: f64| {
        *samples
            .borrow_mut()
            .entry(t.to_bits())
            .or_insert_with(|| evaluate(t))
    };
    // Splitting the absolute tolerance also bounds the norm of the combined
    // error estimate. No relative floor changes the caller's requested bound.
    let integrate = |component: &dyn Fn(f64) -> f64| {
        integrate_gk_breakpoints(
            component,
            0.0,
            upper,
            &breaks,
            GkRule::G7K15,
            tolerance * 0.5,
            0.0,
            max_panels,
        )
    };
    // Each closure has a distinct type. Use a borrowed dynamic scalar
    // callback only at the common quadrature call boundary.
    let re = integrate(&|t| at(t).re)?;
    let im = integrate(&|t| at(t).im)?;
    let error = re.error.hypot(im.error);
    if !error.is_finite() {
        return Err(GkError::NonFiniteEstimate);
    }
    Ok(ContourIntegral {
        value: Complex64::new(re.value, im.value),
        error,
        evaluations: evaluations.get(),
        converged: re.converged && im.converged && error <= tolerance,
    })
}

fn valid_segment(segment: &ContourSegment) -> bool {
    match *segment {
        | ContourSegment::Line { start, end } => {
            start.is_finite() && end.is_finite() && (end - start).is_finite()
        },
        | ContourSegment::Arc {
            center,
            radius,
            theta0,
            theta1,
        } => {
            let sweep = theta1 - theta0;
            center.is_finite()
                && radius.is_finite()
                && radius > 0.0
                && theta0.is_finite()
                && theta1.is_finite()
                && sweep.is_finite()
                && (radius * sweep).is_finite()
                && [center.re, center.im]
                    .into_iter()
                    .all(|axis| (axis + radius).is_finite() && (axis - radius).is_finite())
        },
    }
}

fn active_segment(segment: ContourSegment) -> bool {
    match segment {
        | ContourSegment::Line { start, end } => start != end,
        | ContourSegment::Arc { theta0, theta1, .. } => (theta1 - theta0).abs() > 0.0,
    }
}

fn sample_segment(
    segment: ContourSegment,
    parameter: f64,
) -> (Complex64, Complex64) {
    match segment {
        | ContourSegment::Line { start, end } => {
            let tangent = end - start;
            (start + tangent * parameter, tangent)
        },
        | ContourSegment::Arc {
            center,
            radius,
            theta0,
            theta1,
        } => {
            let sweep = theta1 - theta0;
            let direction = Complex64::from_polar(1.0, theta0 + sweep * parameter);
            (
                center + radius * direction,
                I * (radius * sweep) * direction,
            )
        },
    }
}

/// Numerical estimate for one finite oriented segment, including termination status.
///
/// See [`integrate_contour_checked`] for the tolerance and panel semantics.
///
/// # Errors
/// Rejects malformed geometry or a failed quadrature evaluation.
pub fn integrate_segment_checked<F: Fn(Complex64) -> Complex64>(
    f: F,
    segment: ContourSegment,
    tolerance: f64,
    max_panels: usize,
) -> Result<ContourIntegral, ContourQuadratureError> {
    if max_panels == 0 {
        return Err(GkError::PanelBudget.into());
    }
    if !valid_segment(&segment) {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    let pieces = usize::from(active_segment(segment));
    Ok(complex_partition(
        pieces,
        |_, parameter| {
            let (point, tangent) = sample_segment(segment, parameter);
            (f(point) * tangent, 1)
        },
        tolerance,
        max_panels,
    )?)
}

/// Estimates the integral along a piecewise contour and preserves quadrature status.
///
/// `tolerance` is an absolute error estimate for the complex result. The
/// `max_panels` budget of 1 through 4096 applies to the whole path per scalar
/// component. Embedded quadrature error estimates are heuristic. A finite
/// result can be returned after budget exhaustion with `converged = false`.
///
/// # Errors
/// Rejects malformed segments, invalid budgets or tolerance, nonfinite
/// evaluation results, and cancellation.
pub fn integrate_contour_checked<F: Fn(Complex64) -> Complex64>(
    f: F,
    segments: &[ContourSegment],
    tolerance: f64,
    max_panels: usize,
) -> Result<ContourIntegral, ContourQuadratureError> {
    if let [segment] = segments {
        return integrate_segment_checked(f, *segment, tolerance, max_panels);
    }
    if segments.len() > max_panels {
        return Err(GkError::PanelBudget.into());
    }
    if !segments.iter().all(valid_segment) {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    let active: Vec<_> = segments
        .iter()
        .filter(|segment| active_segment(**segment))
        .copied()
        .collect();
    Ok(complex_partition(
        active.len(),
        |index, t| {
            let Some(segment) = active.get(index) else {
                return (Complex64::new(f64::NAN, f64::NAN), 0);
            };
            let (z, jacobian) = sample_segment(*segment, t);
            (f(z) * jacobian, 1)
        },
        tolerance,
        max_panels,
    )?)
}

/// Estimates the residue at infinity using adaptive integration along a positive circle.
///
/// The caller supplies a positive finite radius enclosing the finite singularities.
/// The error estimate and tolerance are in residue units.
///
/// # Errors
/// Rejects invalid geometry, nonfinite estimates or cancellation.
pub fn residue_at_infinity_checked<F: Fn(Complex64) -> Complex64>(
    f: F,
    radius: f64,
    tolerance: f64,
    max_panels: usize,
) -> Result<ContourIntegral, ContourQuadratureError> {
    let mut result = integrate_segment_checked(
        f,
        ContourSegment::Arc {
            center: Complex64::new(0.0, 0.0),
            radius,
            theta0: 0.0,
            theta1: TAU,
        },
        tolerance * TAU,
        max_panels,
    )?;
    result.value *= I / TAU;
    result.error /= TAU;
    result.converged &= result.error <= tolerance;
    Ok(result)
}

#[derive(Clone, Copy)]
enum Segment {
    Plain(f64, f64),
    Fold { pole: f64, radius: f64 },
}

fn principal_value_pieces(
    a: f64,
    b: f64,
    poles: &[f64],
) -> Option<Vec<Segment>> {
    if !a.is_finite()
        || !b.is_finite()
        || !(b - a).is_finite()
        || poles.iter().any(|&p| !(a < p && p < b))
    {
        return None;
    }
    let mut poles = poles.to_vec();
    poles.sort_by(f64::total_cmp);
    let mut out = Vec::with_capacity(poles.len() * 2 + 1);
    let mut cursor = a;
    for (index, &pole) in poles.iter().enumerate() {
        let right_limit = poles
            .get(index + 1)
            .map_or(b, |&next| f64::midpoint(pole, next));
        let radius = (pole - cursor).min(right_limit - pole);
        if radius <= 0.0 || pole - radius >= pole || pole + radius <= pole {
            return None;
        }
        if cursor < pole - radius {
            out.push(Segment::Plain(cursor, pole - radius));
        }
        out.push(Segment::Fold { pole, radius });
        cursor = pole + radius;
    }
    if cursor < b {
        out.push(Segment::Plain(cursor, b));
    }
    Some(out)
}

/// Estimates a Cauchy principal value by pairing samples symmetrically around each simple pole.
///
/// Poles must be distinct and strictly inside the finite limits. Reversing the
/// limits negates the estimate; neither orientation certifies the simple-pole
/// assumption. `evaluations` includes both evaluations of each paired sample.
///
/// # Errors
/// Rejects invalid geometry, nonfinite estimates or cancellation.
pub fn cauchy_principal_value_checked<F: Fn(f64) -> Complex64>(
    f: F,
    a: f64,
    b: f64,
    poles: &[f64],
    tolerance: f64,
    max_panels: usize,
) -> Result<ContourIntegral, ContourQuadratureError> {
    if poles.len() > max_panels {
        return Err(GkError::PanelBudget.into());
    }
    if !a.is_finite() || !b.is_finite() {
        return Err(ContourQuadratureError::InvalidGeometry);
    }
    let pieces = principal_value_pieces(a.min(b), a.max(b), poles)
        .ok_or(ContourQuadratureError::InvalidGeometry)?;
    let mut result = complex_partition(
        pieces.len(),
        |index, t| {
            match pieces.get(index) {
                | Some(Segment::Plain(a, b)) => {
                    let width = b - a;
                    (f(a + width * t) * width, 1)
                },
                | Some(Segment::Fold { pole, radius }) => {
                    let delta = radius * t;
                    let left = f(pole - delta);
                    if !left.is_finite() {
                        return (left, 1);
                    }
                    ((left + f(pole + delta)) * radius, 2)
                },
                | None => (Complex64::new(f64::NAN, f64::NAN), 0),
            }
        },
        tolerance,
        max_panels,
    )?;
    if a > b {
        result.value = -result.value;
    }
    Ok(result)
}
