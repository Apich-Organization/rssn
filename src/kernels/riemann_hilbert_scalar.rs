//! Status-preserving quadrature for the normalized scalar Riemann–Hilbert arc.

use std::f64::consts::FRAC_1_PI;
use std::f64::consts::FRAC_PI_2;
use std::f64::consts::TAU;

use num_complex::Complex64;

use super::complex_quotient;
use crate::kernels::contour_deform::ContourIntegral;
use crate::kernels::contour_deform::ContourQuadratureError;
use crate::kernels::contour_deform::ContourSegment;
use crate::kernels::contour_deform::integrate_contour_checked;
use crate::kernels::quadrature::GkError;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);
const INVALID: Complex64 = Complex64::new(f64::NAN, f64::NAN);
const MAX_PANELS: usize = 1024;
pub(super) const DEFAULT_PANELS: usize = 300;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalarRhError {
    Quad(ContourQuadratureError),
    OffArcRequired,
    OpenArcRequired,
    PhaseResolution,
    NonFinite,
}

impl From<ContourQuadratureError> for ScalarRhError {
    fn from(error: ContourQuadratureError) -> Self {
        Self::Quad(error)
    }
}

impl From<GkError> for ScalarRhError {
    fn from(error: GkError) -> Self {
        Self::Quad(error.into())
    }
}

impl std::fmt::Display for ScalarRhError {
    fn fmt(
        &self,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            | Self::Quad(error) => error.fmt(f),
            | Self::OffArcRequired => {
                f.write_str("Riemann–Hilbert evaluation requires a finite point off [-1, 1]")
            },
            | Self::OpenArcRequired => {
                f.write_str("boundary values require x strictly inside (-1, 1)")
            },
            | Self::PhaseResolution => {
                f.write_str(
                    "the jump must be nonzero with phase resolved by the supplied sampling grid",
                )
            },
            | Self::NonFinite => {
                f.write_str("Riemann–Hilbert integrand or estimate is outside binary64 range")
            },
        }
    }
}

pub struct BoundaryEstimate {
    pub upper: Complex64,
    pub lower: Complex64,
    pub error: f64,
    pub evaluations: usize,
    pub converged: bool,
}

pub struct FactorEstimate {
    pub estimate: ContourIntegral,
    pub logarithm: Complex64,
    pub log_error: f64,
}

fn validate_options(
    tolerance: f64,
    panels: usize,
    initial_panels: usize,
) -> Result<(), ScalarRhError> {
    if !tolerance.is_finite() || 0.5 * tolerance <= 0.0 {
        return Err(GkError::InvalidTolerance.into());
    }
    if panels < initial_panels || panels > MAX_PANELS {
        return Err(GkError::PanelBudget.into());
    }
    if crate::cancellation::is_cancelled() {
        return Err(GkError::Cancelled.into());
    }
    Ok(())
}

fn validate_off_arc(z: Complex64) -> Result<(), ScalarRhError> {
    if z.is_finite() && (z.im != 0.0 || !(-1.0..=1.0).contains(&z.re)) {
        Ok(())
    } else {
        Err(ScalarRhError::OffArcRequired)
    }
}

const fn line(
    a: f64,
    b: f64,
) -> ContourSegment {
    ContourSegment::Line {
        start: Complex64::new(a, 0.0),
        end: Complex64::new(b, 0.0),
    }
}

fn cauchy_weighted(
    numerator: Complex64,
    reference: Complex64,
    denominator: Complex64,
) -> Complex64 {
    // Apply 1/(2*pi*i) before division so the unscaled integrand need not be representable.
    let unit = FRAC_1_PI * 0.5;
    let scaled = (numerator * unit - reference * unit) * -I;
    complex_quotient(scaled, denominator)
}

fn cauchy_with_reference<F: Fn(f64) -> Complex64>(
    phi: &F,
    z: Complex64,
    tolerance: f64,
    panels: usize,
    cuts: &[ContourSegment],
    reference: Option<Complex64>,
) -> Result<ContourIntegral, ScalarRhError> {
    let at_x = reference.unwrap_or(ZERO);
    let mut result = integrate_contour_checked(
        |t| cauchy_weighted(phi(t.re), at_x, t - z),
        cuts,
        tolerance,
        panels,
    )?;
    if reference.is_some() {
        let log = (Complex64::new(1.0, 0.0) - z).ln() - (Complex64::new(-1.0, 0.0) - z).ln();
        result.value += (at_x * (FRAC_1_PI * 0.5)) * (-I * log);
        result.evaluations += 1;
    }
    if result.value.is_finite() {
        Ok(result)
    } else {
        Err(ScalarRhError::NonFinite)
    }
}

pub fn cauchy_transform_arc_checked<F: Fn(f64) -> Complex64>(
    phi: &F,
    z: Complex64,
    tolerance: f64,
    panels: usize,
) -> Result<ContourIntegral, ScalarRhError> {
    validate_off_arc(z)?;
    let near = z.re.abs() < 1.0 && z.im.abs() < 0.25;
    validate_options(tolerance, panels, if near { 2 } else { 1 })?;
    if near {
        let value = phi(z.re);
        if !value.is_finite() {
            return Err(ScalarRhError::NonFinite);
        }
        cauchy_with_reference(
            phi,
            z,
            tolerance,
            panels,
            &[line(-1.0, z.re), line(z.re, 1.0)],
            Some(value),
        )
    } else {
        cauchy_with_reference(phi, z, tolerance, panels, &[line(-1.0, 1.0)], None)
    }
}

pub fn plemelj_boundary_values_checked<F: Fn(f64) -> Complex64>(
    phi: &F,
    x: f64,
    tolerance: f64,
    panels: usize,
) -> Result<BoundaryEstimate, ScalarRhError> {
    if !(-1.0 < x && x < 1.0) {
        return Err(ScalarRhError::OpenArcRequired);
    }
    validate_options(tolerance, panels, 2)?;
    let value = phi(x);
    if !value.is_finite() {
        return Err(ScalarRhError::NonFinite);
    }
    let pv = integrate_contour_checked(
        |t| {
            let distance = t.re - x;
            if distance == 0.0 {
                INVALID
            } else {
                cauchy_weighted(phi(t.re), value, Complex64::new(distance, 0.0))
            }
        },
        &[line(-1.0, x), line(x, 1.0)],
        tolerance,
        panels,
    )?;
    let logarithm = (-x).ln_1p() - x.ln_1p();
    let average = pv.value + (value * (FRAC_1_PI * 0.5)) * (-I * logarithm);
    let upper = average + value * 0.5;
    let lower = average - value * 0.5;
    if !upper.is_finite() || !lower.is_finite() {
        return Err(ScalarRhError::NonFinite);
    }
    Ok(BoundaryEstimate {
        upper,
        lower,
        error: pv.error,
        evaluations: pv.evaluations + 1,
        converged: pv.converged,
    })
}

fn log_jump_value(value: Complex64) -> Option<Complex64> {
    let scale = value.re.abs().max(value.im.abs());
    if !value.is_finite() || scale == 0.0 {
        return None;
    }
    Some(Complex64::new(
        scale.ln() + (value / scale).norm().ln(),
        value.arg(),
    ))
}

fn arc_node(
    index: usize,
    intervals: usize,
) -> f64 {
    -1.0 + 2.0 * (index as f64) / (intervals as f64)
}

fn lift_phase(
    phase: f64,
    reference: f64,
) -> Option<f64> {
    let lifted = phase + TAU * ((reference - phase) / TAU).round();
    ((lifted - reference).abs() <= FRAC_PI_2).then_some(lifted)
}

pub fn scalar_factorization_checked<F: Fn(f64) -> Complex64>(
    jump: &F,
    z: Complex64,
    tolerance: f64,
    phase_intervals: usize,
    panels: usize,
) -> Result<FactorEstimate, ScalarRhError> {
    validate_off_arc(z)?;
    if !(8..=256).contains(&phase_intervals) {
        return Err(ScalarRhError::PhaseResolution);
    }
    validate_options(tolerance, panels, phase_intervals)?;
    let mut phase = Vec::with_capacity(phase_intervals + 1);
    for i in 0..=phase_intervals {
        if crate::cancellation::is_cancelled() {
            return Err(GkError::Cancelled.into());
        }
        let arg = log_jump_value(jump(arc_node(i, phase_intervals)))
            .ok_or(ScalarRhError::PhaseResolution)?
            .im;
        let lifted = match phase.last() {
            | Some(&last) => lift_phase(arg, last).ok_or(ScalarRhError::PhaseResolution)?,
            | None => arg,
        };
        phase.push(lifted);
    }
    let unresolved = std::cell::Cell::new(false);
    let log_jump = |t: f64| {
        let result = log_jump_value(jump(t)).and_then(|value| {
            let u = f64::midpoint(t, 1.0) * (phase_intervals as f64);
            let index = (u.floor() as usize).min(phase_intervals - 1);
            let fraction = u - (index as f64);
            let reference = phase[index] * (1.0 - fraction) + phase[index + 1] * fraction;
            lift_phase(value.im, reference).map(|angle| Complex64::new(value.re, angle))
        });
        result.unwrap_or_else(|| {
            unresolved.set(true);
            INVALID
        })
    };
    let cuts: Vec<_> = (0..phase_intervals)
        .map(|i| {
            line(
                arc_node(i, phase_intervals),
                arc_node(i + 1, phase_intervals),
            )
        })
        .collect();
    let local = if z.re.abs() < 1.0 && z.im.abs() < 0.25 {
        Some(log_jump(z.re))
    } else {
        None
    };
    if unresolved.get() {
        return Err(ScalarRhError::PhaseResolution);
    }
    let integral = cauchy_with_reference(&log_jump, z, tolerance.min(0.25), panels, &cuts, local);
    if unresolved.get() {
        return Err(ScalarRhError::PhaseResolution);
    }
    let mut estimate = integral?;
    estimate.evaluations += phase.len();
    let logarithm = estimate.value;
    let log_error = estimate.error;
    estimate.value = logarithm.exp();
    estimate.error = if log_error == 0.0 {
        0.0
    } else {
        let log_expm1 = if log_error < 1.0 {
            log_error.exp_m1().ln()
        } else {
            log_error + (-(-log_error).exp()).ln_1p()
        };
        (logarithm.re + log_expm1).exp()
    };
    if !estimate.value.is_finite() || !estimate.error.is_finite() {
        return Err(ScalarRhError::NonFinite);
    }
    estimate.converged &= estimate.error <= tolerance;
    Ok(FactorEstimate {
        estimate,
        logarithm,
        log_error,
    })
}
