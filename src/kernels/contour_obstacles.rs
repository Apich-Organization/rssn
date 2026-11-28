//! Conservative Float validation for pole bypasses and branch-ray detours.

use num_complex::Complex64;
use std::f64::consts::PI;
use std::f64::consts::TAU;

use super::BranchRay;
use super::ContourSegment;

const GEOMETRY_GUARD: f64 = 64.0 * f64::EPSILON;
const ANGLE_GUARD: f64 = 256.0 * f64::EPSILON;
/// Maximum total number of poles and branch rays.
pub const MAX_OBSTACLES: usize = 64;

/// Failure to build and verify a bounded contour deformation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DeformationError {
    /// Input or constructed coordinates cannot represent finite geometry.
    InvalidGeometry,
    /// More obstacles were supplied than the checked construction allows.
    WorkLimit,
    /// The computed path cannot be separated from a pole by the requested distance.
    PoleClearance,
    /// The computed path cannot be separated from a branch ray.
    BranchRayCrossing,
    /// The caller cancelled construction.
    Cancelled,
}

impl std::fmt::Display for DeformationError {
    fn fmt(
        &self,
        output: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        output.write_str(match self {
            | Self::InvalidGeometry => {
                "contour deformation needs finite coordinates and positive representable clearance"
            },
            | Self::WorkLimit => "contour deformation exceeds the obstacle budget",
            | Self::PoleClearance => "contour detour does not have the requested pole clearance",
            | Self::BranchRayCrossing => "contour may cross or touch a branch ray",
            | Self::Cancelled => "contour deformation cancelled",
        })
    }
}

impl std::error::Error for DeformationError {}

const fn extent(z: Complex64) -> f64 {
    z.re.abs().max(z.im.abs())
}

fn vertex(
    center: Complex64,
    radius: f64,
    angle: f64,
) -> Complex64 {
    center + Complex64::from_polar(radius, angle)
}

fn dot(
    a: Complex64,
    b: Complex64,
) -> f64 {
    a.re * b.re + a.im * b.im
}

fn on_sweep(
    theta0: f64,
    theta1: f64,
    angle: f64,
    margin: f64,
) -> bool {
    let sweep = theta1 - theta0;
    let offset = ((angle - theta0) * sweep.signum()).rem_euclid(TAU);
    offset <= sweep.abs() + margin || TAU - offset <= margin
}

#[derive(Clone, Copy)]
struct PreparedRay {
    origin: Complex64,
    unit: Complex64,
    angle: f64,
}

impl From<BranchRay> for PreparedRay {
    fn from(ray: BranchRay) -> Self {
        let unit = Complex64::from_polar(1.0, ray.angle);
        Self {
            origin: ray.origin,
            unit,
            angle: unit.arg(),
        }
    }
}

fn line_clears_ray(
    a: Complex64,
    b: Complex64,
    ray: PreparedRay,
) -> bool {
    let absolute_scale = extent(a).max(extent(b)).max(extent(ray.origin));
    let a = a - ray.origin;
    let b = b - ray.origin;
    if !a.is_finite() || !b.is_finite() {
        return false;
    }
    let scale = extent(a).max(extent(b));
    if scale == 0.0 {
        return false;
    }
    let rotation = ray.unit.conj();
    let a = (a / scale) * rotation;
    let b = (b / scale) * rotation;
    let margin = (absolute_scale / scale).max(1.0) * GEOMETRY_GUARD;
    if a.re.max(b.re) < -margin || a.im.min(b.im) > margin || a.im.max(b.im) < -margin {
        return true;
    }
    if a.im.abs().min(b.im.abs()) <= margin {
        return false;
    }
    // Opposite signed endpoints. Bound the interpolation parameter using
    // the uncertainty on both imaginary coordinates, without a shallow-angle
    // division by a determinant.
    let a_height = a.im.abs();
    let b_height = b.im.abs();
    let total_height = a_height + b_height;
    let lo = (a_height - margin) / total_height;
    let hi = (a_height + margin) / total_height;
    let interp = |t: f64| b.re.mul_add(t, a.re * (1.0 - t));
    interp(lo).max(interp(hi)) < -margin
}

fn angular_interval_misses_sweep(
    theta0: f64,
    theta1: f64,
    lo: f64,
    hi: f64,
) -> bool {
    let start = theta0.min(theta1);
    let width = (theta1 - theta0).abs();
    let offset = (lo - start).rem_euclid(TAU);
    let end = offset + (hi - lo);
    offset > width + ANGLE_GUARD && end < TAU - ANGLE_GUARD
}

fn arc_clears_ray(
    center: Complex64,
    radius: f64,
    theta0: f64,
    theta1: f64,
    ray: PreparedRay,
) -> bool {
    let absolute_scale = extent(center).max(extent(ray.origin)).max(radius);
    let offset = center - ray.origin;
    if !offset.is_finite() {
        return false;
    }
    let scale = extent(offset).max(radius);
    let center = (offset / scale) * ray.unit.conj();
    let radius = radius / scale;
    let margin = (absolute_scale / scale).max(1.0) * GEOMETRY_GUARD;
    if center.im.abs() > radius + margin || center.re + radius < -margin {
        return true;
    }
    if radius <= margin {
        return false;
    }
    // At or near tangency the inverse sine is poorly conditioned. Check
    // the entire conservative root interval and both semicircles.
    let lo = ((-center.im - margin) / radius).clamp(-1.0, 1.0).asin();
    let hi = ((-center.im + margin) / radius).clamp(-1.0, 1.0).asin();
    [(lo, hi), (PI - hi, PI - lo)]
        .into_iter()
        .all(|(lower, upper)| {
            let max_cos = if lower <= 0.0 && upper >= 0.0 {
                1.0
            } else {
                lower.cos().max(upper.cos())
            };
            center.re + radius * max_cos < -margin
                || angular_interval_misses_sweep(
                    theta0,
                    theta1,
                    ray.angle + lower,
                    ray.angle + upper,
                )
        })
}

fn clears_ray(
    segment: ContourSegment,
    ray: PreparedRay,
) -> bool {
    match segment {
        | ContourSegment::Line { start, end } => line_clears_ray(start, end, ray),
        | ContourSegment::Arc {
            center,
            radius,
            theta0,
            theta1,
        } => arc_clears_ray(center, radius, theta0, theta1, ray),
    }
}

fn line_clears_pole(
    a: Complex64,
    b: Complex64,
    pole: Complex64,
    clearance: f64,
) -> bool {
    let absolute_scale = extent(a).max(extent(b)).max(extent(pole));
    let a = a - pole;
    let b = b - pole;
    if !a.is_finite() || !b.is_finite() {
        return false;
    }
    let scale = extent(a).max(extent(b)).max(clearance);
    let a = a / scale;
    let d = b / scale - a;
    let length = d.norm();
    let distance = if length > 0.0 {
        let unit = d / length;
        (a + unit * (-dot(a, unit)).clamp(0.0, length)).norm()
    } else {
        a.norm()
    };
    distance > clearance / scale + (absolute_scale / scale).max(1.0) * GEOMETRY_GUARD
}

fn arc_clears_pole(
    center: Complex64,
    radius: f64,
    theta0: f64,
    theta1: f64,
    pole: Complex64,
    clearance: f64,
) -> bool {
    let absolute_scale = extent(pole).max(extent(center)).max(radius);
    let offset = pole - center;
    if !offset.is_finite() {
        return false;
    }
    let scale = extent(offset).max(radius).max(clearance);
    let offset = offset / scale;
    let radius = radius / scale;
    let margin = (absolute_scale / scale).max(1.0) * GEOMETRY_GUARD;
    let angular_margin = ANGLE_GUARD + (2.0 * margin / offset.norm()).min(PI);
    let distance = if on_sweep(theta0, theta1, offset.arg(), angular_margin) {
        (offset.norm() - radius).abs()
    } else {
        (offset - Complex64::from_polar(radius, theta0))
            .norm()
            .min((offset - Complex64::from_polar(radius, theta1)).norm())
    };
    distance > clearance / scale + margin
}

fn clears_pole(
    segment: ContourSegment,
    pole: Complex64,
    clearance: f64,
) -> bool {
    match segment {
        | ContourSegment::Line { start, end } => line_clears_pole(start, end, pole, clearance),
        | ContourSegment::Arc {
            center,
            radius,
            theta0,
            theta1,
        } => arc_clears_pole(center, radius, theta0, theta1, pole, clearance),
    }
}

fn ray_detour(
    start: Complex64,
    end: Complex64,
    ray: PreparedRay,
    poles: &[Complex64],
    clearance: f64,
) -> Result<Vec<ContourSegment>, DeformationError> {
    let a = start - ray.origin;
    let b = end - ray.origin;
    let radius = a.norm().min(b.norm());
    if !radius.is_finite() || radius <= 0.0 {
        return Err(DeformationError::BranchRayCrossing);
    }
    // Keep both endpoint arguments within the same open branch sheet.
    // Connecting them on that sheet takes the long route when the chord
    // crosses the ray away from its origin.
    let theta0 = ray.angle + (a.arg() - ray.angle).rem_euclid(TAU);
    let theta1 = ray.angle + (b.arg() - ray.angle).rem_euclid(TAU);
    let start_arc = vertex(ray.origin, radius, theta0);
    let end_arc = vertex(ray.origin, radius, theta1);
    let mut path = pole_path(start, start_arc, poles, clearance)?;
    path.push(ContourSegment::Arc {
        center: ray.origin,
        radius,
        theta0,
        theta1,
    });
    path.extend(pole_path(end_arc, end, poles, clearance)?);
    Ok(path)
}

struct Bypass {
    start: f64,
    end: f64,
    pole: Complex64,
}

fn pole_path(
    start: Complex64,
    end: Complex64,
    poles: &[Complex64],
    clearance: f64,
) -> Result<Vec<ContourSegment>, DeformationError> {
    if poles.is_empty() {
        return Ok(vec![ContourSegment::Line { start, end }]);
    }
    let local_scale = extent(start).max(extent(end)).max(clearance);
    let radius = (8.0 * GEOMETRY_GUARD).mul_add(local_scale, clearance);
    if !radius.is_finite() {
        return Err(DeformationError::InvalidGeometry);
    }
    let delta = end - start;
    let length = delta.norm();
    if !length.is_finite() {
        return Err(DeformationError::InvalidGeometry);
    }
    if length == 0.0 {
        return Ok(vec![ContourSegment::Line { start, end }]);
    }
    let along = delta / length;
    let perp = Complex64::new(-along.im, along.re);
    let mut bypasses = Vec::new();
    for &pole in poles {
        let relative = pole - start;
        let normal = dot(relative, perp).abs();
        let center = dot(relative, along);
        if !normal.is_finite() || !center.is_finite() {
            return Err(DeformationError::InvalidGeometry);
        }
        if normal >= radius || center <= 0.0 || center >= length {
            continue;
        }
        let ratio = normal / radius;
        let offset = radius * ((1.0 - ratio) * (1.0 + ratio)).sqrt();
        let enter = center - offset;
        let leave = center + offset;
        if !(0.0 < enter && enter < leave && leave < length) {
            return Err(DeformationError::PoleClearance);
        }
        bypasses.push(Bypass {
            start: enter,
            end: leave,
            pole,
        });
    }
    bypasses.sort_by(|a, b| a.start.total_cmp(&b.start));
    let mut path = Vec::with_capacity(bypasses.len() * 2 + 1);
    let mut cursor = 0.0;
    let mut location = start;
    for bypass in bypasses {
        if bypass.start <= cursor {
            return Err(DeformationError::PoleClearance);
        }
        let entry = start + along * bypass.start;
        let exit = start + along * bypass.end;
        let theta0 = (entry - bypass.pole).arg();
        let forward = ((exit - bypass.pole).arg() - theta0).rem_euclid(TAU);
        let theta1 = theta0
            + if forward <= PI {
                forward
            } else {
                forward - TAU
            };
        path.push(ContourSegment::Line {
            start: location,
            end: vertex(bypass.pole, radius, theta0),
        });
        path.push(ContourSegment::Arc {
            center: bypass.pole,
            radius,
            theta0,
            theta1,
        });
        cursor = bypass.end;
        location = vertex(bypass.pole, radius, theta1);
    }
    path.push(ContourSegment::Line { start: location, end });
    Ok(path)
}

fn representable(segment: &ContourSegment) -> bool {
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
            radius.is_finite()
                && radius > 0.0
                && (theta1 - theta0).abs() > 0.0
                && [center.re, center.im].into_iter().all(|component| {
                    let lower = component - radius;
                    let upper = component + radius;
                    lower.is_finite() && upper.is_finite() && lower < component && component < upper
                })
        },
    }
}

fn verify_path(
    path: Vec<ContourSegment>,
    poles: &[Complex64],
    rays: &[PreparedRay],
    clearance: f64,
) -> Result<Vec<ContourSegment>, DeformationError> {
    for &segment in &path {
        if crate::cancellation::is_cancelled() {
            return Err(DeformationError::Cancelled);
        }
        if !representable(&segment) {
            return Err(DeformationError::InvalidGeometry);
        }
        if rays.iter().any(|&ray| !clears_ray(segment, ray)) {
            return Err(DeformationError::BranchRayCrossing);
        }
        if poles
            .iter()
            .any(|&pole| !clears_pole(segment, pole, clearance))
        {
            return Err(DeformationError::PoleClearance);
        }
    }
    Ok(path)
}

/// Constructs and verifies one open deformation in the supplied coordinate frame.
///
/// # Errors
/// Rejects unsupported geometry, work limits, and cancellation.
pub fn deformation_checked(
    start: Complex64,
    end: Complex64,
    poles: &[Complex64],
    cuts: &[BranchRay],
    clearance: f64,
) -> Result<Vec<ContourSegment>, DeformationError> {
    if poles.len().saturating_add(cuts.len()) > MAX_OBSTACLES {
        return Err(DeformationError::WorkLimit);
    }
    if !start.is_finite()
        || !end.is_finite()
        || !clearance.is_finite()
        || clearance <= 0.0
        || poles.iter().any(|pole| !pole.is_finite())
        || cuts
            .iter()
            .any(|ray| !ray.origin.is_finite() || !ray.angle.is_finite())
    {
        return Err(DeformationError::InvalidGeometry);
    }
    if crate::cancellation::is_cancelled() {
        return Err(DeformationError::Cancelled);
    }
    let rays: Vec<_> = cuts.iter().copied().map(PreparedRay::from).collect();
    let candidate = pole_path(start, end, poles, clearance)?;
    let error = match verify_path(candidate, poles, &rays, clearance) {
        | Ok(path) => return Ok(path),
        | Err(DeformationError::BranchRayCrossing) => DeformationError::BranchRayCrossing,
        | Err(error) => return Err(error),
    };
    for &ray in &rays {
        if line_clears_ray(start, end, ray) {
            continue;
        }
        let candidate = ray_detour(start, end, ray, poles, clearance)
            .and_then(|path| verify_path(path, poles, &rays, clearance));
        match candidate {
            | Ok(path) => return Ok(path),
            | Err(DeformationError::Cancelled) => return Err(DeformationError::Cancelled),
            | Err(_) => {},
        }
    }
    Err(error)
}
