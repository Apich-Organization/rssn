//! Coordinate transformations between Cartesian, polar, cylindrical and spherical systems.
//!
//! # Coordinate transformations
//!
//! Conversions of plain `f64` points between Cartesian, polar, cylindrical
//! and spherical coordinates, and numerical Jacobians of those maps.
//!
//! Conventions:
//!
//! * Cartesian: `(x, y)` or `(x, y, z)`.
//! * Polar: `(r, theta)` with `theta` the angle from the x axis.
//! * Cylindrical: `(r, theta, z)` with `theta` the angle from the x axis.
//! * Spherical: `(r, theta, phi)` with `theta` the polar angle measured from
//!   the z axis (`0..=pi`) and `phi` the azimuth from the x axis.
//!
//! # Limitations and coverage
//!
//! * Angles are radians and the conventions above are the only ones: polar and
//!   cylindrical `theta` is measured from the x axis, spherical `theta` from
//!   the z axis. A caller holding the "math" `(r, phi, theta)` ordering must
//!   permute the point; no `theta`/`phi` swap option is provided.
//! * Conversion between non-Cartesian systems goes through Cartesian. A
//!   round trip can round angles, and 2-D and 3-D systems cannot be mixed.
//!   [`cartesian_to_spherical`] and [`spherical_to_cartesian`] accept and
//!   return three-element arrays without allocating. A negative spherical
//!   radius reflects the point through the origin. The origin maps to
//!   `theta = 0`, with its azimuth given by `atan2(y, x)` on the signed zeros.
//! * [`numerical_jacobian`] and [`jacobian_determinant`] use central
//!   differences with relative step `6e-6`, so they are accurate to roughly
//!   `1e-9` relative for smooth maps and are meaningless at the singular
//!   points of a map (`r = 0` for polar, the origin and the z axis for
//!   spherical). No analytic Jacobian, metric or Christoffel symbol is
//!   provided.

/// A coordinate system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateSystem {
    /// Cartesian `(x, y)` or `(x, y, z)`.
    Cartesian,
    /// Plane polar `(r, theta)`.
    Polar,
    /// Cylindrical `(r, theta, z)`.
    Cylindrical,
    /// Spherical `(r, theta, phi)`, `theta` measured from the z axis.
    Spherical,
}

/// Relative step for the central-difference Jacobian, near the optimum
/// `eps^(1/3)` of a second-order formula.
const JACOBIAN_STEP: f64 = 6.0e-6;
/// 2⁵⁴ brings even the smallest subnormal f64 into the normal range.
const SUBNORMAL_SCALE: f64 = 18_014_398_509_481_984.0;

fn component(
    point: &[f64],
    index: usize,
) -> f64 {
    point.get(index).copied().unwrap_or(0.0)
}

fn check_len(
    point: &[f64],
    system: CoordinateSystem,
) -> Result<(), String> {
    let ok = match system {
        | CoordinateSystem::Cartesian => point.len() == 2 || point.len() == 3,
        | CoordinateSystem::Polar => point.len() == 2,
        | CoordinateSystem::Cylindrical | CoordinateSystem::Spherical => point.len() == 3,
    };
    if ok {
        Ok(())
    } else {
        Err(format!(
            "{system:?} point cannot have {} components",
            point.len()
        ))
    }
}

/// Converts Cartesian `[x, y, z]` to spherical `[r, theta, phi]`.
///
/// `theta` is measured from the positive z axis in `[0, π]` and `phi` is
/// the azimuth from the positive x axis in `[-π, π]` (both in radians).
/// With horizontal length `h = hypot(x, y)`, the polar angle is
/// `atan2(h, z)`. It is zero at the origin. Finite coordinates can produce
/// an infinite radius if the norm overflows, with both angles still defined.
/// An infinite z component produces a NaN polar angle.
#[must_use]
pub fn cartesian_to_spherical([x, y, z]: [f64; 3]) -> [f64; 3] {
    let horizontal = x.hypot(y);
    let azimuth = y.atan2(x);

    if horizontal.is_infinite() && x.is_finite() && y.is_finite() && z.is_finite() {
        // The norm overflows, but taking the angle at half scale is safe.
        let half_horizontal = (0.5 * x).hypot(0.5 * y);
        return [f64::INFINITY, half_horizontal.atan2(0.5 * z), azimuth];
    }
    if horizontal.is_subnormal() && z.abs() <= f64::MAX / SUBNORMAL_SCALE {
        // Restore precision lost by rounding the horizontal norm to a subnormal.
        let horizontal = (x * SUBNORMAL_SCALE).hypot(y * SUBNORMAL_SCALE);
        let z = z * SUBNORMAL_SCALE;
        let radius = horizontal.hypot(z) / SUBNORMAL_SCALE;
        return [radius, horizontal.atan2(z), azimuth];
    }

    let radius = horizontal.hypot(z);
    let theta = if radius == 0.0 {
        0.0
    } else if z.is_infinite() {
        f64::NAN
    } else {
        horizontal.atan2(z)
    };
    [radius, theta, azimuth]
}

/// Converts spherical `[r, theta, phi]` in radians to Cartesian `[x, y, z]`.
///
/// Uses the angles in [`cartesian_to_spherical`]. A negative radius reflects
/// the Cartesian point through the origin. See [DLMF 1.5.16] for the map.
///
/// [DLMF 1.5.16]: https://dlmf.nist.gov/1.5.E16
#[must_use]
pub fn spherical_to_cartesian([r, theta, phi]: [f64; 3]) -> [f64; 3] {
    let (sin_theta, cos_theta) = theta.sin_cos();
    let (sin_phi, cos_phi) = phi.sin_cos();
    if r.is_subnormal() {
        // Form the unit direction before rounding at the scale of the radius.
        return [
            r * (sin_theta * cos_phi),
            r * (sin_theta * sin_phi),
            r * cos_theta,
        ];
    }
    let horizontal = r * sin_theta;
    [horizontal * cos_phi, horizontal * sin_phi, r * cos_theta]
}

/// Converts a point of `from` to Cartesian coordinates, without checking
/// its length: missing components are taken as zero.
#[must_use]
pub fn to_cartesian_pure(
    point: &[f64],
    from: CoordinateSystem,
) -> Vec<f64> {
    let (a, b, c) = (
        component(point, 0),
        component(point, 1),
        component(point, 2),
    );
    match from {
        | CoordinateSystem::Cartesian => point.to_vec(),
        | CoordinateSystem::Polar => vec![a * b.cos(), a * b.sin()],
        | CoordinateSystem::Cylindrical => vec![a * b.cos(), a * b.sin(), c],
        | CoordinateSystem::Spherical => spherical_to_cartesian([a, b, c]).to_vec(),
    }
}

/// Converts a Cartesian point to the system `to`, without checking its
/// length: missing components are taken as zero. The polar angle of the
/// origin in spherical coordinates is 0.
#[must_use]
pub fn from_cartesian_pure(
    point: &[f64],
    to: CoordinateSystem,
) -> Vec<f64> {
    let (x, y, z) = (
        component(point, 0),
        component(point, 1),
        component(point, 2),
    );
    match to {
        | CoordinateSystem::Cartesian => point.to_vec(),
        | CoordinateSystem::Polar => vec![x.hypot(y), y.atan2(x)],
        | CoordinateSystem::Cylindrical => vec![x.hypot(y), y.atan2(x), z],
        | CoordinateSystem::Spherical => cartesian_to_spherical([x, y, z]).to_vec(),
    }
}

/// Inverse map to Cartesian: converts a point of `from` to Cartesian.
///
/// # Errors
///
/// Returns an error if `point` has the wrong number of components.
pub fn to_cartesian(
    point: &[f64],
    from: CoordinateSystem,
) -> Result<Vec<f64>, String> {
    check_len(point, from)?;
    Ok(to_cartesian_pure(point, from))
}

/// Converts a Cartesian point to the system `to`.
///
/// # Errors
///
/// Returns an error if `point` has the wrong number of components (2 for
/// polar, 3 for cylindrical and spherical).
pub fn from_cartesian(
    point: &[f64],
    to: CoordinateSystem,
) -> Result<Vec<f64>, String> {
    check_len(point, CoordinateSystem::Cartesian)?;
    check_len(point, to)?;
    Ok(from_cartesian_pure(point, to))
}

/// Transforms `point` from the system `from` to the system `to`, going
/// through Cartesian coordinates, with validation of the dimensions.
///
/// # Errors
///
/// Returns an error if `point` has the wrong number of components for `from`,
/// or if the dimension of `from` and `to` differ (2D polar against 3D
/// systems).
pub fn transform_point(
    point: &[f64],
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Result<Vec<f64>, String> {
    check_len(point, from)?;
    if from == to {
        return Ok(point.to_vec());
    }
    if from == CoordinateSystem::Cartesian {
        check_len(point, to)?;
        return Ok(from_cartesian_pure(point, to));
    }
    let cartesian = to_cartesian_pure(point, from);
    if to == CoordinateSystem::Cartesian {
        Ok(cartesian)
    } else {
        from_cartesian(&cartesian, to)
    }
}

/// Like [`transform_point`] but without validation: a direct, allocation
/// light evaluation for hot loops and for differentiating the map. Missing
/// components are taken as zero.
#[must_use]
pub fn transform_point_pure(
    point: &[f64],
    from: CoordinateSystem,
    to: CoordinateSystem,
) -> Vec<f64> {
    match (from, to) {
        | (from, to) if from == to => point.to_vec(),
        | (CoordinateSystem::Cartesian, to) => from_cartesian_pure(point, to),
        | (from, CoordinateSystem::Cartesian) => to_cartesian_pure(point, from),
        | _ => from_cartesian_pure(&to_cartesian_pure(point, from), to),
    }
}

/// Jacobian `J[i][j] = d(to_i) / d(from_j)` of the map `from -> to` at
/// `at_point`, by central differences.
///
/// # Errors
///
/// Returns an error if the point is invalid for the transformation, see
/// [`transform_point`].
pub fn numerical_jacobian(
    from: CoordinateSystem,
    to: CoordinateSystem,
    at_point: &[f64],
) -> Result<Vec<Vec<f64>>, String> {
    let center = transform_point(at_point, from, to)?;
    let mut jacobian = vec![vec![0.0; at_point.len()]; center.len()];
    for j in 0..at_point.len() {
        let h = JACOBIAN_STEP * at_point[j].abs().max(1.0);
        let mut shifted = at_point.to_vec();
        shifted[j] = at_point[j] + h;
        let plus = transform_point_pure(&shifted, from, to);
        shifted[j] = at_point[j] - h;
        let minus = transform_point_pure(&shifted, from, to);
        for (i, row) in jacobian.iter_mut().enumerate() {
            row[j] = (plus[i] - minus[i]) / (2.0 * h);
        }
    }
    Ok(jacobian)
}

/// Determinant of a small square matrix via [`super::dense::lu_factor`].
fn determinant(matrix: &[Vec<f64>]) -> f64 {
    let Ok(mat) = super::dense::Mat::from_rows(matrix) else {
        return 0.0;
    };
    super::dense::lu_factor(&mat).map_or(0.0, |lu| lu.det())
}

/// Determinant of the Jacobian of the map `from -> to` at `at_point`.
///
/// This is the local volume-element factor `dV_to = |det J| dV_from` of the
/// two coordinate charts: for example `r` for polar-to-Cartesian and
/// `r^2 sin(theta)` for spherical-to-Cartesian, while the inverse maps carry
/// the reciprocal determinant. It is built from the central-difference
/// Jacobian of [`numerical_jacobian`] and inherits its accuracy and precision
/// caveat (no analytic Jacobian).
///
/// # Errors
/// As [`numerical_jacobian`].
pub fn jacobian_determinant(
    from: CoordinateSystem,
    to: CoordinateSystem,
    at_point: &[f64],
) -> Result<f64, String> {
    let jacobian = numerical_jacobian(from, to, at_point)?;
    let rows = jacobian.len();
    if rows != jacobian.first().map_or(0, Vec::len) {
        return Err("Jacobian is not square".to_string());
    }
    Ok(determinant(&jacobian))
}
