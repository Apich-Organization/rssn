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

/// Converts a point of `from` to Cartesian coordinates, without checking
/// its length: missing components are taken as zero.
#[must_use]
pub fn to_cartesian_pure(
    point: &[f64],
    from: CoordinateSystem,
) -> Vec<f64> {
    let (a, b, c) = (component(point, 0), component(point, 1), component(point, 2));
    match from {
        | CoordinateSystem::Cartesian => point.to_vec(),
        | CoordinateSystem::Polar => vec![a * b.cos(), a * b.sin()],
        | CoordinateSystem::Cylindrical => vec![a * b.cos(), a * b.sin(), c],
        | CoordinateSystem::Spherical => {
            vec![a * b.sin() * c.cos(), a * b.sin() * c.sin(), a * b.cos()]
        },
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
    let (x, y, z) = (component(point, 0), component(point, 1), component(point, 2));
    match to {
        | CoordinateSystem::Cartesian => point.to_vec(),
        | CoordinateSystem::Polar => vec![x.hypot(y), y.atan2(x)],
        | CoordinateSystem::Cylindrical => vec![x.hypot(y), y.atan2(x), z],
        | CoordinateSystem::Spherical => {
            let r = x.hypot(y).hypot(z);
            let theta = if r == 0.0 { 0.0 } else { (z / r).clamp(-1.0, 1.0).acos() };
            vec![r, theta, y.atan2(x)]
        },
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
    let cartesian = to_cartesian_pure(point, from);
    from_cartesian(&cartesian, to)
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
    if from == to {
        return point.to_vec();
    }
    from_cartesian_pure(&to_cartesian_pure(point, from), to)
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
