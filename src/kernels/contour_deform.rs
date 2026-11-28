//! Numerical contour integration and deformation.
//!
//! Includes keyhole, dogbone, indented semicircle, and Hankel contours, Cauchy principal values,
//! residues at infinity, and adaptive deformation around singularities and branch cuts.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;

#[path = "contour_quadrature.rs"]
mod checked;
#[path = "contour_obstacles.rs"]
mod obstacles;
pub use checked::ContourIntegral;
pub use checked::ContourQuadratureError;
pub use checked::cauchy_principal_value_checked;
pub use checked::integrate_contour_checked;
pub use checked::integrate_segment_checked;
pub use checked::residue_at_infinity_checked;
pub(crate) use obstacles::DeformationError;
pub(crate) use obstacles::MAX_OBSTACLES;
pub(crate) use obstacles::deformation_checked;

const I: Complex64 = Complex64::new(0.0, 1.0);

/// A piecewise-smooth contour segment in `C`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ContourSegment {
    /// Straight line segment from `start` to `end`.
    Line {
        /// Starting point.
        start: Complex64,
        /// Ending point.
        end: Complex64,
    },
    /// Circular arc centered at `center` with `radius` from `theta0` to `theta1`.
    Arc {
        /// Circle center.
        center: Complex64,
        /// Circle radius (positive).
        radius: f64,
        /// Starting angle in radians.
        theta0: f64,
        /// Ending angle in radians.
        theta1: f64,
    },
}

/// Branch-cut ray `{ origin + r * exp(i * angle) : r >= 0 }`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BranchRay {
    /// Branch point origin.
    pub origin: Complex64,
    /// Ray angle in radians.
    pub angle: f64,
}

fn legacy_tolerance(tolerance: f64) -> f64 {
    // The historical tolerance applied to each scalar component.
    if tolerance.is_finite() && tolerance >= 0.0 {
        (2.0 * tolerance.min(f64::MAX * 0.5)).max(f64::MIN_POSITIVE)
    } else {
        tolerance
    }
}

/// Integrates `f(z)` along a single [`ContourSegment`] using adaptive Gauss-Kronrod.
pub fn integrate_segment<F>(
    f: &F,
    seg: ContourSegment,
    tol: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    if let ContourSegment::Arc { radius, .. } = seg
        && radius <= 0.0
    {
        return Complex64::ZERO;
    }
    integrate_segment_checked(f, seg, legacy_tolerance(tol), 24)
        .map_or(Complex64::new(f64::NAN, f64::NAN), |result| result.value)
}

/// Integrates `f(z)` along a sequence of [`ContourSegment`]s.
pub fn integrate_contour<F>(
    f: &F,
    segments: &[ContourSegment],
    tol: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let per_seg_tol = (tol / (segments.len().max(1) as f64)).max(1e-14);
    let mut sum = Complex64::new(0.0, 0.0);
    for &seg in segments {
        sum += integrate_segment(f, seg, per_seg_tol);
    }
    sum
}

/// Keyhole contour around a branch cut at `branch_pt` along `cut_angle`.
#[must_use]
pub fn keyhole_contour(
    branch_pt: Complex64,
    r_inner: f64,
    r_outer: f64,
    cut_angle: f64,
    gap: f64,
) -> Vec<ContourSegment> {
    let theta_above = cut_angle + gap;
    let theta_below = cut_angle + TAU - gap;
    let p_out_above = branch_pt + Complex64::from_polar(r_outer, theta_above);
    let p_out_below = branch_pt + Complex64::from_polar(r_outer, theta_below);
    let p_in_above = branch_pt + Complex64::from_polar(r_inner, theta_above);
    let p_in_below = branch_pt + Complex64::from_polar(r_inner, theta_below);
    vec![
        ContourSegment::Arc {
            center: branch_pt,
            radius: r_outer,
            theta0: theta_above,
            theta1: theta_below,
        },
        ContourSegment::Line {
            start: p_out_below,
            end: p_in_below,
        },
        ContourSegment::Arc {
            center: branch_pt,
            radius: r_inner,
            theta0: theta_below,
            theta1: theta_above,
        },
        ContourSegment::Line {
            start: p_in_above,
            end: p_out_above,
        },
    ]
}

/// Dogbone (dumbbell) contour surrounding the finite branch cut `[a, b]`.
#[must_use]
pub fn dogbone_contour(
    a: Complex64,
    b: Complex64,
    radius: f64,
) -> Vec<ContourSegment> {
    let dir = b - a;
    let len = dir.norm();
    if len == 0.0 {
        return vec![ContourSegment::Arc {
            center: a,
            radius,
            theta0: 0.0,
            theta1: TAU,
        }];
    }
    let phi = dir.arg();
    let perp = Complex64::from_polar(radius, phi + 0.5 * PI);
    vec![
        ContourSegment::Line {
            start: a - perp,
            end: b - perp,
        },
        ContourSegment::Arc {
            center: b,
            radius,
            theta0: phi - 0.5 * PI,
            theta1: phi + 0.5 * PI,
        },
        ContourSegment::Line {
            start: b + perp,
            end: a + perp,
        },
        ContourSegment::Arc {
            center: a,
            radius,
            theta0: phi + 0.5 * PI,
            theta1: phi + 1.5 * PI,
        },
    ]
}

/// Indented semicircle contour on `[-radius, radius]` bypassing `real_poles`
/// with small semicircular indents into the upper (or lower) half-plane.
#[must_use]
pub fn indented_semicircle_contour(
    radius: f64,
    real_poles: &[f64],
    indent_r: f64,
    upper_half: bool,
) -> Vec<ContourSegment> {
    let mut poles: Vec<f64> = real_poles
        .iter()
        .copied()
        .filter(|&p| p > -radius + indent_r && p < radius - indent_r)
        .collect();
    poles.sort_by(f64::total_cmp);

    let mut segs = Vec::new();
    let mut cur_x = -radius;
    for (idx, &p) in poles.iter().enumerate() {
        let mut r_use = indent_r;
        if p - r_use < cur_x {
            r_use = (p - cur_x).max(0.0);
        }
        if let Some(&next_p) = poles.get(idx + 1) {
            r_use = r_use.min(0.45 * (next_p - p));
        }
        if p - r_use > cur_x + 1e-14 {
            segs.push(ContourSegment::Line {
                start: Complex64::new(cur_x, 0.0),
                end: Complex64::new(p - r_use, 0.0),
            });
        }
        let (t0, t1) = if upper_half {
            (PI, 0.0)
        } else {
            (-PI, 0.0)
        };
        if r_use > 1e-14 {
            segs.push(ContourSegment::Arc {
                center: Complex64::new(p, 0.0),
                radius: r_use,
                theta0: t0,
                theta1: t1,
            });
        }
        cur_x = p + r_use;
    }
    if radius > cur_x {
        segs.push(ContourSegment::Line {
            start: Complex64::new(cur_x, 0.0),
            end: Complex64::new(radius, 0.0),
        });
    }
    let (arc_t0, arc_t1) = if upper_half {
        (0.0, PI)
    } else {
        (0.0, -PI)
    };
    segs.push(ContourSegment::Arc {
        center: Complex64::new(0.0, 0.0),
        radius,
        theta0: arc_t0,
        theta1: arc_t1,
    });
    segs
}

/// Open Hankel contour around `0` avoiding the branch cut along `cut_angle`.
///
///
/// Enters from `r_outer` along `cut_angle + gap`, wraps `2*pi - 2*gap` around
/// the origin at `r_inner`, and exits to `r_outer` along `cut_angle - gap`.
#[must_use]
pub fn hankel_contour(
    r_inner: f64,
    r_outer: f64,
    cut_angle: f64,
    gap: f64,
) -> Vec<ContourSegment> {
    let theta_above = cut_angle + gap;
    let theta_below_wrap = cut_angle + TAU - gap;
    vec![
        ContourSegment::Line {
            start: Complex64::from_polar(r_outer, theta_above),
            end: Complex64::from_polar(r_inner, theta_above),
        },
        ContourSegment::Arc {
            center: Complex64::new(0.0, 0.0),
            radius: r_inner,
            theta0: theta_above,
            theta1: theta_below_wrap,
        },
        ContourSegment::Line {
            start: Complex64::from_polar(r_inner, theta_below_wrap),
            end: Complex64::from_polar(r_outer, theta_below_wrap),
        },
    ]
}

/// Vertical Bromwich contour `c + i * y` for `y in [-y_max, y_max]`.
#[must_use]
pub fn bromwich_contour(
    c: f64,
    y_max: f64,
) -> Vec<ContourSegment> {
    vec![ContourSegment::Line {
        start: Complex64::new(c, -y_max),
        end: Complex64::new(c, y_max),
    }]
}

/// Counterclockwise rectangular contour `[x0, x1] x [y0, y1]`.
#[must_use]
pub fn rectangular_contour(
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
) -> Vec<ContourSegment> {
    let p00 = Complex64::new(x0, y0);
    let p10 = Complex64::new(x1, y0);
    let p11 = Complex64::new(x1, y1);
    let p01 = Complex64::new(x0, y1);
    vec![
        ContourSegment::Line { start: p00, end: p10 },
        ContourSegment::Line { start: p10, end: p11 },
        ContourSegment::Line { start: p11, end: p01 },
        ContourSegment::Line { start: p01, end: p00 },
    ]
}

/// Closed sector (wedge) contour from `0` to `radius`, along arc `0..angle`, and back to `0`.
#[must_use]
pub fn sector_contour(
    radius: f64,
    angle: f64,
) -> Vec<ContourSegment> {
    let p1 = Complex64::new(radius, 0.0);
    let p2 = Complex64::from_polar(radius, angle);
    vec![
        ContourSegment::Line {
            start: Complex64::new(0.0, 0.0),
            end: p1,
        },
        ContourSegment::Arc {
            center: Complex64::new(0.0, 0.0),
            radius,
            theta0: 0.0,
            theta1: angle,
        },
        ContourSegment::Line {
            start: p2,
            end: Complex64::new(0.0, 0.0),
        },
    ]
}

/// Returns an oriented path avoiding the supplied branch rays and maintaining
/// the given distance from poles.
///
/// A direct path is reused when it is clear. An empty result means that
/// the bounded construction could not verify a path.
#[must_use]
pub fn deform_path_around_obstacles(
    start: Complex64,
    end: Complex64,
    poles: &[Complex64],
    cuts: &[BranchRay],
    clearance: f64,
) -> Vec<ContourSegment> {
    deformation_checked(start, end, poles, cuts, clearance).unwrap_or_default()
}

/// Integrates an automatically deformed contour from `start` to `end`.
///
/// Returns complex NaN if the geometry cannot be verified or quadrature does
/// not converge. The numerical ABI also exposes checked error and work data.
pub fn integrate_deformed<F>(
    f: &F,
    start: Complex64,
    end: Complex64,
    poles: &[Complex64],
    cuts: &[BranchRay],
    clearance: f64,
    tol: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let nan = Complex64::new(f64::NAN, f64::NAN);
    let Ok(path) = deformation_checked(start, end, poles, cuts, clearance) else {
        return nan;
    };
    integrate_contour_checked(f, &path, tol, 4096).map_or(nan, |estimate| {
        if estimate.converged {
            estimate.value
        } else {
            nan
        }
    })
}

/// Numerical residue at infinity: `Res_{z=oo} f(z) = -1/(2*pi*i) * oint_{|z|=R} f(z) dz`.
pub fn residue_at_infinity_numeric<F>(
    f: F,
    radius: f64,
    n_samples: usize,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let integral = crate::kernels::complex::circle_integral(
        f,
        Complex64::new(0.0, 0.0),
        radius,
        n_samples.max(32),
    );
    -integral / (TAU * I)
}

/// Cauchy principal value `P.V. int_a^b f(x) dx` with real simple poles in `(a, b)`
/// via symmetric folding `int_0^delta (f(p + t) + f(p - t)) dt`.
pub fn cauchy_principal_value<F>(
    f: F,
    a: f64,
    b: f64,
    poles: &[f64],
    tol: f64,
) -> Complex64
where
    F: Fn(f64) -> Complex64,
{
    if a >= b {
        return Complex64::new(0.0, 0.0);
    }
    let mut interior: Vec<f64> = poles.iter().copied().filter(|&p| p > a && p < b).collect();
    interior.sort_by(f64::total_cmp);
    interior.dedup_by(|x, y| (*x - *y).abs() < 1e-12);

    let panels = 24_usize
        .saturating_mul(2_usize.saturating_mul(interior.len()).saturating_add(1))
        .min(4096);
    cauchy_principal_value_checked(f, a, b, &interior, legacy_tolerance(tol), panels)
        .map_or(Complex64::new(f64::NAN, f64::NAN), |result| result.value)
}
