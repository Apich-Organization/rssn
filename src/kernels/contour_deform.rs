//! Numerical contour integration, standard complex contours (keyhole,
//! dogbone, indented semicircle, Hankel), Cauchy principal value, residue at
//! infinity, and adaptive singularity/branch-cut contour deformation.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;

use crate::kernels::integrate::gauss_kronrod_any;

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

/// Integrates `f(z)` along a single [`ContourSegment`] using adaptive Gauss-Kronrod.
pub fn integrate_segment<F>(
    f: &F,
    seg: ContourSegment,
    tol: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    match seg {
        | ContourSegment::Line { start, end } => {
            let dz = end - start;
            if dz.norm() == 0.0 {
                return Complex64::new(0.0, 0.0);
            }
            let re = gauss_kronrod_any(|t| (f(start + dz * t) * dz).re, 0.0, 1.0, tol, 24).value;
            let im = gauss_kronrod_any(|t| (f(start + dz * t) * dz).im, 0.0, 1.0, tol, 24).value;
            Complex64::new(re, im)
        },
        | ContourSegment::Arc {
            center,
            radius,
            theta0,
            theta1,
        } => {
            let dtheta = theta1 - theta0;
            if dtheta.abs() == 0.0 || radius <= 0.0 {
                return Complex64::new(0.0, 0.0);
            }
            let eval_integrand = |t: f64| {
                let theta = theta0 + dtheta * t;
                let dir = Complex64::from_polar(1.0, theta);
                let z = center + radius * dir;
                let dz_dt = I * radius * dtheta * dir;
                f(z) * dz_dt
            };
            let re = gauss_kronrod_any(|t| eval_integrand(t).re, 0.0, 1.0, tol, 24).value;
            let im = gauss_kronrod_any(|t| eval_integrand(t).im, 0.0, 1.0, tol, 24).value;
            Complex64::new(re, im)
        },
    }
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
    for p in poles {
        if p - indent_r > cur_x {
            segs.push(ContourSegment::Line {
                start: Complex64::new(cur_x, 0.0),
                end: Complex64::new(p - indent_r, 0.0),
            });
        }
        let (t0, t1) = if upper_half {
            (PI, 0.0)
        } else {
            (-PI, 0.0)
        };
        segs.push(ContourSegment::Arc {
            center: Complex64::new(p, 0.0),
            radius: indent_r,
            theta0: t0,
            theta1: t1,
        });
        cur_x = p + indent_r;
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

/// Open Hankel contour around `0` along `cut_angle`.
#[must_use]
pub fn hankel_contour(
    r_inner: f64,
    r_outer: f64,
    cut_angle: f64,
    gap: f64,
) -> Vec<ContourSegment> {
    let theta_below = cut_angle - gap;
    let theta_above = cut_angle + gap;
    vec![
        ContourSegment::Line {
            start: Complex64::from_polar(r_outer, theta_below),
            end: Complex64::from_polar(r_inner, theta_below),
        },
        ContourSegment::Arc {
            center: Complex64::new(0.0, 0.0),
            radius: r_inner,
            theta0: theta_below,
            theta1: theta_above,
        },
        ContourSegment::Line {
            start: Complex64::from_polar(r_inner, theta_above),
            end: Complex64::from_polar(r_outer, theta_above),
        },
    ]
}

/// Checks if segment `[a, b]` intersects ray `origin + r * exp(i * angle)` (`r >= 0`).
fn segment_hits_ray(
    a: Complex64,
    b: Complex64,
    ray: BranchRay,
) -> bool {
    let d = b - a;
    let dir = Complex64::from_polar(1.0, ray.angle);
    // Solve a + t*d = origin + s*dir => cross(d, dir) * t = cross(origin - a, dir)
    let cross_dd = d.re * dir.im - d.im * dir.re;
    if cross_dd.abs() < 1e-13 {
        return false;
    }
    let diff = ray.origin - a;
    let t = (diff.re * dir.im - diff.im * dir.re) / cross_dd;
    let s = (diff.re * d.im - diff.im * d.re) / cross_dd;
    t > 0.0 && t < 1.0 && s >= 0.0
}

/// Deforms the straight path `[start, end]` around nearby `poles` and crossed
/// `BranchRay`s by inserting circular detour arcs of radius `clearance`.
#[must_use]
pub fn deform_path_around_obstacles(
    start: Complex64,
    end: Complex64,
    poles: &[Complex64],
    cuts: &[BranchRay],
    clearance: f64,
) -> Vec<ContourSegment> {
    let d = end - start;
    let len = d.norm();
    if len == 0.0 || clearance <= 0.0 {
        return vec![ContourSegment::Line { start, end }];
    }
    let unit = d / len;
    let mut obstacles: Vec<(f64, Complex64, bool)> = Vec::new();

    for &p in poles {
        let rel = p - start;
        let proj = rel.re * unit.re + rel.im * unit.im;
        let perp = rel.re * unit.im - rel.im * unit.re;
        if proj > 0.0 && proj < len && perp.abs() < clearance {
            obstacles.push((proj, p, perp >= 0.0));
        }
    }
    for &ray in cuts {
        if segment_hits_ray(start, end, ray) {
            let rel = ray.origin - start;
            let proj = (rel.re * unit.re + rel.im * unit.im)
                .clamp(clearance, (len - clearance).max(clearance));
            // Detour around the branch point origin on the side opposite the ray direction.
            let dir = Complex64::from_polar(1.0, ray.angle);
            let side = (dir.re * unit.im - dir.im * unit.re) < 0.0;
            obstacles.push((proj, ray.origin, side));
        }
    }

    if obstacles.is_empty() {
        return vec![ContourSegment::Line { start, end }];
    }
    obstacles.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut segs = Vec::new();
    let mut cur = start;
    let line_angle = d.arg();
    for (_proj, center, detour_right) in obstacles {
        let entry_pt = center - unit * clearance;
        let exit_pt = center + unit * clearance;
        if (entry_pt - cur).norm() > 1e-12 && (entry_pt - start).norm() < (exit_pt - start).norm() {
            segs.push(ContourSegment::Line {
                start: cur,
                end: entry_pt,
            });
            let t_entry = line_angle + PI;
            let t_exit = if detour_right {
                line_angle + TAU
            } else {
                line_angle
            };
            segs.push(ContourSegment::Arc {
                center,
                radius: clearance,
                theta0: t_entry,
                theta1: t_exit,
            });
            cur = exit_pt;
        }
    }
    if (end - cur).norm() > 1e-12 {
        segs.push(ContourSegment::Line { start: cur, end });
    }
    segs
}

/// Integrates `f(z)` from `start` to `end` along an automatically deformed path
/// avoiding `poles` and `cuts` by at least `clearance`.
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
    let segs = deform_path_around_obstacles(start, end, poles, cuts, clearance);
    integrate_contour(f, &segs, tol)
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

    if interior.is_empty() {
        let re = gauss_kronrod_any(|x| f(x).re, a, b, tol, 24).value;
        let im = gauss_kronrod_any(|x| f(x).im, a, b, tol, 24).value;
        return Complex64::new(re, im);
    }

    let mut total = Complex64::new(0.0, 0.0);
    let mut cur = a;
    for (idx, &p) in interior.iter().enumerate() {
        let right_bound = interior.get(idx + 1).map_or(b, |&next| 0.5 * (p + next));
        let delta = (p - cur).min(right_bound - p);
        if p - delta > cur {
            let re = gauss_kronrod_any(|x| f(x).re, cur, p - delta, tol, 24).value;
            let im = gauss_kronrod_any(|x| f(x).im, cur, p - delta, tol, 24).value;
            total += Complex64::new(re, im);
        }
        // Symmetric folded integral on [0, delta]
        let re_fold = gauss_kronrod_any(|t| (f(p + t) + f(p - t)).re, 0.0, delta, tol, 24).value;
        let im_fold = gauss_kronrod_any(|t| (f(p + t) + f(p - t)).im, 0.0, delta, tol, 24).value;
        total += Complex64::new(re_fold, im_fold);
        cur = p + delta;
    }
    if b > cur {
        let re = gauss_kronrod_any(|x| f(x).re, cur, b, tol, 24).value;
        let im = gauss_kronrod_any(|x| f(x).im, cur, b, tol, 24).value;
        total += Complex64::new(re, im);
    }
    total
}
