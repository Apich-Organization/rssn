//! Complex coordinate systems, Riemann sphere stereographic projection, and topology.
//!
//! Covers hyperbolic disk/upper half-plane metrics, custom branch-cut control,
//! continuous phase unwrapping, and algebraic Riemann surface topology.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;

/// Polar coordinates `(r, theta)` with `r = |z|` and the signed principal angle.
/// On the negative real axis, signed zero distinguishes the two lips at ±π.
#[must_use]
pub fn to_polar(z: Complex64) -> (f64, f64) {
    (z.norm(), z.arg())
}

/// Reconstructs `r * exp(i * theta)`.
#[must_use]
pub fn from_polar(
    r: f64,
    theta: f64,
) -> Complex64 {
    Complex64::from_polar(r, theta)
}

/// `cis(theta)` returns the unit complex number `cos(theta) + i * sin(theta)`.
/// The angle is in radians.
#[must_use]
pub fn cis(theta: f64) -> Complex64 {
    let (s, c) = theta.sin_cos();
    Complex64::new(c, s)
}

/// Squared modulus `|z|^2 = re^2 + im^2`.
#[must_use]
pub const fn abs_sq(z: Complex64) -> f64 {
    z.re * z.re + z.im * z.im
}

/// Complex signum `z / |z|`, or `0` if `z == 0`.
#[must_use]
pub fn sign_complex(z: Complex64) -> Complex64 {
    let magnitude = z.norm();
    if magnitude.is_normal() {
        return z / magnitude;
    }
    let scale = z.re.abs().max(z.im.abs());
    let direction = if z.is_finite() && scale > 0.0 {
        // Normalize before dividing by the norm. The true norm may be more
        // than Float::MAX, or be rounded coarsely in the subnormal range.
        z / scale
    } else if z.re.is_infinite() || z.im.is_infinite() {
        Complex64::new(
            if z.re.is_infinite() {
                z.re.signum()
            } else {
                0.0
            },
            if z.im.is_infinite() {
                z.im.signum()
            } else {
                0.0
            },
        )
    } else {
        return Complex64::new(0.0, 0.0);
    };
    direction / direction.norm()
}

/// Stereographic projection from `C` to the unit Riemann sphere `S^2`.
///
/// Maps `0` to `(0, 0, -1)` and `infinity` to `(0, 0, 1)`.
/// Uses `1/z` scaling when `|z| > 1` to avoid overflow for huge `|z|`.
#[must_use]
pub fn to_riemann_sphere(z: Complex64) -> (f64, f64, f64) {
    if !z.re.is_finite() || !z.im.is_finite() {
        return (0.0, 0.0, 1.0);
    }
    let scale = z.re.abs().max(z.im.abs());
    if scale <= 1.0 {
        let r2 = abs_sq(z);
        let d = 1.0 + r2;
        return (2.0 * z.re / d, 2.0 * z.im / d, (r2 - 1.0) / d);
    }
    // Let w = 1 / z: X = 2 Re(w) / (1 + |w|^2), Y = -2 Im(w) / (1 + |w|^2), Z = (1 - |w|^2) / (1 + |w|^2).
    let reciprocal = crate::kernels::complex_scalar::Complex::new(z.re, z.im).recip();
    let w = Complex64::new(reciprocal.re, reciprocal.im);
    let w2 = abs_sq(w);
    let d = 1.0 + w2;
    (2.0 * w.re / d, -2.0 * w.im / d, (1.0 - w2) / d)
}

/// Inverse stereographic projection from `S^2` to `C`.
/// Returns `None` at the north pole `(0, 0, 1)` (infinity).
#[must_use]
pub fn from_riemann_sphere(
    x: f64,
    y: f64,
    z: f64,
) -> Option<Complex64> {
    if !x.is_finite() || !y.is_finite() || !z.is_finite() {
        return None;
    }
    if z > 0.0 {
        // On the unit sphere, x² + y² = (1 - z)(1 + z). The transverse
        // coordinates retain the point even if its height rounds to 1.
        return divide_chart_coordinates(Complex64::new(1.0 + z, 0.0), Complex64::new(x, -y));
    }
    let denom = 1.0 - z;
    Some(Complex64::new(x / denom, y / denom))
}

/// `chordal_distance(z1, z2)` returns the Euclidean distance of the images on the unit sphere.
/// `2 * |z1 - z2| / sqrt((1 + |z1|^2) * (1 + |z2|^2))`.
/// Overflow-safe for large moduli.
#[must_use]
pub fn chordal_distance(
    z1: Complex64,
    z2: Complex64,
) -> f64 {
    let p1 = to_riemann_sphere(z1);
    let p2 = to_riemann_sphere(z2);
    let dx = p1.0 - p2.0;
    let dy = p1.1 - p2.1;
    let dz = p1.2 - p2.2;
    dx.hypot(dy).hypot(dz)
}

/// `spherical_metric_factor(z)` returns the conformal factor `2 / (1 + |z|^2)`.
#[must_use]
pub fn spherical_metric_factor(z: Complex64) -> f64 {
    if !z.re.is_finite() || !z.im.is_finite() {
        return 0.0;
    }
    let s = z.re.abs().max(z.im.abs());
    if s > 1.0 {
        let inverse_scale = s.recip();
        let relative = z / s;
        // Scale the denominator before the final product so subnormal
        // results round only once at the output magnitude.
        let denominator = inverse_scale * inverse_scale + abs_sq(relative);
        return ((2.0 * inverse_scale) / denominator) * inverse_scale;
    }
    2.0 / (1.0 + abs_sq(z))
}

fn divide_chart_coordinates(
    numerator: Complex64,
    denominator: Complex64,
) -> Option<Complex64> {
    use crate::kernels::complex_scalar::Complex;

    if denominator == Complex64::new(0.0, 0.0) {
        return None;
    }
    let numerator = Complex::new(numerator.re, numerator.im);
    let denominator = Complex::new(denominator.re, denominator.im);
    let quotient = numerator.div(denominator);
    Some(Complex64::new(quotient.re, quotient.im))
}

/// Cayley transform mapping the upper half-plane `{ Im(z) > 0 }` to the unit disk `{ |w| < 1 }`:
/// `w = (z - i) / (z + i)`. Returns `None` at `z = -i`.
#[must_use]
pub fn uhp_to_disk(z: Complex64) -> Option<Complex64> {
    let i = Complex64::new(0.0, 1.0);
    divide_chart_coordinates(z - i, z + i)
}

/// Inverse Cayley transform mapping the unit disk `{ |w| < 1 }` to the upper half-plane:
/// `z = i * (1 + w) / (1 - w)`. Returns `None` at `w = 1`.
#[must_use]
pub fn disk_to_uhp(w: Complex64) -> Option<Complex64> {
    let one = Complex64::new(1.0, 0.0);
    let denominator = one - w;
    let denominator_sq = abs_sq(denominator);
    if denominator_sq.is_normal() {
        // The imaginary numerator is the signed radial gap. Computing it
        // through the chart predicate retains the side of the unit circle.
        return Some(Complex64::new(
            (-2.0 * w.im) / denominator_sq,
            poincare_disk_gap(w) / denominator_sq,
        ));
    }
    let i = Complex64::new(0.0, 1.0);
    divide_chart_coordinates(i * (one + w), denominator)
}

/// Signed radial gap `1 - |z|²` for the open-disk chart.
pub(crate) fn poincare_disk_gap(z: Complex64) -> f64 {
    let re_sq = z.re * z.re;
    let im_sq = z.im * z.im;
    let norm_sq = re_sq + im_sq;
    let gap = 1.0 - norm_sq;
    if gap.abs() < 1.0 / 256.0 {
        return disk_gap_near_circle(z, re_sq, im_sq, norm_sq, gap);
    }
    gap
}

#[cold]
#[allow(
    clippy::single_call_fn,
    reason = "keeps compensated work off the ordinary chart path"
)]
fn disk_gap_near_circle(
    z: Complex64,
    re_sq: f64,
    im_sq: f64,
    norm_sq: f64,
    gap: f64,
) -> f64 {
    use crate::kernels::dd::Dd;

    // A nonzero component next to an exact unit component is exterior,
    // even when its square is too small for binary64.
    let unit_bits = 1.0_f64.to_bits();
    if (re_sq.to_bits() == unit_bits && im_sq == 0.0 && z.im != 0.0)
        || (im_sq.to_bits() == unit_bits && re_sq == 0.0 && z.re != 0.0)
    {
        return -f64::from_bits(1);
    }

    // Near one, 1 - norm_sq is exact. Retain the two product errors and
    // the error of adding the squares as a short nonoverlapping expansion.
    let norm_error = re_sq.min(im_sq) - (norm_sq - re_sq.max(im_sq));
    let errors = [
        -norm_error,
        -z.re.mul_add(z.re, -re_sq),
        -z.im.mul_add(z.im, -im_sq),
    ];
    let mut expansion = [gap, 0.0, 0.0, 0.0];
    for (length, mut carry) in (1..).zip(errors) {
        let (previous, next) = expansion.split_at_mut(length);
        for limb in previous {
            let sum = Dd::from_parts(carry, *limb);
            *limb = sum.lo;
            carry = sum.hi;
        }
        if let Some(limb) = next.first_mut() {
            *limb = carry;
        }
    }
    expansion.into_iter().sum()
}

/// Hyperbolic distance in the Poincaré unit disk `|z| < 1`.
#[must_use]
pub fn poincare_disk_distance(
    z1: Complex64,
    z2: Complex64,
) -> Option<f64> {
    let gap1 = poincare_disk_gap(z1);
    let gap2 = poincare_disk_gap(z2);
    if !(gap1 > 0.0 && gap2 > 0.0) {
        return None;
    }
    // sinh(d/2) retains separation when the pseudohyperbolic ratio
    // would round to 1 for two interior points near opposite edges.
    let ratio = (z1 - z2).norm() / (gap1.sqrt() * gap2.sqrt());
    Some(2.0 * ratio.asinh())
}

/// Hyperbolic distance in the upper half-plane `Im(z) > 0`.
#[must_use]
pub fn uhp_distance(
    z1: Complex64,
    z2: Complex64,
) -> Option<f64> {
    if !z1.is_finite() || !z2.is_finite() || z1.im <= 0.0 || z2.im <= 0.0 {
        return None;
    }
    let separation = (z1 - z2).norm();
    let ratio = separation / (z1.im.sqrt() * z2.im.sqrt());
    if ratio <= f64::EPSILON.sqrt() {
        // The cubic correction to 2*asinh(ratio/2) is below rounding error.
        return Some(ratio);
    }
    if ratio.is_finite() {
        return Some(2.0 * (0.5 * ratio).asinh());
    }
    // Use the logarithm of the ratio if the separation or the ratio
    // itself overflowed. Quartering leaves room for both component
    // differences and for their hypotenuse.
    let log_separation = if separation.is_finite() {
        separation.ln()
    } else {
        (0.25 * z1 - 0.25 * z2).norm().ln() + 2.0 * std::f64::consts::LN_2
    };
    let log_ratio = log_separation - z1.im.ln().midpoint(z2.im.ln());
    Some(if log_ratio > 40.0 {
        // At this scale, 2*asinh(exp(log_ratio)/2) rounds as 2*log_ratio.
        2.0 * log_ratio
    } else {
        2.0 * (0.5 * log_ratio.exp()).asinh()
    })
}

/// Argument of `z` in `(theta0, theta0 + 2*pi]`, placing the branch cut along ray `theta0`.
#[must_use]
pub fn arg_branch(
    z: Complex64,
    theta0: f64,
) -> f64 {
    let a = z.arg();
    let mut shifted = (a - theta0).rem_euclid(TAU);
    if shifted == 0.0 {
        shifted = TAU;
    }
    theta0 + shifted
}

/// Logarithm `ln|z| + i * arg_branch(z, theta0)` with branch cut along ray `theta0`.
#[must_use]
pub fn log_cut(
    z: Complex64,
    theta0: f64,
) -> Complex64 {
    Complex64::new(z.norm().ln(), arg_branch(z, theta0))
}

/// Square root `sqrt(|z|) * cis(arg_branch(z, theta0) / 2)` with branch cut along ray `theta0`.
#[must_use]
pub fn sqrt_cut(
    z: Complex64,
    theta0: f64,
) -> Complex64 {
    Complex64::from_polar(z.norm().sqrt(), 0.5 * arg_branch(z, theta0))
}

/// Complex power `exp(w * log_cut(z, theta0))` with branch cut along ray `theta0`.
#[must_use]
pub fn pow_cut(
    z: Complex64,
    w: Complex64,
    theta0: f64,
) -> Complex64 {
    (w * log_cut(z, theta0)).exp()
}

// Keeping the wrapped angle and the integer number of turns separate avoids
// reducing an ever larger accumulated angle after each sample.
fn phase_sequence(samples: &[Complex64]) -> impl Iterator<Item = f64> + '_ {
    let mut previous = 0.0;
    let mut turns = 0.0_f64;
    samples.iter().enumerate().map(move |(index, point)| {
        let angle = point.arg();
        let delta = angle - previous;
        if index != 0 {
            if delta <= -PI {
                turns += 1.0;
            } else if delta > PI {
                turns -= 1.0;
            }
        }
        if angle.is_nan() {
            turns = f64::NAN;
        }
        previous = angle;
        if turns == 0.0 {
            angle
        } else {
            turns.mul_add(TAU, angle)
        }
    })
}

/// Unwraps the phase `arg(samples[k])` continuously along a sampled path.
/// The increment convention is `(-π, π]`.
#[must_use]
pub fn unwrap_phase(samples: &[Complex64]) -> Vec<f64> {
    phase_sequence(samples).collect()
}

/// Tracks a continuous branch of `z^(1/n)` along `points` starting on `initial_branch`.
#[must_use]
pub fn analytic_continuation_root(
    points: &[Complex64],
    n: u32,
    initial_branch: i32,
) -> Vec<Complex64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return points.to_vec();
    }
    let degree = f64::from(n);
    let inv_degree = degree.recip();
    let branch = i64::from(initial_branch).rem_euclid(i64::from(n));
    let offset = TAU * f64::from(u32::try_from(branch).unwrap_or_default()) * inv_degree;
    points
        .iter()
        .zip(phase_sequence(points))
        .map(|(&point, phase)| {
            let magnitude = point.norm();
            let scale = point.re.abs().max(point.im.abs());
            let radius = if magnitude.is_normal() || !point.is_finite() || scale == 0.0 {
                magnitude.powf(inv_degree)
            } else {
                scale.powf(inv_degree) * (point / scale).norm().powf(inv_degree)
            };
            Complex64::from_polar(radius, phase.mul_add(inv_degree, offset))
        })
        .collect()
}

/// Topological invariants of a branched covering `P(z, w) = 0` over the Riemann sphere.
#[derive(Clone, Debug, PartialEq)]
pub struct RiemannSurfaceInfo {
    /// Degree `n` in `w` (number of sheets).
    pub degree: usize,
    /// Finite branch points in `C`.
    pub branch_points: Vec<Complex64>,
    /// Whether `z = infinity` is a branch point.
    pub has_branch_at_infinity: bool,
    /// Ramification cycle lengths above each finite branch point (and infinity if branched).
    pub ramification_indices: Vec<Vec<usize>>,
    /// Monodromy permutation generators (each a bijection on `0..degree`).
    pub monodromy_generators: Vec<Vec<usize>>,
    /// Topological genus `g` by the Riemann-Hurwitz formula.
    pub genus: usize,
}

/// Topological invariants for a hyperelliptic curve `w^2 = prod_{k=1}^d (z - r_k)`
/// assuming the `roots` are distinct.
#[must_use]
pub fn algebraic_riemann_surface_hyperelliptic(roots: &[Complex64]) -> RiemannSurfaceInfo {
    algebraic_riemann_surface_superelliptic(roots, 2)
}

/// Topological invariants for a superelliptic curve `w^n = prod_{k=1}^d (z - r_k)`.
///
/// Roots within `1e-12` count once when computing the genus.
/// The generators use positive local loops, including a loop in the coordinate
/// `1/z` at infinity. The finite permutations followed by the infinity
/// permutation compose to the identity.
#[must_use]
pub fn algebraic_riemann_surface_superelliptic(
    roots: &[Complex64],
    n: usize,
) -> RiemannSurfaceInfo {
    let n = n.max(2);
    let mut distinct: Vec<Complex64> = Vec::with_capacity(roots.len());
    for &r in roots {
        if !distinct.iter().any(|&u| {
            let difference = (u - r) / 1e-12;
            difference.norm_sqr() < 1.0
        }) {
            distinct.push(r);
        }
    }
    let d = distinct.len();
    if d == 0 {
        return RiemannSurfaceInfo {
            degree: n,
            branch_points: Vec::new(),
            has_branch_at_infinity: false,
            ramification_indices: Vec::new(),
            monodromy_generators: Vec::new(),
            genus: 0,
        };
    }
    let g_nd = num_integer::gcd(n, d);
    let inf_cycle_len = n / g_nd;
    let has_branch_at_infinity = inf_cycle_len > 1;

    let mut cycle_n = Vec::with_capacity(n);
    for k in 0..n {
        cycle_n.push((k + 1) % n);
    }

    let mut ramification_indices = vec![vec![n]; d];
    let mut monodromy_generators = vec![cycle_n; d];
    let mut total_ramification = d.saturating_mul(n.saturating_sub(1));

    if has_branch_at_infinity {
        let inf_cycles = vec![inf_cycle_len; g_nd];
        ramification_indices.push(inf_cycles);
        total_ramification = total_ramification.saturating_add(n.saturating_sub(g_nd));

        let mut inf_perm = vec![0_usize; n];
        let shift = (n - d % n) % n;
        for (idx, slot) in inf_perm.iter_mut().enumerate() {
            *slot = if idx < n - shift {
                idx + shift
            } else {
                idx - (n - shift)
            };
        }
        monodromy_generators.push(inf_perm);
    }

    // Riemann-Hurwitz: 2g - 2 = -2n + sum (e_p - 1) => g = 1 - n + total_ramification / 2.
    let genus = 1_usize
        .saturating_add(total_ramification / 2)
        .saturating_sub(n);

    RiemannSurfaceInfo {
        degree: n,
        branch_points: distinct,
        has_branch_at_infinity,
        ramification_indices,
        monodromy_generators,
        genus,
    }
}

fn fiber_coefficients(
    polynomials: &[Vec<Complex64>],
    z: Complex64,
) -> Option<crate::kernels::reuse::Scratch<'static, Complex64>> {
    let mut coefficients = crate::kernels::reuse::fft_scratch_for(polynomials.len());
    coefficients.extend(polynomials.iter().map(|row| {
        row.iter()
            .rev()
            .fold(Complex64::ZERO, |tail, &term| tail * z + term)
    }));
    let scale = coefficients.iter().fold(0.0_f64, |acc, value| {
        acc.max(value.re.abs()).max(value.im.abs())
    });
    if !scale.is_finite() || scale <= 0.0 || coefficients.iter().any(|value| !value.is_finite()) {
        return None;
    }
    for coefficient in coefficients.iter_mut() {
        *coefficient /= scale;
    }
    if coefficients
        .last()
        .is_none_or(|&value| value == Complex64::ZERO)
    {
        return None;
    }
    Some(coefficients)
}

fn fiber_value(
    coefficients: &[Complex64],
    w: Complex64,
) -> Option<(Complex64, Complex64, f64)> {
    let radius = w.norm();
    let mut value = Complex64::ZERO;
    let mut derivative = Complex64::ZERO;
    let mut bound: f64 = 0.0;
    for &term in coefficients.iter().rev() {
        derivative = derivative * w + value;
        value = value * w + term;
        bound = bound.mul_add(radius, term.norm());
    }
    (value.is_finite() && derivative.is_finite() && bound.is_finite())
        .then_some((value, derivative, bound))
}

fn tracking_divide(
    a: Complex64,
    b: Complex64,
) -> Complex64 {
    use crate::kernels::complex_scalar::Complex;
    let ratio = Complex::new(a.re, a.im).div(Complex::new(b.re, b.im));
    Complex64::new(ratio.re, ratio.im)
}

fn minimum_separation(roots: &[Complex64]) -> f64 {
    roots
        .iter()
        .enumerate()
        .map(|(i, root)| {
            roots
                .iter()
                .skip(i + 1)
                .map(|other| (other - root).norm())
                .fold(f64::INFINITY, f64::min)
        })
        .fold(f64::INFINITY, f64::min)
}

fn continue_fiber(
    coefficients: &[Complex64],
    roots: &[Complex64],
) -> Option<Vec<Complex64>> {
    if let [constant, linear, quadratic] = coefficients
        && *linear == Complex64::ZERO
        && let [old_first, old_second] = roots
    {
        let root = tracking_divide(-constant, *quadratic).sqrt();
        if !root.is_finite() {
            return None;
        }
        let direct = (root - old_first).norm_sqr() + (-root - old_second).norm_sqr();
        let reverse = (-root - old_first).norm_sqr() + (root - old_second).norm_sqr();
        return if direct <= reverse {
            Some(vec![root, -root])
        } else {
            Some(vec![-root, root])
        };
    }
    let mut next = roots.to_vec();
    for root in &mut next {
        for _ in 0..16 {
            let (p, dp, bound) = fiber_value(coefficients, *root)?;
            if p.norm() <= 1e-13 * bound {
                break;
            }
            if dp == Complex64::ZERO {
                return None;
            }
            *root -= tracking_divide(p, dp);
        }
        let (p, _, bound) = fiber_value(coefficients, *root)?;
        if p.norm() > 1e-11 * bound {
            return None;
        }
    }
    Some(next)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MonodromyError {
    Inputs,
    Resolution,
    Budget,
    Cancelled,
}

fn track_fiber_segment(
    poly: &[Vec<Complex64>],
    start: Complex64,
    end: Complex64,
    sheets: &[Complex64],
    depth: u32,
    budget: &mut usize,
) -> Result<Vec<Complex64>, MonodromyError> {
    if crate::cancellation::is_cancelled() {
        return Err(MonodromyError::Cancelled);
    }
    *budget = budget.checked_sub(1).ok_or(MonodromyError::Budget)?;
    if let Some(coefficients) = fiber_coefficients(poly, end)
        && let Some(next) = continue_fiber(&coefficients, sheets)
    {
        let max_move = next
            .iter()
            .zip(sheets)
            .map(|(a, b)| (a - b).norm())
            .fold(0.0_f64, f64::max);
        let min_sep = minimum_separation(&next).min(minimum_separation(sheets));
        let extent = next
            .iter()
            .map(|value| value.norm())
            .fold(1.0_f64, f64::max);
        if min_sep > 128.0 * f64::EPSILON * extent && max_move <= 0.35 * min_sep {
            return Ok(next);
        }
    }
    if depth >= 12 {
        return Err(MonodromyError::Resolution);
    }
    let middle = Complex64::new(start.re.midpoint(end.re), start.im.midpoint(end.im));
    if middle == start || middle == end {
        return Err(MonodromyError::Resolution);
    }
    let halfway = track_fiber_segment(poly, start, middle, sheets, depth + 1, budget)?;
    track_fiber_segment(poly, middle, end, &halfway, depth + 1, budget)
}

/// Returns the numerical monodromy permutation of `P(z,w)=0` along a closed
/// piecewise-linear loop.
///
/// The final vertex is connected to the first one. Unresolved fibers, poles
/// and exhausted subdivision return `None`.
#[must_use]
pub fn monodromy_around_loop(
    poly_coeffs_in_w: &[Vec<Complex64>],
    loop_pts: &[Complex64],
) -> Option<Vec<usize>> {
    let budget = loop_pts.len().saturating_mul(64).min(1_000_000);
    monodromy_checked(poly_coeffs_in_w, loop_pts, budget).ok()
}

pub(crate) fn monodromy_checked(
    polynomial: &[Vec<Complex64>],
    path: &[Complex64],
    mut budget: usize,
) -> Result<Vec<usize>, MonodromyError> {
    if polynomial.len() < 2
        || path.len() < 2
        || path.iter().any(|z| !z.is_finite())
        || polynomial.iter().flatten().any(|a| !a.is_finite())
    {
        return Err(MonodromyError::Inputs);
    }
    let start = path[0];
    let first = fiber_coefficients(polynomial, start).ok_or(MonodromyError::Resolution)?;
    let mut initial = crate::kernels::complex_roots::aberth_ehrlich(&first, 160, 1e-13)
        .map_err(|_| MonodromyError::Resolution)?;
    let extent = initial
        .iter()
        .map(|value| value.norm())
        .fold(1.0_f64, f64::max);
    if initial.len() + 1 != polynomial.len()
        || minimum_separation(&initial) <= 128.0 * f64::EPSILON * extent
    {
        return Err(MonodromyError::Resolution);
    }
    initial.sort_by(|a, b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
    let mut tracked = initial.clone();
    let mut previous = start;
    for next in path.iter().skip(1).copied().chain(std::iter::once(start)) {
        if previous == next {
            continue;
        }
        let middle = Complex64::new(previous.re.midpoint(next.re), previous.im.midpoint(next.im));
        tracked = track_fiber_segment(polynomial, previous, middle, &tracked, 0, &mut budget)?;
        tracked = track_fiber_segment(polynomial, middle, next, &tracked, 0, &mut budget)?;
        previous = next;
    }
    let mut used = vec![false; initial.len()];
    let separation = minimum_separation(&initial);
    tracked
        .into_iter()
        .map(|end| {
            let (index, error) = initial
                .iter()
                .enumerate()
                .map(|(i, &value)| (i, (value - end).norm()))
                .min_by(|(_, a), (_, b)| a.total_cmp(b))
                .ok_or(MonodromyError::Resolution)?;
            if used[index] || error > (1e-9 * extent).min(0.25 * separation) {
                return Err(MonodromyError::Resolution);
            }
            used[index] = true;
            Ok(index)
        })
        .collect()
}

/// Numerical branch-cut jump `lim_{eps -> 0+} (f(p + i*eps) - f(p - i*eps))`
/// using Richardson extrapolation across `eps` and `eps / 2`.
pub fn branch_cut_jump_numeric<F>(
    f: F,
    p: Complex64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    branch_jump_checked(f, p, Complex64::new(0.0, 1e-6))
        .unwrap_or(Complex64::new(f64::NAN, f64::NAN))
}

pub(crate) fn branch_jump_checked(
    f: impl Fn(Complex64) -> Complex64,
    point: Complex64,
    normal_offset: Complex64,
) -> Option<Complex64> {
    let half_offset = 0.5 * normal_offset;
    let points = [
        point + normal_offset,
        point - normal_offset,
        point + half_offset,
        point - half_offset,
    ];
    if !point.is_finite()
        || points.iter().any(|&v| !v.is_finite() || v == point)
        || points
            .iter()
            .enumerate()
            .any(|(i, value)| points.iter().skip(i + 1).any(|other| other == value))
    {
        return None;
    }
    let [above, below, above_half, below_half] = points.map(f);
    let full = above - below;
    let half = above_half - below_half;
    // A general analytic lip has a first-order normal term. Two grids remove
    // that term without assuming an even expansion across the cut.
    let estimate = half + (half - full);
    estimate.is_finite().then_some(estimate)
}
