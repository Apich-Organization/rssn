//! Complex coordinate systems, Riemann sphere stereographic projection, and topology.
//!
//! Covers hyperbolic disk/upper half-plane metrics, custom branch-cut control,
//! continuous phase unwrapping, and algebraic Riemann surface topology.

use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;

/// Polar coordinates `(r, theta)` with `r = |z|` and `theta = arg(z) in (-pi, pi]`.
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

/// Unit complex number `cos(theta) + i * sin(theta)`.
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
    let r = z.norm();
    if r == 0.0 || !r.is_finite() {
        if z.re.is_infinite() || z.im.is_infinite() {
            let re = if z.re.is_infinite() {
                z.re.signum()
            } else {
                0.0
            };
            let im = if z.im.is_infinite() {
                z.im.signum()
            } else {
                0.0
            };
            let dir = Complex64::new(re, im);
            let n = dir.norm();
            if n > 0.0 {
                return dir / n;
            }
        }
        return Complex64::new(0.0, 0.0);
    }
    z / r
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
    let w = Complex64::new(1.0, 0.0) / z;
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
    let denom = 1.0 - z;
    if denom.abs() <= 1e-15 {
        return None;
    }
    Some(Complex64::new(x / denom, y / denom))
}

/// Spherical chordal distance on the Riemann sphere:
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
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Spherical metric conformal factor `2 / (1 + |z|^2)`.
#[must_use]
pub fn spherical_metric_factor(z: Complex64) -> f64 {
    if !z.re.is_finite() || !z.im.is_finite() {
        return 0.0;
    }
    let s = z.re.abs().max(z.im.abs());
    if s > 1e150 {
        return 0.0;
    }
    2.0 / (1.0 + abs_sq(z))
}

/// Cayley transform mapping the upper half-plane `{ Im(z) > 0 }` to the unit disk `{ |w| < 1 }`:
/// `w = (z - i) / (z + i)`. Returns `None` at `z = -i`.
#[must_use]
pub fn uhp_to_disk(z: Complex64) -> Option<Complex64> {
    let i = Complex64::new(0.0, 1.0);
    let denom = z + i;
    if denom.norm() == 0.0 {
        return None;
    }
    Some((z - i) / denom)
}

/// Inverse Cayley transform mapping the unit disk `{ |w| < 1 }` to the upper half-plane:
/// `z = i * (1 + w) / (1 - w)`. Returns `None` at `w = 1`.
#[must_use]
pub fn disk_to_uhp(w: Complex64) -> Option<Complex64> {
    let i = Complex64::new(0.0, 1.0);
    let denom = Complex64::new(1.0, 0.0) - w;
    if denom.norm() == 0.0 {
        return None;
    }
    Some(i * (Complex64::new(1.0, 0.0) + w) / denom)
}

/// Hyperbolic distance in the Poincaré unit disk `|z| < 1`.
#[must_use]
pub fn poincare_disk_distance(
    z1: Complex64,
    z2: Complex64,
) -> Option<f64> {
    if abs_sq(z1) >= 1.0 || abs_sq(z2) >= 1.0 {
        return None;
    }
    let denom = Complex64::new(1.0, 0.0) - z1.conj() * z2;
    let rho = ((z1 - z2) / denom).norm();
    if rho >= 1.0 {
        return None;
    }
    Some(2.0 * rho.atanh())
}

/// Hyperbolic distance in the upper half-plane `Im(z) > 0`.
#[must_use]
pub fn uhp_distance(
    z1: Complex64,
    z2: Complex64,
) -> Option<f64> {
    if z1.im <= 0.0 || z2.im <= 0.0 {
        return None;
    }
    let diff = z1 - z2;
    let arg = 1.0 + abs_sq(diff) / (2.0 * z1.im * z2.im);
    Some(arg.max(1.0).acosh())
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

/// Unwraps the phase `arg(samples[k])` continuously along a sampled path.
#[must_use]
pub fn unwrap_phase(samples: &[Complex64]) -> Vec<f64> {
    let mut out = Vec::with_capacity(samples.len());
    let Some(&first) = samples.first() else {
        return out;
    };
    let mut prev = first.arg();
    out.push(prev);
    for &z in samples.iter().skip(1) {
        let raw = z.arg();
        let mut delta = raw - (prev.rem_euclid(TAU));
        while delta <= -PI {
            delta += TAU;
        }
        while delta > PI {
            delta -= TAU;
        }
        prev += delta;
        out.push(prev);
    }
    out
}

/// Tracks a continuous branch of `z^(1/n)` along `points` starting on `initial_branch`.
#[must_use]
pub fn analytic_continuation_root(
    points: &[Complex64],
    n: u32,
    initial_branch: i32,
) -> Vec<Complex64> {
    if n == 0 || points.is_empty() {
        return Vec::new();
    }
    let phases = unwrap_phase(points);
    let nf = f64::from(n);
    let offset = TAU * f64::from(initial_branch);
    points
        .iter()
        .zip(phases.iter())
        .map(|(z, &phi)| Complex64::from_polar(z.norm().powf(1.0 / nf), (phi + offset) / nf))
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

const fn gcd_usize(
    mut a: usize,
    mut b: usize,
) -> usize {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Topological invariants for a hyperelliptic curve `w^2 = prod_{k=1}^d (z - r_k)`
/// assuming the `roots` are distinct.
#[must_use]
pub fn algebraic_riemann_surface_hyperelliptic(roots: &[Complex64]) -> RiemannSurfaceInfo {
    algebraic_riemann_surface_superelliptic(roots, 2)
}

/// Topological invariants for a superelliptic curve `w^n = prod_{k=1}^d (z - r_k)`
/// assuming `n >= 2` and the `roots` are distinct.
#[must_use]
pub fn algebraic_riemann_surface_superelliptic(
    roots: &[Complex64],
    n: usize,
) -> RiemannSurfaceInfo {
    let n = n.max(2);
    let d = roots.len();
    let g_nd = gcd_usize(n, d);
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
        for (idx, slot) in inf_perm.iter_mut().enumerate() {
            *slot = (idx + d) % n;
        }
        monodromy_generators.push(inf_perm);
    }

    // Riemann-Hurwitz: 2g - 2 = -2n + sum (e_p - 1) => g = 1 - n + total_ramification / 2.
    let genus = 1_usize
        .saturating_add(total_ramification / 2)
        .saturating_sub(n);

    RiemannSurfaceInfo {
        degree: n,
        branch_points: roots.to_vec(),
        has_branch_at_infinity,
        ramification_indices,
        monodromy_generators,
        genus,
    }
}

/// Evaluates `P(z, w) = sum_{j=0}^n c_j(z) w^j` and `dP/dw(z, w)`, where
/// `poly_coeffs_in_w[j]` is `[a_0, a_1, ...]` for `c_j(z) = sum_k a_k z^k`.
fn eval_bivariate_and_dw(
    poly_coeffs_in_w: &[Vec<Complex64>],
    z: Complex64,
    w: Complex64,
) -> (Complex64, Complex64) {
    let eval_c = |coeffs: &[Complex64]| -> Complex64 {
        let mut acc = Complex64::new(0.0, 0.0);
        for &c in coeffs.iter().rev() {
            acc = acc * z + c;
        }
        acc
    };
    let mut p_val = Complex64::new(0.0, 0.0);
    let mut dp_val = Complex64::new(0.0, 0.0);
    for (j, coeffs) in poly_coeffs_in_w.iter().enumerate().rev() {
        let cj = eval_c(coeffs);
        if j < poly_coeffs_in_w.len().saturating_sub(1) {
            dp_val = dp_val * w + p_val;
        }
        p_val = p_val * w + cj;
    }
    (p_val, dp_val)
}

/// Solves univariate `P(z, w) = 0` for `w` via Aberth-Ehrlich.
fn aberth_solve_fiber(
    poly_coeffs_in_w: &[Vec<Complex64>],
    z: Complex64,
) -> Option<Vec<Complex64>> {
    let deg = poly_coeffs_in_w.len().checked_sub(1)?;
    if deg == 0 {
        return None;
    }
    let mut roots = Vec::with_capacity(deg);
    for k in 0..deg {
        let angle = TAU * (k as f64 + 0.13) / (deg as f64);
        roots.push(Complex64::from_polar(1.5 + 0.1 * (k as f64), angle));
    }
    for _ in 0..120 {
        let mut max_step = 0.0_f64;
        for i in 0..deg {
            let &wi = roots.get(i)?;
            let (p, dp) = eval_bivariate_and_dw(poly_coeffs_in_w, z, wi);
            if p.norm() < 1e-13 {
                continue;
            }
            if dp.norm() < 1e-30 {
                continue;
            }
            let ratio = p / dp;
            let mut sum_inv = Complex64::new(0.0, 0.0);
            for (j, &wj) in roots.iter().enumerate() {
                if i != j {
                    let diff = wi - wj;
                    if diff.norm() > 1e-30 {
                        sum_inv += Complex64::new(1.0, 0.0) / diff;
                    }
                }
            }
            let denom = Complex64::new(1.0, 0.0) - ratio * sum_inv;
            if denom.norm() < 1e-30 {
                continue;
            }
            let delta = ratio / denom;
            if let Some(slot) = roots.get_mut(i) {
                *slot -= delta;
            }
            max_step = max_step.max(delta.norm());
        }
        if max_step < 1e-13 {
            break;
        }
    }
    Some(roots)
}

/// Tracks the `n` sheets of `P(z, w) = 0` along a closed loop `loop_pts`
/// (from `loop_pts[0]` around and back to `loop_pts[0]`) and returns the
/// induced monodromy permutation on `0..n`.
#[must_use]
pub fn monodromy_around_loop(
    poly_coeffs_in_w: &[Vec<Complex64>],
    loop_pts: &[Complex64],
) -> Option<Vec<usize>> {
    let &z0 = loop_pts.first()?;
    let start_sheets = aberth_solve_fiber(poly_coeffs_in_w, z0)?;
    let n = start_sheets.len();
    let mut tracked = start_sheets.clone();

    for &zk in loop_pts.iter().skip(1).chain(std::iter::once(&z0)) {
        for w in &mut tracked {
            for _ in 0..16 {
                let (p, dp) = eval_bivariate_and_dw(poly_coeffs_in_w, zk, *w);
                if p.norm() < 1e-12 {
                    break;
                }
                if dp.norm() < 1e-30 {
                    return None;
                }
                *w -= p / dp;
            }
        }
    }

    // Match each final tracked sheet to the closest starting sheet.
    let mut perm = vec![0_usize; n];
    let mut used = vec![false; n];
    for (i, &wt) in tracked.iter().enumerate() {
        let mut best_j = 0;
        let mut best_d = f64::INFINITY;
        for (j, &ws) in start_sheets.iter().enumerate() {
            if !used.get(j).copied().unwrap_or(true) {
                let d = (wt - ws).norm();
                if d < best_d {
                    best_d = d;
                    best_j = j;
                }
            }
        }
        if best_d > 1e-3 {
            return None;
        }
        if let Some(u) = used.get_mut(best_j) {
            *u = true;
        }
        if let Some(p) = perm.get_mut(i) {
            *p = best_j;
        }
    }
    Some(perm)
}
