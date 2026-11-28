//! Conformal mappings, Möbius geometry, Joukowsky/Kármán–Trefftz/Koebe maps,
//! Schwarz–Christoffel upper half-plane solver, Szegő/Kerzman–Stein disk mapping,
//! and Bergman kernel utilities.

use crate::kernels::complex::Mobius;
use crate::kernels::complex::complex_derivative;
use crate::kernels::complex::differentiation;
use crate::kernels::complex_scalar::Complex as ScalarComplex;
use crate::kernels::contour_deform::ContourIntegral;
use crate::kernels::contour_deform::ContourQuadratureError;
use crate::kernels::contour_deform::ContourSegment;
use crate::kernels::contour_deform::integrate_segment_checked;
use crate::kernels::dense::solve_linear_rows;
use crate::kernels::geometry::segment_intersection_2d;
use crate::kernels::integrate::gauss_kronrod_any;
use crate::kernels::quadrature::GkError;
use crate::kernels::quadrature::GkRule;
use crate::kernels::quadrature::gk_panel_complex;
use num_complex::Complex;
use std::f64::consts::PI;

/// `cross_ratio(z1, z2, z3, z4)` computes `((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))`.
#[must_use]
pub fn cross_ratio(
    z1: Complex<f64>,
    z2: Complex<f64>,
    z3: Complex<f64>,
    z4: Complex<f64>,
) -> Complex<f64> {
    ((z1 - z3) * (z2 - z4)) / ((z1 - z4) * (z2 - z3))
}

fn mobius_to_01inf(
    z1: Complex<f64>,
    z2: Complex<f64>,
    z3: Complex<f64>,
) -> Mobius {
    // T(z) = ((z - z1)(z2 - z3)) / ((z - z3)(z2 - z1))
    let k1 = z2 - z3;
    let k2 = z2 - z1;
    Mobius {
        a: k1,
        b: -z1 * k1,
        c: k2,
        d: -z3 * k2,
    }
}

// A common binary scale preserves exact Float values until the result crosses
// a binary64 endpoint. Splitting the exponent keeps every multiplier finite.
fn scale_power_two(
    mut value: Complex<f64>,
    mut exponent: i32,
) -> Complex<f64> {
    while exponent != 0 {
        let shift = exponent.clamp(-512, 512);
        value *= 2.0_f64.powi(shift);
        exponent -= shift;
    }
    value
}

fn projective_exponent(value: Complex<f64>) -> i32 {
    value.re.abs().max(value.im.abs()).log2().floor() as i32
}

fn projective_chart(points: [Complex<f64>; 3]) -> Result<([Complex<f64>; 3], i32), String> {
    if points.iter().any(|v| !v.is_finite()) {
        return Err("Möbius points must be finite".into());
    }
    let extent = points
        .iter()
        .fold(0.0_f64, |acc, v| acc.max(v.re.abs()).max(v.im.abs()));
    if extent <= 0.0 {
        return Err("Möbius points must be distinct".into());
    }
    let exponent = projective_exponent(Complex::from(extent));
    let [a, b, c] = points.map(|v| scale_power_two(v, -exponent));
    if a == b || b == c || a == c {
        return Err("Möbius points must remain distinct in the normalized chart".into());
    }
    Ok(([a, b, c], exponent))
}

/// Constructs the unique Möbius transformation mapping `(z1, z2, z3) -> (w1, w2, w3)`.
///
/// # Errors
///
/// Rejects nonfinite or repeated points and coefficients outside the representable projective range.
pub fn mobius_from_three_points(
    z1: Complex<f64>,
    z2: Complex<f64>,
    z3: Complex<f64>,
    w1: Complex<f64>,
    w2: Complex<f64>,
    w3: Complex<f64>,
) -> Result<Mobius, String> {
    let ([a, b, c], source) = projective_chart([z1, z2, z3])?;
    let ([x, y, z], target) = projective_chart([w1, w2, w3])?;
    let chart = mobius_to_01inf(x, y, z)
        .inverse()
        .compose(&mobius_to_01inf(a, b, c));
    let coeffs = [chart.a, chart.b, chart.c, chart.d];
    let exponents = [target, source + target, 0, source];
    let (minimum, maximum) = coeffs
        .iter()
        .zip(exponents)
        .filter(|(v, _)| **v != Complex::default())
        .map(|(&v, exponent)| projective_exponent(v) + exponent)
        .fold((i32::MAX, i32::MIN), |(lo, hi), exp| {
            (lo.min(exp), hi.max(exp))
        });
    // Choose a common power-of-two projective factor that leaves headroom
    // for the subsequent complex products. Intersect all four exponent ranges.
    let lower = maximum.saturating_sub(1016);
    let upper = minimum.saturating_add(1070);
    if maximum == i32::MIN || lower > upper {
        return Err("Möbius coefficients exceed the binary64 projective range".into());
    }
    let common = 0_i32.clamp(lower, upper);
    let [(a, ea), (b, eb), (c, ec), (d, ed)] =
        std::array::from_fn(|i| (coeffs[i], exponents[i] - common));
    let result = Mobius::new(
        scale_power_two(a, ea),
        scale_power_two(b, eb),
        scale_power_two(c, ec),
        scale_power_two(d, ed),
    );
    Ok(result)
}

/// Classification of a non-singular Möbius transformation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MobiusClass {
    /// Identity map.
    Identity,
    /// Conjugated to a rotation ($\text{tr}^2 \in [0, 4)$).
    Elliptic,
    /// Single fixed point ($\text{tr}^2 = 4$).
    Parabolic,
    /// Conjugated to a pure scaling ($\text{tr}^2 \in (4, \infty)$).
    Hyperbolic,
    /// General spiral scaling ($\text{tr}^2 \notin [0, \infty)$).
    Loxodromic,
}

#[allow(
    clippy::arithmetic_side_effects,
    reason = "normalizes the projective representative"
)]
fn normalize_mobius(m: &Mobius) -> Option<Mobius> {
    let coefficients = [m.a, m.b, m.c, m.d];
    if coefficients.iter().any(|z| !z.is_finite()) {
        return None;
    }
    let scale = coefficients
        .iter()
        .fold(0.0_f64, |max, z| max.max(z.re.abs()).max(z.im.abs()));
    (scale > 0.0).then(|| Mobius::new(m.a / scale, m.b / scale, m.c / scale, m.d / scale))
}

fn divide_projective(
    a: Complex<f64>,
    b: Complex<f64>,
) -> Complex<f64> {
    let value = ScalarComplex::new(a.re, a.im).div(ScalarComplex::new(b.re, b.im));
    Complex::new(value.re, value.im)
}

/// Classifies a Möbius transformation by its normalized trace squared $\tau = (a+d)^2 / (ad-bc)$.
///
/// # Errors
///
/// Returns an error for a degenerate transformation.
#[allow(
    clippy::arithmetic_side_effects,
    reason = "forms the projective trace invariant"
)]
pub fn classify_mobius(m: &Mobius) -> Result<MobiusClass, String> {
    let m = normalize_mobius(m).ok_or_else(|| "Invalid Möbius coefficients".to_owned())?;
    let det = m.a * m.d - m.b * m.c;
    if det.norm() < 1e-15 {
        return Err("Degenerate Möbius transformation".to_string());
    }
    let zero = Complex::new(0.0, 0.0);
    if m.b == zero && m.c == zero && m.a == m.d {
        return Ok(MobiusClass::Identity);
    }
    let tr = m.a + m.d;
    let tau = (tr * tr) / det;
    if (tau - Complex::new(4.0, 0.0)).norm() < 1e-9 {
        return Ok(MobiusClass::Parabolic);
    }
    if tau.im.abs() < 1e-9 {
        if tau.re >= 0.0 && tau.re < 4.0 {
            return Ok(MobiusClass::Elliptic);
        }
        if tau.re > 4.0 {
            return Ok(MobiusClass::Hyperbolic);
        }
    }
    Ok(MobiusClass::Loxodromic)
}

/// Returns the finite fixed points of $m(z) = z$, i.e. roots of $c z^2 + (d - a) z - b = 0$.
#[must_use]
#[allow(
    clippy::arithmetic_side_effects,
    reason = "solves the complex fixed-point equation"
)]
pub fn mobius_fixed_points(m: &Mobius) -> Vec<Complex<f64>> {
    let Some(m) = normalize_mobius(m) else {
        return Vec::new();
    };
    let zero = Complex::new(0.0, 0.0);
    let d_minus_a = m.d - m.a;
    if m.c == zero {
        return if d_minus_a == zero {
            Vec::new()
        } else {
            vec![divide_projective(m.b, d_minus_a)]
        };
    }
    let half_difference = 0.5 * d_minus_a;
    let discriminant = (half_difference * half_difference + m.c * m.b).sqrt();
    if discriminant == zero {
        return vec![divide_projective(-half_difference, m.c)];
    }
    let plus = -half_difference + discriminant;
    let minus = -half_difference - discriminant;
    let numerator = if plus.norm() >= minus.norm() {
        plus
    } else {
        minus
    };
    // The product of the two roots is -b/c. Divide by the stable
    // numerator directly so a tiny c does not overflow an intermediate.
    vec![
        divide_projective(numerator, m.c),
        divide_projective(-m.b, numerator),
    ]
}

/// Unit disk automorphism $f(z) = e^{i\theta} \frac{z - \alpha}{1 - \bar{\alpha} z}$ for $|\alpha| < 1$.
///
/// # Errors
///
/// Returns an error if an argument is nonfinite or the translation lies on or outside the unit circle.
pub fn disk_automorphism(
    alpha: Complex<f64>,
    theta: f64,
) -> Result<Mobius, String> {
    if !alpha.is_finite() || !theta.is_finite() {
        return Err("Disk automorphism requires finite arguments".into());
    }
    if alpha.norm() >= 1.0 {
        return Err("Disk automorphism requires |alpha| < 1".to_string());
    }
    let rot = Complex::from_polar(1.0, theta);
    Ok(Mobius {
        a: rot,
        b: -rot * alpha,
        c: -alpha.conj(),
        d: Complex::new(1.0, 0.0),
    })
}

/// `joukowsky(z)` returns the Joukowsky airfoil transformation `(z + 1/z)/2`.
#[must_use]
pub fn joukowsky(z: Complex<f64>) -> Complex<f64> {
    let reciprocal_half = ScalarComplex::new(0.5, 0.0).div(ScalarComplex::new(z.re, z.im));
    z * 0.5 + Complex::new(reciprocal_half.re, reciprocal_half.im)
}

/// Two branches of the inverse Joukowsky map $w \pm \sqrt{w^2 - 1}$.
#[must_use]
pub fn joukowsky_inverse(w: Complex<f64>) -> (Complex<f64>, Complex<f64>) {
    let scale = w.re.abs().max(w.im.abs());
    // Form the unscaled discriminant near the branch points at ±1.
    if !scale.is_finite() || scale <= 8.0 {
        let rad = (w * w - Complex::new(1.0, 0.0)).sqrt();
        return (w + rad, w - rad);
    }
    let normalized = w / scale;
    let reciprocal_scale = scale.recip();
    let rad = (normalized * normalized - reciprocal_scale * reciprocal_scale).sqrt();
    let plus = normalized + rad;
    let minus = normalized - rad;
    // The roots have product one. Take the reciprocal of the larger
    // root before undoing the scaling so the smaller one survives even
    // if the larger one overflows.
    let plus_is_larger = plus.re.abs().max(plus.im.abs()) >= minus.re.abs().max(minus.im.abs());
    let large = if plus_is_larger {
        plus
    } else {
        minus
    };
    let inverse = ScalarComplex::new(large.re, large.im).recip();
    let small = Complex::new(inverse.re, inverse.im) * reciprocal_scale;
    let large = large * scale;
    if plus_is_larger {
        (large, small)
    } else {
        (small, large)
    }
}

/// `karman_trefftz(z, exponent)` evaluates the Kármán–Trefftz map with principal complex powers.
/// For exponent `n`, the map is `n * ((z+1)^n + (z-1)^n) / ((z+1)^n - (z-1)^n)`.
#[must_use]
pub fn karman_trefftz(
    z: Complex<f64>,
    n: f64,
) -> Complex<f64> {
    let zp = (z + Complex::new(1.0, 0.0)).powf(n);
    let zm = (z - Complex::new(1.0, 0.0)).powf(n);
    Complex::new(n, 0.0) * (zp + zm) / (zp - zm)
}

/// `koebe_map(z)` returns the Koebe function `z / (1 - z)^2`.
#[must_use]
pub fn koebe_map(z: Complex<f64>) -> Complex<f64> {
    let one_m_z = Complex::new(1.0, 0.0) - z;
    z / (one_m_z * one_m_z)
}

/// Numerically verifies conformality at `z`: Cauchy–Riemann holds and $|f'(z)| > \text{tol}$.
pub fn is_conformal_numeric<F>(
    f: &F,
    z: Complex<f64>,
    h: f64,
    tol: f64,
) -> bool
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    if !tol.is_finite() || tol < 0.0 {
        return false;
    }
    differentiation::wirtinger(f, z, h)
        .is_some_and(|value| value.residual(1.0) <= tol && value.dz.norm() > tol)
}

/// `cayley_transform(z)` returns the Möbius map `(z - i) / (z + i)`.
/// Its pole is at `-i`. On the upper half-plane its image is the open unit disk.
#[must_use]
pub fn cayley_transform(z: Complex<f64>) -> Complex<f64> {
    crate::kernels::complex_coords::uhp_to_disk(z).unwrap_or(Complex::new(f64::NAN, f64::NAN))
}

/// Returns the local conformal magnification $|f'(z)|$ and rotation angle $\arg f'(z)$.
pub fn conformal_scale_and_rotation<F>(
    f: &F,
    z: Complex<f64>,
    h: f64,
) -> (f64, f64)
where
    F: Fn(Complex<f64>) -> Complex<f64>,
{
    let df = if h.is_finite() && h > 0.0 {
        crate::kernels::complex::wirtinger_dz_numeric(f, z, h)
    } else {
        complex_derivative(f, z)
    };
    (df.norm(), df.arg())
}

fn integrate_pure_segment(
    evaluator: impl Fn(Complex<f64>) -> Complex<f64>,
    start: Complex<f64>,
    end: Complex<f64>,
    tolerance: f64,
    panels: usize,
) -> Result<ContourIntegral, ContourQuadratureError> {
    let segment = ContourSegment::Line { start, end };
    let delta = end - start;
    if !delta.is_finite()
        || delta == Complex::default()
        || !(1..=4096).contains(&panels)
        || !tolerance.is_finite()
        || tolerance <= 0.0
        || crate::cancellation::is_cancelled()
    {
        return integrate_segment_checked(evaluator, segment, tolerance, panels);
    }
    let (value, [re_error, im_error], calls) = gk_panel_complex(0.0, 1.0, GkRule::G7K15, |t| {
        evaluator(start + delta * t) * delta
    });
    if crate::cancellation::is_cancelled() {
        return Err(GkError::Cancelled.into());
    }
    let error = re_error.hypot(im_error);
    if !value.is_finite() || !error.is_finite() {
        return Err(GkError::NonFiniteEstimate.into());
    }
    let converged = re_error <= tolerance * 0.5 && im_error <= tolerance * 0.5;
    if converged || panels == 1 {
        return Ok(ContourIntegral {
            value,
            error,
            evaluations: calls,
            converged,
        });
    }
    let mut refined = integrate_segment_checked(evaluator, segment, tolerance, panels - 1)?;
    refined.evaluations += calls;
    Ok(refined)
}

fn integrate_complex_01<F>(
    f: F,
    tol: f64,
    max_evals: usize,
) -> Complex<f64>
where
    F: Fn(f64) -> Complex<f64>,
{
    let re = gauss_kronrod_any(|t| f(t).re, 0.0, 1.0, tol, max_evals).value;
    let im = gauss_kronrod_any(|t| f(t).im, 0.0, 1.0, tol, max_evals).value;
    Complex::new(re, im)
}

/// Schwarz–Christoffel map from the upper half-plane $\mathbb{H}^+$ to a polygon:
/// $f(z) = \text{shift} + \text{scale} \int_{`z_0}^z` \prod_{k=1}^m (\zeta - `x_k)^{\alpha_k` - 1} d\zeta$.
#[derive(Debug, Clone)]
pub struct SchwarzChristoffelMap {
    /// Prevertices $`x_1` < `x_2` < \dots < `x_m`$ on the real axis.
    pub prevertices: Vec<f64>,
    /// Interior angle factors `\alpha_k \in (0, 2)` (so exponent is `\alpha_k - 1`).
    pub interior_angles_pi: Vec<f64>,
    /// Complex scaling factor $C$.
    pub scale: Complex<f64>,
    /// Additive offset $`f(z_0)`$ at basepoint $`z_0` = i$.
    pub shift: Complex<f64>,
}

impl SchwarzChristoffelMap {
    fn integrand(
        &self,
        zeta: Complex<f64>,
    ) -> Complex<f64> {
        let mut prod = Complex::new(1.0, 0.0);
        for (&xk, &alpha) in self.prevertices.iter().zip(self.interior_angles_pi.iter()) {
            let diff = zeta - Complex::new(xk, 0.0);
            // Principal branch for Im(zeta) > 0 keeps diff in upper half-plane
            prod *= diff.powf(alpha - 1.0);
        }
        prod
    }

    /// Evaluates the Schwarz–Christoffel integral from `z₀ = i` to a point
    /// in the open upper half-plane using a common complex sample grid.
    /// # Errors
    /// Rejects malformed parameters, nonfinite results, cancellation, or a
    /// quadrature estimate that misses the default tolerance.
    pub fn eval(
        &self,
        z: Complex<f64>,
    ) -> Result<Complex<f64>, String> {
        let coefficient_scale = self.scale.re.abs().max(self.scale.im.abs()).max(1e-100);
        let result = self.eval_checked(z, 1e-10 * coefficient_scale, 128)?;
        if result.converged {
            Ok(result.value)
        } else {
            Err("Schwarz–Christoffel integration did not meet the default tolerance".into())
        }
    }

    pub(crate) fn eval_checked(
        &self,
        z: Complex<f64>,
        tolerance: f64,
        panels: usize,
    ) -> Result<ContourIntegral, String> {
        if !z.is_finite()
            || z.im <= 0.0
            || !self.scale.is_finite()
            || self.scale == Complex::default()
            || !self.shift.is_finite()
            || self.prevertices.len() != self.interior_angles_pi.len()
            || self.prevertices.iter().any(|x| !x.is_finite())
            || self.prevertices.windows(2).any(|pair| pair[0] >= pair[1])
            || self
                .interior_angles_pi
                .iter()
                .any(|&angle| !(angle > 0.0 && angle < 2.0))
        {
            return Err("Schwarz–Christoffel evaluation requires ordered finite prevertices, interior angle factors in (0,2), a finite nonzero scale, and an interior point".into());
        }

        let start = Complex::new(0.0, 1.0);
        let mut estimate = integrate_pure_segment(
            |point| self.scale * self.integrand(point),
            start,
            z,
            tolerance,
            panels,
        )
        .map_err(|error| error.to_string())?;
        estimate.value += self.shift;
        if !estimate.value.is_finite() {
            return Err("Schwarz–Christoffel value exceeds Float range".into());
        }
        Ok(estimate)
    }

    fn segment_integral(
        &self,
        a: f64,
        b: f64,
    ) -> Complex<f64> {
        // Compound tanh-sinh-like / double-sin^2 substitution for acute endpoint singularities:
        // x(u) = a + span * sin^2(pi/2 * sin^2(pi/2 * u)) clusters doubly near endpoints.
        let span = b - a;
        integrate_complex_01(
            |u| {
                let s1 = (0.5 * PI * u).sin();
                let c1 = (0.5 * PI * u).cos();
                let v = s1 * s1;
                let dv = PI * s1 * c1;
                let s2 = (0.5 * PI * v).sin();
                let c2 = (0.5 * PI * v).cos();
                let x = a + span * s2 * s2;
                let dx = span * PI * s2 * c2 * dv;
                self.integrand(Complex::new(x, 1e-14)) * dx
            },
            1e-10,
            400,
        )
    }
}

fn polygon_interior_angles_pi(vertices: &[Complex<f64>]) -> Result<Vec<f64>, String> {
    let n = vertices.len();
    let mut angles = Vec::with_capacity(n);
    let edge = |i: usize| {
        let a = vertices[i];
        let b = vertices[(i + 1) % n];
        ([a.re, a.im], [b.re, b.im])
    };
    let mut area_twice = 0.0;
    for i in 0..n {
        crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
        let prev = vertices[(i + n - 1) % n];
        let curr = vertices[i];
        let next = vertices[(i + 1) % n];
        let incoming = curr - prev;
        let outgoing = next - curr;
        if incoming == Complex::default() || outgoing == Complex::default() {
            return Err("Schwarz–Christoffel polygon has a repeated edge endpoint".into());
        }
        area_twice += curr.re * next.im - curr.im * next.re;
        let angle = 1.0 - divide_projective(outgoing, incoming).arg() / PI;
        if !angle.is_finite() || angle <= 0.0 || angle >= 2.0 {
            return Err("Schwarz–Christoffel polygon has an unresolved interior angle".into());
        }
        angles.push(angle);
        let (a, b) = edge(i);
        for j in (i + 2)..n {
            if i == 0 && j + 1 == n {
                continue;
            }
            let (c, d) = edge(j);
            if segment_intersection_2d(a, b, c, d).is_some() {
                return Err("Schwarz–Christoffel polygon edges intersect".into());
            }
        }
    }
    if !(area_twice.is_finite() && area_twice > 0.0) {
        return Err("Schwarz–Christoffel polygon needs positive counterclockwise area".into());
    }
    Ok(angles)
}

/// Solves for the prevertices of a Schwarz–Christoffel map.
///
/// Maps the upper half-plane to the simple counterclockwise polygon `polygon_vertices`. Fixes
/// the first, second and last prevertices at `(-1, 0, 1)` and finds the remaining `n - 3` by damped Newton
/// iteration.
///
/// # Errors
///
/// Returns an error for an invalid polygon or when the quadrature or Newton solve fails.
pub fn solve_sc_uhp(
    polygon_vertices: &[Complex<f64>],
    max_iter: usize,
    tol: f64,
) -> Result<SchwarzChristoffelMap, String> {
    let n = polygon_vertices.len();
    if n < 3 || !tol.is_finite() || tol <= 0.0 || polygon_vertices.iter().any(|v| !v.is_finite()) {
        return Err(
            "Schwarz–Christoffel solve needs a finite polygon and positive finite tolerance".into(),
        );
    }
    let origin = polygon_vertices[0];
    let extent = polygon_vertices
        .iter()
        .fold(0.0_f64, |acc, v| acc.max(v.re.abs()).max(v.im.abs()));
    if extent <= 0.0 {
        return Err("Schwarz–Christoffel polygon has zero extent".into());
    }
    let local_origin = origin / extent;
    let polygon_vertices: Vec<_> = polygon_vertices
        .iter()
        .map(|&v| v / extent - local_origin)
        .collect();
    let alphas = polygon_interior_angles_pi(&polygon_vertices)?;
    let mut prevertices = Vec::with_capacity(n);
    prevertices.push(-1.0);
    prevertices.push(0.0);
    for k in 2..n.saturating_sub(1) {
        prevertices.push(((k - 1) as f64) / ((n - 2) as f64));
    }
    prevertices.push(1.0);

    let mut target_lengths = Vec::with_capacity(n.saturating_sub(1));
    for i in 0..n.saturating_sub(1) {
        let a = *polygon_vertices.get(i).ok_or_else(|| "Index".to_string())?;
        let b = *polygon_vertices
            .get(i + 1)
            .ok_or_else(|| "Index".to_string())?;
        target_lengths.push((b - a).norm());
    }
    let len0 = *target_lengths
        .first()
        .ok_or_else(|| "Missing edge".to_string())?;
    let target_ratios: Vec<f64> = target_lengths.iter().skip(1).map(|&l| l / len0).collect();

    // Unconstrained coordinates u_k parameterizing 0 < x_2 < ... < x_{n-2} < 1 via softmax gaps
    let m_free = n.saturating_sub(3);
    let mut u = vec![0.0_f64; m_free];

    let build_prevertices = |u_vec: &[f64]| -> Vec<f64> {
        let mut xs = vec![-1.0, 0.0];
        let mut gaps = Vec::with_capacity(u_vec.len() + 1);
        let mut total = 1.0_f64; // last gap has exp(0) = 1
        for &uk in u_vec {
            let g = uk.clamp(-12.0, 12.0).exp();
            total += g;
            gaps.push(g);
        }
        let mut acc = 0.0_f64;
        for g in gaps {
            acc += g / total;
            xs.push(acc);
        }
        xs.push(1.0);
        xs
    };

    let eval_ratios = |xs: &[f64]| -> Result<(Vec<f64>, Complex<f64>), String> {
        let sc = SchwarzChristoffelMap {
            prevertices: xs.to_vec(),
            interior_angles_pi: alphas.clone(),
            scale: Complex::new(1.0, 0.0),
            shift: Complex::new(0.0, 0.0),
        };
        let mut segs = Vec::with_capacity(n.saturating_sub(1));
        for i in 0..n.saturating_sub(1) {
            let a = *xs.get(i).ok_or_else(|| "Index".to_string())?;
            let b = *xs.get(i + 1).ok_or_else(|| "Index".to_string())?;
            segs.push(sc.segment_integral(a, b));
        }
        let s0 = *segs.first().ok_or_else(|| "Missing seg".to_string())?;
        let s0_norm = s0.norm();
        if s0_norm <= 0.0 || !s0_norm.is_finite() {
            return Err("Schwarz–Christoffel reference edge integral is degenerate".into());
        }
        let ratios: Vec<f64> = segs.iter().skip(1).map(|s| s.norm() / s0_norm).collect();
        if ratios.iter().any(|ratio| !ratio.is_finite()) {
            return Err("Schwarz–Christoffel side length ratio exceeds Float range".into());
        }
        Ok((ratios, s0))
    };

    let mut converged = m_free == 0;
    if m_free > 0 {
        for _ in 0..max_iter {
            crate::cancellation::checkpoint().map_err(|error| error.to_string())?;
            let xs = build_prevertices(&u);
            let (ratios, _) = eval_ratios(&xs)?;
            let mut res = Vec::with_capacity(m_free);
            let mut norm2 = 0.0_f64;
            for k in 0..m_free {
                let r_k = *ratios.get(k).ok_or_else(|| "Ratio index".to_string())?;
                let t_k = *target_ratios
                    .get(k)
                    .ok_or_else(|| "Target index".to_string())?;
                let diff = r_k - t_k;
                norm2 += diff * diff;
                res.push(diff);
            }
            if norm2.sqrt() < tol {
                converged = true;
                break;
            }
            // Finite-difference Jacobian
            let eps = 1e-5;
            let mut jac = vec![vec![0.0_f64; m_free]; m_free];
            for col in 0..m_free {
                let mut u_pert = u.clone();
                if let Some(slot) = u_pert.get_mut(col) {
                    *slot += eps;
                }
                let xs_p = build_prevertices(&u_pert);
                let (ratios_p, _) = eval_ratios(&xs_p)?;
                for row in 0..m_free {
                    let r_p = *ratios_p.get(row).ok_or_else(|| "Ratio p".to_string())?;
                    let r_0 = *ratios.get(row).ok_or_else(|| "Ratio 0".to_string())?;
                    if let Some(cell) = jac.get_mut(row).and_then(|r| r.get_mut(col)) {
                        *cell = (r_p - r_0) / eps;
                    }
                }
            }
            let delta = solve_linear_rows(&jac, &res).map_err(|error| error.to_string())?;
            let mut alpha_step = 1.0_f64;
            while alpha_step > 1e-4 {
                let u_try: Vec<f64> = u
                    .iter()
                    .zip(delta.iter())
                    .map(|(&uk, &dk)| uk - alpha_step * dk)
                    .collect();
                let xs_try = build_prevertices(&u_try);
                let (ratios_try, _) = eval_ratios(&xs_try)?;
                let mut try_norm2 = 0.0_f64;
                for k in 0..m_free {
                    let rk = *ratios_try.get(k).unwrap_or(&0.0);
                    let tk = *target_ratios.get(k).unwrap_or(&0.0);
                    try_norm2 += (rk - tk) * (rk - tk);
                }
                if try_norm2 < norm2 {
                    u = u_try;
                    break;
                }
                alpha_step *= 0.5;
            }
        }
        prevertices = build_prevertices(&u);
    }
    if !converged {
        return Err("Schwarz-Christoffel solver failed to converge within max_iter".to_string());
    }

    let (_, s0) = eval_ratios(&prevertices)?;
    let v0 = *polygon_vertices
        .first()
        .ok_or_else(|| "Index".to_string())?;
    let v1 = *polygon_vertices.get(1).ok_or_else(|| "Index".to_string())?;
    let scale = divide_projective(v1 - v0, s0);
    let mut sc = SchwarzChristoffelMap {
        prevertices,
        interior_angles_pi: alphas,
        scale,
        shift: Complex::new(0.0, 0.0),
    };
    // Integrate from z0 = i to a point near x0 = -1 along an interior approach
    let f_at_x0 = sc.eval(Complex::new(-1.0, 1e-6))?;
    sc.scale *= extent;
    sc.shift = origin + (v0 - f_at_x0) * extent;
    if sc.scale.is_finite() && sc.shift.is_finite() && sc.scale != Complex::default() {
        Ok(sc)
    } else {
        Err("Schwarz–Christoffel physical parameters exceed Float range".into())
    }
}

#[path = "conformal_szego.rs"]
mod szego;
pub use szego::SzegoKerzmanSteinSolver;
pub use szego::szego_kerzman_stein_map;

/// Exact Bergman kernel of the unit disk: $K_{\mathbb{D}}(z, w) = \frac{1}{\pi (1 - z \bar{w})^2}$.
#[must_use]
pub fn bergman_kernel_disk(
    z: Complex<f64>,
    w: Complex<f64>,
) -> Complex<f64> {
    let denom = Complex::new(1.0, 0.0) - z * w.conj();
    Complex::new(1.0 / PI, 0.0) / (denom * denom)
}

/// Evaluates a Bergman kernel expansion `K(z, w) = \sum_{k=0}^M \phi_k(z) \overline{\phi_k(w)}`
/// from a list of orthonormal polynomial basis vectors (each in ascending power coefficients).
#[must_use]
pub fn bergman_kernel_polynomial_basis(
    orthonormal_polys: &[Vec<Complex<f64>>],
    z: Complex<f64>,
    w: Complex<f64>,
) -> Complex<f64> {
    if z == w {
        let sum = orthonormal_polys
            .iter()
            .map(|coefficients| {
                coefficients
                    .iter()
                    .rev()
                    .fold(Complex::<f64>::default(), |tail, &term| tail * z + term)
                    .norm_sqr()
            })
            .sum();
        return Complex::new(sum, 0.0);
    }
    orthonormal_polys
        .iter()
        .map(|coefficients| {
            let (a, b) = coefficients.iter().rev().fold(
                (Complex::<f64>::default(), Complex::<f64>::default()),
                |(left, right), &term| (left * z + term, right * w + term),
            );
            a * b.conj()
        })
        .sum()
}
