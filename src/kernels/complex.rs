//! Numerical complex analysis: contour integrals, residues, the argument
//! principle, complex derivatives and Möbius transformations.
//!
//! Closed circular contours use the trapezoidal rule, which converges
//! geometrically for functions analytic in an annulus around the circle;
//! open polygonal paths use composite Simpson on each segment.

use std::cmp::Ordering;
use std::f64::consts::PI;
use std::f64::consts::TAU;

use num_complex::Complex64;
use serde::Deserialize;
use serde::Serialize;

use crate::kernels::fft;
use crate::kernels::reuse::Scratch;
use crate::kernels::reuse::fft_scratch_for;
use crate::kernels::reuse::global_pool;

const I: Complex64 = Complex64::new(0.0, 1.0);

/// The point at angle `2 pi k / n` on the circle `|z - center| = radius`.
fn on_circle(
    center: Complex64,
    radius: f64,
    k: usize,
    n: usize,
) -> (Complex64, Complex64) {
    #[allow(clippy::cast_precision_loss)]
    let theta = 2.0 * PI * k as f64 / n as f64;
    let direction = Complex64::from_polar(1.0, theta);
    (center + radius * direction, I * radius * direction)
}

/// `∮ f(z) dz` over the circle `|z - center| = radius`, counter-clockwise,
/// with `n` equally spaced nodes.
pub fn circle_integral<F>(
    f: F,
    center: Complex64,
    radius: f64,
    n: usize,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let n = n.max(8);
    let mut sum = Complex64::new(0.0, 0.0);
    for k in 0..n {
        let (z, dz_dtheta) = on_circle(center, radius, k, n);
        sum += f(z) * dz_dtheta;
    }
    #[allow(clippy::cast_precision_loss)]
    let step = 2.0 * PI / n as f64;
    sum * step
}

/// `∫ f(z) dz` along the polygon through `path`, by Simpson's rule on each
/// segment subdivided `per_segment` times.
pub fn contour_integral<F>(
    f: F,
    path: &[Complex64],
    per_segment: usize,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let m = per_segment.max(1);
    let mut total = Complex64::new(0.0, 0.0);
    for pair in path.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        #[allow(clippy::cast_precision_loss)]
        let h = (b - a) / m as f64;
        for j in 0..m {
            #[allow(clippy::cast_precision_loss)]
            let z0 = a + h * j as f64;
            let z1 = z0 + h;
            total += (f(z0) + 4.0 * f((z0 + z1) / 2.0) + f(z1)) / 6.0 * h;
        }
    }
    total
}

/// The residue of `f` at `z0`: `(1 / 2 pi i) ∮ f dz` over a small circle
/// that must enclose no other singularity.
pub fn residue<F>(
    f: F,
    z0: Complex64,
    radius: f64,
    n: usize,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    circle_integral(f, z0, radius, n) / (2.0 * PI * I)
}

// Compare a sum of products with zero. A nonzero rounded sum already has the
// exact sign by monotonic rounding of the two product magnitudes. Zero or
// NaN can also come from cancellation or a lost product, so resolve those
// cases from the exact integer significands.
fn product_sum_order(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
) -> Ordering {
    let sum = a * b + c * d;
    if sum > 0.0 {
        Ordering::Greater
    } else if sum < 0.0 {
        Ordering::Less
    } else {
        exact_product_sum_order(a, b, c, d)
    }
}

fn dyadic_parts(value: f64) -> (u64, i32) {
    let magnitude = value.abs().to_bits();
    let biased = (magnitude >> 52) as i32;
    let fraction = magnitude & 0x000f_ffff_ffff_ffff;
    let significand = fraction | (u64::from(biased != 0) << 52);
    (significand, biased.max(1) - 1075)
}

// Normalize the 106-bit product to compare by exponent and then significand.
fn product_key(
    a: f64,
    b: f64,
) -> (i32, u128) {
    let (ma, ea) = dyadic_parts(a);
    let (mb, eb) = dyadic_parts(b);
    let product = u128::from(ma) * u128::from(mb);
    if product == 0 {
        return (i32::MIN, 0);
    }
    let shift = product.leading_zeros();
    (ea + eb - shift as i32, product << shift)
}

#[cold]
fn exact_product_sum_order(
    a: f64,
    b: f64,
    c: f64,
    d: f64,
) -> Ordering {
    let left_negative = a.is_sign_negative() != b.is_sign_negative();
    let right_negative = c.is_sign_negative() != d.is_sign_negative();
    if left_negative != right_negative {
        let order = product_key(a, b).cmp(&product_key(c, d));
        return if left_negative {
            order.reverse()
        } else {
            order
        };
    }
    if (a == 0.0 || b == 0.0) && (c == 0.0 || d == 0.0) {
        Ordering::Equal
    } else if left_negative {
        Ordering::Less
    } else {
        Ordering::Greater
    }
}

// Count signed crossings of the positive horizontal ray. The half-open
// convention counts a vertex on the ray in just one adjoining edge.
fn winding_edge(
    from: Complex64,
    to: Complex64,
    limit_quarter_turn: bool,
) -> Option<i64> {
    if limit_quarter_turn && product_sum_order(from.re, to.re, from.im, to.im).is_lt() {
        return None;
    }
    let above_from = from.im > 0.0;
    let above_to = to.im > 0.0;
    if above_from == above_to {
        if from.im == 0.0 && to.im == 0.0 && from.re.is_sign_negative() != to.re.is_sign_negative()
        {
            return None;
        }
        return Some(0);
    }
    match product_sum_order(from.re, to.im, -from.im, to.re) {
        | Ordering::Equal => None,
        | Ordering::Greater if above_to => Some(1),
        | Ordering::Less if above_from => Some(-1),
        | _ => Some(0),
    }
}

fn winding(
    points: impl IntoIterator<Item = Complex64>,
    limit_quarter_turn: bool,
) -> Option<i64> {
    let mut points = points.into_iter();
    let first = points.next()?;
    if !first.is_finite() || first == Complex64::ZERO {
        return None;
    }
    let mut previous = first;
    let mut total = 0;
    for point in points {
        if !point.is_finite() || point == Complex64::ZERO {
            return None;
        }
        total += winding_edge(previous, point, limit_quarter_turn)?;
        previous = point;
    }
    Some(total + winding_edge(previous, first, limit_quarter_turn)?)
}

/// Winding number about the origin of the closed polygonal path through `points`.
///
/// Joins the last vertex to the first. Counterclockwise turns count as
/// positive. Returns `None` for an empty path, a nonfinite vertex or an edge
/// through the origin. A single nonzero finite vertex has winding number zero.
#[must_use]
pub fn winding_number(points: &[Complex64]) -> Option<i64> {
    winding(points.iter().copied(), false)
}

/// The number of zeros minus the number of poles of `f` inside a circle.
///
/// The circle is `|z - center| = radius` (argument principle), and the count
/// is the winding number of the polygon through the sampled values of `f`.
/// Samples `f` once per node, in counterclockwise order, with at least 64
/// nodes. Returns `None` if a sampled value is zero or nonfinite or if a
/// sampled change of argument exceeds `pi / 2`. Undersampling can still
/// alias a winding whose changes between nodes are not observed.
#[must_use]
pub fn count_zeros_poles<F>(
    f: F,
    center: Complex64,
    radius: f64,
    n: usize,
) -> Option<i64>
where
    F: Fn(Complex64) -> Complex64,
{
    let n = n.max(64);
    winding((0..n).map(|k| f(on_circle(center, radius, k, n).0)), true)
}

/// `f'(z)` by the central difference along the real axis (valid for
/// analytic `f`).
pub fn complex_derivative<F>(
    f: &F,
    z: Complex64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let h = 1e-6 * (1.0 + z.norm());
    (f(z + h) - f(z - h)) / (2.0 * h)
}

/// `f^(n)(z) = n! / (2 pi i) ∮ f(w) / (w - z)^(n+1) dw` over the circle of
/// `radius` around `z` (Cauchy's integral formula); accurate for high
/// orders where finite differences are not.
pub fn cauchy_derivative<F>(
    f: F,
    z: Complex64,
    order: u32,
    radius: f64,
    n: usize,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    if !radius.is_finite() || radius <= 0.0 || n == 0 {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    // Compute n! / radius^n in log-space when order is large to avoid f64 overflow.
    let n_f = n.max(1) as f64;
    let mut sum = Complex64::new(0.0, 0.0);
    let mut log_r_pow = 0.0_f64;
    let mut log_fact = 0.0_f64;
    for k in 1..=order {
        log_fact += f64::from(k).ln();
        log_r_pow += radius.ln();
    }
    let scale = (log_fact - log_r_pow).exp();
    if !scale.is_finite() {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    let dtheta = 2.0 * PI / n_f;
    for j in 0..n {
        let theta = dtheta * (j as f64);
        let (s, c) = theta.sin_cos();
        let offset = Complex64::new(radius * c, radius * s);
        let phase = Complex64::from_polar(1.0, -(f64::from(order) * theta));
        sum += f(z + offset) * phase;
    }
    sum * (scale / n_f)
}

/// Newton's method for an analytic `f` with derivative `df`, from `z0`.
/// Returns the root once a step is below `tolerance`, or `None` after
/// `max_iterations` or at a vanishing derivative.
pub fn newton<F, D>(
    f: F,
    df: D,
    z0: Complex64,
    tolerance: f64,
    max_iterations: usize,
) -> Option<Complex64>
where
    F: Fn(Complex64) -> Complex64,
    D: Fn(Complex64) -> Complex64,
{
    let mut z = z0;
    for _ in 0..max_iterations {
        let slope = df(z);
        if slope.norm() == 0.0 || !slope.is_finite() {
            return None;
        }
        let step = f(z) / slope;
        z -= step;
        if step.norm() <= tolerance * (1.0 + z.norm()) {
            return Some(z);
        }
    }
    None
}

/// A Möbius transformation `z -> (a z + b) / (c z + d)`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Mobius {
    /// Coefficient `a`.
    pub a: Complex64,
    /// Coefficient `b`.
    pub b: Complex64,
    /// Coefficient `c`.
    pub c: Complex64,
    /// Coefficient `d`.
    pub d: Complex64,
}

impl Mobius {
    /// `(a z + b) / (c z + d)`.
    #[must_use]
    pub const fn new(
        a: Complex64,
        b: Complex64,
        c: Complex64,
        d: Complex64,
    ) -> Self {
        Self { a, b, c, d }
    }

    /// The identity map.
    #[must_use]
    pub const fn identity() -> Self {
        let (zero, one) = (Complex64::new(0.0, 0.0), Complex64::new(1.0, 0.0));
        Self {
            a: one,
            b: zero,
            c: zero,
            d: one,
        }
    }

    /// The image of `z` (infinite at the pole `-d / c`).
    #[must_use]
    pub fn apply(
        &self,
        z: Complex64,
    ) -> Complex64 {
        (self.a * z + self.b) / (self.c * z + self.d)
    }

    /// `self ∘ other`.
    #[must_use]
    #[allow(clippy::suspicious_operation_groupings)] // 2x2 matrix product, indices are intentional
    pub fn compose(
        &self,
        other: &Self,
    ) -> Self {
        Self {
            a: self.a * other.a + self.b * other.c,
            b: self.a * other.b + self.b * other.d,
            c: self.c * other.a + self.d * other.c,
            d: self.c * other.b + self.d * other.d,
        }
    }

    /// The inverse map (up to the projective scale).
    #[must_use]
    pub fn inverse(&self) -> Self {
        Self {
            a: self.d,
            b: -self.b,
            c: -self.c,
            d: self.a,
        }
    }
}

/// Numerical Wirtinger derivative `d/dz = 0.5 * (d/dx - i * d/dy)` by central
/// finite differences with step `h`.
pub fn wirtinger_dz_numeric<F>(
    f: F,
    z: Complex64,
    h: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let df_dx = (f(Complex64::new(z.re + h, z.im)) - f(Complex64::new(z.re - h, z.im))) / (2.0 * h);
    let df_dy = (f(Complex64::new(z.re, z.im + h)) - f(Complex64::new(z.re, z.im - h))) / (2.0 * h);
    0.5 * (df_dx - I * df_dy)
}

/// Numerical Wirtinger derivative `d/d(conj(z)) = 0.5 * (d/dx + i * d/dy)` by
/// central finite differences with step `h`.
pub fn wirtinger_dzbar_numeric<F>(
    f: F,
    z: Complex64,
    h: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    let df_dx = (f(Complex64::new(z.re + h, z.im)) - f(Complex64::new(z.re - h, z.im))) / (2.0 * h);
    let df_dy = (f(Complex64::new(z.re, z.im + h)) - f(Complex64::new(z.re, z.im - h))) / (2.0 * h);
    0.5 * (df_dx + I * df_dy)
}

/// Checks the Cauchy-Riemann equations numerically at `z` via
/// `|d/d(conj(z)) f(z)| <= tol * (1 + |d/dz f(z)|)`.
pub fn check_cauchy_riemann_numeric<F>(
    f: F,
    z: Complex64,
    h: f64,
    tol: f64,
) -> bool
where
    F: Fn(Complex64) -> Complex64,
{
    if !h.is_finite() || h <= 0.0 || !tol.is_finite() || tol < 0.0 {
        return false;
    }
    let dzbar = wirtinger_dzbar_numeric(&f, z, h);
    let dz = wirtinger_dz_numeric(&f, z, h);
    dzbar.is_finite() && dz.is_finite() && dzbar.norm() <= tol * (1.0 + dz.norm())
}

/// Numerical complex Laplacian `4 * d^2 f / (dz d(conj(z))) = f_xx + f_yy`
/// using central 2nd-order finite differences with step `h`.
pub fn complex_laplacian_numeric<F>(
    f: F,
    z: Complex64,
    h: f64,
) -> Complex64
where
    F: Fn(Complex64) -> Complex64,
{
    if !h.is_finite() || h <= 0.0 {
        return Complex64::new(f64::NAN, f64::NAN);
    }
    let f0 = f(z);
    let fxx = (f(Complex64::new(z.re + h, z.im)) - 2.0 * f0 + f(Complex64::new(z.re - h, z.im)))
        / (h * h);
    let fyy = (f(Complex64::new(z.re, z.im + h)) - 2.0 * f0 + f(Complex64::new(z.re, z.im - h)))
        / (h * h);
    fxx + fyy
}

/// Numerical polar Cauchy-Riemann residual `(r * u_r - v_theta, u_theta + r * v_r)`
/// at `(r, theta)`.
pub fn cauchy_riemann_polar_numeric<F>(
    f: F,
    r: f64,
    theta: f64,
    h: f64,
) -> (f64, f64)
where
    F: Fn(f64, f64) -> Complex64,
{
    if !h.is_finite() || h <= 0.0 {
        return (f64::NAN, f64::NAN);
    }
    let fr_p = f(r + h, theta);
    let fr_m = f(r - h, theta);
    let ft_p = f(r, theta + h);
    let ft_m = f(r, theta - h);
    let ur = (fr_p.re - fr_m.re) / (2.0 * h);
    let vr = (fr_p.im - fr_m.im) / (2.0 * h);
    let ut = (ft_p.re - ft_m.re) / (2.0 * h);
    let vt = (ft_p.im - ft_m.im) / (2.0 * h);
    (r * ur - vt, ut + r * vr)
}

/// Complex-step differentiation `Im(f(x + i*h)) / h` for real analytic
/// functions, avoiding subtractive cancellation down to machine precision.
pub fn complex_step_derivative<F>(
    f: F,
    x: f64,
    h: f64,
) -> f64
where
    F: Fn(Complex64) -> Complex64,
{
    f(Complex64::new(x, h)).im / h
}

/// Samples a circle counter-clockwise and returns its unnormalised spectrum.
pub(crate) fn circle_spectrum(
    f: impl Fn(Complex64) -> Complex64,
    center: Complex64,
    radius: f64,
    n: usize,
) -> Scratch<'static, Complex64> {
    let mut samples = fft_scratch_for(n);
    // A power-of-two length permits the cache's division to move ahead of the
    // multiplication without rounding the sample angle differently.
    if n.is_power_of_two() && n <= 65_536 {
        let roots = global_pool().unit_roots(n);
        samples.extend(roots.iter().map(|root| f(center + radius * root.conj())));
    } else {
        samples.extend((0..n).map(|k| {
            let angle = TAU * (k as f64) / (n as f64);
            let (sin, cos) = angle.sin_cos();
            f(center + Complex64::new(radius * cos, radius * sin))
        }));
    }
    fft::transform_in_place(&mut samples, false, fft::Norm::Backward);
    samples
}

/// Computes `f^(0)(z0), f^(1)(z0), ..., f^(max_order)(z0)` in one pass using
/// Cauchy's integral formula discretized with an FFT (Lyness-Moler method) on
/// `|z - z0| = radius` with `n_samples` equispaced nodes.
#[must_use]
pub fn cauchy_fft_derivatives<F>(
    f: F,
    z0: Complex64,
    radius: f64,
    max_order: usize,
    n_samples: usize,
) -> Vec<Complex64>
where
    F: Fn(Complex64) -> Complex64,
{
    let n = n_samples.max(max_order + 1).max(4);
    if !radius.is_finite() || radius <= 0.0 {
        return vec![Complex64::new(f64::NAN, f64::NAN); max_order + 1];
    }
    let buffer = circle_spectrum(f, z0, radius, n);

    let mut out = Vec::with_capacity(max_order + 1);
    let nf = n as f64;
    let ln_r = radius.ln();
    let mut ln_fact = 0.0_f64;
    for k in 0..=max_order {
        if k > 0 {
            ln_fact += (k as f64).ln();
        }
        let scale = (ln_fact - (k as f64) * ln_r).exp() / nf;
        let ck = buffer.get(k).copied().unwrap_or_default();
        out.push(ck * scale);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circular_fft_samples_in_order_and_recovers_polynomial_derivatives() {
        use std::cell::RefCell;

        let center = Complex64::new(0.25, -0.125);
        let radius = 0.75;
        for count in [31, 32] {
            let points = RefCell::new(Vec::new());
            let derivatives = cauchy_fft_derivatives(
                |z| {
                    points.borrow_mut().push(z);
                    z * z + 3.0 * z + 5.0
                },
                center,
                radius,
                3,
                count,
            );
            let points = points.into_inner();
            assert_eq!(points.len(), count);
            for (k, actual) in points.into_iter().enumerate() {
                let angle = TAU * k as f64 / count as f64;
                let expected = center + Complex64::from_polar(radius, angle);
                assert!((actual - expected).norm() < 1e-15);
            }
            for (actual, expected) in derivatives.into_iter().zip([
                center * center + 3.0 * center + 5.0,
                2.0 * center + 3.0,
                Complex64::new(2.0, 0.0),
                Complex64::ZERO,
            ]) {
                assert!((actual - expected).norm() < 1e-12);
            }
        }
    }

    fn close(
        a: Complex64,
        b: Complex64,
        tol: f64,
    ) -> bool {
        (a - b).norm() <= tol
    }

    #[test]
    fn residues_and_contours() {
        let one = Complex64::new(1.0, 0.0);
        // 1/z around 0: 2 pi i.
        let value = circle_integral(|z| one / z, Complex64::new(0.0, 0.0), 1.0, 64);
        assert!(close(value, 2.0 * PI * I, 1e-12));
        // exp(z)/z^3 at 0: residue 1/2.
        let r = residue(|z| z.exp() / z.powi(3), Complex64::new(0.0, 0.0), 0.5, 64);
        assert!(close(r, Complex64::new(0.5, 0.0), 1e-12), "{r}");
        // 1/(z^2 + 1) at i: residue -i/2.
        let r = residue(|z| one / (z * z + 1.0), I, 0.5, 128);
        assert!(close(r, -I / 2.0, 1e-12), "{r}");
        // Square path around the origin.
        let square = [
            Complex64::new(1.0, -1.0),
            Complex64::new(1.0, 1.0),
            Complex64::new(-1.0, 1.0),
            Complex64::new(-1.0, -1.0),
            Complex64::new(1.0, -1.0),
        ];
        let value = contour_integral(|z| one / z, &square, 200);
        assert!(close(value, 2.0 * PI * I, 1e-8), "{value}");
        // Analytic integrand: zero.
        let value = contour_integral(|z| z * z, &square, 10);
        assert!(close(value, Complex64::new(0.0, 0.0), 1e-12));
    }

    #[test]
    fn argument_principle() {
        let f = |z: Complex64| (z - 0.5) * (z + 0.25 * I) * (z - 3.0) / (z - 0.1).powi(2);
        assert_eq!(
            count_zeros_poles(f, Complex64::new(0.0, 0.0), 1.0, 512),
            Some(0)
        );
        assert_eq!(
            count_zeros_poles(f, Complex64::new(0.0, 0.0), 4.0, 2048),
            Some(1)
        );
        let g = |z: Complex64| z.powi(5) - 1.0;
        assert_eq!(
            count_zeros_poles(g, Complex64::new(0.0, 0.0), 2.0, 512),
            Some(5)
        );
    }

    #[test]
    fn derivatives() {
        let z = Complex64::new(0.3, 0.7);
        let d = complex_derivative(&|w: Complex64| w.sin(), z);
        assert!(close(d, z.cos(), 1e-8));
        let d3 = cauchy_derivative(num_complex::Complex::exp, z, 3, 1.0, 64);
        assert!(close(d3, z.exp(), 1e-12));
    }

    #[test]
    fn newton_finds_complex_roots() {
        let root = newton(
            |z| z * z + 1.0,
            |z| 2.0 * z,
            Complex64::new(0.3, 0.8),
            1e-14,
            50,
        );
        assert!(root.is_some_and(|r| close(r, I, 1e-12)));
        let cube = newton(
            |z| z.powi(3) - 1.0,
            |z| 3.0 * z * z,
            Complex64::new(-1.0, 1.0),
            1e-14,
            100,
        );
        assert!(cube.is_some_and(|r| close(r.powi(3), Complex64::new(1.0, 0.0), 1e-12)));
        assert!(newton(|z| z * z + 1.0, |_| Complex64::new(0.0, 0.0), I, 1e-14, 5).is_none());
    }

    #[test]
    fn mobius() {
        let m = Mobius::new(
            Complex64::new(1.0, 0.0),
            I,
            Complex64::new(2.0, 0.0),
            Complex64::new(1.0, 1.0),
        );
        let z = Complex64::new(0.4, -0.2);
        assert!(close(m.inverse().apply(m.apply(z)), z, 1e-12));
        assert!(close(m.compose(&m).apply(z), m.apply(m.apply(z)), 1e-12));
        assert!(close(Mobius::identity().apply(z), z, 0.0));
    }
}
