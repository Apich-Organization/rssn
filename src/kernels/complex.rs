//! Numerical complex analysis: contour integrals, residues, the argument
//! principle, complex derivatives and Möbius transformations.
//!
//! Closed circular contours use the trapezoidal rule, which converges
//! geometrically for functions analytic in an annulus around the circle;
//! open polygonal paths use composite Simpson on each segment.

use std::f64::consts::PI;

use num_complex::Complex64;
use serde::Deserialize;
use serde::Serialize;

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

/// The number of zeros minus the number of poles of `f` inside a circle.
///
/// The circle is `|z - center| = radius` (argument principle), and the count
/// is the winding number of `f` around the origin. Robust: it follows the continuous argument of
/// `f` instead of integrating `f'/f`. `None` if `f` vanishes or is not
/// finite on the contour, or the sampling is too coarse to follow it.
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
    let mut previous = f(on_circle(center, radius, 0, n).0);
    let mut winding = 0.0;
    for k in 1..=n {
        let value = f(on_circle(center, radius, k % n, n).0);
        if !value.is_finite() || value.norm() == 0.0 {
            return None;
        }
        let step = (value / previous).arg();
        if step.abs() > PI / 2.0 {
            return None;
        }
        winding += step;
        previous = value;
    }
    #[allow(clippy::cast_possible_truncation)]
    Some((winding / (2.0 * PI)).round() as i64)
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
    let factorial: f64 = (1..=order).map(f64::from).product();
    let power = i32::try_from(order).unwrap_or(i32::MAX).saturating_add(1);
    circle_integral(|w| f(w) / (w - z).powi(power), z, radius, n) * factorial / (2.0 * PI * I)
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
        Self { a: one, b: zero, c: zero, d: one }
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
        Self { a: self.d, b: -self.b, c: -self.c, d: self.a }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        assert_eq!(count_zeros_poles(f, Complex64::new(0.0, 0.0), 1.0, 512), Some(0));
        assert_eq!(count_zeros_poles(f, Complex64::new(0.0, 0.0), 4.0, 2048), Some(1));
        let g = |z: Complex64| z.powi(5) - 1.0;
        assert_eq!(count_zeros_poles(g, Complex64::new(0.0, 0.0), 2.0, 512), Some(5));
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
        let root = newton(|z| z * z + 1.0, |z| 2.0 * z, Complex64::new(0.3, 0.8), 1e-14, 50);
        assert!(root.is_some_and(|r| close(r, I, 1e-12)));
        let cube = newton(|z| z.powi(3) - 1.0, |z| 3.0 * z * z, Complex64::new(-1.0, 1.0), 1e-14, 100);
        assert!(cube.is_some_and(|r| close(r.powi(3), Complex64::new(1.0, 0.0), 1e-12)));
        assert!(newton(|z| z * z + 1.0, |_| Complex64::new(0.0, 0.0), I, 1e-14, 5).is_none());
    }

    #[test]
    fn mobius() {
        let m = Mobius::new(Complex64::new(1.0, 0.0), I, Complex64::new(2.0, 0.0), Complex64::new(1.0, 1.0));
        let z = Complex64::new(0.4, -0.2);
        assert!(close(m.inverse().apply(m.apply(z)), z, 1e-12));
        assert!(close(m.compose(&m).apply(z), m.apply(m.apply(z)), 1e-12));
        assert!(close(Mobius::identity().apply(z), z, 0.0));
    }
}
