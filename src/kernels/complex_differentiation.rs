//! Shared finite-difference samples for numerical complex analysis.
use num_complex::Complex64;

/// Nonfinite estimate used by the legacy unchecked entrypoints.
pub const NAN: Complex64 = Complex64::new(f64::NAN, f64::NAN);

#[derive(Clone, Copy)]
struct Axis {
    left: f64,
    right: f64,
    left_gap: f64,
    right_gap: f64,
}

impl Axis {
    fn new(
        center: f64,
        step: f64,
    ) -> Option<Self> {
        let left = center - step;
        let right = center + step;
        let left_gap = center - left;
        let right_gap = right - center;
        (left < center && center < right && left_gap.is_finite() && right_gap.is_finite())
            .then_some(Self {
                left,
                right,
                left_gap,
                right_gap,
            })
    }

    fn slope(
        self,
        left: Complex64,
        right: Complex64,
    ) -> Complex64 {
        let span = self.left_gap + self.right_gap;
        if span.is_finite() {
            divided_difference(left, right, span)
        } else {
            let half_span = 0.5_f64.mul_add(self.left_gap, 0.5 * self.right_gap);
            0.5 * divided_difference(left, right, half_span)
        }
    }

    fn curvature(
        self,
        left: Complex64,
        mid: Complex64,
        right: Complex64,
    ) -> Complex64 {
        let left_slope = divided_difference(left, mid, self.left_gap);
        let right_slope = divided_difference(mid, right, self.right_gap);
        let span = self.left_gap + self.right_gap;
        // Divide by half the span except when halving a subnormal step would
        // coarsen it. The latter case uses two times the full-span quotient.
        if span >= f64::MIN_POSITIVE {
            let half_span = 0.5_f64.mul_add(self.left_gap, 0.5 * self.right_gap);
            divided_difference(left_slope, right_slope, half_span)
        } else {
            2.0 * divided_difference(left_slope, right_slope, span)
        }
    }
}

fn divided_difference(
    left: Complex64,
    right: Complex64,
    span: f64,
) -> Complex64 {
    let divide = |left: f64, right: f64| {
        let difference = right - left;
        if difference.is_finite() {
            difference / span
        } else {
            right / span - left / span
        }
    };
    Complex64::new(divide(left.re, right.re), divide(left.im, right.im))
}

fn samples(
    f: impl Fn(Complex64) -> Complex64,
    z: Complex64,
    h: f64,
) -> Option<(Axis, Axis, [Complex64; 4])> {
    let x = Axis::new(z.re, h)?;
    let y = Axis::new(z.im, h)?;
    let points = [
        Complex64::new(x.left, z.im),
        Complex64::new(x.right, z.im),
        Complex64::new(z.re, y.left),
        Complex64::new(z.re, y.right),
    ];
    let mut values = [Complex64::default(); 4];
    for (value, point) in values.iter_mut().zip(points) {
        *value = f(point);
        if !value.is_finite() {
            return None;
        }
    }
    Some((x, y, values))
}

/// Both first Wirtinger derivatives from one Cartesian callback stencil.
#[derive(Clone, Copy)]
pub struct Wirtinger {
    /// Derivative with respect to z.
    pub dz: Complex64,
    /// Derivative with respect to conjugate z.
    pub dzbar: Complex64,
}

impl Wirtinger {
    /// Relative Cauchy–Riemann residual in a frame with a positive coordinate unit.
    pub fn residual(
        self,
        coordinate_scale: f64,
    ) -> f64 {
        let (dz, dzbar, reference) = if coordinate_scale <= 1.0 {
            (
                self.dz * coordinate_scale,
                self.dzbar * coordinate_scale,
                1.0,
            )
        } else {
            (self.dz, self.dzbar, coordinate_scale.recip())
        };
        let scale = reference
            .max(dz.re.abs())
            .max(dz.im.abs())
            .max(dzbar.re.abs())
            .max(dzbar.im.abs());
        (dzbar / scale).norm() / (reference / scale + (dz / scale).norm())
    }
}

/// Calculates both Wirtinger estimates after verifying the rounded sample locations.
pub fn wirtinger(
    f: impl Fn(Complex64) -> Complex64,
    z: Complex64,
    h: f64,
) -> Option<Wirtinger> {
    let (x, y, [xl, xr, yl, yr]) = samples(f, z, h)?;
    let dx = x.slope(xl, xr);
    let dy = y.slope(yl, yr);
    if !dx.is_finite() || !dy.is_finite() {
        return None;
    }
    Some(Wirtinger {
        dz: Complex64::new(dx.re.midpoint(dy.im), dx.im.midpoint(-dy.re)),
        dzbar: Complex64::new(dx.re.midpoint(-dy.im), dx.im.midpoint(dy.re)),
    })
}

/// Calculates a Cartesian Laplacian with independently rounded axis samples.
pub fn laplacian(
    f: impl Fn(Complex64) -> Complex64,
    z: Complex64,
    h: f64,
) -> Option<Complex64> {
    let (x, y, [xl, xr, yl, yr]) = samples(&f, z, h)?;
    let center = f(z);
    let estimate = x.curvature(xl, center, xr) + y.curvature(yl, center, yr);
    estimate.is_finite().then_some(estimate)
}

/// Calculates the polar Cauchy–Riemann residual using radial and angular steps.
pub fn polar(
    f: impl Fn(f64, f64) -> Complex64,
    radius: f64,
    angle: f64,
    radial_step: f64,
    angular_step: f64,
) -> Option<(f64, f64)> {
    let r = Axis::new(radius, radial_step)?;
    let theta = Axis::new(angle, angular_step)?;
    if r.left < 0.0 {
        return None;
    }
    let radial = r.slope(f(r.left, angle), f(r.right, angle));
    let angular = theta.slope(f(radius, theta.left), f(radius, theta.right));
    let first = radius.mul_add(radial.re, -angular.im);
    let second = radius.mul_add(radial.im, angular.re);
    (first.is_finite() && second.is_finite()).then_some((first, second))
}

/// Estimates a real analytic derivative using a single complex-step sample.
pub fn complex_step(
    f: impl Fn(Complex64) -> Complex64,
    x: f64,
    step: f64,
) -> Option<f64> {
    if !x.is_finite() || !step.is_finite() || step == 0.0 {
        return None;
    }
    let sample = f(Complex64::new(x, step));
    let derivative = sample.im / step;
    (sample.is_finite() && derivative.is_finite()).then_some(derivative)
}
