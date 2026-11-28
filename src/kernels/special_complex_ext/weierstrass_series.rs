//! Scaled theta series for the Weierstrass sigma and zeta functions.
#![allow(
    clippy::arithmetic_side_effects,
    reason = "these complex recurrences use binary64"
)]

use std::f64::consts::FRAC_PI_2;
use std::f64::consts::PI;
use std::sync::Arc;

use num_complex::Complex64;

use crate::kernels::reuse::global_pool;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);
const NAN: Complex64 = Complex64::new(f64::NAN, f64::NAN);
const MAX_EXACT_SHIFT: f64 = 4_503_599_627_370_496.0;
const TERMS: usize = 8;

fn divide(
    a: Complex64,
    b: Complex64,
) -> Complex64 {
    let scale = b.re.abs().max(b.im.abs());
    (a / scale) / (b / scale)
}

// Reduce the basis by integral operations until |q| <= exp(-pi sqrt(3)/2).
#[allow(clippy::single_call_fn)]
fn reduce_periods(
    mut first: Complex64,
    mut second: Complex64,
) -> Option<(Complex64, Complex64)> {
    for _ in 0..64 {
        let mut ratio = divide(second, first);
        if !ratio.is_finite() || ratio.im == 0.0 {
            return None;
        }
        if ratio.im < 0.0 {
            second = -second;
            ratio = -ratio;
        }
        let shift = ratio.re.round();
        second = Complex64::new(
            (-shift).mul_add(first.re, second.re),
            (-shift).mul_add(first.im, second.im),
        );
        let ratio = divide(second, first);
        if !ratio.is_finite() || ratio.im <= 0.0 {
            return None;
        }
        if ratio.norm() >= 1.0 - 4.0 * f64::EPSILON {
            return Some((first, ratio));
        }
        (first, second) = (second, -first);
    }
    None
}

struct ThetaPlan {
    first: Complex64,
    tau: Complex64,
    scale: Complex64,
    phases: [Complex64; TERMS],
    derivative_zero: Complex64,
    quadratic: Complex64,
}

struct Centered {
    argument: Complex64,
    shift: Complex64,
    vertical: f64,
    odd: bool,
}

impl ThetaPlan {
    #[allow(clippy::single_call_fn, reason = "constructs the immutable cache plan")]
    fn new(
        first: Complex64,
        second: Complex64,
    ) -> Option<Self> {
        let (first, tau) = reduce_periods(first, second)?;
        let scale = divide(Complex64::new(FRAC_PI_2, 0.0), first);
        if !scale.is_finite() || scale == ZERO {
            return None;
        }
        let mut derivative_zero = ONE;
        let mut triple_derivative_zero = -ONE;
        let phases = std::array::from_fn(|index| {
            let n = index as f64 + 1.0;
            let sign = if index.is_multiple_of(2) {
                -1.0
            } else {
                1.0
            };
            let phase = Complex64::from_polar(sign, PI * tau.re * n * (n + 1.0));
            let weight = phase * (-PI * tau.im * n * (n + 1.0)).exp();
            let frequency = 2.0 * n + 1.0;
            derivative_zero += frequency * weight;
            triple_derivative_zero -= frequency.powi(3) * weight;
            phase
        });

        // Cancel the cubic theta term, as sigma(z)/z has no z^2 term.
        // These coefficients use the theta argument pi z / (2 omega1).
        let quadratic = -triple_derivative_zero / (6.0 * derivative_zero);
        Some(Self {
            first,
            tau,
            scale,
            phases,
            derivative_zero,
            quadratic,
        })
    }

    fn center(
        &self,
        z: Complex64,
    ) -> Option<Centered> {
        let relative = 0.5 * divide(z, self.first);
        let vertical = (relative.im / self.tau.im).round();
        let horizontal = (-vertical).mul_add(self.tau.re, relative.re).round();
        if !horizontal.is_finite()
            || !vertical.is_finite()
            || horizontal.abs().max(vertical.abs()) >= MAX_EXACT_SHIFT
        {
            return None;
        }
        let argument = PI
            * Complex64::new(
                (-vertical).mul_add(self.tau.re, relative.re - horizontal),
                (-vertical).mul_add(self.tau.im, relative.im),
            );
        let shift = PI * (horizontal + vertical * self.tau);
        let odd = horizontal.rem_euclid(2.0) != 0.0 || vertical.rem_euclid(2.0) != 0.0;
        (argument.is_finite() && shift.is_finite()).then_some(Centered {
            argument,
            shift,
            vertical,
            odd,
        })
    }

    fn theta(
        &self,
        position: Complex64,
    ) -> (Complex64, Complex64) {
        let (mut theta, mut derivative) = scaled_trigonometry(position, 1.0);
        let relative_height = (position.im.abs() / self.tau.im) / PI;
        for (index, phase) in self.phases.iter().enumerate() {
            let n = index as f64 + 1.0;
            // Combine the two exponentials before either overflows.
            let exponent = -PI * n * (n + 1.0 - 2.0 * relative_height) * self.tau.im;
            let amplitude = exponent.exp();
            if amplitude < 1e-24 {
                break;
            }
            let frequency = 2.0 * n + 1.0;
            let factor = amplitude * phase;
            let (sine, cosine) = scaled_trigonometry(position, frequency);
            theta += factor * sine;
            derivative += factor * frequency * cosine;
        }
        (theta, derivative)
    }

    fn log_shift(
        &self,
        point: &Centered,
    ) -> Complex64 {
        2.0 * (self.quadratic * point.shift - I * point.vertical)
    }
}

fn scaled_trigonometry(
    w: Complex64,
    frequency: f64,
) -> (Complex64, Complex64) {
    let (sin, cos) = (frequency * w.re).sin_cos();
    let tail = (-2.0 * frequency * w.im.abs()).exp_m1();
    let cosh = 1.0 + 0.5 * tail;
    let sinh = -0.5 * tail * w.im.signum();
    (
        Complex64::new(sin * cosh, cos * sinh),
        Complex64::new(cos * cosh, -sin * sinh),
    )
}

fn plan(
    first: Complex64,
    second: Complex64,
) -> Arc<Option<ThetaPlan>> {
    let key = [
        first.re.to_bits(),
        first.im.to_bits(),
        second.re.to_bits(),
        second.im.to_bits(),
    ];
    global_pool()
        .named_plan_cache::<[u64; 4], Option<ThetaPlan>>("special.weierstrass_theta", 16)
        .get_or_insert_with(key, || Arc::new(ThetaPlan::new(first, second)))
}

#[allow(
    clippy::single_call_fn,
    reason = "keeps quasiperiodic series internals private"
)]
pub(super) fn zeta(
    z: Complex64,
    first: Complex64,
    second: Complex64,
) -> Complex64 {
    let retained = plan(first, second);
    let Some(plan) = retained.as_ref() else {
        return NAN;
    };
    let Some(point) = plan.center(z) else {
        return NAN;
    };
    if point.argument == ZERO {
        return NAN;
    }
    let (theta, derivative) = plan.theta(point.argument);
    plan.scale
        * (2.0 * plan.quadratic * point.argument
            + divide(derivative, theta)
            + plan.log_shift(&point))
}

#[allow(
    clippy::single_call_fn,
    reason = "keeps quasiperiodic series internals private"
)]
pub(super) fn sigma(
    z: Complex64,
    first: Complex64,
    second: Complex64,
) -> Complex64 {
    let retained = plan(first, second);
    let Some(plan) = retained.as_ref() else {
        return NAN;
    };
    let Some(point) = plan.center(z) else {
        return NAN;
    };
    let w = point.argument;
    if w == ZERO {
        return ZERO;
    }
    let (theta, _) = plan.theta(w);
    let exponent =
        plan.quadratic * w * w + plan.log_shift(&point) * (w + 0.5 * point.shift) + w.im.abs();
    let value = exponent.exp() * divide(divide(theta, plan.derivative_zero), plan.scale);
    if point.odd { -value } else { value }
}
