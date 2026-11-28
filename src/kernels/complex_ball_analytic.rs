//! Outward enclosures for the elementary functions used by complex disk jets.

use num_complex::Complex64;

use super::ComplexBall;
use super::distance_lower;
use super::down_mul;
use super::norm_lower;
use super::up_add;

const WHOLE: ComplexBall = ComplexBall {
    mid: Complex64::new(0.0, 0.0),
    rad: f64::INFINITY,
};

fn reciprocal_integer(n: u32) -> ComplexBall {
    let rounded = 1.0 / f64::from(n);
    // One upper ulp bounds the rounding in 1/n in either direction.
    ComplexBall {
        mid: Complex64::new(rounded, 0.0),
        rad: rounded.next_up() - rounded,
    }
}

pub(super) fn exp(z: ComplexBall) -> ComplexBall {
    let mut bound = z.modulus_upper();
    // The fixed work budget covers at most 22 squarings. Large or
    // unrepresentable enclosures fall back to the whole plane.
    if !bound.is_finite() || bound > f64::from(1 << 20) {
        return WHOLE;
    }
    let mut scale = 1.0;
    let mut squarings = 0;
    while bound > 0.25 {
        bound *= 0.5;
        scale *= 0.5;
        squarings += 1;
    }
    let scaled = z.scale_real(scale);
    if scaled.modulus_upper() > 0.5 {
        return WHOLE;
    }
    let mut term = ComplexBall::one();
    let mut sum = term;
    for n in 1..=18 {
        term = term * scaled * reciprocal_integer(n);
        sum = sum + term;
    }
    // For |z| <= 1/2 the omitted tail is at most 2*(1/2)^19 / 19!.
    sum.rad = up_add(sum.rad, 1e-22);
    for _ in 0..squarings {
        sum = sum.sqr();
    }
    sum
}

const fn rotate_i(z: ComplexBall) -> ComplexBall {
    ComplexBall {
        mid: Complex64::new(-z.mid.im, z.mid.re),
        rad: z.rad,
    }
}

pub(super) fn sin_cos(z: ComplexBall) -> (ComplexBall, ComplexBall) {
    let iz = rotate_i(z);
    let (difference, sum) = exp_pair(iz);
    (-rotate_i(difference), sum)
}

pub(super) fn sinh_cosh(z: ComplexBall) -> (ComplexBall, ComplexBall) {
    exp_pair(z)
}

fn exp_pair(z: ComplexBall) -> (ComplexBall, ComplexBall) {
    let pos = exp(z);
    let neg = if z.mid.re == 0.0 {
        // An imaginary-axis disk satisfies -D = conjugate(D).
        pos.conj()
    } else {
        // A verified reciprocal reuses the exponential when its enclosure
        // excludes zero. Wide enclosures use a second exponential instead.
        pos.recip()
            .filter(ComplexBall::is_finite)
            .unwrap_or_else(|| exp(-z))
    };
    ((pos - neg).scale_real(0.5), (pos + neg).scale_real(0.5))
}

fn misses_cut(
    z: ComplexBall,
    endpoint: f64,
) -> bool {
    let distance = if z.mid.re <= endpoint {
        z.mid.im.abs()
    } else {
        distance_lower(z.mid, Complex64::new(endpoint, 0.0))
    };
    z.is_finite() && distance > z.rad
}

pub(super) fn sqrt(z: ComplexBall) -> Option<ComplexBall> {
    if !z.is_finite() {
        return None;
    }
    if z.contains_zero() {
        return ComplexBall::new(Complex64::new(0.0, 0.0), z.modulus_upper().sqrt().next_up());
    }
    if !misses_cut(z, 0.0) {
        return None;
    }
    let guess = z.mid.sqrt();
    // The exact principal root and this guess lie in the same quadrant.
    // Hence |guess + root| >= |guess|, and the polynomial residual certifies
    // the center independently of the rounding of the initial complex sqrt.
    let mid = Complex64::new(guess.re.abs(), guess.im.abs().copysign(z.mid.im));
    let center = ComplexBall::of(mid);
    let mid_lower = norm_lower(mid);
    let denom = down_mul(2.0, z.modulus_lower().sqrt().next_down());
    if mid_lower <= 0.0 || denom <= 0.0 {
        return None;
    }
    let residue = (center.sqr() - ComplexBall::of(z.mid)).modulus_upper();
    let center_error = (residue / mid_lower).next_up();
    // Along each segment in the disk, |sqrt'(w)| <= 1/(2*sqrt(min |w|)).
    let radius = up_add(center_error, (z.rad / denom).next_up());
    radius
        .is_finite()
        .then_some(ComplexBall { mid, rad: radius })
}

pub(super) fn ln(z: ComplexBall) -> Option<ComplexBall> {
    if !misses_cut(z, 0.0) {
        return None;
    }
    let one = ComplexBall::one();
    let mut root = z;
    let mut scale = 1.0;
    for _ in 0..=16 {
        if let Some(ratio) = (root - one).div_ball(&(root + one))
            && ratio.modulus_upper() <= 0.5
        {
            return Some(log_series(ratio).scale_real(scale));
        }
        root = sqrt(root)?;
        scale *= 2.0;
    }
    None
}

fn log_series(ratio: ComplexBall) -> ComplexBall {
    // ln(w) = 2*atanh((w-1)/(w+1)). The ratio disk lies within |ratio|<=1/2.
    let squared = ratio.sqr();
    let mut power = ratio;
    let mut sum = power;
    for n in 1..=23 {
        power = power * squared;
        sum = sum + power * reciprocal_integer(2 * n + 1);
    }
    let mut enclosure = sum.scale_real(2.0);
    // The tail is at most 2*(1/2)^49 / (49*(1-1/4)), below 2^-53.
    enclosure.rad = up_add(enclosure.rad, f64::EPSILON * 0.5);
    enclosure
}

pub(super) fn acosh(z: ComplexBall) -> Option<(ComplexBall, ComplexBall)> {
    if !misses_cut(z, 1.0) {
        return None;
    }
    // Both factors use the principal square root separately. Replacing
    // their product with sqrt(z*z-1) can switch its sign.
    let one = ComplexBall::one();
    let denom = sqrt(z - one)? * sqrt(z + one)?;
    let value = ln(z + denom)?;
    Some((value, denom))
}

fn misses_real_outer_cuts(z: ComplexBall) -> bool {
    let reflected = ComplexBall {
        mid: Complex64::new(-z.mid.re.abs(), z.mid.im),
        rad: z.rad,
    };
    misses_cut(reflected, -1.0)
}

pub(super) fn asin(z: ComplexBall) -> Option<(ComplexBall, ComplexBall)> {
    if !misses_real_outer_cuts(z) {
        return None;
    }
    let denom = sqrt(ComplexBall::one() - z.sqr())?;
    let value = -rotate_i(ln(rotate_i(z) + denom)?);
    Some((value, denom))
}

pub(super) fn acos(z: ComplexBall) -> Option<(ComplexBall, ComplexBall)> {
    use std::f64::consts::FRAC_PI_2;
    let (asin, denom) = asin(z)?;
    let half_pi = ComplexBall {
        mid: Complex64::new(FRAC_PI_2, 0.0),
        rad: FRAC_PI_2.next_up() - FRAC_PI_2,
    };
    Some((half_pi - asin, denom))
}

pub(super) fn atan(z: ComplexBall) -> Option<(ComplexBall, ComplexBall)> {
    let iz = rotate_i(z);
    if !misses_real_outer_cuts(iz) {
        return None;
    }
    let one = ComplexBall::one();
    let ratio = (one - iz).div_ball(&(one + iz))?;
    let value = rotate_i(ln(ratio)?).scale_real(0.5);
    Some((value, one + z.sqr()))
}

pub(super) fn asinh(z: ComplexBall) -> Option<(ComplexBall, ComplexBall)> {
    if !misses_real_outer_cuts(rotate_i(z)) {
        return None;
    }
    let denom = sqrt(z.sqr() + ComplexBall::one())?;
    let value = ln(z + denom)?;
    Some((value, denom))
}

pub(super) fn atanh(z: ComplexBall) -> Option<(ComplexBall, ComplexBall)> {
    if !misses_real_outer_cuts(z) {
        return None;
    }
    let one = ComplexBall::one();
    let ratio = (one + z).div_ball(&(one - z))?;
    let value = ln(ratio)?.scale_real(0.5);
    Some((value, one - z.sqr()))
}
