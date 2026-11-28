//! Special-function disks from explicit series and remainders, without midpoint libm error guesses.
#![allow(
    clippy::arithmetic_side_effects,
    reason = "positive bounds are rounded outwards and fixed loop counts fit the integer arithmetic"
)]

use num_complex::Complex64;

use super::ComplexBall as B;
use super::analytic;
use super::up_add;
use super::up_mul;

const TERMS: u32 = 256;

fn rounded(value: f64) -> B {
    B {
        mid: Complex64::new(value, 0.0),
        rad: (value.next_up() - value).max(value - value.next_down()),
    }
}

const fn radius_disk(radius: f64) -> B {
    B {
        mid: Complex64::new(0.0, 0.0),
        rad: radius,
    }
}

fn finite(disk: B) -> Option<B> {
    disk.is_finite().then_some(disk)
}

fn add_tail(
    mut sum: B,
    term: B,
    ratio: f64,
) -> Option<B> {
    if !(0.0..=0.5).contains(&ratio) {
        return None;
    }
    // The next absolute term is at most |term| * ratio, and all later
    // ratios are no greater. 2 covers 1/(1-ratio).
    let tail = up_mul(up_mul(term.modulus_upper(), ratio), 2.0);
    if tail > f64::EPSILON * 0.5 {
        return None;
    }
    sum.rad = up_add(sum.rad, tail);
    finite(sum)
}

#[allow(
    clippy::single_call_fn,
    reason = "keep the enclosure algorithm separate from the public disk method"
)]
pub(super) fn erf(z: B) -> Option<B> {
    if !z.is_finite() {
        return None;
    }
    let constant = rounded(std::f64::consts::FRAC_2_SQRT_PI);
    if z.mid.re == 0.0 && z.mid.im == 0.0 && z.rad == 0.0 {
        return Some(B::zero());
    }
    let z_squared = z.sqr();
    let magnitude = z.modulus_upper();
    if magnitude <= 8.0 {
        // DLMF 7.6.1. Maintain z*(-z²)^n/n! while summing the additional
        // coefficient 1/(2n+1). Dropping that denominator bounds the tail.
        let minus_z_squared = -z_squared;
        let numerator_bound = minus_z_squared.modulus_upper();
        let mut term = z;
        let mut sum = z;
        for n in 1..=TERMS {
            term = term * minus_z_squared * rounded(1.0 / f64::from(n));
            sum = sum + term * rounded(1.0 / f64::from(2 * n + 1));
            let ratio = (numerator_bound / f64::from(n + 1)).next_up();
            if let Some(result) = add_tail(sum, term, ratio) {
                return finite(constant * result);
            }
        }
    }
    // Integrate exp(-u² z²) for u in [0,1]. The real exponent is at
    // most max(0,-Re(z²)) for every point in the disk.
    let exponent = (-z_squared.mid.re + z_squared.rad).next_up().max(0.0);
    let exponential = analytic::exp(B::point(exponent, 0.0)).modulus_upper();
    finite(radius_disk(up_mul(
        up_mul(constant.modulus_upper(), magnitude),
        exponential,
    )))
}

pub(super) fn ln_gamma(z: B) -> Option<B> {
    if !z.is_finite() {
        return None;
    }
    let lower = (z.mid.re - z.rad).next_down();
    let offset = (16.0 - lower).ceil().max(0.0);
    if !offset.is_finite() || offset > 256.0 {
        return None;
    }
    // Aim for Re(current) >= 16, with a verified lower bound of 15. Applying the recurrence
    // horizontally retains the logarithmic branch off the negative real axis.
    let current = z + B::point(offset, 0.0);
    if (current.mid.re - current.rad).next_down() < 15.0 {
        return None;
    }
    let inv = current.recip()?;
    let inv_squared = inv.sqr();
    let correction = [
        (1.0, 12.0),
        (-1.0, 360.0),
        (1.0, 1260.0),
        (-1.0, 1680.0),
        (1.0, 1188.0),
        (-691.0, 360360.0),
        (1.0, 156.0),
    ]
    .into_iter()
    .rev()
    .fold(B::zero(), |value, (numer, denom)| {
        rounded(numer / denom) + inv_squared * value
    }) * inv;
    // DLMF 5.11.1, 5.11(ii), n=8. The first omitted coefficient
    // is |B16|/(16*15). Re z >= 15 limits the remaining factor
    // by 2^8 / 15^15 < 6e-16, giving total remainder < 2e-17.
    let mut result = (current - B::point(0.5, 0.0)) * analytic::ln(current)? - current
        + rounded(0.918_938_533_204_672_7)
        + correction;
    result.rad = up_add(result.rad, 2e-17);
    for n in 0..(offset as u32) {
        result = result - analytic::ln(z + B::point(f64::from(n), 0.0))?;
    }
    finite(result)
}

#[allow(
    clippy::single_call_fn,
    reason = "keep the enclosure algorithm separate from the public disk method"
)]
pub(super) fn bessel_j(
    nu: f64,
    z: B,
) -> Option<B> {
    if !z.is_finite() || !nu.is_finite() || nu < 0.0 && nu.fract() != 0.0 {
        return None;
    }
    let order = nu.abs();
    let integral_order = order.fract() == 0.0;
    if integral_order && order > f64::from(TERMS) {
        // DLMF 10.9.2 bounds all integer orders by exp(|Im(z)|).
        let bound = up_add(z.mid.im.abs(), z.rad);
        return finite(radius_disk(
            analytic::exp(B::point(bound, 0.0)).modulus_upper(),
        ));
    }
    // At z=0 the principal function has the continuous value 0 for nu>0.
    if z.mid.re == 0.0 && z.mid.im == 0.0 && z.rad == 0.0 {
        return Some(if order == 0.0 {
            B::one()
        } else {
            B::zero()
        });
    }
    let half = z.scale_real(0.5);
    let mut prefactor = if integral_order {
        (1..=order as u32).fold(B::one(), |term, k| {
            term * half * rounded(1.0 / f64::from(k))
        })
    } else {
        let gamma_argument = B::point(order, 0.0) + B::one();
        let log = analytic::ln(half)?.scale_real(order) - ln_gamma(gamma_argument)?;
        analytic::exp(log)
    };
    if nu < 0.0 && order % 2.0 > 0.0 {
        prefactor = -prefactor;
    }
    let square = -half.sqr();
    let bound = square.modulus_upper();
    if bound <= 64.0 {
        // DLMF 10.2.2. At each n the absolute ratios from the last
        // included term onwards are bounded by |z/2|²/(n+1)².
        let mut term = B::one();
        let mut sum = term;
        for n in 1..=TERMS {
            let divisor = B::point(order, 0.0) + B::point(f64::from(n), 0.0);
            term = (term * square * rounded(1.0 / f64::from(n))).div_ball(&divisor)?;
            sum = sum + term;
            let ratio = (bound / f64::from(n + 1).powi(2)).next_up();
            if let Some(result) = add_tail(sum, term, ratio) {
                return finite(prefactor * result);
            }
        }
    }
    // DLMF 10.9.4. After dividing out the prefactor the integrand
    // exp(i z t) is averaged on [-1,1] with a positive normalized weight.
    let bound = up_add(z.mid.im.abs(), z.rad);
    let exponential = analytic::exp(B::point(bound, 0.0)).modulus_upper();
    finite(radius_disk(up_mul(prefactor.modulus_upper(), exponential)))
}
