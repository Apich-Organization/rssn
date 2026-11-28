//! Certified binary approximations to mathematical constants.
//!
//! The shared cache stores outward fixed-point endpoints. A value is rounded
//! only when both endpoints round to the same dyadic. The working precision
//! can exceed the precision limit of a returned dyadic.

use std::sync::Arc;
use std::sync::Mutex;
use std::sync::OnceLock;

use num_bigint::BigInt;

use super::ball::Ball;
use super::fixed::Fx;
use super::real::MAX_BITS;
use super::real::Real;
use super::real::Round;

const MAX_WORKING_BITS: u32 = MAX_BITS + 512;
const GUARD_BITS: u32 = 32;

struct Bounds {
    low: BigInt,
    high: BigInt,
    scale: u32,
}

impl Bounds {
    fn rounded(
        &self,
        prec: u32,
    ) -> Result<Option<Real>, String> {
        let exponent = -i64::from(self.scale);
        let low = Real::rounded(self.low.clone(), exponent, prec, Round::NearestEven)?;
        let high = Real::rounded(self.high.clone(), exponent, prec, Round::NearestEven)?;
        Ok((low.cmp_value(&high) == 0).then_some(low))
    }

    fn enclose(
        &self,
        prec: u32,
    ) -> Result<Option<Ball>, String> {
        let Some(mid) = self.rounded(prec)? else {
            return Ok(None);
        };
        // The working scale has more fractional bits than the midpoint.
        let exponent = -i64::from(self.scale);
        let shift = (mid.exponent() - exponent).unsigned_abs();
        let centre = mid.mantissa() << shift;
        let radius = (&centre - &self.low).max(&self.high - &centre);
        let radius = Real::rounded(radius, exponent, prec, Round::TowardPositive)?;
        Ball::new(mid, radius).map(Some)
    }

    fn scaled_floor(
        &self,
        bits: u32,
    ) -> Option<BigInt> {
        let shift = self.scale - bits;
        let low = &self.low >> shift;
        let high = &self.high >> shift;
        (low == high).then_some(low)
    }
}

type Compute = fn(u32) -> Result<Bounds, String>;

/// One lazily extended constant.
struct Cache {
    slot: OnceLock<Mutex<Option<Arc<Bounds>>>>,
}

impl Cache {
    const fn new() -> Self {
        Self {
            slot: OnceLock::new(),
        }
    }

    fn bounds(
        &self,
        scale: u32,
        compute: &dyn Fn(u32) -> Result<Bounds, String>,
    ) -> Result<Arc<Bounds>, String> {
        let mutex = self.slot.get_or_init(|| Mutex::new(None));
        let mut guard = mutex
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let scale = match guard.as_ref() {
            | Some(cached) if cached.scale >= scale => return Ok(Arc::clone(cached)),
            | Some(cached) => scale.max(cached.scale.saturating_mul(2)),
            | None => scale,
        }
        .min(MAX_WORKING_BITS);
        let bounds = Arc::new(compute(scale)?);
        *guard = Some(Arc::clone(&bounds));
        drop(guard);
        Ok(bounds)
    }

    fn choose<T>(
        &self,
        prec: u32,
        compute: &dyn Fn(u32) -> Result<Bounds, String>,
        choose: impl Fn(&Bounds) -> Result<Option<T>, String>,
    ) -> Result<T, String> {
        let mut scale = prec + GUARD_BITS;
        loop {
            let bounds = self.bounds(scale, compute)?;
            if let Some(value) = choose(&bounds)? {
                return Ok(value);
            }
            if bounds.scale == MAX_WORKING_BITS {
                return Err("constant bounds do not determine the requested rounding".into());
            }
            scale = (bounds.scale + GUARD_BITS).min(MAX_WORKING_BITS);
        }
    }

    fn get(
        &self,
        prec: u32,
        compute: &dyn Fn(u32) -> Result<Bounds, String>,
    ) -> Result<Real, String> {
        let prec = prec.clamp(1, MAX_BITS);
        self.choose(prec, compute, |bounds| bounds.rounded(prec))
    }
}

static PI: Cache = Cache::new();
static LN2: Cache = Cache::new();
static E: Cache = Cache::new();
static SQRT2: Cache = Cache::new();
static PHI: Cache = Cache::new();

/// Correctly rounded `pi` at `prec` bits.
///
/// # Errors
/// Fails if the working precision limit cannot certify rounding.
pub fn pi(prec: u32) -> Result<Real, String> {
    PI.get(prec, &machin)
}

/// Correctly rounded `ln 2` at `prec` bits.
///
/// # Errors
/// Fails if the working precision limit cannot certify rounding.
pub fn ln2(prec: u32) -> Result<Real, String> {
    LN2.get(prec, &ln2_at)
}

/// Outward cached bounds for ln(2) at the requested fixed-point scale.
pub(super) fn ln2_at_scale(scale: u32) -> Result<(BigInt, BigInt), String> {
    at_scale(&LN2, &ln2_at, scale)
}

/// Outward cached bounds for pi at the requested fixed-point scale.
pub(super) fn pi_at_scale(scale: u32) -> Result<(BigInt, BigInt), String> {
    at_scale(&PI, &machin, scale)
}

fn at_scale(
    cache: &Cache,
    compute: &dyn Fn(u32) -> Result<Bounds, String>,
    scale: u32,
) -> Result<(BigInt, BigInt), String> {
    if scale > MAX_WORKING_BITS {
        return Err("constant scale out of range".into());
    }
    let bounds = cache.bounds(scale, compute)?;
    let shift = bounds.scale - scale;
    Ok((&bounds.low >> shift, ((&bounds.high - 1_i32) >> shift) + 1))
}

/// Correctly rounded `e` at `prec` bits.
///
/// # Errors
/// Fails if the working precision limit cannot certify rounding.
pub fn e(prec: u32) -> Result<Real, String> {
    E.get(prec, &e_at)
}

/// Correctly rounded `sqrt 2` at `prec` bits.
///
/// # Errors
/// Fails if the working precision limit cannot certify rounding.
pub fn sqrt2(prec: u32) -> Result<Real, String> {
    SQRT2.get(prec, &sqrt2_at)
}

/// Correctly rounded golden ratio `(1 + sqrt 5) / 2` at `prec` bits.
///
/// # Errors
/// Fails if the working precision limit cannot certify rounding.
pub fn phi(prec: u32) -> Result<Real, String> {
    PHI.get(prec, &phi_at)
}

/// Bounds the odd-power series for atan(1/n) or atanh(1/n).
fn odd_series(
    n: i64,
    alternating: bool,
    scale: u32,
) -> Result<Bounds, String> {
    let mut power = Fx::one(scale).div_i64(n);
    let mut sum = Fx::zero(scale);
    let step = n * n;
    for k in 0..=i64::from(scale) {
        // For positive integers, successive floors in these divisions are
        // the same as one floor of the exact ratio for this term.
        let term = power.div_i64(2 * k + 1);
        if term.is_zero() {
            // Each of the k computed terms has error below one unit. The
            // alternating tail is below one unit. For the positive series,
            // its ratio is at most 1/9, so the tail is below two units.
            return if alternating {
                Ok(Bounds {
                    low: &sum.v - (k + 1),
                    high: sum.v + (k + 1),
                    scale,
                })
            } else {
                Ok(Bounds {
                    low: sum.v.clone(),
                    high: sum.v + (k + 2),
                    scale,
                })
            };
        }
        sum = if alternating && k % 2 == 1 {
            sum.sub(&term)
        } else {
            sum.add(&term)
        };
        power = power.div_i64(step);
    }
    Err("constant series did not converge".into())
}

/// Uses Machin's identity with directed combination of both series.
fn machin(scale: u32) -> Result<Bounds, String> {
    let a = odd_series(5, true, scale)?;
    let b = odd_series(239, true, scale)?;
    Ok(Bounds {
        low: (a.low << 4_u32) - (b.high << 2_u32),
        high: (a.high << 4_u32) - (b.low << 2_u32),
        scale,
    })
}

/// Uses `ln(2) = 2 atanh(1/3)`.
fn ln2_at(scale: u32) -> Result<Bounds, String> {
    let Bounds { low, high, .. } = odd_series(3, false, scale)?;
    Ok(Bounds {
        low: low << 1_u32,
        high: high << 1_u32,
        scale,
    })
}

/// Bounds e via the positive factorial series.
fn e_at(scale: u32) -> Result<Bounds, String> {
    let mut term = Fx::one(scale);
    let mut sum = term.clone();
    for k in 1..=i64::from(scale) + 2 {
        term = term.div_i64(k);
        if term.is_zero() {
            // Positive successive floors compute floor(2^scale / k!).
            // After the first omitted term the ratio is at most 1/2.
            return Ok(Bounds {
                low: sum.v.clone(),
                high: sum.v + (k + 2),
                scale,
            });
        }
        sum = sum.add(&term);
    }
    Err("constant factorial series did not converge".into())
}

fn sqrt2_at(scale: u32) -> Result<Bounds, String> {
    let shift = scale.checked_mul(2).ok_or("constant scale out of range")?;
    let low = (BigInt::from(2) << shift).sqrt();
    Ok(Bounds {
        low: low.clone(),
        high: low + 1,
        scale,
    })
}

fn phi_at(scale: u32) -> Result<Bounds, String> {
    let half_scale = scale.checked_sub(1).ok_or("constant scale out of range")?;
    let low = (BigInt::from(1) << half_scale) + (BigInt::from(5) << (2 * half_scale)).sqrt();
    Ok(Bounds {
        low: low.clone(),
        high: low + 1,
        scale,
    })
}

fn lookup(name: &str) -> Result<(&'static Cache, Compute, i64), String> {
    match name {
        | "pi" => Ok((&PI, machin, 0)),
        | "e" => Ok((&E, e_at, 0)),
        | "tau" => Ok((&PI, machin, 1)),
        | "ln2" => Ok((&LN2, ln2_at, 0)),
        | "sqrt2" => Ok((&SQRT2, sqrt2_at, 0)),
        | "phi" => Ok((&PHI, phi_at, 0)),
        | _ => Err(format!("unknown constant `{name}`")),
    }
}

/// The correctly rounded constant `name` at `prec` bits.
///
/// # Errors
/// Fails for an unknown name or when the working precision limit cannot
/// certify rounding.
pub fn named(
    name: &str,
    prec: u32,
) -> Result<Real, String> {
    let (cache, compute, shift) = lookup(name)?;
    cache.get(prec, &compute)?.scale_pow2(shift)
}

/// A certified dyadic enclosure of the named constant.
///
/// # Errors
/// Fails for an unknown name or when the working precision limit cannot
/// certify the midpoint.
pub fn named_ball(
    name: &str,
    prec: u32,
) -> Result<Ball, String> {
    let prec = prec.clamp(1, MAX_BITS);
    let (cache, compute, shift) = lookup(name)?;
    let enclosure = cache.choose(prec, &compute, |bounds| bounds.enclose(prec))?;
    Ball::new(
        enclosure.mid.scale_pow2(shift)?,
        enclosure.rad.scale_pow2(shift)?,
    )
}

/// The names accepted by [`named`] and [`named_ball`].
#[must_use]
pub const fn names() -> &'static [&'static str] {
    &["pi", "e", "tau", "ln2", "sqrt2", "phi"]
}

/// The exact integer floor of `pi * 2^bits`.
///
/// # Errors
/// Fails when `bits` exceeds [`MAX_BITS`] or the working precision limit
/// cannot certify the result.
pub fn pi_scaled(bits: u32) -> Result<BigInt, String> {
    if bits > MAX_BITS {
        return Err(format!("pi fractional bits must be within 0..={MAX_BITS}"));
    }
    PI.choose(bits, &machin, |bounds| Ok(bounds.scaled_floor(bits)))
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::Cell;

    #[test]
    fn warm_cache_returns_the_requested_precision_without_recomputing() -> Result<(), String> {
        let cache = Cache {
            slot: OnceLock::new(),
        };
        let computations = Cell::new(0_u32);
        let compute = |precision| {
            computations.set(computations.get().saturating_add(1));
            machin(precision)
        };

        let cold = cache.get(8, &compute)?;
        assert_eq!(cold.precision(), 8);
        assert_eq!(computations.get(), 1);

        let wide = cache.get(256, &compute)?;
        assert_eq!(wide.precision(), 256);
        assert_eq!(computations.get(), 2);

        let warm = cache.get(8, &compute)?;
        assert_eq!(warm.precision(), cold.precision());
        assert_eq!(warm.mantissa(), cold.mantissa());
        assert_eq!(warm.exponent(), cold.exponent());
        assert_eq!(cache.get(256, &compute)?.mantissa(), wide.mantissa());
        assert_eq!(computations.get(), 2);
        Ok(())
    }

    #[test]
    fn pi_scaled_uses_the_requested_fractional_bits() -> Result<(), String> {
        assert_eq!(pi_scaled(8)?, BigInt::from(804));
        Ok(())
    }
}
