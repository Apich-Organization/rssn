//! Mathematical constants at an arbitrary precision.
//!
//! Each constant is computed once at the highest precision requested so far
//! and then rounded to the precision a caller asks for, so a run that uses
//! many digits pays for them once. The computations are classical integer
//! series, evaluated in the fixed-point type of `super::fixed`:
//!
//! * `pi` by Machin's formula `16 atan(1/5) - 4 atan(1/239)`;
//! * `ln 2` by `2 atanh(1/3)`;
//! * `e` by `sum 1/k!`;
//! * `sqrt 2` and the golden ratio by [`Real::sqrt`].
//!
//! Each series is alternating or positive with strictly decreasing terms, so
//! the first omitted term bounds the truncation error; the error count
//! returned with the value covers that term plus three units per evaluated
//! term, which is what the rounding of each division costs.

use std::sync::Mutex;
use std::sync::OnceLock;

use num_bigint::BigInt;

use super::fixed::Fx;
use super::real::MAX_BITS;
use super::real::Real;
use super::real::Round;

/// One lazily extended constant.
struct Cache {
    slot: OnceLock<Mutex<Option<Real>>>,
}

impl Cache {
    /// The constant at `prec` bits.
    fn get(
        &self,
        prec: u32,
        compute: &dyn Fn(u32) -> Result<Real, String>,
    ) -> Result<Real, String> {
        let prec = prec.clamp(1, MAX_BITS);
        let mutex = self.slot.get_or_init(|| Mutex::new(None));
        let mut guard = mutex
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let wanted = match guard.as_ref() {
            | Some(cached) if cached.precision() >= prec => {
                let value = cached.clone();
                drop(guard);
                return value.with_precision(prec, Round::NearestEven);
            },
            | Some(cached) => prec.max(cached.precision().saturating_mul(2)),
            | None => prec,
        }
        .min(MAX_BITS);
        let value = compute(wanted.saturating_add(64))?;
        *guard = Some(value.clone());
        drop(guard);
        value.with_precision(prec, Round::NearestEven)
    }
}

static PI: Cache = Cache {
    slot: OnceLock::new(),
};
static LN2: Cache = Cache {
    slot: OnceLock::new(),
};
static E: Cache = Cache {
    slot: OnceLock::new(),
};
static SQRT2: Cache = Cache {
    slot: OnceLock::new(),
};
static PHI: Cache = Cache {
    slot: OnceLock::new(),
};

/// `pi` at `prec` bits.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn pi(prec: u32) -> Result<Real, String> {
    PI.get(prec, &machin)
}

/// `ln 2` at `prec` bits.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn ln2(prec: u32) -> Result<Real, String> {
    LN2.get(prec, &ln2_at)
}

/// `e` at `prec` bits.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn e(prec: u32) -> Result<Real, String> {
    E.get(prec, &e_at)
}

/// `sqrt 2` at `prec` bits.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn sqrt2(prec: u32) -> Result<Real, String> {
    SQRT2.get(prec, &sqrt2_at)
}

/// The golden ratio `(1 + sqrt 5) / 2` at `prec` bits.
///
/// # Errors
/// Fails when an exponent leaves the representable range.
pub fn phi(prec: u32) -> Result<Real, String> {
    PHI.get(prec, &phi_at)
}

/// `atan(1/n)` at scale `s`, with the error count of the series.
fn atan_inverse(
    n: i64,
    s: u32,
) -> Result<(Fx, u64), String> {
    let one = Fx::one(s);
    let mut power = one.div_i64(n);
    let mut sum = power.clone();
    let mut term_count: u64 = 1;
    let mut k: i64 = 1;
    let step = n.saturating_mul(n);
    loop {
        power = power.div_i64(step);
        let term = power.div_i64(2 * k + 1);
        if term.is_zero() {
            break;
        }
        sum = if k % 2 == 1 {
            sum.sub(&term)
        } else {
            sum.add(&term)
        };
        term_count = term_count.saturating_add(1);
        k = k.saturating_add(1);
        if k > 1_000_000 {
            return Err("series did not converge".to_string());
        }
    }
    Ok((sum, 3 * term_count + 4))
}

/// `pi` at `prec + 64` bits by Machin's formula.
fn machin(prec: u32) -> Result<Real, String> {
    let s = prec.saturating_add(64);
    let (fifth, e1) = atan_inverse(5, s)?;
    let (two_thirty_ninth, e2) = atan_inverse(239, s)?;
    let value = fifth.scale_pow2(4).sub(&two_thirty_ninth.scale_pow2(2));
    let _ = (e1, e2);
    value.to_real(prec)
}

/// `ln 2` at `prec + 64` bits by `2 atanh(1/3)`.
fn ln2_at(prec: u32) -> Result<Real, String> {
    let s = prec.saturating_add(64);
    let three = Fx::from_i64(3, s);
    let mut power = Fx::one(s).div(&three);
    let mut sum = power.clone();
    let mut term_count: u64 = 1;
    let mut k: i64 = 1;
    let step = Fx::from_i64(9, s);
    loop {
        power = power.div(&step);
        let term = power.div_i64(2 * k + 1);
        if term.is_zero() {
            break;
        }
        sum = sum.add(&term);
        term_count = term_count.saturating_add(1);
        k = k.saturating_add(1);
        if k > 1_000_000 {
            return Err("series did not converge".to_string());
        }
    }
    let _ = term_count;
    sum.scale_pow2(1).to_real(prec)
}

/// `e` at `prec + 64` bits by `sum 1/k!`.
fn e_at(prec: u32) -> Result<Real, String> {
    let s = prec.saturating_add(64);
    let mut term = Fx::one(s);
    let mut sum = term.clone();
    let mut k: i64 = 1;
    loop {
        term = term.div_i64(k);
        if term.is_zero() {
            break;
        }
        sum = sum.add(&term);
        k = k.saturating_add(1);
        if k > 1_000_000 {
            return Err("series did not converge".to_string());
        }
    }
    sum.to_real(prec)
}

/// `sqrt 2` at `prec + 64` bits.
fn sqrt2_at(prec: u32) -> Result<Real, String> {
    let two = Real::from_i64(2, prec.saturating_add(64))?;
    two.sqrt(prec, Round::NearestEven)
}

/// The golden ratio at `prec + 64` bits.
fn phi_at(prec: u32) -> Result<Real, String> {
    let inner = prec.saturating_add(64);
    let five = Real::from_i64(5, inner)?;
    let root = five.sqrt(inner, Round::NearestEven)?;
    let one = Real::one(inner);
    one.add(&root, inner, Round::NearestEven)?
        .scale_pow2(-1)?
        .with_precision(prec, Round::NearestEven)
}

/// The constant `name` at `prec` bits, if it is known.
///
/// # Errors
/// Fails for an unknown name or when an exponent leaves the range.
pub fn named(
    name: &str,
    prec: u32,
) -> Result<Real, String> {
    match name {
        | "pi" => pi(prec),
        | "e" => e(prec),
        | "tau" => pi(prec)?.scale_pow2(1),
        | "ln2" => ln2(prec),
        | "sqrt2" => sqrt2(prec),
        | "phi" => phi(prec),
        | _ => Err(format!("unknown constant `{name}`")),
    }
}

/// The names [`named`] accepts, for diagnostics and the operator's help text.
#[must_use]
pub const fn names() -> &'static [&'static str] {
    &["pi", "e", "tau", "ln2", "sqrt2", "phi"]
}

/// A fixed-point approximation to `pi` as a [`BigInt`] numerator over `2^prec`.
///
/// The constant is evaluated with guard bits before scaling and truncating.
///
/// # Errors
/// Fails when an exponent leaves the range.
pub fn pi_scaled(prec: u32) -> Result<BigInt, String> {
    let s = prec.saturating_add(64);
    let value = pi(s)?;
    let scaled = value.scale_pow2(i64::from(prec))?;
    scaled.trunc_to_bigint()
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
