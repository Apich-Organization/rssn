//! Session-owned modular reducers for the checked integer ABI.
use super::native_object_api as api;
use crate::graph::RuleError;
use crate::graph::dim::DimRule;
use crate::graph::rule::Installer;
use crate::kernel::OpError;
use crate::kernels::modular::Barrett;
use crate::kernels::modular::Montgomery;
use crate::kernels::number_theory::is_prime_miller_rabin;
use api::Handle;
use api::NativeType;
use api::NewObject;
use api::error;
use api::object_binding;
use std::mem::size_of;

impl NativeType for Montgomery {
    const KIND: &'static str = "rssn/native/modular/montgomery/v1";
    const MAX_BYTES: usize = size_of::<Self>();

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

impl NativeType for Barrett {
    const KIND: &'static str = "rssn/native/modular/barrett/v1";
    const MAX_BYTES: usize = size_of::<Self>();

    fn validate(&self) -> Result<(), OpError> {
        Ok(())
    }
}

fn require_residue(
    value: u64,
    modulus: u64,
) -> Result<(), OpError> {
    if value >= modulus {
        return Err(error("encoded residue must be smaller than the modulus"));
    }
    Ok(())
}

// Declare the method and its ABI metadata together. The object scope and exact
// integer conversions come from the existing checked native binding.
macro_rules! reducer_op {
    ($plan:ident, $name:ident, $op:literal, $source:literal, $doc:literal,
     $method:ident ($($arg:ident : $ty:ty),* $(,)?) -> $ret:ty) => {
        object_binding!(
            $name, $op, $source, $doc, true,
            (plan: Handle<$plan> $(, $arg : $ty)*) -> $ret,
            |ctx, plan $(, $arg)*| plan.with(ctx, |plan| Ok(plan.$method($($arg),*)))
        );
    };
}

macro_rules! checked_montgomery {
    ($name:ident, $op:literal, $source:literal, $doc:literal,
     $method:ident ($($arg:ident : $ty:ty),+), check ($($checked:ident),+)) => {
        object_binding!(
            $name, $op, $source, $doc, true,
            (plan: Handle<Montgomery>, $($arg : $ty),+) -> u64,
            |ctx, plan, $($arg),+| plan.with(ctx, |plan| {
                $(require_residue($checked, plan.modulus())?;)+
                Ok(plan.$method($($arg),+))
            })
        );
    };
}

object_binding!(
    montgomery_new,
    "num_modular_montgomery_new",
    "crate::kernels::modular::Montgomery::new",
    "`num_modular_montgomery_new(modulus)` prepares reusable reduction constants for an unsigned 64-bit odd modulus greater than one.",
    true,
    (modulus: u64) -> NewObject<Montgomery>,
    |_, modulus| Montgomery::new(modulus)
        .map(NewObject::normalized)
        .ok_or_else(|| error("Montgomery requires an odd modulus greater than one"))
);
object_binding!(
    montgomery_new_prime,
    "num_modular_montgomery_new_prime",
    "crate::kernels::modular::Montgomery::new_prime",
    "`num_modular_montgomery_new_prime(modulus)` certifies an odd prime before preparing reduction constants.",
    true,
    (modulus: u64) -> NewObject<Montgomery>,
    |_, modulus| {
        if !(modulus > 2 && is_prime_miller_rabin(modulus)) {
            return Err(error("Montgomery prime plan requires an odd prime"));
        }
        Ok(NewObject::normalized(Montgomery::new_prime(modulus)))
    }
);

reducer_op!(
    Montgomery, montgomery_modulus, "num_modular_montgomery_modulus",
    "crate::kernels::modular::Montgomery::modulus",
    "`num_modular_montgomery_modulus(plan)` Returns the prepared odd modulus.",
    modulus() -> u64
);
reducer_op!(
    Montgomery, montgomery_n_prime, "num_modular_montgomery_n_prime",
    "crate::kernels::modular::Montgomery::n_prime",
    "`num_modular_montgomery_n_prime(plan)` Returns the Montgomery reduction multiplier -p^(-1) mod 2^64.",
    n_prime() -> u64
);
reducer_op!(
    Montgomery, montgomery_r_mod, "num_modular_montgomery_r_mod",
    "crate::kernels::modular::Montgomery::r_mod",
    "`num_modular_montgomery_r_mod(plan)` Returns 2^64 reduced modulo the prepared modulus.",
    r_mod() -> u64
);
reducer_op!(
    Montgomery, montgomery_r2_mod, "num_modular_montgomery_r2_mod",
    "crate::kernels::modular::Montgomery::r2_mod",
    "`num_modular_montgomery_r2_mod(plan)` Returns 2^128 reduced modulo the prepared modulus.",
    r2_mod() -> u64
);
reducer_op!(
    Montgomery, montgomery_one, "num_modular_montgomery_one",
    "crate::kernels::modular::Montgomery::one",
    "`num_modular_montgomery_one(plan)` Returns the multiplicative identity in standard representation.",
    one() -> u64
);
reducer_op!(
    Montgomery, montgomery_one_mont, "num_modular_montgomery_one_mont",
    "crate::kernels::modular::Montgomery::one_mont",
    "`num_modular_montgomery_one_mont(plan)` Returns the Montgomery representation of the multiplicative identity.",
    one_mont() -> u64
);
object_binding!(
    montgomery_redc,
    "num_modular_montgomery_redc",
    "crate::kernels::modular::Montgomery::redc",
    "`num_modular_montgomery_redc(plan, value)` Returns value times the inverse Montgomery radix modulo p. Requires value < p * 2^64.",
    true,
    (plan: Handle<Montgomery>, value: u128) -> u64,
    |ctx, plan, value| plan.with(ctx, |plan| {
        if (value >> 64) >= u128::from(plan.modulus()) {
            return Err(error("Montgomery reduction requires value < modulus * 2^64"));
        }
        Ok(plan.redc(value))
    })
);
reducer_op!(
    Montgomery,
    montgomery_to_mont,
    "num_modular_montgomery_to_mont",
    "crate::kernels::modular::Montgomery::to_mont",
    "`num_modular_montgomery_to_mont(plan, value)` Converts an unsigned 64-bit word to a reduced Montgomery residue.",
    to_mont(value: u64) -> u64
);
reducer_op!(
    Montgomery,
    montgomery_to_mont_u128,
    "num_modular_montgomery_to_mont_u128",
    "crate::kernels::modular::Montgomery::to_mont_u128",
    "`num_modular_montgomery_to_mont_u128(plan, value)` Converts an unsigned 128-bit word to a reduced Montgomery residue.",
    to_mont_u128(value: u128) -> u64
);
checked_montgomery!(
    montgomery_from_mont,
    "num_modular_montgomery_from_mont",
    "crate::kernels::modular::Montgomery::from_mont",
    "`num_modular_montgomery_from_mont(plan, value)` Converts a canonical Montgomery residue to standard representation.",
    from_mont(value: u64),
    check(value)
);
checked_montgomery!(
    montgomery_add_mont,
    "num_modular_montgomery_add_mont",
    "crate::kernels::modular::Montgomery::add_mont",
    "`num_modular_montgomery_add_mont(plan, a, b)` Adds two canonical Montgomery residues.",
    add_mont(a: u64, b: u64),
    check(a, b)
);
checked_montgomery!(
    montgomery_sub_mont,
    "num_modular_montgomery_sub_mont",
    "crate::kernels::modular::Montgomery::sub_mont",
    "`num_modular_montgomery_sub_mont(plan, a, b)` Subtracts two canonical Montgomery residues.",
    sub_mont(a: u64, b: u64),
    check(a, b)
);
checked_montgomery!(
    montgomery_neg_mont,
    "num_modular_montgomery_neg_mont",
    "crate::kernels::modular::Montgomery::neg_mont",
    "`num_modular_montgomery_neg_mont(plan, value)` Negates a canonical Montgomery residue.",
    neg_mont(value: u64),
    check(value)
);
checked_montgomery!(
    montgomery_mul_mont,
    "num_modular_montgomery_mul_mont",
    "crate::kernels::modular::Montgomery::mul_mont",
    "`num_modular_montgomery_mul_mont(plan, a, b)` Multiplies two canonical Montgomery residues.",
    mul_mont(a: u64, b: u64),
    check(a, b)
);
checked_montgomery!(
    montgomery_fma_mont,
    "num_modular_montgomery_fma_mont",
    "crate::kernels::modular::Montgomery::fma_mont",
    "`num_modular_montgomery_fma_mont(plan, a, b, accumulator)` Adds the product to the accumulator in Montgomery representation. All three inputs must be canonical.",
    fma_mont(a: u64, b: u64, accumulator: u64),
    check(a, b, accumulator)
);
checked_montgomery!(
    montgomery_half,
    "num_modular_montgomery_half",
    "crate::kernels::modular::Montgomery::half",
    "`num_modular_montgomery_half(plan, value)` Halves a canonical residue for the odd modulus, in either representation.",
    half(value: u64),
    check(value)
);
checked_montgomery!(
    montgomery_pow_mont,
    "num_modular_montgomery_pow_mont",
    "crate::kernels::modular::Montgomery::pow_mont",
    "`num_modular_montgomery_pow_mont(plan, base, exponent)` Raises a canonical Montgomery residue to an unsigned 64-bit exponent.",
    pow_mont(base: u64, exponent: u64),
    check(base)
);
reducer_op!(
    Montgomery,
    montgomery_add,
    "num_modular_montgomery_add",
    "crate::kernels::modular::Montgomery::add",
    "`num_modular_montgomery_add(plan, a, b)` Adds two unsigned words modulo the prepared modulus.",
    add(a: u64, b: u64) -> u64
);
reducer_op!(
    Montgomery,
    montgomery_sub,
    "num_modular_montgomery_sub",
    "crate::kernels::modular::Montgomery::sub",
    "`num_modular_montgomery_sub(plan, a, b)` Subtracts two unsigned words modulo the prepared modulus.",
    sub(a: u64, b: u64) -> u64
);
reducer_op!(
    Montgomery,
    montgomery_neg,
    "num_modular_montgomery_neg",
    "crate::kernels::modular::Montgomery::neg",
    "`num_modular_montgomery_neg(plan, value)` Negates an unsigned word modulo the prepared modulus.",
    neg(value: u64) -> u64
);
reducer_op!(
    Montgomery,
    montgomery_mul,
    "num_modular_montgomery_mul",
    "crate::kernels::modular::Montgomery::mul",
    "`num_modular_montgomery_mul(plan, a, b)` Multiplies two unsigned words modulo the prepared modulus.",
    mul(a: u64, b: u64) -> u64
);
reducer_op!(
    Montgomery,
    montgomery_pow,
    "num_modular_montgomery_pow",
    "crate::kernels::modular::Montgomery::pow",
    "`num_modular_montgomery_pow(plan, base, exponent)` Raises an unsigned word to an unsigned 64-bit exponent modulo the prepared modulus.",
    pow(base: u64, exponent: u64) -> u64
);
reducer_op!(
    Montgomery,
    montgomery_inv,
    "num_modular_montgomery_inv",
    "crate::kernels::modular::Montgomery::inv",
    "`num_modular_montgomery_inv(plan, value)` Returns the modular inverse, or None when the value is not coprime to the modulus.",
    inv(value: u64) -> Option<u64>
);
reducer_op!(
    Montgomery,
    montgomery_div,
    "num_modular_montgomery_div",
    "crate::kernels::modular::Montgomery::div",
    "`num_modular_montgomery_div(plan, numerator, denominator)` Returns the modular quotient, or None when the denominator is not invertible.",
    div(numerator: u64, denominator: u64) -> Option<u64>
);
object_binding!(
    barrett_new,
    "num_modular_barrett_new",
    "crate::kernels::modular::Barrett::new",
    "`num_modular_barrett_new(modulus)` prepares reusable reduction constants for an unsigned 64-bit modulus of at least two.",
    true,
    (modulus: u64) -> NewObject<Barrett>,
    |_, modulus| Barrett::new(modulus)
        .map(NewObject::normalized)
        .ok_or_else(|| error("Barrett requires a modulus of at least two"))
);
object_binding!(
    barrett_new_prime,
    "num_modular_barrett_new_prime",
    "crate::kernels::modular::Barrett::new_prime",
    "`num_modular_barrett_new_prime(modulus)` certifies a prime before preparing reduction constants.",
    true,
    (modulus: u64) -> NewObject<Barrett>,
    |_, modulus| {
        if !(is_prime_miller_rabin(modulus)) {
            return Err(error("Barrett prime plan requires a prime"));
        }
        Ok(NewObject::normalized(Barrett::new_prime(modulus)))
    }
);

reducer_op!(
    Barrett, barrett_modulus, "num_modular_barrett_modulus",
    "crate::kernels::modular::Barrett::modulus",
    "`num_modular_barrett_modulus(plan)` Returns the prepared modulus.",
    modulus() -> u64
);
reducer_op!(
    Barrett, barrett_mu, "num_modular_barrett_mu",
    "crate::kernels::modular::Barrett::mu",
    "`num_modular_barrett_mu(plan)` Returns floor(2^64 / modulus), the high word of the reciprocal.",
    mu() -> u64
);
reducer_op!(
    Barrett,
    barrett_reduce_u64,
    "num_modular_barrett_reduce_u64",
    "crate::kernels::modular::Barrett::reduce_u64",
    "`num_modular_barrett_reduce_u64(plan, value)` Reduces an unsigned 64-bit word modulo the prepared modulus.",
    reduce_u64(value: u64) -> u64
);
reducer_op!(
    Barrett,
    barrett_reduce_u128,
    "num_modular_barrett_reduce_u128",
    "crate::kernels::modular::Barrett::reduce_u128",
    "`num_modular_barrett_reduce_u128(plan, value)` Reduces an unsigned 128-bit word modulo the prepared modulus.",
    reduce_u128(value: u128) -> u64
);
reducer_op!(
    Barrett,
    barrett_add,
    "num_modular_barrett_add",
    "crate::kernels::modular::Barrett::add",
    "`num_modular_barrett_add(plan, a, b)` Adds two unsigned words modulo the prepared modulus.",
    add(a: u64, b: u64) -> u64
);
reducer_op!(
    Barrett,
    barrett_sub,
    "num_modular_barrett_sub",
    "crate::kernels::modular::Barrett::sub",
    "`num_modular_barrett_sub(plan, a, b)` Subtracts two unsigned words modulo the prepared modulus.",
    sub(a: u64, b: u64) -> u64
);
reducer_op!(
    Barrett,
    barrett_neg,
    "num_modular_barrett_neg",
    "crate::kernels::modular::Barrett::neg",
    "`num_modular_barrett_neg(plan, value)` Negates an unsigned word modulo the prepared modulus.",
    neg(value: u64) -> u64
);
reducer_op!(
    Barrett,
    barrett_mul,
    "num_modular_barrett_mul",
    "crate::kernels::modular::Barrett::mul",
    "`num_modular_barrett_mul(plan, a, b)` Multiplies two unsigned words modulo the prepared modulus.",
    mul(a: u64, b: u64) -> u64
);
reducer_op!(
    Barrett,
    barrett_pow,
    "num_modular_barrett_pow",
    "crate::kernels::modular::Barrett::pow",
    "`num_modular_barrett_pow(plan, base, exponent)` Raises an unsigned word to an unsigned 64-bit exponent modulo the prepared modulus.",
    pow(base: u64, exponent: u64) -> u64
);
reducer_op!(
    Barrett,
    barrett_inv,
    "num_modular_barrett_inv",
    "crate::kernels::modular::Barrett::inv",
    "`num_modular_barrett_inv(plan, value)` Returns the modular inverse, or None when the value is not coprime to the modulus.",
    inv(value: u64) -> Option<u64>
);
reducer_op!(
    Barrett,
    barrett_div,
    "num_modular_barrett_div",
    "crate::kernels::modular::Barrett::div",
    "`num_modular_barrett_div(plan, numerator, denominator)` Returns the modular quotient, or None when the denominator is not invertible.",
    div(numerator: u64, denominator: u64) -> Option<u64>
);
#[allow(
    clippy::single_call_fn,
    reason = "Keeps reducer registrations and their dimension policy with the native adapter"
)]
pub(super) fn install(i: &mut Installer<'_>) -> Result<(), RuleError> {
    api::install_with_dimension(
        i,
        &[
            &montgomery_new::BINDING,
            &montgomery_new_prime::BINDING,
            &montgomery_modulus::BINDING,
            &montgomery_n_prime::BINDING,
            &montgomery_r_mod::BINDING,
            &montgomery_r2_mod::BINDING,
            &montgomery_one::BINDING,
            &montgomery_one_mont::BINDING,
            &montgomery_redc::BINDING,
            &montgomery_to_mont::BINDING,
            &montgomery_to_mont_u128::BINDING,
            &montgomery_from_mont::BINDING,
            &montgomery_add_mont::BINDING,
            &montgomery_sub_mont::BINDING,
            &montgomery_neg_mont::BINDING,
            &montgomery_mul_mont::BINDING,
            &montgomery_fma_mont::BINDING,
            &montgomery_half::BINDING,
            &montgomery_pow_mont::BINDING,
            &montgomery_add::BINDING,
            &montgomery_sub::BINDING,
            &montgomery_neg::BINDING,
            &montgomery_mul::BINDING,
            &montgomery_pow::BINDING,
            &montgomery_inv::BINDING,
            &montgomery_div::BINDING,
            &barrett_new::BINDING,
            &barrett_new_prime::BINDING,
            &barrett_modulus::BINDING,
            &barrett_mu::BINDING,
            &barrett_reduce_u64::BINDING,
            &barrett_reduce_u128::BINDING,
            &barrett_add::BINDING,
            &barrett_sub::BINDING,
            &barrett_neg::BINDING,
            &barrett_mul::BINDING,
            &barrett_pow::BINDING,
            &barrett_inv::BINDING,
            &barrett_div::BINDING,
        ],
        "numerical_special",
        DimRule::Dimensionless,
    )
}
