//! Handle-based FFI API for numerical number theory operations.

use crate::numerical::number_theory as nt;
use std::ptr;

/// Computes the greatest common divisor (GCD).
#[no_mangle]

pub extern "C" fn rssn_num_nt_gcd(
    a: u64,
    b: u64,
) -> u64 {

    nt::gcd(a, b)
}

/// Computes the least common multiple (LCM).
#[no_mangle]

pub extern "C" fn rssn_num_nt_lcm(
    a: u64,
    b: u64,
) -> u64 {

    nt::lcm(a, b)
}

/// Computes (base^exp) % modulus.
#[no_mangle]

pub extern "C" fn rssn_num_nt_mod_pow(
    base: u128,
    exp: u64,
    modulus: u64,
) -> u64 {

    nt::mod_pow(base, exp, modulus)
}

/// Finds the modular multiplicative inverse.
/// Returns 0 if no inverse exists (modulus cannot be 0).
#[no_mangle]

pub extern "C" fn rssn_num_nt_mod_inverse(
    a: i64,
    m: i64,
) -> i64 {

    nt::mod_inverse(a, m).unwrap_or(0)
}

/// Tests if a number is prime using Miller-Rabin.
#[no_mangle]

pub extern "C" fn rssn_num_nt_is_prime(
    n: u64
) -> bool {

    nt::is_prime_miller_rabin(n)
}

/// Computes Euler's totient function φ(n).
#[no_mangle]

pub extern "C" fn rssn_num_nt_phi(
    n: u64
) -> u64 {

    nt::phi(n)
}

/// Returns the number of prime factors and writes them to out_factors.
/// out_factors must be large enough.
#[no_mangle]

pub unsafe extern "C" fn rssn_num_nt_factorize(
    n: u64,
    out_factors: *mut u64,
) -> usize {

    let factors = nt::factorize(n);

    if !out_factors.is_null() {

        ptr::copy_nonoverlapping(
            factors.as_ptr(),
            out_factors,
            factors.len(),
        );
    }

    factors.len()
}
