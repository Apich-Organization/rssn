use crate::symbolic::core::Expr;
use crate::symbolic::number_theory::{
    chinese_remainder, extended_gcd, is_prime, solve_diophantine,
};
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};
use std::slice;

/// Solves a Diophantine equation.
///
/// # Safety
/// `equation` must be a valid pointer to an `Expr`.
/// `vars_ptr` must be a valid pointer to an array of C strings of length `vars_len`.
#[no_mangle]
pub extern "C" fn rssn_solve_diophantine_handle(
    equation: *const Expr,
    vars_ptr: *const *const c_char,
    vars_len: c_int,
) -> *mut Expr {
    let equation_ref = unsafe {
        if equation.is_null() {
            return std::ptr::null_mut();
        }
        &*equation
    };

    let vars: Vec<String> = unsafe {
        slice::from_raw_parts(vars_ptr, vars_len as usize)
            .iter()
            .map(|&p| CStr::from_ptr(p).to_string_lossy().into_owned())
            .collect()
    };
    let vars_str: Vec<&str> = vars.iter().map(|s| s.as_str()).collect();

    match solve_diophantine(equation_ref, &vars_str) {
        Ok(solutions) => {
            let result = Expr::new_vector(solutions);
            Box::into_raw(Box::new(result))
        }
        Err(_) => std::ptr::null_mut(),
    }
}

/// Computes the Extended GCD of two expressions.
///
/// # Safety
/// `a` and `b` must be valid pointers to `Expr`.
#[no_mangle]
pub extern "C" fn rssn_extended_gcd_handle(a: *const Expr, b: *const Expr) -> *mut Expr {
    let a_ref = unsafe {
        if a.is_null() {
            return std::ptr::null_mut();
        }
        &*a
    };
    let b_ref = unsafe {
        if b.is_null() {
            return std::ptr::null_mut();
        }
        &*b
    };

    let (g, x, y) = extended_gcd(a_ref, b_ref);
    let result = Expr::new_tuple(vec![g, x, y]);
    Box::into_raw(Box::new(result))
}

/// Checks if an expression is prime.
///
/// # Safety
/// `n` must be a valid pointer to an `Expr`.
#[no_mangle]
pub extern "C" fn rssn_is_prime_handle(n: *const Expr) -> *mut Expr {
    let n_ref = unsafe {
        if n.is_null() {
            return std::ptr::null_mut();
        }
        &*n
    };

    let result = is_prime(n_ref);
    Box::into_raw(Box::new(result))
}

/// Solves a system of congruences using the Chinese Remainder Theorem.
///
/// # Safety
/// `remainders` and `moduli` must be valid pointers to arrays of `Expr` pointers of length `len`.
#[no_mangle]
pub extern "C" fn rssn_chinese_remainder_handle(
    remainders: *const *const Expr,
    moduli: *const *const Expr,
    len: c_int,
) -> *mut Expr {
    let remainders_slice = unsafe { slice::from_raw_parts(remainders, len as usize) };
    let moduli_slice = unsafe { slice::from_raw_parts(moduli, len as usize) };

    let mut congruences = Vec::new();
    for i in 0..len as usize {
        let r_ptr = remainders_slice[i];
        let m_ptr = moduli_slice[i];
        if r_ptr.is_null() || m_ptr.is_null() {
            return std::ptr::null_mut();
        }
        let r = unsafe { &*r_ptr };
        let m = unsafe { &*m_ptr };
        congruences.push((r.clone(), m.clone()));
    }

    match chinese_remainder(&congruences) {
        Some(result) => Box::into_raw(Box::new(result)),
        None => std::ptr::null_mut(),
    }
}
