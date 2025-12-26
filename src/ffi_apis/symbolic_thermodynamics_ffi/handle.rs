//! Handle-based FFI API for thermodynamics functions.

use crate::symbolic::core::Expr;
use crate::symbolic::thermodynamics;
use std::os::raw::c_char;

/// Calculates ideal gas Law expression: PV - nRT.
#[no_mangle]

pub unsafe extern "C" fn rssn_ideal_gas_law(
    p: *const Expr,
    v: *const Expr,
    n: *const Expr,
    r: *const Expr,
    t: *const Expr,
) -> *mut Expr {

    if p.is_null() || v.is_null() || n.is_null() || r.is_null() || t.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(thermodynamics::ideal_gas_law(
        &*p, &*v, &*n, &*r, &*t,
    )))
}

/// Calculates enthalpy: U + PV.
#[no_mangle]

pub unsafe extern "C" fn rssn_enthalpy(
    u: *const Expr,
    p: *const Expr,
    v: *const Expr,
) -> *mut Expr {

    if u.is_null() || p.is_null() || v.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(thermodynamics::enthalpy(&*u, &*p, &*v)))
}

/// Calculates Gibbs Free Energy: H - TS.
#[no_mangle]

pub unsafe extern "C" fn rssn_gibbs_free_energy(
    h: *const Expr,
    t: *const Expr,
    s: *const Expr,
) -> *mut Expr {

    if h.is_null() || t.is_null() || s.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(thermodynamics::gibbs_free_energy(&*h, &*t, &*s)))
}

/// Calculates Carnot Efficiency: 1 - Tc/Th.
#[no_mangle]

pub unsafe extern "C" fn rssn_carnot_efficiency(
    tc: *const Expr,
    th: *const Expr,
) -> *mut Expr {

    if tc.is_null() || th.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(thermodynamics::carnot_efficiency(&*tc, &*th)))
}

/// Calculates Boltzmann Distribution.
#[no_mangle]

pub unsafe extern "C" fn rssn_boltzmann_distribution(
    energy: *const Expr,
    temperature: *const Expr,
    partition_function: *const Expr,
) -> *mut Expr {

    if energy.is_null() || temperature.is_null() || partition_function.is_null() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(thermodynamics::boltzmann_distribution(
        &*energy,
        &*temperature,
        &*partition_function,
    )))
}
