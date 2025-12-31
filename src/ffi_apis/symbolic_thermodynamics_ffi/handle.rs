//! Handle-based FFI API for thermodynamics functions.

use crate::symbolic::core::Expr;
use crate::symbolic::thermodynamics;

/// Calculates ideal gas Law expression: PV - nRT.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_ideal_gas_law(
    p: *const Expr,
    v: *const Expr,
    n: *const Expr,
    r: *const Expr,
    t: *const Expr,
) -> *mut Expr {

    unsafe {

        if p.is_null()
            || v.is_null()
            || n.is_null()
            || r.is_null()
            || t.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        thermodynamics::ideal_gas_law(
            &*p, &*v, &*n, &*r, &*t,
        ),
    ))
    }
}

/// Calculates enthalpy: U + PV.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_enthalpy(
    u: *const Expr,
    p: *const Expr,
    v: *const Expr,
) -> *mut Expr {

    unsafe {

        if u.is_null()
            || p.is_null()
            || v.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
            thermodynamics::enthalpy(
                &*u, &*p, &*v,
            ),
        ))
    }
}

/// Calculates Gibbs Free Energy: H - TS.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_gibbs_free_energy(
    h: *const Expr,
    t: *const Expr,
    s: *const Expr,
) -> *mut Expr {

    unsafe {

        if h.is_null()
            || t.is_null()
            || s.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        thermodynamics::gibbs_free_energy(&*h, &*t, &*s),
    ))
    }
}

/// Calculates Carnot Efficiency: 1 - Tc/Th.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_carnot_efficiency(
    tc: *const Expr,
    th: *const Expr,
) -> *mut Expr {

    unsafe {

        if tc.is_null() || th.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        thermodynamics::carnot_efficiency(&*tc, &*th),
    ))
    }
}

/// Calculates Boltzmann Distribution.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_boltzmann_distribution(
    energy: *const Expr,
    temperature: *const Expr,
    partition_function: *const Expr,
) -> *mut Expr {

    unsafe {

        if energy.is_null()
            || temperature.is_null()
            || partition_function
                .is_null()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
        thermodynamics::boltzmann_distribution(
            &*energy,
            &*temperature,
            &*partition_function,
        ),
    ))
    }
}
