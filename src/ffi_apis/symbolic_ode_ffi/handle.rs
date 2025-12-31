//! Handle-based FFI API for symbolic ODE functions.

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::symbolic::core::Expr;
use crate::symbolic::ode;

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

unsafe fn c_str_to_str<'a>(
    s: *const c_char
) -> Option<&'a str> {

    unsafe {

        if s.is_null() {

            None
        } else {

            CStr::from_ptr(s)
                .to_str()
                .ok()
        }
    }
}

/// Solves an ordinary differential equation.
///
/// # Safety
/// The caller must ensure `ode_expr` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_ode(
    ode_expr: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    unsafe {

        if ode_expr.is_null()
            || func.is_null()
            || var.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let ode_ref = &*ode_expr;

        let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        let var_str = match c_str_to_str(
        var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        Box::into_raw(Box::new(
            ode::solve_ode(
                ode_ref,
                func_str,
                var_str,
                None,
            ),
        ))
    }
}

/// Solves a separable ODE.
///
/// # Safety
/// The caller must ensure `equation` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_separable_ode(
    equation: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    unsafe {

        if equation.is_null()
            || func.is_null()
            || var.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let eq_ref = &*equation;

        let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        let var_str = match c_str_to_str(
        var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        match ode::solve_separable_ode(
            eq_ref,
            func_str,
            var_str,
        ) {
            | Some(solution) => {
                Box::into_raw(Box::new(
                    solution,
                ))
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Solves a first-order linear ODE.
///
/// # Safety
/// The caller must ensure `equation` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_first_order_linear_ode(
    equation: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    unsafe {

        if equation.is_null()
            || func.is_null()
            || var.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let eq_ref = &*equation;

        let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        let var_str = match c_str_to_str(
        var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        match ode::solve_first_order_linear_ode(
        eq_ref,
        func_str,
        var_str,
    ) {
        | Some(solution) => Box::into_raw(Box::new(solution)),
        | None => std::ptr::null_mut(),
    }
    }
}

/// Solves a Bernoulli ODE.
///
/// # Safety
/// The caller must ensure `equation` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_bernoulli_ode(
    equation: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    unsafe {

        if equation.is_null()
            || func.is_null()
            || var.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let eq_ref = &*equation;

        let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        let var_str = match c_str_to_str(
        var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        match ode::solve_bernoulli_ode(
            eq_ref,
            func_str,
            var_str,
        ) {
            | Some(solution) => {
                Box::into_raw(Box::new(
                    solution,
                ))
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Solves a Riccati ODE with a known particular solution.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_riccati_ode(
    equation: *const Expr,
    func: *const c_char,
    var: *const c_char,
    y1: *const Expr,
) -> *mut Expr {

    unsafe {

        if equation.is_null()
            || func.is_null()
            || var.is_null()
            || y1.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let eq_ref = &*equation;

        let y1_ref = &*y1;

        let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        let var_str = match c_str_to_str(
        var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        match ode::solve_riccati_ode(
            eq_ref,
            func_str,
            var_str,
            y1_ref,
        ) {
            | Some(solution) => {
                Box::into_raw(Box::new(
                    solution,
                ))
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Solves a Cauchy-Euler ODE.
///
/// # Safety
/// The caller must ensure `equation` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_cauchy_euler_ode(
    equation: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    unsafe {

        if equation.is_null()
            || func.is_null()
            || var.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let eq_ref = &*equation;

        let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        let var_str = match c_str_to_str(
        var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        match ode::solve_cauchy_euler_ode(
        eq_ref,
        func_str,
        var_str,
    ) {
        | Some(solution) => {
            Box::into_raw(Box::new(
                solution,
            ))
        },
        | None => std::ptr::null_mut(),
    }
    }
}

/// Solves an exact ODE.
///
/// # Safety
/// The caller must ensure `equation` is a valid Expr pointer, and `func` and `var` are valid C strings.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_exact_ode(
    equation: *const Expr,
    func: *const c_char,
    var: *const c_char,
) -> *mut Expr {

    unsafe {

        if equation.is_null()
            || func.is_null()
            || var.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let eq_ref = &*equation;

        let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        let var_str = match c_str_to_str(
        var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        match ode::solve_exact_ode(
            eq_ref,
            func_str,
            var_str,
        ) {
            | Some(solution) => {
                Box::into_raw(Box::new(
                    solution,
                ))
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Solves a second-order ODE by reduction of order with a known solution.
///
/// # Safety
/// The caller must ensure all pointers are valid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_solve_by_reduction_of_order(
    equation: *const Expr,
    func: *const c_char,
    var: *const c_char,
    y1: *const Expr,
) -> *mut Expr {

    unsafe {

        if equation.is_null()
            || func.is_null()
            || var.is_null()
            || y1.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let eq_ref = &*equation;

        let y1_ref = &*y1;

        let func_str = match c_str_to_str(
        func,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        let var_str = match c_str_to_str(
        var,
    ) {
        | Some(s) => s,
        | None => {
            return std::ptr::null_mut()
        },
    };

        match ode::solve_by_reduction_of_order(
        eq_ref,
        func_str,
        var_str,
        y1_ref,
    ) {
        | Some(solution) => Box::into_raw(Box::new(solution)),
        | None => std::ptr::null_mut(),
    }
    }
}
