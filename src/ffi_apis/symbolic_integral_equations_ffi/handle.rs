//! Handle-based FFI API for symbolic integral equations.

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::symbolic::core::Expr;
use crate::symbolic::integral_equations::solve_airfoil_equation;
use crate::symbolic::integral_equations::FredholmEquation;
use crate::symbolic::integral_equations::VolterraEquation;

/// Creates a new Fredholm integral equation.
#[no_mangle]

pub extern "C" fn rssn_fredholm_new(
    y_x: *const Expr,
    f_x: *const Expr,
    lambda: *const Expr,
    kernel: *const Expr,
    lower_bound: *const Expr,
    upper_bound: *const Expr,
    var_x: *const c_char,
    var_t: *const c_char,
) -> *mut FredholmEquation {

    unsafe {

        if y_x.is_null()
            || f_x.is_null()
            || lambda.is_null()
            || kernel.is_null()
            || lower_bound.is_null()
            || upper_bound.is_null()
            || var_x.is_null()
            || var_t.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let var_x_str = match CStr::from_ptr(var_x).to_str() {
            | Ok(s) => s.to_string(),
            | Err(_) => return std::ptr::null_mut(),
        };

        let var_t_str = match CStr::from_ptr(var_t).to_str() {
            | Ok(s) => s.to_string(),
            | Err(_) => return std::ptr::null_mut(),
        };

        let eq = FredholmEquation::new(
            (*y_x).clone(),
            (*f_x).clone(),
            (*lambda).clone(),
            (*kernel).clone(),
            (*lower_bound).clone(),
            (*upper_bound).clone(),
            var_x_str,
            var_t_str,
        );

        Box::into_raw(Box::new(eq))
    }
}

/// Frees a Fredholm integral equation.
#[no_mangle]

pub extern "C" fn rssn_fredholm_free(
    ptr: *mut FredholmEquation
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Solves a Fredholm equation using the Neumann series method.
#[no_mangle]

pub extern "C" fn rssn_fredholm_solve_neumann(
    eq: *const FredholmEquation,
    iterations: usize,
) -> *mut Expr {

    unsafe {

        if eq.is_null() {

            return std::ptr::null_mut(
            );
        }

        let result = (*eq)
            .solve_neumann_series(
                iterations,
            );

        Box::into_raw(Box::new(result))
    }
}

/// Solves a Fredholm equation with a separable kernel.
#[no_mangle]

pub extern "C" fn rssn_fredholm_solve_separable(
    eq: *const FredholmEquation,
    a_funcs: *const *const Expr,
    a_len: usize,
    b_funcs: *const *const Expr,
    b_len: usize,
) -> *mut Expr {

    unsafe {

        if eq.is_null()
            || a_funcs.is_null()
            || b_funcs.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let a_slice =
            std::slice::from_raw_parts(
                a_funcs,
                a_len,
            );

        let b_slice =
            std::slice::from_raw_parts(
                b_funcs,
                b_len,
            );

        let mut a_vec =
            Vec::with_capacity(a_len);

        for &ptr in a_slice {

            if ptr.is_null() {

                return std::ptr::null_mut();
            }

            a_vec.push((*ptr).clone());
        }

        let mut b_vec =
            Vec::with_capacity(b_len);

        for &ptr in b_slice {

            if ptr.is_null() {

                return std::ptr::null_mut();
            }

            b_vec.push((*ptr).clone());
        }

        match (*eq)
            .solve_separable_kernel(
                a_vec, b_vec,
            ) {
            | Ok(result) => {
                Box::into_raw(Box::new(
                    result,
                ))
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Creates a new Volterra integral equation.
#[no_mangle]

pub extern "C" fn rssn_volterra_new(
    y_x: *const Expr,
    f_x: *const Expr,
    lambda: *const Expr,
    kernel: *const Expr,
    lower_bound: *const Expr,
    var_x: *const c_char,
    var_t: *const c_char,
) -> *mut VolterraEquation {

    unsafe {

        if y_x.is_null()
            || f_x.is_null()
            || lambda.is_null()
            || kernel.is_null()
            || lower_bound.is_null()
            || var_x.is_null()
            || var_t.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let var_x_str = match CStr::from_ptr(var_x).to_str() {
            | Ok(s) => s.to_string(),
            | Err(_) => return std::ptr::null_mut(),
        };

        let var_t_str = match CStr::from_ptr(var_t).to_str() {
            | Ok(s) => s.to_string(),
            | Err(_) => return std::ptr::null_mut(),
        };

        let eq = VolterraEquation::new(
            (*y_x).clone(),
            (*f_x).clone(),
            (*lambda).clone(),
            (*kernel).clone(),
            (*lower_bound).clone(),
            var_x_str,
            var_t_str,
        );

        Box::into_raw(Box::new(eq))
    }
}

/// Frees a Volterra integral equation.
#[no_mangle]

pub extern "C" fn rssn_volterra_free(
    ptr: *mut VolterraEquation
) {

    if !ptr.is_null() {

        unsafe {

            let _ = Box::from_raw(ptr);
        }
    }
}

/// Solves a Volterra equation using successive approximations.
#[no_mangle]

pub extern "C" fn rssn_volterra_solve_successive(
    eq: *const VolterraEquation,
    iterations: usize,
) -> *mut Expr {

    unsafe {

        if eq.is_null() {

            return std::ptr::null_mut(
            );
        }

        let result = (*eq).solve_successive_approximations(iterations);

        Box::into_raw(Box::new(result))
    }
}

/// Solves a Volterra equation by differentiation.
#[no_mangle]

pub extern "C" fn rssn_volterra_solve_by_differentiation(
    eq: *const VolterraEquation
) -> *mut Expr {

    unsafe {

        if eq.is_null() {

            return std::ptr::null_mut(
            );
        }

        match (*eq)
            .solve_by_differentiation()
        {
            | Ok(result) => {
                Box::into_raw(Box::new(
                    result,
                ))
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

/// Solves the airfoil singular integral equation.
#[no_mangle]

pub extern "C" fn rssn_solve_airfoil_equation(
    f_x: *const Expr,
    var_x: *const c_char,
    var_t: *const c_char,
) -> *mut Expr {

    unsafe {

        if f_x.is_null()
            || var_x.is_null()
            || var_t.is_null()
        {

            return std::ptr::null_mut(
            );
        }

        let var_x_str = match CStr::from_ptr(var_x).to_str() {
            | Ok(s) => s,
            | Err(_) => return std::ptr::null_mut(),
        };

        let var_t_str = match CStr::from_ptr(var_t).to_str() {
            | Ok(s) => s,
            | Err(_) => return std::ptr::null_mut(),
        };

        let result =
            solve_airfoil_equation(
                &(*f_x),
                var_x_str,
                var_t_str,
            );

        Box::into_raw(Box::new(result))
    }
}
