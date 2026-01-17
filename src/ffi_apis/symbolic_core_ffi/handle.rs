//! Handle-based FFI API for symbolic core operations.

use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;
use std::sync::Arc;

use crate::symbolic::core::Expr;

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


// --- Binary Operations ---
handle_binary_ffi!(rssn_core_create_add, new_add);
handle_binary_ffi!(rssn_core_create_sub, new_sub);
handle_binary_ffi!(rssn_core_create_mul, new_mul);
handle_binary_ffi!(rssn_core_create_div, new_div);
handle_binary_ffi!(rssn_core_create_pow, new_pow);
handle_binary_ffi!(rssn_core_create_log_base, new_log_base);
handle_binary_ffi!(rssn_core_create_atan2, new_atan2);

handle_binary_ffi!(rssn_core_create_complex, new_complex);
handle_binary_ffi!(rssn_core_create_matrix_mul, new_matrix_mul);
handle_binary_ffi!(rssn_core_create_matrix_vec_mul, new_matrix_vec_mul);
handle_binary_ffi!(rssn_core_create_xor, new_xor);
handle_binary_ffi!(rssn_core_create_implies, new_implies);
handle_binary_ffi!(rssn_core_create_equivalent, new_equivalent);
handle_binary_ffi!(rssn_core_create_beta, new_beta);
handle_binary_ffi!(rssn_core_create_bessel_j, new_bessel_j);
handle_binary_ffi!(rssn_core_create_bessel_y, new_bessel_y);
handle_binary_ffi!(rssn_core_create_legendre_p, new_legendre_p);
handle_binary_ffi!(rssn_core_create_laguerre_l, new_laguerre_l);
handle_binary_ffi!(rssn_core_create_hermite_h, new_hermite_h);
handle_binary_ffi!(rssn_core_create_kronecker_delta, new_kronecker_delta);
handle_binary_ffi!(rssn_core_create_apply, new_apply);

// --- Unary Operations ---
handle_unary_ffi!(rssn_core_create_sin, new_sin);
handle_unary_ffi!(rssn_core_create_cos, new_cos);
handle_unary_ffi!(rssn_core_create_tan, new_tan);
handle_unary_ffi!(rssn_core_create_exp, new_exp);
handle_unary_ffi!(rssn_core_create_log, new_log);
handle_unary_ffi!(rssn_core_create_neg, new_neg);
handle_unary_ffi!(rssn_core_create_abs, new_abs);
handle_unary_ffi!(rssn_core_create_sqrt, new_sqrt);
handle_unary_ffi!(rssn_core_create_asin, new_arcsin);
handle_unary_ffi!(rssn_core_create_acos, new_arccos);
handle_unary_ffi!(rssn_core_create_atan, new_arctan);
handle_unary_ffi!(rssn_core_create_sinh, new_sinh);
handle_unary_ffi!(rssn_core_create_cosh, new_cosh);
handle_unary_ffi!(rssn_core_create_tanh, new_tanh);
handle_unary_ffi!(rssn_core_create_transpose, new_transpose);
handle_unary_ffi!(rssn_core_create_inverse, new_inverse);
handle_unary_ffi!(rssn_core_create_sec, new_sec);
handle_unary_ffi!(rssn_core_create_csc, new_csc);
handle_unary_ffi!(rssn_core_create_cot, new_cot);
handle_unary_ffi!(rssn_core_create_asec, new_arcsec);
handle_unary_ffi!(rssn_core_create_acsc, new_arccsc);
handle_unary_ffi!(rssn_core_create_acot, new_arccot);
handle_unary_ffi!(rssn_core_create_sech, new_sech);
handle_unary_ffi!(rssn_core_create_csch, new_csch);
handle_unary_ffi!(rssn_core_create_coth, new_coth);
handle_unary_ffi!(rssn_core_create_asinh, new_arcsinh);
handle_unary_ffi!(rssn_core_create_acosh, new_arccosh);
handle_unary_ffi!(rssn_core_create_atanh, new_arctanh);
handle_unary_ffi!(rssn_core_create_asech, new_arcsech);
handle_unary_ffi!(rssn_core_create_acsch, new_arccsch);
handle_unary_ffi!(rssn_core_create_acoth, new_arccoth);
handle_unary_ffi!(rssn_core_create_not, new_not);
handle_unary_ffi!(rssn_core_create_floor, new_floor);
handle_unary_ffi!(rssn_core_create_gamma, new_gamma);
handle_unary_ffi!(rssn_core_create_erf, new_erf);
handle_unary_ffi!(rssn_core_create_erfc, new_erfc);
handle_unary_ffi!(rssn_core_create_erfi, new_erfi);
handle_unary_ffi!(rssn_core_create_zeta, new_zeta);
handle_unary_ffi!(rssn_core_create_digamma, new_digamma);

// --- N-ary Operations ---
handle_nary_ffi!(rssn_core_create_vector, new_vector);
handle_nary_ffi!(rssn_core_create_and, new_and);
handle_nary_ffi!(rssn_core_create_or, new_or);
handle_nary_ffi!(rssn_core_create_union, new_union);
handle_nary_ffi!(rssn_core_create_polynomial, new_polynomial);
handle_nary_ffi!(rssn_core_create_tuple, new_tuple);


/// Creates a new derivative expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
pub unsafe extern "C" fn rssn_core_create_derivative(
    expr: *const Expr,
    var: *const c_char
) -> *mut Expr {
    unsafe {
        if expr.is_null() { return std::ptr::null_mut(); }
        let var_str = match c_str_to_str(var) {
            Some(s) => s,
            None => return std::ptr::null_mut(),
        };
        let expr_ref = &*expr;
        Box::into_raw(Box::new(Expr::new_derivative(expr_ref, var_str.to_string())))
    }
}

/// Creates a new ForAll quantifier expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
pub unsafe extern "C" fn rssn_core_create_forall(
    var: *const c_char,
    expr: *const Expr
) -> *mut Expr {
    unsafe {
        if expr.is_null() { return std::ptr::null_mut(); }
        let var_str = match c_str_to_str(var) {
            Some(s) => s,
            None => return std::ptr::null_mut(),
        };
        let expr_ref = &*expr;
        Box::into_raw(Box::new(Expr::new_forall(var_str, expr_ref)))
    }
}

/// Creates a new Exists quantifier expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
pub unsafe extern "C" fn rssn_core_create_exists(
    var: *const c_char,
    expr: *const Expr
) -> *mut Expr {
    unsafe {
        if expr.is_null() { return std::ptr::null_mut(); }
        let var_str = match c_str_to_str(var) {
            Some(s) => s,
            None => return std::ptr::null_mut(),
        };
        let expr_ref = &*expr;
        Box::into_raw(Box::new(Expr::new_exists(var_str, expr_ref)))
    }
}

/// Creates a new higher-order derivative expression (derivativen).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub unsafe extern "C" fn rssn_core_create_derivativen(
    function: *const Expr,
    variable: *const c_char,
    grades: *const Expr
) -> *mut Expr {
    unsafe {
        if function.is_null() || grades.is_null() { return std::ptr::null_mut(); }
        let var_str = match c_str_to_str(variable) {
            Some(s) => s,
            None => return std::ptr::null_mut(),
        };
        Box::into_raw(Box::new(Expr::new_derivativen(
            &*function,
            var_str.to_string(),
            &*grades
        )))
    }
}

/// Creates a new Interval expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub unsafe extern "C" fn rssn_core_create_interval(
    lower: *const Expr,
    upper: *const Expr,
    incl_lower: bool,
    incl_upper: bool
) -> *mut Expr {
    unsafe {
        if lower.is_null() || upper.is_null() { return std::ptr::null_mut(); }
        Box::into_raw(Box::new(Expr::new_interval(
            &*lower,
            &*upper,
            incl_lower,
            incl_upper
        )))
    }
}

/// Creates a new Predicate expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub unsafe extern "C" fn rssn_core_create_predicate(
    name: *const c_char,
    args: *const *const Expr,
    len: usize
) -> *mut Expr {
    unsafe {
        let name_str = match c_str_to_str(name) {
            Some(s) => s,
            None => return std::ptr::null_mut(),
        };
        
        if args.is_null() { return std::ptr::null_mut(); }
        let mut exprs = Vec::with_capacity(len);
        for i in 0..len {
            let ptr = *args.add(i);
            if ptr.is_null() { return std::ptr::null_mut(); }
            exprs.push((&*ptr).clone());
        }
        
        Box::into_raw(Box::new(Expr::new_predicate(name_str, exprs)))
    }
}

/// Creates a new Matrix expression from a flat array of elements.
/// Elements should be in row-major order.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub unsafe extern "C" fn rssn_core_create_matrix(
    elements: *const *const Expr,
    rows: usize,
    cols: usize
) -> *mut Expr {
    unsafe {
        if elements.is_null() { return std::ptr::null_mut(); }
        
        let len = rows * cols;
        let mut grid = Vec::with_capacity(rows);
        
        for r in 0..rows {
            let mut row_vec = Vec::with_capacity(cols);
            for c in 0..cols {
                let idx = r * cols + c;
                let ptr = *elements.add(idx);
                if ptr.is_null() { return std::ptr::null_mut(); }
                row_vec.push((&*ptr).clone());
            }
            grid.push(row_vec);
        }
        
        Box::into_raw(Box::new(Expr::new_matrix(grid)))
    }
}

/// Creates a new Pi expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub extern "C" fn rssn_core_create_pi() -> *mut Expr {
    Box::into_raw(Box::new(Expr::new_pi()))
}

/// Creates a new E expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub extern "C" fn rssn_core_create_e() -> *mut Expr {
    Box::into_raw(Box::new(Expr::new_e()))
}

/// Creates a new Infinity expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub extern "C" fn rssn_core_create_infinity() -> *mut Expr {
    Box::into_raw(Box::new(Expr::new_infinity()))
}

/// Creates a new NegativeInfinity expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub extern "C" fn rssn_core_create_negative_infinity() -> *mut Expr {
    Box::into_raw(Box::new(Expr::new_negative_infinity()))
}

/// Creates a new variable expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_core_create_variable(
    name: *const c_char
) -> *mut Expr {

    unsafe {

        let name_str = match c_str_to_str(
            name,
        ) {
            | Some(s) => s,
            | None => {
                return std::ptr::null_mut()
            },
        };

        Box::into_raw(Box::new(
            Expr::new_variable(
                name_str,
            ),
        ))
    }
}

/// Creates a new symbolic constant (f64).
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.
pub extern "C" fn rssn_core_create_constant(
    val: f64
) -> *mut Expr {

    Box::into_raw(Box::new(
        Expr::new_constant(val),
    ))
}

/// Frees a symbolic expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_core_free_expr(
    expr: *mut Expr
) {

    unsafe {

        if !expr.is_null() {

            let _ = Box::from_raw(expr);
        }
    }
}

/// Returns the string representation of an expression.
/// The caller is responsible for freeing the returned string using `rssn_core_free_string`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_core_to_string(
    expr: *const Expr
) -> *mut c_char {

    unsafe {

        if expr.is_null() {

            return std::ptr::null_mut();
        }

        let expr_ref = &*expr;
        let s = format!("{expr_ref:?}");
        
        match CString::new(s) {
            Ok(c_str) => c_str.into_raw(),
            Err(_) => std::ptr::null_mut(),
        }
    }
}

/// Frees a C string returned by `rssn_core_to_string`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_core_free_string(
    s: *mut c_char
) {

    unsafe {

        if !s.is_null() {

            let _ = CString::from_raw(s);
        }
    }
}
