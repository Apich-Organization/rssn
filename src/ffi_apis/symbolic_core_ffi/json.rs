//! JSON-based FFI API for symbolic core operations.

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::core::Expr;
use crate::symbolic::core::SparsePolynomial;

/// Creates a variable expression using JSON.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_json_create_variable(
    name: *const c_char
) -> *mut c_char {

    let name_str = unsafe {

        if name.is_null() {

            return std::ptr::null_mut(
            );
        }

        match CStr::from_ptr(name).to_str() {
            Ok(s) => s,
            Err(_) => return std::ptr::null_mut(),
        }
    };

    let expr =
        Expr::new_variable(name_str);

    to_json_string(&expr)
}

/// Creates a new SparsePolynomial expression using JSON.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_json_create_sparse_polynomial(
    poly_json: *const c_char
) -> *mut c_char {

    let poly: Option<SparsePolynomial> =
        from_json_string(poly_json);

    if let Some(p) = poly {

        let expr =
            Expr::new_sparse_polynomial(
                p,
            );

        to_json_string(&expr)
    } else {

        std::ptr::null_mut()
    }
}

/// Creates a constant expression using JSON.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_json_create_constant(
    val: f64
) -> *mut c_char {

    let expr = Expr::new_constant(val);

    to_json_string(&expr)
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

pub extern "C" fn rssn_core_json_create_derivativen(
    function_json: *const c_char,
    variable: *const c_char,
    grades_json: *const c_char,
) -> *mut c_char {

    let function: Option<Expr> =
        from_json_string(function_json);

    let grades: Option<Expr> =
        from_json_string(grades_json);

    let var_str = unsafe {

        if variable.is_null() {

            None
        } else {

            CStr::from_ptr(variable)
                .to_str()
                .ok()
        }
    };

    if let (Some(f), Some(g), Some(v)) = (
        function,
        grades,
        var_str,
    ) {

        let res = Expr::new_derivativen(
            &f,
            v.to_string(),
            &g,
        );

        to_json_string(&res)
    } else {

        std::ptr::null_mut()
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

pub extern "C" fn rssn_core_json_create_interval(
    lower_json: *const c_char,
    upper_json: *const c_char,
    incl_lower: bool,
    incl_upper: bool,
) -> *mut c_char {

    let lower: Option<Expr> =
        from_json_string(lower_json);

    let upper: Option<Expr> =
        from_json_string(upper_json);

    if let (Some(l), Some(u)) =
        (lower, upper)
    {

        let res = Expr::new_interval(
            &l,
            &u,
            incl_lower,
            incl_upper,
        );

        to_json_string(&res)
    } else {

        std::ptr::null_mut()
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

pub extern "C" fn rssn_core_json_create_predicate(
    name: *const c_char,
    args_json: *const c_char,
) -> *mut c_char {

    let args: Option<Vec<Expr>> =
        from_json_string(args_json);

    let name_str = unsafe {

        if name.is_null() {

            None
        } else {

            CStr::from_ptr(name)
                .to_str()
                .ok()
        }
    };

    if let (Some(a), Some(n)) =
        (args, name_str)
    {

        let res =
            Expr::new_predicate(n, a);

        to_json_string(&res)
    } else {

        std::ptr::null_mut()
    }
}

/// Creates a new Matrix expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_json_create_matrix(
    elements_json: *const c_char
) -> *mut c_char {

    let elements: Option<
        Vec<Vec<Expr>>,
    > = from_json_string(elements_json);

    if let Some(e) = elements {

        let res =
            std::panic::catch_unwind(
                || Expr::new_matrix(e),
            );

        match res {
            | Ok(expr) => {
                to_json_string(&expr)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}

/// Creates a Pi expression using JSON.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_json_create_pi()
-> *mut c_char {

    let expr = Expr::new_pi();

    to_json_string(&expr)
}

/// Creates an E expression using JSON.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_json_create_e()
-> *mut c_char {

    let expr = Expr::new_e();

    to_json_string(&expr)
}

/// Creates an Infinity expression using JSON.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_json_create_infinity()
-> *mut c_char {

    let expr = Expr::new_infinity();

    to_json_string(&expr)
}

/// Creates a NegativeInfinity expression using JSON.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_json_create_negative_infinity()
-> *mut c_char {

    let expr =
        Expr::new_negative_infinity();

    to_json_string(&expr)
}


// --- Binary Operations ---
json_binary_ffi!(
    rssn_core_json_add,
    new_add
);

json_binary_ffi!(
    rssn_core_json_sub,
    new_sub
);

json_binary_ffi!(
    rssn_core_json_mul,
    new_mul
);

json_binary_ffi!(
    rssn_core_json_div,
    new_div
);

json_binary_ffi!(
    rssn_core_json_pow,
    new_pow
);

json_binary_ffi!(
    rssn_core_json_log_base,
    new_log_base
);

json_binary_ffi!(
    rssn_core_json_atan2,
    new_atan2
);

json_binary_ffi!(
    rssn_core_json_complex,
    new_complex
);

json_binary_ffi!(
    rssn_core_json_matrix_mul,
    new_matrix_mul
);

json_binary_ffi!(
    rssn_core_json_matrix_vec_mul,
    new_matrix_vec_mul
);

json_binary_ffi!(
    rssn_core_json_xor,
    new_xor
);

json_binary_ffi!(
    rssn_core_json_implies,
    new_implies
);

json_binary_ffi!(
    rssn_core_json_equivalent,
    new_equivalent
);

json_binary_ffi!(
    rssn_core_json_beta,
    new_beta
);

json_binary_ffi!(
    rssn_core_json_bessel_j,
    new_bessel_j
);

json_binary_ffi!(
    rssn_core_json_bessel_y,
    new_bessel_y
);

json_binary_ffi!(
    rssn_core_json_legendre_p,
    new_legendre_p
);

json_binary_ffi!(
    rssn_core_json_laguerre_l,
    new_laguerre_l
);

json_binary_ffi!(
    rssn_core_json_hermite_h,
    new_hermite_h
);

json_binary_ffi!(
    rssn_core_json_kronecker_delta,
    new_kronecker_delta
);

json_binary_ffi!(
    rssn_core_json_apply,
    new_apply
);

// --- Unary Operations ---
json_unary_ffi!(
    rssn_core_json_sin,
    new_sin
);

json_unary_ffi!(
    rssn_core_json_cos,
    new_cos
);

json_unary_ffi!(
    rssn_core_json_tan,
    new_tan
);

json_unary_ffi!(
    rssn_core_json_exp,
    new_exp
);

json_unary_ffi!(
    rssn_core_json_log,
    new_log
);

json_unary_ffi!(
    rssn_core_json_neg,
    new_neg
);

json_unary_ffi!(
    rssn_core_json_abs,
    new_abs
);

json_unary_ffi!(
    rssn_core_json_sqrt,
    new_sqrt
);

json_unary_ffi!(
    rssn_core_json_asin,
    new_arcsin
);

json_unary_ffi!(
    rssn_core_json_acos,
    new_arccos
);

json_unary_ffi!(
    rssn_core_json_atan,
    new_arctan
);

json_unary_ffi!(
    rssn_core_json_sinh,
    new_sinh
);

json_unary_ffi!(
    rssn_core_json_cosh,
    new_cosh
);

json_unary_ffi!(
    rssn_core_json_tanh,
    new_tanh
);

json_unary_ffi!(
    rssn_core_json_transpose,
    new_transpose
);

json_unary_ffi!(
    rssn_core_json_inverse,
    new_inverse
);

json_unary_ffi!(
    rssn_core_json_sec,
    new_sec
);

json_unary_ffi!(
    rssn_core_json_csc,
    new_csc
);

json_unary_ffi!(
    rssn_core_json_cot,
    new_cot
);

json_unary_ffi!(
    rssn_core_json_asec,
    new_arcsec
);

json_unary_ffi!(
    rssn_core_json_acsc,
    new_arccsc
);

json_unary_ffi!(
    rssn_core_json_acot,
    new_arccot
);

json_unary_ffi!(
    rssn_core_json_sech,
    new_sech
);

json_unary_ffi!(
    rssn_core_json_csch,
    new_csch
);

json_unary_ffi!(
    rssn_core_json_coth,
    new_coth
);

json_unary_ffi!(
    rssn_core_json_asinh,
    new_arcsinh
);

json_unary_ffi!(
    rssn_core_json_acosh,
    new_arccosh
);

json_unary_ffi!(
    rssn_core_json_atanh,
    new_arctanh
);

json_unary_ffi!(
    rssn_core_json_asech,
    new_arcsech
);

json_unary_ffi!(
    rssn_core_json_acsch,
    new_arccsch
);

json_unary_ffi!(
    rssn_core_json_acoth,
    new_arccoth
);

json_unary_ffi!(
    rssn_core_json_not,
    new_not
);

json_unary_ffi!(
    rssn_core_json_floor,
    new_floor
);

json_unary_ffi!(
    rssn_core_json_gamma,
    new_gamma
);

json_unary_ffi!(
    rssn_core_json_erf,
    new_erf
);

json_unary_ffi!(
    rssn_core_json_erfc,
    new_erfc
);

json_unary_ffi!(
    rssn_core_json_erfi,
    new_erfi
);

json_unary_ffi!(
    rssn_core_json_zeta,
    new_zeta
);

json_unary_ffi!(
    rssn_core_json_digamma,
    new_digamma
);

// --- N-ary Operations ---
json_nary_ffi!(
    rssn_core_json_vector,
    new_vector
);

json_nary_ffi!(
    rssn_core_json_and,
    new_and
);

json_nary_ffi!(
    rssn_core_json_or,
    new_or
);

json_nary_ffi!(
    rssn_core_json_union,
    new_union
);

json_nary_ffi!(
    rssn_core_json_polynomial,
    new_polynomial
);

json_nary_ffi!(
    rssn_core_json_tuple,
    new_tuple
);
