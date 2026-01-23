//! Bincode-based FFI API for symbolic core operations.

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::core::SparsePolynomial;

/// Creates a variable expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_variable(
    name: *const c_char
) -> BincodeBuffer {

    let name_str = unsafe {

        if name.is_null() {

            return BincodeBuffer::empty();
        }

        match CStr::from_ptr(name).to_str() {
            Ok(s) => s,
            Err(_) => return BincodeBuffer::empty(),
        }
    };

    let expr =
        Expr::new_variable(name_str);

    to_bincode_buffer(&expr)
}

/// Creates a constant expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_constant(
    val: f64
) -> BincodeBuffer {

    let expr = Expr::new_constant(val);

    to_bincode_buffer(&expr)
}

/// Creates a new derivative expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_derivative(
    expr_buf: BincodeBuffer,
    var: *const c_char,
) -> BincodeBuffer {

    let expr_opt: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(expr), Some(v)) =
        (expr_opt, var_str)
    {

        let res = Expr::new_derivative(
            &expr,
            v.to_string(),
        );

        to_bincode_buffer(&res)
    } else {

        BincodeBuffer::empty()
    }
}

/// Creates a new ForAll quantifier expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_forall(
    var: *const c_char,
    expr_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr_opt: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(expr), Some(v)) =
        (expr_opt, var_str)
    {

        let res =
            Expr::new_forall(v, &expr);

        to_bincode_buffer(&res)
    } else {

        BincodeBuffer::empty()
    }
}

/// Creates a new Exists quantifier expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_exists(
    var: *const c_char,
    expr_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr_opt: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var_str = unsafe {

        if var.is_null() {

            None
        } else {

            CStr::from_ptr(var)
                .to_str()
                .ok()
        }
    };

    if let (Some(expr), Some(v)) =
        (expr_opt, var_str)
    {

        let res =
            Expr::new_exists(v, &expr);

        to_bincode_buffer(&res)
    } else {

        BincodeBuffer::empty()
    }
}


/// Creates a new SparsePolynomial expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_sparse_polynomial(
    poly_buf: BincodeBuffer
) -> BincodeBuffer {

    let poly: Option<SparsePolynomial> =
        from_bincode_buffer(&poly_buf);

    if let Some(p) = poly {

        let expr =
            Expr::new_sparse_polynomial(
                p,
            );

        to_bincode_buffer(&expr)
    } else {

        BincodeBuffer::empty()
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

pub extern "C" fn rssn_core_bincode_create_derivativen(
    function_buf: BincodeBuffer,
    variable: *const c_char,
    grades_buf: BincodeBuffer,
) -> BincodeBuffer {

    let function: Option<Expr> =
        from_bincode_buffer(
            &function_buf,
        );

    let grades: Option<Expr> =
        from_bincode_buffer(
            &grades_buf,
        );

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

        to_bincode_buffer(&res)
    } else {

        BincodeBuffer::empty()
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

pub extern "C" fn rssn_core_bincode_create_interval(
    lower_buf: BincodeBuffer,
    upper_buf: BincodeBuffer,
    incl_lower: bool,
    incl_upper: bool,
) -> BincodeBuffer {

    let lower: Option<Expr> =
        from_bincode_buffer(&lower_buf);

    let upper: Option<Expr> =
        from_bincode_buffer(&upper_buf);

    if let (Some(l), Some(u)) =
        (lower, upper)
    {

        let res = Expr::new_interval(
            &l,
            &u,
            incl_lower,
            incl_upper,
        );

        to_bincode_buffer(&res)
    } else {

        BincodeBuffer::empty()
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

pub extern "C" fn rssn_core_bincode_create_predicate(
    name: *const c_char,
    args_buf: BincodeBuffer,
) -> BincodeBuffer {

    // Expects serialized Vec<Expr> for args
    let args: Option<Vec<Expr>> =
        from_bincode_buffer(&args_buf);

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

        to_bincode_buffer(&res)
    } else {

        BincodeBuffer::empty()
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

pub extern "C" fn rssn_core_bincode_create_matrix(
    elements_buf: BincodeBuffer /* Serialized Vec<Vec<Expr>> */
) -> BincodeBuffer {

    let elements: Option<
        Vec<Vec<Expr>>,
    > = from_bincode_buffer(
        &elements_buf,
    );

    if let Some(e) = elements {

        let res =
            std::panic::catch_unwind(
                || Expr::new_matrix(e),
            );

        match res {
            | Ok(expr) => {
                to_bincode_buffer(&expr)
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Creates a Pi expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_pi()
-> BincodeBuffer {

    let expr = Expr::new_pi();

    to_bincode_buffer(&expr)
}

/// Creates an E expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_e()
-> BincodeBuffer {

    let expr = Expr::new_e();

    to_bincode_buffer(&expr)
}

/// Creates an Infinity expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_infinity()
-> BincodeBuffer {

    let expr = Expr::new_infinity();

    to_bincode_buffer(&expr)
}

/// Creates a NegativeInfinity expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub extern "C" fn rssn_core_bincode_create_negative_infinity()
-> BincodeBuffer {

    let expr =
        Expr::new_negative_infinity();

    to_bincode_buffer(&expr)
}


// --- Binary Operations ---
bincode_binary_ffi!(
    rssn_core_bincode_add,
    new_add
);

bincode_binary_ffi!(
    rssn_core_bincode_sub,
    new_sub
);

bincode_binary_ffi!(
    rssn_core_bincode_mul,
    new_mul
);

bincode_binary_ffi!(
    rssn_core_bincode_div,
    new_div
);

bincode_binary_ffi!(
    rssn_core_bincode_pow,
    new_pow
);

bincode_binary_ffi!(
    rssn_core_bincode_log_base,
    new_log_base
);

bincode_binary_ffi!(
    rssn_core_bincode_atan2,
    new_atan2
);

bincode_binary_ffi!(
    rssn_core_bincode_complex,
    new_complex
);

bincode_binary_ffi!(
    rssn_core_bincode_matrix_mul,
    new_matrix_mul
);

bincode_binary_ffi!(
    rssn_core_bincode_matrix_vec_mul,
    new_matrix_vec_mul
);

bincode_binary_ffi!(
    rssn_core_bincode_xor,
    new_xor
);

bincode_binary_ffi!(
    rssn_core_bincode_implies,
    new_implies
);

bincode_binary_ffi!(
    rssn_core_bincode_equivalent,
    new_equivalent
);

bincode_binary_ffi!(
    rssn_core_bincode_beta,
    new_beta
);

bincode_binary_ffi!(
    rssn_core_bincode_bessel_j,
    new_bessel_j
);

bincode_binary_ffi!(
    rssn_core_bincode_bessel_y,
    new_bessel_y
);

bincode_binary_ffi!(
    rssn_core_bincode_legendre_p,
    new_legendre_p
);

bincode_binary_ffi!(
    rssn_core_bincode_laguerre_l,
    new_laguerre_l
);

bincode_binary_ffi!(
    rssn_core_bincode_hermite_h,
    new_hermite_h
);

bincode_binary_ffi!(
    rssn_core_bincode_kronecker_delta,
    new_kronecker_delta
);

bincode_binary_ffi!(
    rssn_core_bincode_apply,
    new_apply
);

// --- Unary Operations ---
bincode_unary_ffi!(
    rssn_core_bincode_sin,
    new_sin
);

bincode_unary_ffi!(
    rssn_core_bincode_cos,
    new_cos
);

bincode_unary_ffi!(
    rssn_core_bincode_tan,
    new_tan
);

bincode_unary_ffi!(
    rssn_core_bincode_exp,
    new_exp
);

bincode_unary_ffi!(
    rssn_core_bincode_log,
    new_log
);

bincode_unary_ffi!(
    rssn_core_bincode_neg,
    new_neg
);

bincode_unary_ffi!(
    rssn_core_bincode_abs,
    new_abs
);

bincode_unary_ffi!(
    rssn_core_bincode_sqrt,
    new_sqrt
);

bincode_unary_ffi!(
    rssn_core_bincode_asin,
    new_arcsin
);

bincode_unary_ffi!(
    rssn_core_bincode_acos,
    new_arccos
);

bincode_unary_ffi!(
    rssn_core_bincode_atan,
    new_arctan
);

bincode_unary_ffi!(
    rssn_core_bincode_sinh,
    new_sinh
);

bincode_unary_ffi!(
    rssn_core_bincode_cosh,
    new_cosh
);

bincode_unary_ffi!(
    rssn_core_bincode_tanh,
    new_tanh
);

bincode_unary_ffi!(
    rssn_core_bincode_transpose,
    new_transpose
);

bincode_unary_ffi!(
    rssn_core_bincode_inverse,
    new_inverse
);

bincode_unary_ffi!(
    rssn_core_bincode_sec,
    new_sec
);

bincode_unary_ffi!(
    rssn_core_bincode_csc,
    new_csc
);

bincode_unary_ffi!(
    rssn_core_bincode_cot,
    new_cot
);

bincode_unary_ffi!(
    rssn_core_bincode_asec,
    new_arcsec
);

bincode_unary_ffi!(
    rssn_core_bincode_acsc,
    new_arccsc
);

bincode_unary_ffi!(
    rssn_core_bincode_acot,
    new_arccot
);

bincode_unary_ffi!(
    rssn_core_bincode_sech,
    new_sech
);

bincode_unary_ffi!(
    rssn_core_bincode_csch,
    new_csch
);

bincode_unary_ffi!(
    rssn_core_bincode_coth,
    new_coth
);

bincode_unary_ffi!(
    rssn_core_bincode_asinh,
    new_arcsinh
);

bincode_unary_ffi!(
    rssn_core_bincode_acosh,
    new_arccosh
);

bincode_unary_ffi!(
    rssn_core_bincode_atanh,
    new_arctanh
);

bincode_unary_ffi!(
    rssn_core_bincode_asech,
    new_arcsech
);

bincode_unary_ffi!(
    rssn_core_bincode_acsch,
    new_arccsch
);

bincode_unary_ffi!(
    rssn_core_bincode_acoth,
    new_arccoth
);

bincode_unary_ffi!(
    rssn_core_bincode_not,
    new_not
);

bincode_unary_ffi!(
    rssn_core_bincode_floor,
    new_floor
);

bincode_unary_ffi!(
    rssn_core_bincode_gamma,
    new_gamma
);

bincode_unary_ffi!(
    rssn_core_bincode_erf,
    new_erf
);

bincode_unary_ffi!(
    rssn_core_bincode_erfc,
    new_erfc
);

bincode_unary_ffi!(
    rssn_core_bincode_erfi,
    new_erfi
);

bincode_unary_ffi!(
    rssn_core_bincode_zeta,
    new_zeta
);

bincode_unary_ffi!(
    rssn_core_bincode_digamma,
    new_digamma
);

// --- N-ary Operations ---
bincode_nary_ffi!(
    rssn_core_bincode_vector,
    new_vector
);

bincode_nary_ffi!(
    rssn_core_bincode_and,
    new_and
);

bincode_nary_ffi!(
    rssn_core_bincode_or,
    new_or
);

bincode_nary_ffi!(
    rssn_core_bincode_union,
    new_union
);

bincode_nary_ffi!(
    rssn_core_bincode_polynomial,
    new_polynomial
);

bincode_nary_ffi!(
    rssn_core_bincode_tuple,
    new_tuple
);
