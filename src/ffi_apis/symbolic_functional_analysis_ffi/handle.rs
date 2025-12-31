use std::ffi::CStr;
use std::os::raw::c_char;

use crate::symbolic::core::Expr;
use crate::symbolic::functional_analysis::{HilbertSpace, BanachSpace, LinearOperator, inner_product, norm, banach_norm, are_orthogonal, project, gram_schmidt};

// --- HilbertSpace ---

/// Creates a new Hilbert space over a specified interval.
///
/// # Arguments
/// * `var` - The name of the independent variable defining the space.
/// * `lower_bound` - Symbolic expression for the lower bound of the interval.
/// * `upper_bound` - Symbolic expression for the upper bound of the interval.
///
/// # Returns
/// A raw pointer to the newly created `HilbertSpace`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_hilbert_space_create(
    var: *const c_char,
    lower_bound: *const Expr,
    upper_bound: *const Expr,
) -> *mut HilbertSpace {

    unsafe {

        let var_str =
            CStr::from_ptr(var)
                .to_str()
                .unwrap();

        let space = HilbertSpace::new(
            var_str,
            (*lower_bound).clone(),
            (*upper_bound).clone(),
        );

        Box::into_raw(Box::new(space))
    }
}

/// Frees a Hilbert space handle.
///
/// # Arguments
/// * `ptr` - Pointer to the `HilbertSpace` to free.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_hilbert_space_free(
    ptr: *mut HilbertSpace
) {

    unsafe {

        if !ptr.is_null() {

            let _ = Box::from_raw(ptr);
        }
    }
}

// --- BanachSpace ---

/// Creates a new Banach space $L^p$ over a specified interval.
///
/// # Arguments
/// * `var` - The name of the independent variable.
/// * `lower_bound` - Symbolic expression for the lower bound.
/// * `upper_bound` - Symbolic expression for the upper bound.
/// * `p` - Symbolic expression representing the $p$ parameter of the $L^p$ norm.
///
/// # Returns
/// A raw pointer to the newly created `BanachSpace`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_banach_space_create(
    var: *const c_char,
    lower_bound: *const Expr,
    upper_bound: *const Expr,
    p: *const Expr,
) -> *mut BanachSpace {

    unsafe {

        let var_str =
            CStr::from_ptr(var)
                .to_str()
                .unwrap();

        let space = BanachSpace::new(
            var_str,
            (*lower_bound).clone(),
            (*upper_bound).clone(),
            (*p).clone(),
        );

        Box::into_raw(Box::new(space))
    }
}

/// Frees a Banach space handle.
///
/// # Arguments
/// * `ptr` - Pointer to the `BanachSpace` to free.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_banach_space_free(
    ptr: *mut BanachSpace
) {

    unsafe {

        if !ptr.is_null() {

            let _ = Box::from_raw(ptr);
        }
    }
}

// --- LinearOperator ---

/// Creates a linear derivative operator handle.
///
/// # Arguments
/// * `var` - The name of the variable with respect to which differentiation is performed.
///
/// # Returns
/// A raw pointer to the `LinearOperator`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_linear_operator_derivative_create(
    var: *const c_char
) -> *mut LinearOperator {

    unsafe {

        let var_str =
            CStr::from_ptr(var)
                .to_str()
                .unwrap();

        let op =
            LinearOperator::Derivative(
                var_str.to_string(),
            );

        Box::into_raw(Box::new(op))
    }
}

/// Creates a linear integral operator handle.
///
/// # Arguments
/// * `lower_bound` - Symbolic expression for the lower integration limit.
/// * `var` - The integration variable.
///
/// # Returns
/// A raw pointer to the `LinearOperator`.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

/// # Panics
///
/// This function may panic if the FFI input is malformed, null where not expected,
/// or if internal state synchronization fails (e.g., poisoned locks).

pub unsafe extern "C" fn rssn_linear_operator_integral_create(
    lower_bound: *const Expr,
    var: *const c_char,
) -> *mut LinearOperator {

    unsafe {

        let var_str =
            CStr::from_ptr(var)
                .to_str()
                .unwrap();

        let op =
            LinearOperator::Integral(
                (*lower_bound).clone(),
                var_str.to_string(),
            );

        Box::into_raw(Box::new(op))
    }
}

/// Applies a linear operator to a symbolic expression.
///
/// # Arguments
/// * `op` - Handle to the linear operator.
/// * `expr` - Handle to the expression to operate on.
///
/// # Returns
/// A raw pointer to the resulting symbolic expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_linear_operator_apply(
    op: *const LinearOperator,
    expr: *const Expr,
) -> *mut Expr {

    unsafe {

        let result =
            (*op).apply(&*expr);

        Box::into_raw(Box::new(result))
    }
}

/// Frees a linear operator handle.
///
/// # Arguments
/// * `ptr` - Pointer to the `LinearOperator` to free.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_linear_operator_free(
    ptr: *mut LinearOperator
) {

    unsafe {

        if !ptr.is_null() {

            let _ = Box::from_raw(ptr);
        }
    }
}

// --- Functions ---

/// Computes the inner product of two functions in a Hilbert space.
///
/// # Arguments
/// * `space` - Handle to the Hilbert space.
/// * `f` - Handle to the first expression.
/// * `g` - Handle to the second expression.
///
/// # Returns
/// A raw pointer to the symbolic expression representing the inner product.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_inner_product(
    space: *const HilbertSpace,
    f: *const Expr,
    g: *const Expr,
) -> *mut Expr {

    unsafe {

        let result = inner_product(
            &*space,
            &*f,
            &*g,
        );

        Box::into_raw(Box::new(result))
    }
}

/// Computes the $L^2$ norm of a function in a Hilbert space.
///
/// # Arguments
/// * `space` - Handle to the Hilbert space.
/// * `f` - Handle to the expression.
///
/// # Returns
/// A raw pointer to the symbolic expression representing the norm.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_norm(
    space: *const HilbertSpace,
    f: *const Expr,
) -> *mut Expr {

    unsafe {

        let result = norm(&*space, &*f);

        Box::into_raw(Box::new(result))
    }
}

/// Computes the $L^p$ norm of a function in a Banach space.
///
/// # Arguments
/// * `space` - Handle to the Banach space.
/// * `f` - Handle to the expression.
///
/// # Returns
/// A raw pointer to the symbolic expression representing the norm.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_banach_norm(
    space: *const BanachSpace,
    f: *const Expr,
) -> *mut Expr {

    unsafe {

        let result =
            banach_norm(&*space, &*f);

        Box::into_raw(Box::new(result))
    }
}

/// Checks if two functions are orthogonal in a Hilbert space.
///
/// # Arguments
/// * `space` - Handle to the Hilbert space.
/// * `f` - Handle to the first expression.
/// * `g` - Handle to the second expression.
///
/// # Returns
/// `true` if orthogonal, `false` otherwise.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_are_orthogonal(
    space: *const HilbertSpace,
    f: *const Expr,
    g: *const Expr,
) -> bool {

    unsafe {

        are_orthogonal(
            &*space,
            &*f,
            &*g,
        )
    }
}

/// Projects one function onto another within a Hilbert space.
///
/// # Arguments
/// * `space` - Handle to the Hilbert space.
/// * `f` - Handle to the function to project.
/// * `g` - Handle to the function onto which $f$ is projected.
///
/// # Returns
/// A raw pointer to the resulting projection expression.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_project(
    space: *const HilbertSpace,
    f: *const Expr,
    g: *const Expr,
) -> *mut Expr {

    unsafe {

        let result =
            project(&*space, &*f, &*g);

        Box::into_raw(Box::new(result))
    }
}

/// Performs the Gram-Schmidt process to produce an orthogonal basis.
///
/// # Arguments
/// * `space` - Handle to the Hilbert space.
/// * `basis_ptr` - Array of handles to the input basis expressions.
/// * `basis_len` - Number of input expressions.
/// * `out_len` - Pointer to store the number of orthogonal expressions produced.
///
/// # Returns
/// A raw pointer to an array of handles to the orthogonalized basis expressions.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_gram_schmidt(
    space: *const HilbertSpace,
    basis_ptr: *const *const Expr,
    basis_len: usize,
    out_len: *mut usize,
) -> *mut *mut Expr {

    unsafe {

        let basis_slice =
            std::slice::from_raw_parts(
                basis_ptr,
                basis_len,
            );

        let basis: Vec<Expr> =
            basis_slice
                .iter()
                .map(|&p| (*p).clone())
                .collect();

        let orthogonal_basis =
            gram_schmidt(
                &*space,
                &basis,
            );

        *out_len =
            orthogonal_basis.len();

        let mut out_ptrs =
            Vec::with_capacity(
                orthogonal_basis.len(),
            );

        for expr in orthogonal_basis {

            out_ptrs.push(
                Box::into_raw(
                    Box::new(expr),
                ),
            );
        }

        let ptr = out_ptrs.as_mut_ptr();

        std::mem::forget(out_ptrs);

        ptr
    }
}
