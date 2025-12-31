use std::os::raw::c_char;

use crate::ffi_apis::common::{to_json_string, from_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::lie_groups_and_algebras::{so3, su2, lie_bracket, exponential_map, adjoint_representation_group, adjoint_representation_algebra, LieAlgebra, commutator_table, check_jacobi_identity, so3_generators, su2_generators};

// --- LieAlgebra Creation ---

#[unsafe(no_mangle)]

/// Constructs the Lie algebra \(\mathfrak{so}(3)\) and returns it as JSON.
///
/// The Lie algebra \(\mathfrak{so}(3)\) consists of skew-symmetric \(3\times3\) matrices
/// associated with the rotation group SO(3).
///
/// # Arguments
///
/// This function takes no arguments.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded [`LieAlgebra`] for \(\mathfrak{so}(3)\),
/// or null if serialization fails.
///
/// # Safety
///
/// This function is unsafe because it returns ownership of a heap-allocated C
/// string that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_lie_algebra_so3()
-> *mut c_char {

    let algebra = so3();

    to_json_string(&algebra)
}

#[unsafe(no_mangle)]

/// Constructs the Lie algebra \(\mathfrak{su}(2)\) and returns it as JSON.
///
/// The Lie algebra \(\mathfrak{su}(2)\) consists of traceless skew-Hermitian
/// \(2\times2\) matrices associated with the special unitary group SU(2).
///
/// # Arguments
///
/// This function takes no arguments.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded [`LieAlgebra`] for \(\mathfrak{su}(2)\),
/// or null if serialization fails.
///
/// # Safety
///
/// This function is unsafe because it returns ownership of a heap-allocated C
/// string that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_lie_algebra_su2()
-> *mut c_char {

    let algebra = su2();

    to_json_string(&algebra)
}

// --- Lie Bracket ---

#[unsafe(no_mangle)]

/// Computes the Lie bracket of two elements of a Lie algebra using JSON serialization.
///
/// # Arguments
///
/// * `x_json` - C string pointer with JSON-encoded `Expr` representing \(x\).
/// * `y_json` - C string pointer with JSON-encoded `Expr` representing \(y\).
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the bracket \([x, y]\), or
/// null if deserialization fails or the computation encounters an error.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_lie_bracket(
    x_json: *const c_char,
    y_json: *const c_char,
) -> *mut c_char {

    let x: Expr = match from_json_string(
        x_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    let y: Expr = match from_json_string(
        y_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    match lie_bracket(&x, &y) {
        | Ok(result) => {
            to_json_string(&result)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

// --- Exponential Map ---

#[unsafe(no_mangle)]

/// Applies the exponential map to a Lie algebra element using JSON serialization.
///
/// The exponential map is approximated by a truncated series of the given order.
///
/// # Arguments
///
/// * `x_json` - C string pointer with JSON-encoded `Expr` representing the Lie algebra element.
/// * `order` - Truncation order for the series expansion of the exponential map.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the corresponding group
/// element, or null if deserialization fails or the computation encounters an error.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_exponential_map(
    x_json: *const c_char,
    order: usize,
) -> *mut c_char {

    let x: Expr = match from_json_string(
        x_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    match exponential_map(&x, order) {
        | Ok(result) => {
            to_json_string(&result)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

// --- Adjoint Representations ---

#[unsafe(no_mangle)]

/// Applies the adjoint representation of a Lie group element to a Lie algebra element using JSON.
///
/// This computes \(\mathrm{Ad}_g(x)\), describing how the group element \(g\)
/// conjugates the Lie algebra element \(x\).
///
/// # Arguments
///
/// * `g_json` - C string pointer with JSON-encoded `Expr` for the group element \(g\).
/// * `x_json` - C string pointer with JSON-encoded `Expr` for the Lie algebra element \(x\).
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the result of the group
/// adjoint action, or null if deserialization fails or the computation encounters
/// an error.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_adjoint_representation_group(
    g_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let g: Expr = match from_json_string(
        g_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    let x: Expr = match from_json_string(
        x_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    match adjoint_representation_group(
        &g, &x,
    ) {
        | Ok(result) => {
            to_json_string(&result)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

#[unsafe(no_mangle)]

/// Applies the adjoint representation of a Lie algebra element to another element using JSON.
///
/// This computes \(\mathrm{ad}_x(y) = [x, y]\), the derivation induced by \(x\)
/// on the Lie algebra.
///
/// # Arguments
///
/// * `x_json` - C string pointer with JSON-encoded `Expr` for the Lie algebra element \(x\).
/// * `y_json` - C string pointer with JSON-encoded `Expr` for the Lie algebra element \(y\).
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the result of the algebra
/// adjoint action, or null if deserialization fails or the computation encounters
/// an error.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and
/// returns ownership of a heap-allocated C string that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_adjoint_representation_algebra(
    x_json: *const c_char,
    y_json: *const c_char,
) -> *mut c_char {

    let x: Expr = match from_json_string(
        x_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    let y: Expr = match from_json_string(
        y_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    match adjoint_representation_algebra(
        &x, &y,
    ) {
        | Ok(result) => {
            to_json_string(&result)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

// --- Commutator Table ---

#[unsafe(no_mangle)]

/// Computes the commutator table of a Lie algebra using JSON serialization.
///
/// # Arguments
///
/// * `algebra_json` - C string pointer with JSON-encoded [`LieAlgebra`].
///
/// # Returns
///
/// A C string pointer containing JSON-encoded commutator table (typically a nested
/// collection of brackets of basis elements), or null if deserialization fails or
/// the computation encounters an error.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer and
/// returns ownership of a heap-allocated C string that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_commutator_table(
    algebra_json: *const c_char
) -> *mut c_char {

    let algebra : LieAlgebra = match from_json_string(algebra_json) {
        | Some(a) => a,
        | None => return std::ptr::null_mut(),
    };

    match commutator_table(&algebra) {
        | Ok(table) => {
            to_json_string(&table)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

// --- Jacobi Identity Check ---

#[unsafe(no_mangle)]

/// Checks whether a Lie algebra satisfies the Jacobi identity using JSON serialization.
///
/// The Jacobi identity is a fundamental property of Lie algebras,
/// \([x,[y,z]] + [y,[z,x]] + [z,[x,y]] = 0\).
///
/// # Arguments
///
/// * `algebra_json` - C string pointer with JSON-encoded [`LieAlgebra`].
///
/// # Returns
///
/// `true` if the Jacobi identity holds, `false` otherwise or if deserialization or
/// the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw C string pointer; the
/// caller must ensure it points to a valid JSON string.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_check_jacobi_identity(
    algebra_json: *const c_char
) -> bool {

    let algebra: LieAlgebra =
        match from_json_string(
            algebra_json,
        ) {
            | Some(a) => a,
            | None => return false,
        };

    match check_jacobi_identity(
        &algebra,
    ) {
        | Ok(result) => result,
        | Err(_) => false,
    }
}

// --- Generators ---

#[unsafe(no_mangle)]

/// Returns the standard generators of \(\mathfrak{so}(3)\) as a JSON-encoded list of expressions.
///
/// # Arguments
///
/// This function takes no arguments.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Vec<Expr>` for the generators,
/// or null if serialization fails.
///
/// # Safety
///
/// This function is unsafe because it returns ownership of a heap-allocated C
/// string that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_so3_generators()
-> *mut c_char {

    let generators = so3_generators();

    let exprs: Vec<Expr> = generators
        .into_iter()
        .map(|g| g.0)
        .collect();

    to_json_string(&exprs)
}

/// Returns the SU(2) Lie algebra generators as a JSON string.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_su2_generators()
-> *mut c_char {

    let generators = su2_generators();

    let exprs: Vec<Expr> = generators
        .into_iter()
        .map(|g| g.0)
        .collect();

    to_json_string(&exprs)
}
