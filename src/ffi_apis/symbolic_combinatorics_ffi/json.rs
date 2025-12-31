use std::os::raw::c_char;

use crate::ffi_apis::common::from_json_string;
use crate::ffi_apis::common::to_json_string;
use crate::symbolic::combinatorics::bell_number;
use crate::symbolic::combinatorics::catalan_number;
use crate::symbolic::combinatorics::combinations;
use crate::symbolic::combinatorics::permutations;
use crate::symbolic::combinatorics::stirling_number_second_kind;
use crate::symbolic::core::Expr;


/// Computes the number of permutations symbolically using JSON-encoded `Expr` arguments.
///
/// This corresponds to the falling factorial \( P(n,k) = n! / (n-k)! \) when `n` and `k`
/// are integers, but operates on general symbolic expressions.
///
/// # Arguments
///
/// * `n_json` - C string pointer to JSON encoding an `Expr` for the population size `n`.
/// * `k_json` - C string pointer to JSON encoding an `Expr` for the selection size `k`.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the symbolic permutation count,
/// or null on deserialization failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and returns
/// ownership of a heap-allocated C string.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_permutations(
    n_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let n: Expr = match from_json_string(
        n_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    let k: Expr = match from_json_string(
        k_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    let result = permutations(n, k);

    to_json_string(&result)
}

/// Computes the number of combinations symbolically using JSON-encoded `Expr` arguments.
///
/// This corresponds to the binomial coefficient \( C(n,k) = n! / (k!(n-k)!) \) when
/// `n` and `k` are integers, but operates on general symbolic expressions.
///
/// # Arguments
///
/// * `n_json` - C string pointer to JSON encoding an `Expr` for the population size `n`.
/// * `k_json` - C string pointer to JSON encoding an `Expr` for the selection size `k`.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the symbolic combination count,
/// or null on deserialization failure.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw C string pointers and returns
/// ownership of a heap-allocated C string.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_combinations(
    n_json: *const c_char,
    k_json: *const c_char,
) -> *mut c_char {

    let n: Expr = match from_json_string(
        n_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    let k: Expr = match from_json_string(
        k_json,
    ) {
        | Some(e) => e,
        | None => {
            return std::ptr::null_mut();
        },
    };

    let result = combinations(&n, k);

    to_json_string(&result)
}

/// Computes the \(n\)-th Catalan number symbolically and returns it as JSON-encoded `Expr`.
///
/// Catalan numbers count many combinatorial structures, such as binary trees, Dyck paths,
/// and non-crossing partitions.
///
/// # Arguments
///
/// * `n` - Index of the Catalan number to compute.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the \(n\)-th Catalan number.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point and returns
/// ownership of a heap-allocated C string.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_catalan_number(
    n: usize
) -> *mut c_char {

    let result = catalan_number(n);

    to_json_string(&result)
}

/// Computes a Stirling number of the second kind symbolically and returns it as JSON-encoded `Expr`.
///
/// Stirling numbers of the second kind \( S(n,k) \) count partitions of an \(n\)-element
/// set into \(k\) non-empty unlabeled blocks.
///
/// # Arguments
///
/// * `n` - Total number of elements.
/// * `k` - Number of non-empty blocks.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for \( S(n,k) \).
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point and returns
/// ownership of a heap-allocated C string.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_stirling_number_second_kind(
    n: usize,
    k: usize,
) -> *mut c_char {

    let result =
        stirling_number_second_kind(
            n, k,
        );

    to_json_string(&result)
}


/// Computes the \(n\)-th Bell number symbolically and returns it as JSON-encoded `Expr`.
///
/// Bell numbers count the total number of set partitions of an \(n\)-element set.
///
/// # Arguments
///
/// * `n` - Index of the Bell number to compute.
///
/// # Returns
///
/// A C string pointer containing JSON-encoded `Expr` for the \(n\)-th Bell number.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point and returns
/// ownership of a heap-allocated C string.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_bell_number(
    n: usize
) -> *mut c_char {

    let result = bell_number(n);

    to_json_string(&result)
}
