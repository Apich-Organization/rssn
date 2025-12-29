use crate::symbolic::combinatorics::bell_number;
use crate::symbolic::combinatorics::catalan_number;
use crate::symbolic::combinatorics::combinations;
use crate::symbolic::combinatorics::permutations;
use crate::symbolic::combinatorics::stirling_number_second_kind;
use crate::symbolic::core::Expr;

/// Computes the number of permutations symbolically as an `Expr`.
///
/// This corresponds to \( P(n,k) = n! / (n-k)! \) for integer `n` and `k`, but also
/// supports symbolic `Expr` arguments.
///
/// # Arguments
///
/// * `n` - Pointer to an `Expr` representing the population size.
/// * `k` - Pointer to an `Expr` representing the selection size.
///
/// # Returns
///
/// A newly allocated `Expr` pointer representing the symbolic permutation count.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw `Expr` pointers and returns
/// ownership of a heap-allocated `Expr` to the caller.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_permutations(
    n: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let result = permutations(
        (*n).clone(),
        (*k).clone(),
    );

    Box::into_raw(Box::new(result))
}


/// Computes the number of combinations symbolically as an `Expr`.
///
/// This corresponds to the binomial coefficient \( C(n,k) = n! / (k!(n-k)!) \) for
/// integer `n` and `k`, but also supports symbolic `Expr` arguments.
///
/// # Arguments
///
/// * `n` - Pointer to an `Expr` representing the population size.
/// * `k` - Pointer to an `Expr` representing the selection size.
///
/// # Returns
///
/// A newly allocated `Expr` pointer representing the symbolic combination count.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw `Expr` pointers and returns
/// ownership of a heap-allocated `Expr` to the caller.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_combinations(
    n: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let result = combinations(
        &(*n),
        (*k).clone(),
    );

    Box::into_raw(Box::new(result))
}


/// Computes the \(n\)-th Catalan number symbolically and returns it as an `Expr` pointer.
///
/// Catalan numbers count many combinatorial structures, such as full binary trees,
/// Dyck paths, and non-crossing partitions.
///
/// # Arguments
///
/// * `n` - Index of the Catalan number to compute.
///
/// # Returns
///
/// A newly allocated `Expr` pointer representing the \(n\)-th Catalan number.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point and returns
/// ownership of a heap-allocated `Expr` to the caller.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_catalan_number(
    n: usize
) -> *mut Expr {

    let result = catalan_number(n);

    Box::into_raw(Box::new(result))
}


/// Computes a Stirling number of the second kind symbolically and returns it as an `Expr` pointer.
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
/// A newly allocated `Expr` pointer representing \( S(n,k) \).
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point and returns
/// ownership of a heap-allocated `Expr` to the caller.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_stirling_number_second_kind(
    n: usize,
    k: usize,
) -> *mut Expr {

    let result =
        stirling_number_second_kind(
            n, k,
        );

    Box::into_raw(Box::new(result))
}

/// Computes the \(n\)-th Bell number symbolically and returns it as an `Expr` pointer.
///
/// Bell numbers count the total number of set partitions of an \(n\)-element set.
///
/// # Arguments
///
/// * `n` - Index of the Bell number to compute.
///
/// # Returns
///
/// A newly allocated `Expr` pointer representing the \(n\)-th Bell number.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point and returns
/// ownership of a heap-allocated `Expr` to the caller.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bell_number(
    n: usize
) -> *mut Expr {

    let result = bell_number(n);

    Box::into_raw(Box::new(result))
}
