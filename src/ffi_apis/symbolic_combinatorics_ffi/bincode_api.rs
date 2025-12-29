use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::combinatorics::bell_number;
use crate::symbolic::combinatorics::catalan_number;
use crate::symbolic::combinatorics::combinations;
use crate::symbolic::combinatorics::permutations;
use crate::symbolic::combinatorics::stirling_number_second_kind;
use crate::symbolic::core::Expr;

/// Computes the number of permutations symbolically using bincode-encoded `Expr` arguments.
///
/// This corresponds to the falling factorial \( P(n,k) = n! / (n-k)! \) when `n` and `k`
/// are specialized to integers, but operates on general symbolic expressions.
///
/// # Arguments
///
/// * `n_buf` - Bincode buffer encoding an `Expr` representing the population size `n`.
/// * `k_buf` - Bincode buffer encoding an `Expr` representing the selection size `k`.
///
/// # Returns
///
/// A bincode buffer encoding an `Expr` for the symbolic permutation count.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw bincode buffers that must
/// contain valid serialized `Expr` values.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_permutations(
    n_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n : Expr = match from_bincode_buffer(&n_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let k : Expr = match from_bincode_buffer(&k_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = permutations(n, k);

    to_bincode_buffer(&result)
}

/// Computes the number of combinations symbolically using bincode-encoded `Expr` arguments.
///
/// This corresponds to the binomial coefficient \( C(n,k) = n! / (k!(n-k)!) \) when
/// `n` and `k` are integers, but operates on general symbolic expressions.
///
/// # Arguments
///
/// * `n_buf` - Bincode buffer encoding an `Expr` representing the population size `n`.
/// * `k_buf` - Bincode buffer encoding an `Expr` representing the selection size `k`.
///
/// # Returns
///
/// A bincode buffer encoding an `Expr` for the symbolic combination count.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw bincode buffers that must
/// contain valid serialized `Expr` values.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_combinations(
    n_buf: BincodeBuffer,
    k_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n : Expr = match from_bincode_buffer(&n_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let k : Expr = match from_bincode_buffer(&k_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = combinations(&n, k);

    to_bincode_buffer(&result)
}

/// Computes the \(n\)-th Catalan number symbolically and returns it as a bincode-encoded `Expr`.
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
/// A bincode buffer encoding an `Expr` for the \(n\)-th Catalan number.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point, but it does
/// not dereference raw pointers.

#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_catalan_number(
    n: usize
) -> BincodeBuffer {

    let result = catalan_number(n);

    to_bincode_buffer(&result)
}

#[no_mangle]

/// Computes a Stirling number of the second kind symbolically and returns it via bincode.
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
/// A bincode buffer encoding an `Expr` for \( S(n,k) \).
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point, but it does
/// not dereference raw pointers.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_stirling_number_second_kind(
    n: usize,
    k: usize,
) -> BincodeBuffer {

    let result =
        stirling_number_second_kind(
            n, k,
        );

    to_bincode_buffer(&result)
}

/// Computes the \(n\)-th Bell number symbolically and returns it via bincode.
///
/// Bell numbers count the total number of set partitions of an \(n\)-element set.
///
/// # Arguments
///
/// * `n` - Index of the Bell number to compute.
///
/// # Returns
///
/// A bincode buffer encoding an `Expr` for the \(n\)-th Bell number.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point, but it does
/// not dereference raw pointers.
#[no_mangle]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_bell_number(
    n: usize
) -> BincodeBuffer {

    let result = bell_number(n);

    to_bincode_buffer(&result)
}
