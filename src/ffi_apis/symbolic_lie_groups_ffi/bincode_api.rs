use crate::ffi_apis::common::{BincodeBuffer, to_bincode_buffer, from_bincode_buffer};
use crate::symbolic::core::Expr;
use crate::symbolic::lie_groups_and_algebras::{so3, su2, lie_bracket, exponential_map, adjoint_representation_group, adjoint_representation_algebra, LieAlgebra, commutator_table, check_jacobi_identity, so3_generators, su2_generators};

// --- LieAlgebra Creation ---

#[unsafe(no_mangle)]

/// Constructs the Lie algebra \(\mathfrak{so}(3)\) and returns it via bincode serialization.
///
/// The Lie algebra \(\mathfrak{so}(3)\) consists of skew-symmetric \(3\times3\) matrices and
/// is associated with the rotation group SO(3).
///
/// # Arguments
///
/// This function takes no arguments.
///
/// # Returns
///
/// A `BincodeBuffer` encoding a [`LieAlgebra`] representing \(\mathfrak{so}(3)\).
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_lie_algebra_so3()
-> BincodeBuffer {

    let algebra = so3();

    to_bincode_buffer(&algebra)
}

#[unsafe(no_mangle)]

/// Constructs the Lie algebra \(\mathfrak{su}(2)\) and returns it via bincode serialization.
///
/// The Lie algebra \(\mathfrak{su}(2)\) consists of traceless skew-Hermitian \(2\times2\)
/// matrices and is associated with the special unitary group SU(2).
///
/// # Arguments
///
/// This function takes no arguments.
///
/// # Returns
///
/// A `BincodeBuffer` encoding a [`LieAlgebra`] representing \(\mathfrak{su}(2)\).
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_lie_algebra_su2()
-> BincodeBuffer {

    let algebra = su2();

    to_bincode_buffer(&algebra)
}

// --- Lie Bracket ---

#[unsafe(no_mangle)]

/// Computes the Lie bracket \([x, y]\) of two elements of a Lie algebra using bincode serialization.
///
/// # Arguments
///
/// * `x_buf` - `BincodeBuffer` encoding an `Expr` representing \(x\).
/// * `y_buf` - `BincodeBuffer` encoding an `Expr` representing \(y\).
///
/// # Returns
///
/// A `BincodeBuffer` encoding the Lie bracket \([x, y]\) as an `Expr`, or an empty
/// buffer if deserialization fails or the bracket cannot be computed.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_lie_bracket(
    x_buf: BincodeBuffer,
    y_buf: BincodeBuffer,
) -> BincodeBuffer {

    let x : Expr = match from_bincode_buffer(&x_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let y : Expr = match from_bincode_buffer(&y_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    match lie_bracket(&x, &y) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

// --- Exponential Map ---

#[unsafe(no_mangle)]

/// Computes the exponential map from a Lie algebra element to the corresponding Lie group element.
///
/// The exponential map \(\exp(x)\) is approximated by a truncated series of the
/// given order.
///
/// # Arguments
///
/// * `x_buf` - `BincodeBuffer` encoding an `Expr` representing the Lie algebra element \(x\).
/// * `order` - Truncation order for the series expansion of the exponential map.
///
/// # Returns
///
/// A `BincodeBuffer` encoding the resulting group element as an `Expr`, or an empty
/// buffer if deserialization fails or the computation encounters an error.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_exponential_map(
    x_buf: BincodeBuffer,
    order: usize,
) -> BincodeBuffer {

    let x : Expr = match from_bincode_buffer(&x_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    match exponential_map(&x, order) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

// --- Adjoint Representations ---

#[unsafe(no_mangle)]

/// Applies the adjoint representation of a Lie group element to a Lie algebra element.
///
/// This computes \(\mathrm{Ad}_g(x)\), describing how the group element \(g\)
/// conjugates the Lie algebra element \(x\).
///
/// # Arguments
///
/// * `g_buf` - `BincodeBuffer` encoding an `Expr` for the group element \(g\).
/// * `x_buf` - `BincodeBuffer` encoding an `Expr` for the Lie algebra element \(x\).
///
/// # Returns
///
/// A `BincodeBuffer` encoding the result of the group adjoint action as an `Expr`,
/// or an empty buffer if deserialization fails or the computation encounters an error.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_adjoint_representation_group(
    g_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let g : Expr = match from_bincode_buffer(&g_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let x : Expr = match from_bincode_buffer(&x_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    match adjoint_representation_group(
        &g, &x,
    ) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

#[unsafe(no_mangle)]

/// Applies the adjoint representation of a Lie algebra element to another element.
///
/// This computes \(\mathrm{ad}_x(y) = [x, y]\), the derivation induced by \(x\)
/// on the Lie algebra.
///
/// # Arguments
///
/// * `x_buf` - `BincodeBuffer` encoding an `Expr` for the Lie algebra element \(x\).
/// * `y_buf` - `BincodeBuffer` encoding an `Expr` for the Lie algebra element \(y\).
///
/// # Returns
///
/// A `BincodeBuffer` encoding the result of the algebra adjoint action as an
/// `Expr`, or an empty buffer if deserialization fails or the computation
/// encounters an error.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_adjoint_representation_algebra(
    x_buf: BincodeBuffer,
    y_buf: BincodeBuffer,
) -> BincodeBuffer {

    let x : Expr = match from_bincode_buffer(&x_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let y : Expr = match from_bincode_buffer(&y_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    match adjoint_representation_algebra(
        &x, &y,
    ) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

// --- Commutator Table ---

#[unsafe(no_mangle)]

/// Computes the commutator table (structure constants) of a Lie algebra.
///
/// # Arguments
///
/// * `algebra_buf` - `BincodeBuffer` encoding a [`LieAlgebra`].
///
/// # Returns
///
/// A `BincodeBuffer` encoding the commutator table (typically as a collection of
/// brackets of basis elements), or an empty buffer if deserialization fails or the
/// computation encounters an error.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_commutator_table(
    algebra_buf: BincodeBuffer
) -> BincodeBuffer {

    let algebra : LieAlgebra = match from_bincode_buffer(&algebra_buf) {
        | Some(a) => a,
        | None => return BincodeBuffer::empty(),
    };

    match commutator_table(&algebra) {
        | Ok(table) => {
            to_bincode_buffer(&table)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

// --- Jacobi Identity Check ---

#[unsafe(no_mangle)]

/// Checks whether a Lie algebra satisfies the Jacobi identity.
///
/// The Jacobi identity is a fundamental property of Lie algebras,
/// \([x,[y,z]] + [y,[z,x]] + [z,[x,y]] = 0\).
///
/// # Arguments
///
/// * `algebra_buf` - `BincodeBuffer` encoding a [`LieAlgebra`].
///
/// # Returns
///
/// `true` if the Jacobi identity holds, `false` if it does not or if
/// deserialization or the computation fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must ensure the buffer encodes a valid `LieAlgebra`.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_check_jacobi_identity(
    algebra_buf: BincodeBuffer
) -> bool {

    let algebra: LieAlgebra =
        match from_bincode_buffer(
            &algebra_buf,
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

/// Returns the standard generators of \(\mathfrak{so}(3)\) via bincode serialization.
///
/// # Arguments
///
/// This function takes no arguments.
///
/// # Returns
///
/// A `BincodeBuffer` encoding a `Vec<Expr>` whose entries represent a basis of
/// generators for \(\mathfrak{so}(3)\).
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_so3_generators()
-> BincodeBuffer {

    let generators = so3_generators();

    let exprs: Vec<Expr> = generators
        .into_iter()
        .map(|g| g.0)
        .collect();

    to_bincode_buffer(&exprs)
}

#[unsafe(no_mangle)]

/// Returns the standard generators of \(\mathfrak{su}(2)\) via bincode serialization.
///
/// # Arguments
///
/// This function takes no arguments.
///
/// # Returns
///
/// A `BincodeBuffer` encoding a `Vec<Expr>` whose entries represent a basis of
/// generators for \(\mathfrak{su}(2)\).
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; callers
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_su2_generators()
-> BincodeBuffer {

    let generators = su2_generators();

    let exprs: Vec<Expr> = generators
        .into_iter()
        .map(|g| g.0)
        .collect();

    to_bincode_buffer(&exprs)
}
