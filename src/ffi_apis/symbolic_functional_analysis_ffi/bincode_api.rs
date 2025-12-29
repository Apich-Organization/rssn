use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::functional_analysis::*;

/// Constructs a Hilbert space from a bincode-encoded description.
///
/// The input buffer encodes a [`HilbertSpace`] specification (e.g., underlying
/// function space, inner product, and measure), which is deserialized and
/// returned in canonical internal form.
///
/// # Arguments
///
/// * `buf` - `BincodeBuffer` containing a serialized `HilbertSpace` description.
///
/// # Returns
///
/// A `BincodeBuffer` containing the canonicalized `HilbertSpace`, or an empty
/// buffer if deserialization fails.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; the caller
/// must treat the returned buffer as opaque and only pass it to compatible APIs.
#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_hilbert_space_create(
    buf: BincodeBuffer
) -> BincodeBuffer {

    let space : HilbertSpace = match from_bincode_buffer(&buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&space)
}

#[no_mangle]

/// Computes the inner product of two functions in a Hilbert space using bincode serialization.
///
/// Given a Hilbert space and two symbolic functions \(f\) and \(g\), this evaluates
/// the inner product \(\langle f, g \rangle\) according to the space's inner
/// product structure.
///
/// # Arguments
///
/// * `space_buf` - `BincodeBuffer` encoding a [`HilbertSpace`].
/// * `f_buf` - `BincodeBuffer` encoding an `Expr` for \(f\).
/// * `g_buf` - `BincodeBuffer` encoding an `Expr` for \(g\).
///
/// # Returns
///
/// A `BincodeBuffer` containing the symbolic inner product value (typically an `Expr`).
/// Returns an empty buffer if any input fails to deserialize.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; the caller
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

pub unsafe extern "C" fn rssn_bincode_inner_product(
    space_buf: BincodeBuffer,
    f_buf: BincodeBuffer,
    g_buf: BincodeBuffer,
) -> BincodeBuffer {

    let space : HilbertSpace = match from_bincode_buffer(&space_buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    let f : Expr = match from_bincode_buffer(&f_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let g : Expr = match from_bincode_buffer(&g_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result =
        inner_product(&space, &f, &g);

    to_bincode_buffer(&result)
}

#[no_mangle]

/// Computes the norm of a function in a Hilbert space using bincode serialization.
///
/// Given a Hilbert space and a symbolic function \(f\), this evaluates the norm
/// \(\|f\| = \sqrt{\langle f, f \rangle}\) induced by the inner product.
///
/// # Arguments
///
/// * `space_buf` - `BincodeBuffer` encoding a [`HilbertSpace`].
/// * `f_buf` - `BincodeBuffer` encoding an `Expr` for \(f\).
///
/// # Returns
///
/// A `BincodeBuffer` containing the symbolic norm value (typically an `Expr`).
/// Returns an empty buffer if any input fails to deserialize.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; the caller
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

pub unsafe extern "C" fn rssn_bincode_norm(
    space_buf: BincodeBuffer,
    f_buf: BincodeBuffer,
) -> BincodeBuffer {

    let space : HilbertSpace = match from_bincode_buffer(&space_buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    let f : Expr = match from_bincode_buffer(&f_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = norm(&space, &f);

    to_bincode_buffer(&result)
}

#[no_mangle]

/// Applies the Gram–Schmidt process to produce an orthonormal basis in a Hilbert space.
///
/// Given a Hilbert space and a list of symbolic basis vectors, this performs the
/// Gram–Schmidt orthonormalization procedure to obtain an orthonormal basis.
///
/// # Arguments
///
/// * `space_buf` - `BincodeBuffer` encoding a [`HilbertSpace`].
/// * `basis_buf` - `BincodeBuffer` encoding a `Vec<Expr>` of basis vectors.
///
/// # Returns
///
/// A `BincodeBuffer` containing a `Vec<Expr>` for the orthonormal basis. Returns an
/// empty buffer if any input fails to deserialize.
///
/// # Safety
///
/// This function is unsafe because it is exposed as an FFI entry point; the caller
/// must treat the returned buffer as opaque and only pass it to compatible APIs.

pub unsafe extern "C" fn rssn_bincode_gram_schmidt(
    space_buf: BincodeBuffer,
    basis_buf: BincodeBuffer,
) -> BincodeBuffer {

    let space : HilbertSpace = match from_bincode_buffer(&space_buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    let basis : Vec<Expr> = match from_bincode_buffer(&basis_buf) {
        | Some(b) => b,
        | None => return BincodeBuffer::empty(),
    };

    let result =
        gram_schmidt(&space, &basis);

    to_bincode_buffer(&result)
}
