use std::os::raw::c_char;

use crate::symbolic::core::Expr;
use crate::symbolic::lie_groups_and_algebras::{LieAlgebra, so3, su2, lie_bracket, exponential_map, adjoint_representation_group, adjoint_representation_algebra, commutator_table, check_jacobi_identity, so3_generators, su2_generators};

// --- LieAlgebra ---

#[unsafe(no_mangle)]

/// Constructs the Lie algebra \(\mathfrak{so}(3)\) and returns a heap-allocated handle.
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
/// A pointer to a heap-allocated [`LieAlgebra`] representing \(\mathfrak{so}(3)\).
///
/// # Safety
///
/// This function is unsafe because it returns ownership of a heap-allocated
/// `LieAlgebra` that must later be freed with [`rssn_lie_algebra_free`].

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_lie_algebra_so3_create()
-> *mut LieAlgebra {

    let algebra = so3();

    Box::into_raw(Box::new(algebra))
}

#[unsafe(no_mangle)]

/// Constructs the Lie algebra \(\mathfrak{su}(2)\) and returns a heap-allocated handle.
///
/// The Lie algebra \(\mathfrak{su}(2)\) consists of traceless skew-Hermitian \(2\times2\)
/// matrices associated with the special unitary group SU(2).
///
/// # Arguments
///
/// This function takes no arguments.
///
/// # Returns
///
/// A pointer to a heap-allocated [`LieAlgebra`] representing \(\mathfrak{su}(2)\).
///
/// # Safety
///
/// This function is unsafe because it returns ownership of a heap-allocated
/// `LieAlgebra` that must later be freed with [`rssn_lie_algebra_free`].

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_lie_algebra_su2_create()
-> *mut LieAlgebra {

    let algebra = su2();

    Box::into_raw(Box::new(algebra))
}

#[unsafe(no_mangle)]

/// Frees a Lie algebra previously created by one of the constructor functions.
///
/// # Arguments
///
/// * `ptr` - Pointer to a heap-allocated [`LieAlgebra`] returned by
///   `rssn_lie_algebra_so3_create` or `rssn_lie_algebra_su2_create`.
///
/// # Returns
///
/// This function does not return a value.
///
/// # Safety
///
/// This function is unsafe because it takes ownership of a raw pointer. The pointer
/// must either be null or have been allocated by the corresponding constructor, and
/// must not be used after this call.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_lie_algebra_free(
    ptr: *mut LieAlgebra
) {

    unsafe {

        if !ptr.is_null() {

            let _ = Box::from_raw(ptr);
        }
    }
}

#[unsafe(no_mangle)]

/// Returns the dimension of a Lie algebra.
///
/// # Arguments
///
/// * `ptr` - Pointer to a [`LieAlgebra`] value.
///
/// # Returns
///
/// The dimension of the Lie algebra as a `usize`.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer; the caller must
/// ensure `ptr` points to a valid `LieAlgebra`.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub const unsafe extern "C" fn rssn_lie_algebra_get_dimension(
    ptr: *const LieAlgebra
) -> usize {

    unsafe {

        (*ptr).dimension
    }
}

#[unsafe(no_mangle)]

/// Returns the name of a Lie algebra as a newly allocated C string.
///
/// # Arguments
///
/// * `ptr` - Pointer to a [`LieAlgebra`] value.
///
/// # Returns
///
/// A pointer to a heap-allocated, null-terminated C string containing the algebra
/// name. The caller is responsible for freeing this string using the appropriate
/// deallocation routine.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer and returns
/// ownership of a heap-allocated C string.

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

pub unsafe extern "C" fn rssn_lie_algebra_get_name(
    ptr: *const LieAlgebra
) -> *mut c_char {

    unsafe {

        let name = &(*ptr).name;

        std::ffi::CString::new(
            name.as_str(),
        )
        .unwrap()
        .into_raw()
    }
}

#[unsafe(no_mangle)]

/// Returns a basis element of a Lie algebra as a symbolic expression.
///
/// # Arguments
///
/// * `ptr` - Pointer to a [`LieAlgebra`] value.
/// * `index` - Zero-based index into the algebra's basis.
///
/// # Returns
///
/// A pointer to a newly allocated `Expr` representing the requested basis element,
/// or null if `index` is out of bounds.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer and returns
/// ownership of a heap-allocated `Expr` that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_lie_algebra_get_basis_element(
    ptr: *const LieAlgebra,
    index: usize,
) -> *mut Expr {

    unsafe {

        let algebra = &*ptr;

        if index >= algebra.basis.len()
        {

            return std::ptr::null_mut(
            );
        }

        Box::into_raw(Box::new(
            algebra.basis[index]
                .0
                .clone(),
        ))
    }
}

// --- Lie Bracket ---

#[unsafe(no_mangle)]

/// Computes the Lie bracket of two elements of a Lie algebra.
///
/// # Arguments
///
/// * `x` - Pointer to an `Expr` representing the first element.
/// * `y` - Pointer to an `Expr` representing the second element.
///
/// # Returns
///
/// A pointer to a newly allocated `Expr` representing the Lie bracket \([x, y]\),
/// or null if the bracket cannot be computed.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers and returns
/// ownership of a heap-allocated `Expr` that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_lie_bracket(
    x: *const Expr,
    y: *const Expr,
) -> *mut Expr {

    unsafe {

        match lie_bracket(&*x, &*y) {
            | Ok(result) => {
                Box::into_raw(Box::new(
                    result,
                ))
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

// --- Exponential Map ---

#[unsafe(no_mangle)]

/// Applies the exponential map to a Lie algebra element, returning the corresponding group element.
///
/// The exponential map is approximated by a truncated series of the specified order.
///
/// # Arguments
///
/// * `x` - Pointer to an `Expr` representing the Lie algebra element.
/// * `order` - Truncation order for the series expansion of the exponential map.
///
/// # Returns
///
/// A pointer to a newly allocated `Expr` representing the group element, or null
/// if the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer and returns
/// ownership of a heap-allocated `Expr` that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_exponential_map(
    x: *const Expr,
    order: usize,
) -> *mut Expr {

    unsafe {

        match exponential_map(
            &*x, order,
        ) {
            | Ok(result) => {
                Box::into_raw(Box::new(
                    result,
                ))
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    }
}

// --- Adjoint Representations ---

#[unsafe(no_mangle)]

/// Applies the adjoint representation of a Lie group element to a Lie algebra element.
///
/// This computes \(\mathrm{Ad}_g(x)\), describing how the group element \(g\)
/// conjugates the algebra element \(x\).
///
/// # Arguments
///
/// * `g` - Pointer to an `Expr` representing the group element.
/// * `x` - Pointer to an `Expr` representing the Lie algebra element.
///
/// # Returns
///
/// A pointer to a newly allocated `Expr` representing the result of the group
/// adjoint action, or null if the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers and returns
/// ownership of a heap-allocated `Expr` that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_adjoint_representation_group(
    g: *const Expr,
    x: *const Expr,
) -> *mut Expr {

    unsafe {

        match adjoint_representation_group(
        &*g, &*x,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
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
/// * `x` - Pointer to an `Expr` representing the first Lie algebra element.
/// * `y` - Pointer to an `Expr` representing the second Lie algebra element.
///
/// # Returns
///
/// A pointer to a newly allocated `Expr` representing the result of the algebra
/// adjoint action, or null if the computation fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers and returns
/// ownership of a heap-allocated `Expr` that must be freed by the caller.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_adjoint_representation_algebra(
    x: *const Expr,
    y: *const Expr,
) -> *mut Expr {

    unsafe {

        match adjoint_representation_algebra(
        &*x, &*y,
    ) {
        | Ok(result) => {
            Box::into_raw(Box::new(
                result,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
    }
}

// --- Commutator Table ---

#[unsafe(no_mangle)]

/// Computes the commutator table of a Lie algebra and returns it as a flattened array of expressions.
///
/// # Arguments
///
/// * `algebra` - Pointer to a [`LieAlgebra`] value.
/// * `out_rows` - Pointer to a `usize` that will be filled with the number of rows.
/// * `out_cols` - Pointer to a `usize` that will be filled with the number of columns.
///
/// # Returns
///
/// A pointer to a heap-allocated array of `*mut Expr` representing the commutator
/// table entries in row-major order, or null if the computation fails. On error,
/// `out_rows` and `out_cols` are set to 0.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers and returns
/// ownership of heap-allocated memory. The caller must ensure all pointers are
/// valid and must correctly free the returned expressions and array.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_commutator_table(
    algebra: *const LieAlgebra,
    out_rows: *mut usize,
    out_cols: *mut usize,
) -> *mut *mut Expr {

    unsafe {

        match commutator_table(
            &*algebra,
        ) {
            | Ok(table) => {

                let rows = table.len();

                let cols = if rows > 0 {

                    table[0].len()
                } else {

                    0
                };

                *out_rows = rows;

                *out_cols = cols;

                let mut flat_ptrs =
                    Vec::with_capacity(
                        rows * cols,
                    );

                for row in table {

                    for elem in row {

                        flat_ptrs.push(
                        Box::into_raw(
                            Box::new(
                                elem,
                            ),
                        ),
                    );
                    }
                }

                let ptr = flat_ptrs
                    .as_mut_ptr();

                std::mem::forget(
                    flat_ptrs,
                );

                ptr
            },
            | Err(_) => {

                *out_rows = 0;

                *out_cols = 0;

                std::ptr::null_mut()
            },
        }
    }
}

// --- Jacobi Identity Check ---

#[unsafe(no_mangle)]

/// Checks whether a Lie algebra satisfies the Jacobi identity.
///
/// The Jacobi identity is a fundamental property of Lie algebras,
/// \(\[x,\[y,z\]\] + \[y,\[z,x\]\] + \[z,\[x,y\]\] = 0\).
///
/// # Arguments
///
/// * `algebra` - Pointer to a [`LieAlgebra`] value.
///
/// # Returns
///
/// `true` if the Jacobi identity holds, `false` otherwise or if the computation
/// fails.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer; the caller must
/// ensure `algebra` points to a valid `LieAlgebra`.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_check_jacobi_identity(
    algebra: *const LieAlgebra
) -> bool {

    unsafe {

        match check_jacobi_identity(
            &*algebra,
        ) {
            | Ok(result) => result,
            | Err(_) => false,
        }
    }
}

// --- Generators ---

#[unsafe(no_mangle)]

/// Returns the standard generators of \(\mathfrak{so}(3)\) as a dynamically allocated array of expressions.
///
/// # Arguments
///
/// * `out_len` - Pointer to a `usize` that will be filled with the number of generators.
///
/// # Returns
///
/// A pointer to a heap-allocated array of `*mut Expr` representing the generators.
/// The caller is responsible for freeing each `Expr` and the array itself.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer and returns
/// ownership of heap-allocated memory. The caller must ensure `out_len` is valid
/// and correctly free the returned memory.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_so3_generators(
    out_len: *mut usize
) -> *mut *mut Expr {

    unsafe {

        let generators =
            so3_generators();

        *out_len = generators.len();

        let mut ptrs =
            Vec::with_capacity(
                generators.len(),
            );

        for r#gen in generators {

            ptrs.push(Box::into_raw(
                Box::new(r#gen.0),
            ));
        }

        let ptr = ptrs.as_mut_ptr();

        std::mem::forget(ptrs);

        ptr
    }
}

#[unsafe(no_mangle)]

/// Returns the standard generators of \(\mathfrak{su}(2)\) as a dynamically allocated array of expressions.
///
/// # Arguments
///
/// * `out_len` - Pointer to a `usize` that will be filled with the number of generators.
///
/// # Returns
///
/// A pointer to a heap-allocated array of `*mut Expr` representing the generators.
/// The caller is responsible for freeing each `Expr` and the array itself.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer and returns
/// ownership of heap-allocated memory. The caller must ensure `out_len` is valid
/// and correctly free the returned memory.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_su2_generators(
    out_len: *mut usize
) -> *mut *mut Expr {

    unsafe {

        let generators =
            su2_generators();

        *out_len = generators.len();

        let mut ptrs =
            Vec::with_capacity(
                generators.len(),
            );

        for r#gen in generators {

            ptrs.push(Box::into_raw(
                Box::new(r#gen.0),
            ));
        }

        let ptr = ptrs.as_mut_ptr();

        std::mem::forget(ptrs);

        ptr
    }
}
