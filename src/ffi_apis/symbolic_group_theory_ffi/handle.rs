use std::collections::HashMap;

use crate::symbolic::core::Expr;
use crate::symbolic::group_theory::Group;
use crate::symbolic::group_theory::GroupElement;
use crate::symbolic::group_theory::Representation;
use crate::symbolic::group_theory::character;

// --- Group ---

#[unsafe(no_mangle)]

/// Constructs a group from raw pointers describing its elements, multiplication table, and identity.
///
/// The group operation is encoded as a Cayley table over the given elements, represented by
/// parallel arrays of left factors, right factors, and products.
///
/// # Arguments
///
/// * `elements_ptr` - Pointer to an array of pointers to `Expr` giving the group elements.
/// * `elements_len` - Number of elements referenced by `elements_ptr`.
/// * `keys_a_ptr` - Pointer to an array of pointers to `Expr` for the left factors in the multiplication table.
/// * `keys_b_ptr` - Pointer to an array of pointers to `Expr` for the right factors in the multiplication table.
/// * `values_ptr` - Pointer to an array of pointers to `Expr` for the products in the multiplication table.
/// * `table_len` - Number of entries in the multiplication table arrays.
/// * `identity_ptr` - Pointer to an `Expr` representing the identity element of the group.
///
/// # Returns
///
/// A pointer to a heap-allocated [`Group`] constructed from the supplied data.
///
/// # Safety
///
/// This function is unsafe because it dereferences multiple raw pointers and assumes
/// they form consistent arrays of valid `Expr` objects. The returned `Group` must be
/// freed with [`rssn_group_free`].

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_group_create(
    elements_ptr: *const *const Expr,
    elements_len: usize,
    keys_a_ptr: *const *const Expr,
    keys_b_ptr: *const *const Expr,
    values_ptr: *const *const Expr,
    table_len: usize,
    identity_ptr: *const Expr,
) -> *mut Group {

    unsafe {

        let elements_slice =
            std::slice::from_raw_parts(
                elements_ptr,
                elements_len,
            );

        let elements: Vec<
            GroupElement,
        > = elements_slice
            .iter()
            .map(|&p| {

                GroupElement(
                    (*p).clone(),
                )
            })
            .collect();

        let keys_a_slice =
            std::slice::from_raw_parts(
                keys_a_ptr,
                table_len,
            );

        let keys_b_slice =
            std::slice::from_raw_parts(
                keys_b_ptr,
                table_len,
            );

        let values_slice =
            std::slice::from_raw_parts(
                values_ptr,
                table_len,
            );

        let mut multiplication_table =
            HashMap::new();

        for i in 0 .. table_len {

            let a = GroupElement(
                (*keys_a_slice[i])
                    .clone(),
            );

            let b = GroupElement(
                (*keys_b_slice[i])
                    .clone(),
            );

            let val = GroupElement(
                (*values_slice[i])
                    .clone(),
            );

            multiplication_table
                .insert((a, b), val);
        }

        let identity = GroupElement(
            (*identity_ptr).clone(),
        );

        let group = Group::new(
            elements,
            multiplication_table,
            identity,
        );

        Box::into_raw(Box::new(group))
    }
}

#[unsafe(no_mangle)]

/// Frees a group previously created by [`rssn_group_create`].
///
/// # Arguments
///
/// * `ptr` - Pointer to a heap-allocated [`Group`] returned by `rssn_group_create`.
///
/// # Returns
///
/// This function does not return a value.
///
/// # Safety
///
/// This function is unsafe because it takes ownership of a raw pointer. The pointer
/// must either be null or have been allocated by `rssn_group_create`, and must not
/// be used after this call.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_group_free(
    ptr: *mut Group
) {

    unsafe {

        if !ptr.is_null() {

            let _ = Box::from_raw(ptr);
        }
    }
}

#[unsafe(no_mangle)]

/// Multiplies two group elements using a group handle and returns the product expression.
///
/// # Arguments
///
/// * `group` - Pointer to a [`Group`] value.
/// * `a` - Pointer to an `Expr` representing the left factor.
/// * `b` - Pointer to an `Expr` representing the right factor.
///
/// # Returns
///
/// A pointer to a newly allocated `Expr` representing the product, or null if the
/// elements are not in the group or the multiplication is undefined.
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

pub unsafe extern "C" fn rssn_group_multiply(
    group: *const Group,
    a: *const Expr,
    b: *const Expr,
) -> *mut Expr {

    unsafe {

        let ga =
            GroupElement((*a).clone());

        let gb =
            GroupElement((*b).clone());

        match (*group)
            .multiply(&ga, &gb)
        {
            | Some(result) => {
                Box::into_raw(Box::new(
                    result.0,
                ))
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    }
}

#[unsafe(no_mangle)]

/// Computes the inverse of a group element using a group handle.
///
/// # Arguments
///
/// * `group` - Pointer to a [`Group`] value.
/// * `a` - Pointer to an `Expr` representing the element whose inverse is sought.
///
/// # Returns
///
/// A pointer to a newly allocated `Expr` representing the inverse element, or null
/// if the element has no inverse in the group.
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

pub unsafe extern "C" fn rssn_group_inverse(
    group: *const Group,
    a: *const Expr,
) -> *mut Expr {

    unsafe {

        let ga =
            GroupElement((*a).clone());

        match (*group).inverse(&ga) {
            | Some(result) => {
                Box::into_raw(Box::new(
                    result.0,
                ))
            },
            | None => {
                std::ptr::null_mut()
            },
        }
    }
}

#[unsafe(no_mangle)]

/// Tests whether a group is Abelian using a group handle.
///
/// A group is Abelian if its operation is commutative for all pairs of elements.
///
/// # Arguments
///
/// * `group` - Pointer to a [`Group`] value.
///
/// # Returns
///
/// `true` if the group is Abelian, `false` otherwise.
///
/// # Safety
///
/// This function is unsafe because it dereferences a raw pointer; the caller must
/// ensure `group` points to a valid [`Group`].

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_group_is_abelian(
    group: *const Group
) -> bool {

    unsafe {

        (*group).is_abelian()
    }
}

#[unsafe(no_mangle)]

/// Computes the order of a group element using a group handle.
///
/// The order is the smallest positive integer `n` such that `a^n` is the identity,
/// if such an integer exists.
///
/// # Arguments
///
/// * `group` - Pointer to a [`Group`] value.
/// * `a` - Pointer to an `Expr` representing the element.
///
/// # Returns
///
/// The order of the element as a `usize`, or `0` if the order is undefined.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers; the caller must
/// ensure they point to a valid group and element.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_group_element_order(
    group: *const Group,
    a: *const Expr,
) -> usize {

    unsafe {

        let ga =
            GroupElement((*a).clone());

        (*group)
            .element_order(&ga)
            .unwrap_or(0)
    }
}

#[unsafe(no_mangle)]

/// Computes the center of a group and returns its elements as a dynamically allocated array of expressions.
///
/// The center consists of elements that commute with every element of the group.
///
/// # Arguments
///
/// * `group` - Pointer to a [`Group`] value.
/// * `out_len` - Pointer to a `usize` that will be filled with the number of elements in the center.
///
/// # Returns
///
/// A pointer to a heap-allocated array of `*mut Expr` representing the center
/// elements. The caller is responsible for freeing each `Expr` and the array
/// itself.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers and returns
/// ownership of heap-allocated memory. The caller must ensure `group` and
/// `out_len` are valid pointers and must correctly manage the returned memory.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_group_center(
    group: *const Group,
    out_len: *mut usize,
) -> *mut *mut Expr {

    unsafe {

        let center = (*group).center();

        *out_len = center.len();

        let mut out_ptrs =
            Vec::with_capacity(
                center.len(),
            );

        for elem in center {

            out_ptrs.push(
                Box::into_raw(
                    Box::new(elem.0),
                ),
            );
        }

        let ptr = out_ptrs.as_mut_ptr();

        std::mem::forget(out_ptrs);

        ptr
    }
}

// --- Representation ---

#[unsafe(no_mangle)]

/// Constructs a group representation from raw pointers describing its carrier set and action matrices.
///
/// The representation is defined by a set of group elements and a mapping that
/// assigns a matrix (or linear operator) to each element.
///
/// # Arguments
///
/// * `elements_ptr` - Pointer to an array of pointers to `Expr` giving the group elements.
/// * `elements_len` - Number of elements referenced by `elements_ptr`.
/// * `keys_ptr` - Pointer to an array of pointers to `Expr` for the group elements used as keys.
/// * `values_ptr` - Pointer to an array of pointers to `Expr` for the corresponding matrices/operators.
/// * `map_len` - Number of entries in the representation map.
///
/// # Returns
///
/// A pointer to a heap-allocated [`Representation`] constructed from the supplied
/// data.
///
/// # Safety
///
/// This function is unsafe because it dereferences multiple raw pointers and
/// assumes they form consistent arrays of valid `Expr` objects. The returned
/// `Representation` must be freed with [`rssn_representation_free`].

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_representation_create(
    elements_ptr: *const *const Expr,
    elements_len: usize,
    keys_ptr: *const *const Expr,
    values_ptr: *const *const Expr,
    map_len: usize,
) -> *mut Representation {

    unsafe {

        let elements_slice =
            std::slice::from_raw_parts(
                elements_ptr,
                elements_len,
            );

        let elements: Vec<
            GroupElement,
        > = elements_slice
            .iter()
            .map(|&p| {

                GroupElement(
                    (*p).clone(),
                )
            })
            .collect();

        let keys_slice =
            std::slice::from_raw_parts(
                keys_ptr,
                map_len,
            );

        let values_slice =
            std::slice::from_raw_parts(
                values_ptr,
                map_len,
            );

        let mut matrices =
            HashMap::new();

        for i in 0 .. map_len {

            let key = GroupElement(
                (*keys_slice[i])
                    .clone(),
            );

            let val = (*values_slice
                [i])
                .clone();

            matrices.insert(key, val);
        }

        let rep = Representation::new(
            elements,
            matrices,
        );

        Box::into_raw(Box::new(rep))
    }
}

#[unsafe(no_mangle)]

/// Frees a representation previously created by [`rssn_representation_create`].
///
/// # Arguments
///
/// * `ptr` - Pointer to a heap-allocated [`Representation`] returned by `rssn_representation_create`.
///
/// # Returns
///
/// This function does not return a value.
///
/// # Safety
///
/// This function is unsafe because it takes ownership of a raw pointer. The pointer
/// must either be null or have been allocated by `rssn_representation_create`, and
/// must not be used after this call.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_representation_free(
    ptr: *mut Representation
) {

    unsafe {

        if !ptr.is_null() {

            let _ = Box::from_raw(ptr);
        }
    }
}

#[unsafe(no_mangle)]

/// Checks whether a representation is valid for a given group using handle-based APIs.
///
/// A representation is valid if it defines a group homomorphism from the group to
/// the general linear group on the representation space.
///
/// # Arguments
///
/// * `rep` - Pointer to a [`Representation`] value.
/// * `group` - Pointer to a [`Group`] value.
///
/// # Returns
///
/// `true` if the representation is valid for the group, `false` otherwise.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers; the caller must
/// ensure they point to valid objects.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_representation_is_valid(
    rep: *const Representation,
    group: *const Group,
) -> bool {

    unsafe {

        (*rep).is_valid(&*group)
    }
}

#[unsafe(no_mangle)]

/// Computes the character of a representation and returns its values via raw output buffers.
///
/// The character is the trace of each representation matrix and is returned as a
/// mapping from group elements to scalar values.
///
/// # Arguments
///
/// * `rep` - Pointer to a [`Representation`] value.
/// * `out_len` - Pointer to a `usize` that will be filled with the number of character entries.
/// * `out_keys` - Output pointer that will receive a pointer to an array of `*mut Expr`
///   representing the group elements.
/// * `out_values` - Output pointer that will receive a pointer to an array of `*mut Expr`
///   representing the corresponding character values.
///
/// # Returns
///
/// This function does not return a direct value; results are written through the
/// output pointers.
///
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers and returns
/// ownership of heap-allocated arrays of `Expr`. The caller must ensure all input
/// pointers are valid and is responsible for freeing the returned memory.

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_character(
    rep: *const Representation,
    out_len: *mut usize,
    out_keys: *mut *mut *mut Expr,
    out_values: *mut *mut *mut Expr,
) {

    unsafe {

        let chars = character(&*rep);

        *out_len = chars.len();

        let mut keys_vec =
            Vec::with_capacity(
                chars.len(),
            );

        let mut values_vec =
            Vec::with_capacity(
                chars.len(),
            );

        for (k, v) in chars {

            keys_vec.push(
                Box::into_raw(
                    Box::new(k.0),
                ),
            );

            values_vec.push(
                Box::into_raw(
                    Box::new(v),
                ),
            );
        }

        *out_keys =
            keys_vec.as_mut_ptr();

        *out_values =
            values_vec.as_mut_ptr();

        std::mem::forget(keys_vec);

        std::mem::forget(values_vec);
    }
}
