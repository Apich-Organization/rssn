use std::os::raw::c_char;

use crate::ffi_apis::common::to_json_string;
use crate::symbolic::discrete_groups::cyclic_group;
use crate::symbolic::discrete_groups::dihedral_group;
use crate::symbolic::discrete_groups::klein_four_group;
use crate::symbolic::discrete_groups::symmetric_group;

/// Creates a cyclic group of order `n` and returns it as a JSON string.
///
/// # Arguments
/// * `n` - The order of the cyclic group.
///
/// # Returns
/// A raw pointer to a JSON string representing the group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_cyclic_group_create(
    n: usize
) -> *mut c_char {

    let group = cyclic_group(n);

    to_json_string(&group)
}

/// Creates a dihedral group of order `2n` and returns it as a JSON string.
///
/// # Arguments
/// * `n` - The parameter defining the dihedral group $`D_n`$.
///
/// # Returns
/// A raw pointer to a JSON string representing the group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_dihedral_group_create(
    n: usize
) -> *mut c_char {

    let group = dihedral_group(n);

    to_json_string(&group)
}

/// Creates a symmetric group of degree `n` and returns it as a JSON string.
///
/// # Arguments
/// * `n` - The number of symbols the group acts on.
///
/// # Returns
/// A raw pointer to a JSON string representing the group, or NULL if `n` is invalid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_symmetric_group_create(
    n: usize
) -> *mut c_char {

    match symmetric_group(n) {
        | Ok(group) => {
            to_json_string(&group)
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Creates a Klein four-group and returns it as a JSON string.
///
/// # Returns
/// A raw pointer to a JSON string representing the group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_json_klein_four_group_create()
-> *mut c_char {

    let group = klein_four_group();

    to_json_string(&group)
}
