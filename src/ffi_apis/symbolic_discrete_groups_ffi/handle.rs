use crate::symbolic::discrete_groups::cyclic_group;
use crate::symbolic::discrete_groups::dihedral_group;
use crate::symbolic::discrete_groups::klein_four_group;
use crate::symbolic::discrete_groups::symmetric_group;
use crate::symbolic::group_theory::Group;

/// Creates a cyclic group of order `n` and returns a raw pointer to it.
///
/// # Arguments
/// * `n` - The order of the cyclic group.
///
/// # Returns
/// A raw pointer (`*mut Group`) to the newly created group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_cyclic_group_create(
    n: usize
) -> *mut Group {

    let group = cyclic_group(n);

    Box::into_raw(Box::new(group))
}

/// Creates a dihedral group of order `2n` and returns a raw pointer to it.
///
/// # Arguments
/// * `n` - The parameter defining the dihedral group $`D_n`$.
///
/// # Returns
/// A raw pointer (`*mut Group`) to the newly created group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dihedral_group_create(
    n: usize
) -> *mut Group {

    let group = dihedral_group(n);

    Box::into_raw(Box::new(group))
}

/// Creates a symmetric group of degree `n` and returns a raw pointer to it.
///
/// # Arguments
/// * `n` - The number of symbols the group acts on.
///
/// # Returns
/// A raw pointer (`*mut Group`) to the newly created group, or NULL if `n` is invalid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_symmetric_group_create(
    n: usize
) -> *mut Group {

    match symmetric_group(n) {
        | Ok(group) => {
            Box::into_raw(Box::new(
                group,
            ))
        },
        | Err(_) => {
            std::ptr::null_mut()
        },
    }
}

/// Creates a Klein four-group and returns a raw pointer to it.
///
/// # Returns
/// A raw pointer (`*mut Group`) to the newly created group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_klein_four_group_create()
-> *mut Group {

    let group = klein_four_group();

    Box::into_raw(Box::new(group))
}
