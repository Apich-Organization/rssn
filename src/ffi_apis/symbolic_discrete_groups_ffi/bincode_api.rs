use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::discrete_groups::cyclic_group;
use crate::symbolic::discrete_groups::dihedral_group;
use crate::symbolic::discrete_groups::klein_four_group;
use crate::symbolic::discrete_groups::symmetric_group;

/// Creates a cyclic group of order `n` and returns it as a Bincode buffer.
///
/// # Arguments
/// * `n` - The order of the cyclic group.
///
/// # Returns
/// A `BincodeBuffer` containing the serialized representation of the group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_cyclic_group_create(
    n: usize
) -> BincodeBuffer {

    let group = cyclic_group(n);

    to_bincode_buffer(&group)
}

/// Creates a dihedral group of order `2n` and returns it as a Bincode buffer.
///
/// # Arguments
/// * `n` - The parameter defining the dihedral group $`D_n`$.
///
/// # Returns
/// A `BincodeBuffer` containing the serialized representation of the group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_dihedral_group_create(
    n: usize
) -> BincodeBuffer {

    let group = dihedral_group(n);

    to_bincode_buffer(&group)
}

/// Creates a symmetric group of degree `n` and returns it as a Bincode buffer.
///
/// # Arguments
/// * `n` - The number of symbols the group acts on.
///
/// # Returns
/// A `BincodeBuffer` containing the serialized representation of the group, or an empty buffer if `n` is invalid.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_symmetric_group_create(
    n: usize
) -> BincodeBuffer {

    match symmetric_group(n) {
        | Ok(group) => {
            to_bincode_buffer(&group)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

/// Creates a Klein four-group and returns it as a Bincode buffer.
///
/// # Returns
/// A `BincodeBuffer` containing the serialized representation of the group.
#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_bincode_klein_four_group_create(
) -> BincodeBuffer {

    let group = klein_four_group();

    to_bincode_buffer(&group)
}
