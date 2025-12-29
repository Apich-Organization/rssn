use crate::ffi_apis::common::*;
use crate::symbolic::discrete_groups::*;

/// Creates a cyclic group of order `n` and returns it as a Bincode buffer.
///
/// # Arguments
/// * `n` - The order of the cyclic group.
///
/// # Returns
/// A `BincodeBuffer` containing the serialized representation of the group.
#[no_mangle]
pub unsafe extern "C" fn rssn_bincode_cyclic_group_create(
    n: usize
) -> BincodeBuffer {

    let group = cyclic_group(n);

    to_bincode_buffer(&group)
}

/// Creates a dihedral group of order `2n` and returns it as a Bincode buffer.
///
/// # Arguments
/// * `n` - The parameter defining the dihedral group $D_n$.
///
/// # Returns
/// A `BincodeBuffer` containing the serialized representation of the group.
#[no_mangle]
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
#[no_mangle]
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
#[no_mangle]
pub unsafe extern "C" fn rssn_bincode_klein_four_group_create(
) -> BincodeBuffer {

    let group = klein_four_group();

    to_bincode_buffer(&group)
}
