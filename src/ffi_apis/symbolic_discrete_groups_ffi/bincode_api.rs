use crate::ffi_apis::common::*;
use crate::symbolic::discrete_groups::*;

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_cyclic_group_create(n : usize) -> BincodeBuffer {

    let group = cyclic_group(n);

    to_bincode_buffer(&group)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_dihedral_group_create(n : usize) -> BincodeBuffer {

    let group = dihedral_group(n);

    to_bincode_buffer(&group)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_symmetric_group_create(n : usize) -> BincodeBuffer {

    match symmetric_group(n) {
        | Ok(group) => to_bincode_buffer(&group),
        | Err(_) => BincodeBuffer::empty(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_klein_four_group_create() -> BincodeBuffer {

    let group = klein_four_group();

    to_bincode_buffer(&group)
}
