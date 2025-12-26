use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::discrete_groups::*;

#[no_mangle]

pub unsafe extern "C" fn rssn_json_cyclic_group_create(n : usize) -> *mut c_char {

    let group = cyclic_group(n);

    to_json_string(&group)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dihedral_group_create(n : usize) -> *mut c_char {

    let group = dihedral_group(n);

    to_json_string(&group)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_symmetric_group_create(n : usize) -> *mut c_char {

    match symmetric_group(n) {
        | Ok(group) => to_json_string(&group),
        | Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_json_klein_four_group_create() -> *mut c_char {

    let group = klein_four_group();

    to_json_string(&group)
}
