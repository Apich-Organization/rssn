use crate::symbolic::discrete_groups::*;
use crate::symbolic::group_theory::Group;

#[no_mangle]
pub unsafe extern "C" fn rssn_cyclic_group_create(n: usize) -> *mut Group {
    let group = cyclic_group(n);
    Box::into_raw(Box::new(group))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_dihedral_group_create(n: usize) -> *mut Group {
    let group = dihedral_group(n);
    Box::into_raw(Box::new(group))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_symmetric_group_create(n: usize) -> *mut Group {
    match symmetric_group(n) {
        Ok(group) => Box::into_raw(Box::new(group)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_klein_four_group_create() -> *mut Group {
    let group = klein_four_group();
    Box::into_raw(Box::new(group))
}
