use crate::symbolic::discrete_groups::{cyclic_group, dihedral_group, symmetric_group, klein_four_group};
use crate::symbolic::group_theory::Group;

/// Creates a cyclic group of order `n` and returns a raw pointer to it.
///
/// # Arguments
/// * `n` - The order of the cyclic group.
///
/// # Returns
/// A raw pointer (`*mut Group`) to the newly created group.
#[no_mangle]

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
#[no_mangle]

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
#[no_mangle]

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
#[no_mangle]

pub unsafe extern "C" fn rssn_klein_four_group_create(
) -> *mut Group {

    let group = klein_four_group();

    Box::into_raw(Box::new(group))
}
