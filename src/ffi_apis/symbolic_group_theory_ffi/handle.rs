use crate::symbolic::core::Expr;
use crate::symbolic::group_theory::*;
use std::collections::HashMap;

// --- Group ---

#[no_mangle]

pub unsafe extern "C" fn rssn_group_create(
    elements_ptr: *const *const Expr,
    elements_len: usize,
    keys_a_ptr: *const *const Expr,
    keys_b_ptr: *const *const Expr,
    values_ptr: *const *const Expr,
    table_len: usize,
    identity_ptr: *const Expr,
) -> *mut Group {

    let elements_slice =
        std::slice::from_raw_parts(
            elements_ptr,
            elements_len,
        );

    let elements: Vec<GroupElement> =
        elements_slice
            .iter()
            .map(|&p| {
                GroupElement(
                    (*p).clone(),
                )
            })
            .collect();

    let keys_a_slice =
        std::slice::from_raw_parts(
            keys_a_ptr, table_len,
        );

    let keys_b_slice =
        std::slice::from_raw_parts(
            keys_b_ptr, table_len,
        );

    let values_slice =
        std::slice::from_raw_parts(
            values_ptr, table_len,
        );

    let mut multiplication_table =
        HashMap::new();

    for i in 0..table_len {

        let a = GroupElement(
            (*keys_a_slice[i]).clone(),
        );

        let b = GroupElement(
            (*keys_b_slice[i]).clone(),
        );

        let val = GroupElement(
            (*values_slice[i]).clone(),
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

#[no_mangle]

pub unsafe extern "C" fn rssn_group_free(
    ptr: *mut Group
) {

    if !ptr.is_null() {

        let _ = Box::from_raw(ptr);
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_group_multiply(
    group: *const Group,
    a: *const Expr,
    b: *const Expr,
) -> *mut Expr {

    let ga = GroupElement((*a).clone());

    let gb = GroupElement((*b).clone());

    match (*group).multiply(&ga, &gb) {
        | Some(result) => {
            Box::into_raw(Box::new(
                result.0,
            ))
        },
        | None => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_group_inverse(
    group: *const Group,
    a: *const Expr,
) -> *mut Expr {

    let ga = GroupElement((*a).clone());

    match (*group).inverse(&ga) {
        | Some(result) => {
            Box::into_raw(Box::new(
                result.0,
            ))
        },
        | None => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_group_is_abelian(
    group: *const Group
) -> bool {

    (*group).is_abelian()
}

#[no_mangle]

pub unsafe extern "C" fn rssn_group_element_order(
    group: *const Group,
    a: *const Expr,
) -> usize {

    let ga = GroupElement((*a).clone());

    (*group)
        .element_order(&ga)
        .unwrap_or(0)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_group_center(
    group: *const Group,
    out_len: *mut usize,
) -> *mut *mut Expr {

    let center = (*group).center();

    *out_len = center.len();

    let mut out_ptrs =
        Vec::with_capacity(
            center.len(),
        );

    for elem in center {

        out_ptrs.push(Box::into_raw(
            Box::new(elem.0),
        ));
    }

    let ptr = out_ptrs.as_mut_ptr();

    std::mem::forget(out_ptrs);

    ptr
}

// --- Representation ---

#[no_mangle]

pub unsafe extern "C" fn rssn_representation_create(
    elements_ptr: *const *const Expr,
    elements_len: usize,
    keys_ptr: *const *const Expr,
    values_ptr: *const *const Expr,
    map_len: usize,
) -> *mut Representation {

    let elements_slice =
        std::slice::from_raw_parts(
            elements_ptr,
            elements_len,
        );

    let elements: Vec<GroupElement> =
        elements_slice
            .iter()
            .map(|&p| {
                GroupElement(
                    (*p).clone(),
                )
            })
            .collect();

    let keys_slice =
        std::slice::from_raw_parts(
            keys_ptr, map_len,
        );

    let values_slice =
        std::slice::from_raw_parts(
            values_ptr, map_len,
        );

    let mut matrices = HashMap::new();

    for i in 0..map_len {

        let key = GroupElement(
            (*keys_slice[i]).clone(),
        );

        let val =
            (*values_slice[i]).clone();

        matrices.insert(key, val);
    }

    let rep = Representation::new(
        elements, matrices,
    );

    Box::into_raw(Box::new(rep))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_representation_free(
    ptr: *mut Representation
) {

    if !ptr.is_null() {

        let _ = Box::from_raw(ptr);
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_representation_is_valid(
    rep: *const Representation,
    group: *const Group,
) -> bool {

    (*rep).is_valid(&*group)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_character(
    rep: *const Representation,
    out_len: *mut usize,
    out_keys: *mut *mut *mut Expr,
    out_values: *mut *mut *mut Expr,
) {

    let chars = character(&*rep);

    *out_len = chars.len();

    let mut keys_vec =
        Vec::with_capacity(chars.len());

    let mut values_vec =
        Vec::with_capacity(chars.len());

    for (k, v) in chars {

        keys_vec.push(Box::into_raw(
            Box::new(k.0),
        ));

        values_vec.push(Box::into_raw(
            Box::new(v),
        ));
    }

    *out_keys = keys_vec.as_mut_ptr();

    *out_values =
        values_vec.as_mut_ptr();

    std::mem::forget(keys_vec);

    std::mem::forget(values_vec);
}
