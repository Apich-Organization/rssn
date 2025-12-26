use crate::ffi_apis::common::*;
use crate::symbolic::group_theory::*;

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_group_create(
    buf : BincodeBuffer
) -> BincodeBuffer {

    let group: Group = match from_bincode_buffer(&buf) {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&group)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_group_multiply(
    group_buf : BincodeBuffer,
    a_buf : BincodeBuffer,
    b_buf : BincodeBuffer,
) -> BincodeBuffer {

    let group: Group = match from_bincode_buffer(&group_buf)
    {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let a: GroupElement = match from_bincode_buffer(&a_buf)
    {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let b: GroupElement = match from_bincode_buffer(&b_buf)
    {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = group.multiply(&a, &b);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_group_inverse(
    group_buf : BincodeBuffer,
    a_buf : BincodeBuffer,
) -> BincodeBuffer {

    let group: Group = match from_bincode_buffer(&group_buf)
    {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let a: GroupElement = match from_bincode_buffer(&a_buf)
    {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = group.inverse(&a);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_group_is_abelian(
    group_buf : BincodeBuffer
) -> bool {

    let group : Group =
        match from_bincode_buffer(
            &group_buf,
        ) {
            | Some(g) => g,
            | None => return false,
        };

    group.is_abelian()
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_group_element_order(
    group_buf : BincodeBuffer,
    a_buf : BincodeBuffer,
) -> usize {

    let group : Group =
        match from_bincode_buffer(
            &group_buf,
        ) {
            | Some(g) => g,
            | None => return 0,
        };

    let a : GroupElement =
        match from_bincode_buffer(
            &a_buf,
        ) {
            | Some(e) => e,
            | None => return 0,
        };

    group
        .element_order(&a)
        .unwrap_or(0)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_group_conjugacy_classes(
    group_buf : BincodeBuffer
) -> BincodeBuffer {

    let group: Group = match from_bincode_buffer(&group_buf)
    {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let classes =
        group.conjugacy_classes();

    to_bincode_buffer(&classes)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_group_center(
    group_buf : BincodeBuffer
) -> BincodeBuffer {

    let group: Group = match from_bincode_buffer(&group_buf)
    {
        | Some(g) => g,
        | None => return BincodeBuffer::empty(),
    };

    let center = group.center();

    to_bincode_buffer(&center)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_representation_create(
    buf : BincodeBuffer
) -> BincodeBuffer {

    let rep: Representation =
        match from_bincode_buffer(&buf) {
            | Some(r) => r,
            | None => return BincodeBuffer::empty(),
        };

    to_bincode_buffer(&rep)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_representation_is_valid(
    rep_buf : BincodeBuffer,
    group_buf : BincodeBuffer,
) -> bool {

    let rep : Representation =
        match from_bincode_buffer(
            &rep_buf,
        ) {
            | Some(r) => r,
            | None => return false,
        };

    let group : Group =
        match from_bincode_buffer(
            &group_buf,
        ) {
            | Some(g) => g,
            | None => return false,
        };

    rep.is_valid(&group)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_character(
    rep_buf : BincodeBuffer
) -> BincodeBuffer {

    let rep: Representation =
        match from_bincode_buffer(&rep_buf) {
            | Some(r) => r,
            | None => return BincodeBuffer::empty(),
        };

    let chars = character(&rep);

    to_bincode_buffer(&chars)
}
