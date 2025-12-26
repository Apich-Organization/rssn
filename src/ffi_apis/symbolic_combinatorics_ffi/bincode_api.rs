use crate::ffi_apis::common::*;
use crate::symbolic::combinatorics::*;
use crate::symbolic::core::Expr;

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_permutations(
    n_buf : BincodeBuffer,
    k_buf : BincodeBuffer,
) -> BincodeBuffer {

    let n : Expr = match from_bincode_buffer(&n_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let k : Expr = match from_bincode_buffer(&k_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = permutations(n, k);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_combinations(
    n_buf : BincodeBuffer,
    k_buf : BincodeBuffer,
) -> BincodeBuffer {

    let n : Expr = match from_bincode_buffer(&n_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let k : Expr = match from_bincode_buffer(&k_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = combinations(&n, k);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_catalan_number(n : usize) -> BincodeBuffer {

    let result = catalan_number(n);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_stirling_number_second_kind(
    n : usize,
    k : usize,
) -> BincodeBuffer {

    let result = stirling_number_second_kind(n, k);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_bell_number(n : usize) -> BincodeBuffer {

    let result = bell_number(n);

    to_bincode_buffer(&result)
}
