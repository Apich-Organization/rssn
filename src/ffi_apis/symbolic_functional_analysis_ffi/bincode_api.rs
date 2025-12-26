use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::functional_analysis::*;

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_hilbert_space_create(
    buf: BincodeBuffer
) -> BincodeBuffer {

    let space : HilbertSpace = match from_bincode_buffer(&buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    to_bincode_buffer(&space)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_inner_product(
    space_buf: BincodeBuffer,
    f_buf: BincodeBuffer,
    g_buf: BincodeBuffer,
) -> BincodeBuffer {

    let space : HilbertSpace = match from_bincode_buffer(&space_buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    let f : Expr = match from_bincode_buffer(&f_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let g : Expr = match from_bincode_buffer(&g_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result =
        inner_product(&space, &f, &g);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_norm(
    space_buf: BincodeBuffer,
    f_buf: BincodeBuffer,
) -> BincodeBuffer {

    let space : HilbertSpace = match from_bincode_buffer(&space_buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    let f : Expr = match from_bincode_buffer(&f_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let result = norm(&space, &f);

    to_bincode_buffer(&result)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_gram_schmidt(
    space_buf: BincodeBuffer,
    basis_buf: BincodeBuffer,
) -> BincodeBuffer {

    let space : HilbertSpace = match from_bincode_buffer(&space_buf) {
        | Some(s) => s,
        | None => return BincodeBuffer::empty(),
    };

    let basis : Vec<Expr> = match from_bincode_buffer(&basis_buf) {
        | Some(b) => b,
        | None => return BincodeBuffer::empty(),
    };

    let result =
        gram_schmidt(&space, &basis);

    to_bincode_buffer(&result)
}
