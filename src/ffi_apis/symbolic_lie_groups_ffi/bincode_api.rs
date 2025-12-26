use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::lie_groups_and_algebras::*;

// --- LieAlgebra Creation ---

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_lie_algebra_so3(
) -> BincodeBuffer {

    let algebra = so3();

    to_bincode_buffer(&algebra)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_lie_algebra_su2(
) -> BincodeBuffer {

    let algebra = su2();

    to_bincode_buffer(&algebra)
}

// --- Lie Bracket ---

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_lie_bracket(
    x_buf: BincodeBuffer,
    y_buf: BincodeBuffer,
) -> BincodeBuffer {

    let x : Expr = match from_bincode_buffer(&x_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let y : Expr = match from_bincode_buffer(&y_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    match lie_bracket(&x, &y) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

// --- Exponential Map ---

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_exponential_map(
    x_buf: BincodeBuffer,
    order: usize,
) -> BincodeBuffer {

    let x : Expr = match from_bincode_buffer(&x_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    match exponential_map(&x, order) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

// --- Adjoint Representations ---

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_adjoint_representation_group(
    g_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let g : Expr = match from_bincode_buffer(&g_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let x : Expr = match from_bincode_buffer(&x_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    match adjoint_representation_group(
        &g, &x,
    ) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_adjoint_representation_algebra(
    x_buf: BincodeBuffer,
    y_buf: BincodeBuffer,
) -> BincodeBuffer {

    let x : Expr = match from_bincode_buffer(&x_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    let y : Expr = match from_bincode_buffer(&y_buf) {
        | Some(e) => e,
        | None => return BincodeBuffer::empty(),
    };

    match adjoint_representation_algebra(
        &x, &y,
    ) {
        | Ok(result) => {
            to_bincode_buffer(&result)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

// --- Commutator Table ---

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_commutator_table(
    algebra_buf: BincodeBuffer
) -> BincodeBuffer {

    let algebra : LieAlgebra = match from_bincode_buffer(&algebra_buf) {
        | Some(a) => a,
        | None => return BincodeBuffer::empty(),
    };

    match commutator_table(&algebra) {
        | Ok(table) => {
            to_bincode_buffer(&table)
        },
        | Err(_) => {
            BincodeBuffer::empty()
        },
    }
}

// --- Jacobi Identity Check ---

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_check_jacobi_identity(
    algebra_buf: BincodeBuffer
) -> bool {

    let algebra: LieAlgebra =
        match from_bincode_buffer(
            &algebra_buf,
        ) {
            | Some(a) => a,
            | None => return false,
        };

    match check_jacobi_identity(
        &algebra,
    ) {
        | Ok(result) => result,
        | Err(_) => false,
    }
}

// --- Generators ---

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_so3_generators(
) -> BincodeBuffer {

    let generators = so3_generators();

    let exprs: Vec<Expr> = generators
        .into_iter()
        .map(|g| g.0)
        .collect();

    to_bincode_buffer(&exprs)
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bincode_su2_generators(
) -> BincodeBuffer {

    let generators = su2_generators();

    let exprs: Vec<Expr> = generators
        .into_iter()
        .map(|g| g.0)
        .collect();

    to_bincode_buffer(&exprs)
}
