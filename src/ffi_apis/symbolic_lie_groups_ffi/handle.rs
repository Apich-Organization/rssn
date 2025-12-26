use crate::symbolic::core::Expr;
use crate::symbolic::lie_groups_and_algebras::*;
use std::ffi::CStr;
use std::os::raw::c_char;

// --- LieAlgebra ---

#[no_mangle]

pub unsafe extern "C" fn rssn_lie_algebra_so3_create() -> *mut LieAlgebra {

    let algebra = so3();

    Box::into_raw(Box::new(algebra))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_lie_algebra_su2_create() -> *mut LieAlgebra {

    let algebra = su2();

    Box::into_raw(Box::new(algebra))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_lie_algebra_free(ptr: *mut LieAlgebra) {

    if !ptr.is_null() {

        let _ = Box::from_raw(ptr);
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_lie_algebra_get_dimension(ptr: *const LieAlgebra) -> usize {

    (*ptr).dimension
}

#[no_mangle]

pub unsafe extern "C" fn rssn_lie_algebra_get_name(ptr: *const LieAlgebra) -> *mut c_char {

    let name = &(*ptr).name;

    std::ffi::CString::new(name.as_str())
        .unwrap()
        .into_raw()
}

#[no_mangle]

pub unsafe extern "C" fn rssn_lie_algebra_get_basis_element(
    ptr: *const LieAlgebra,
    index: usize,
) -> *mut Expr {

    let algebra = &*ptr;

    if index >= algebra.basis.len() {

        return std::ptr::null_mut();
    }

    Box::into_raw(Box::new(
        algebra.basis[index]
            .0
            .clone(),
    ))
}

// --- Lie Bracket ---

#[no_mangle]

pub unsafe extern "C" fn rssn_lie_bracket(
    x: *const Expr,
    y: *const Expr,
) -> *mut Expr {

    match lie_bracket(&*x, &*y) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

// --- Exponential Map ---

#[no_mangle]

pub unsafe extern "C" fn rssn_exponential_map(
    x: *const Expr,
    order: usize,
) -> *mut Expr {

    match exponential_map(&*x, order) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

// --- Adjoint Representations ---

#[no_mangle]

pub unsafe extern "C" fn rssn_adjoint_representation_group(
    g: *const Expr,
    x: *const Expr,
) -> *mut Expr {

    match adjoint_representation_group(&*g, &*x) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_adjoint_representation_algebra(
    x: *const Expr,
    y: *const Expr,
) -> *mut Expr {

    match adjoint_representation_algebra(&*x, &*y) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}

// --- Commutator Table ---

#[no_mangle]

pub unsafe extern "C" fn rssn_commutator_table(
    algebra: *const LieAlgebra,
    out_rows: *mut usize,
    out_cols: *mut usize,
) -> *mut *mut Expr {

    match commutator_table(&*algebra) {
        Ok(table) => {

            let rows = table.len();

            let cols = if rows > 0 { table[0].len() } else { 0 };

            *out_rows = rows;

            *out_cols = cols;

            let mut flat_ptrs = Vec::with_capacity(rows * cols);

            for row in table {

                for elem in row {

                    flat_ptrs.push(Box::into_raw(Box::new(elem)));
                }
            }

            let ptr = flat_ptrs.as_mut_ptr();

            std::mem::forget(flat_ptrs);

            ptr
        }
        Err(_) => {

            *out_rows = 0;

            *out_cols = 0;

            std::ptr::null_mut()
        }
    }
}

// --- Jacobi Identity Check ---

#[no_mangle]

pub unsafe extern "C" fn rssn_check_jacobi_identity(algebra: *const LieAlgebra) -> bool {

    match check_jacobi_identity(&*algebra) {
        Ok(result) => result,
        Err(_) => false,
    }
}

// --- Generators ---

#[no_mangle]

pub unsafe extern "C" fn rssn_so3_generators(out_len: *mut usize) -> *mut *mut Expr {

    let generators = so3_generators();

    *out_len = generators.len();

    let mut ptrs = Vec::with_capacity(generators.len());

    for gen in generators {

        ptrs.push(Box::into_raw(Box::new(gen.0)));
    }

    let ptr = ptrs.as_mut_ptr();

    std::mem::forget(ptrs);

    ptr
}

#[no_mangle]

pub unsafe extern "C" fn rssn_su2_generators(out_len: *mut usize) -> *mut *mut Expr {

    let generators = su2_generators();

    *out_len = generators.len();

    let mut ptrs = Vec::with_capacity(generators.len());

    for gen in generators {

        ptrs.push(Box::into_raw(Box::new(gen.0)));
    }

    let ptr = ptrs.as_mut_ptr();

    std::mem::forget(ptrs);

    ptr
}
