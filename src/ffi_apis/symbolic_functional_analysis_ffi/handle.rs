use crate::symbolic::core::Expr;
use crate::symbolic::functional_analysis::*;
use std::ffi::CStr;
use std::os::raw::c_char;

// --- HilbertSpace ---

#[no_mangle]
pub unsafe extern "C" fn rssn_hilbert_space_create(
    var: *const c_char,
    lower_bound: *const Expr,
    upper_bound: *const Expr,
) -> *mut HilbertSpace {
    let var_str = CStr::from_ptr(var).to_str().unwrap();
    let space = HilbertSpace::new(var_str, (*lower_bound).clone(), (*upper_bound).clone());
    Box::into_raw(Box::new(space))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_hilbert_space_free(ptr: *mut HilbertSpace) {
    if !ptr.is_null() {
        let _ = Box::from_raw(ptr);
    }
}

// --- BanachSpace ---

#[no_mangle]
pub unsafe extern "C" fn rssn_banach_space_create(
    var: *const c_char,
    lower_bound: *const Expr,
    upper_bound: *const Expr,
    p: *const Expr,
) -> *mut BanachSpace {
    let var_str = CStr::from_ptr(var).to_str().unwrap();
    let space = BanachSpace::new(
        var_str,
        (*lower_bound).clone(),
        (*upper_bound).clone(),
        (*p).clone(),
    );
    Box::into_raw(Box::new(space))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_banach_space_free(ptr: *mut BanachSpace) {
    if !ptr.is_null() {
        let _ = Box::from_raw(ptr);
    }
}

// --- LinearOperator ---

#[no_mangle]
pub unsafe extern "C" fn rssn_linear_operator_derivative_create(
    var: *const c_char,
) -> *mut LinearOperator {
    let var_str = CStr::from_ptr(var).to_str().unwrap();
    let op = LinearOperator::Derivative(var_str.to_string());
    Box::into_raw(Box::new(op))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_linear_operator_integral_create(
    lower_bound: *const Expr,
    var: *const c_char,
) -> *mut LinearOperator {
    let var_str = CStr::from_ptr(var).to_str().unwrap();
    let op = LinearOperator::Integral((*lower_bound).clone(), var_str.to_string());
    Box::into_raw(Box::new(op))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_linear_operator_apply(
    op: *const LinearOperator,
    expr: *const Expr,
) -> *mut Expr {
    let result = (*op).apply(&*expr);
    Box::into_raw(Box::new(result))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_linear_operator_free(ptr: *mut LinearOperator) {
    if !ptr.is_null() {
        let _ = Box::from_raw(ptr);
    }
}

// --- Functions ---

#[no_mangle]
pub unsafe extern "C" fn rssn_inner_product(
    space: *const HilbertSpace,
    f: *const Expr,
    g: *const Expr,
) -> *mut Expr {
    let result = inner_product(&*space, &*f, &*g);
    Box::into_raw(Box::new(result))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_norm(space: *const HilbertSpace, f: *const Expr) -> *mut Expr {
    let result = norm(&*space, &*f);
    Box::into_raw(Box::new(result))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_banach_norm(space: *const BanachSpace, f: *const Expr) -> *mut Expr {
    let result = banach_norm(&*space, &*f);
    Box::into_raw(Box::new(result))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_are_orthogonal(
    space: *const HilbertSpace,
    f: *const Expr,
    g: *const Expr,
) -> bool {
    are_orthogonal(&*space, &*f, &*g)
}

#[no_mangle]
pub unsafe extern "C" fn rssn_project(
    space: *const HilbertSpace,
    f: *const Expr,
    g: *const Expr,
) -> *mut Expr {
    let result = project(&*space, &*f, &*g);
    Box::into_raw(Box::new(result))
}

#[no_mangle]
pub unsafe extern "C" fn rssn_gram_schmidt(
    space: *const HilbertSpace,
    basis_ptr: *const *const Expr,
    basis_len: usize,
    out_len: *mut usize,
) -> *mut *mut Expr {
    let basis_slice = std::slice::from_raw_parts(basis_ptr, basis_len);
    let basis: Vec<Expr> = basis_slice.iter().map(|&p| (*p).clone()).collect();

    let orthogonal_basis = gram_schmidt(&*space, &basis);

    *out_len = orthogonal_basis.len();
    let mut out_ptrs = Vec::with_capacity(orthogonal_basis.len());
    for expr in orthogonal_basis {
        out_ptrs.push(Box::into_raw(Box::new(expr)));
    }

    let ptr = out_ptrs.as_mut_ptr();
    std::mem::forget(out_ptrs);
    ptr
}
