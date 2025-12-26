use crate::symbolic::combinatorics::*;
use crate::symbolic::core::Expr;

#[no_mangle]

pub unsafe extern "C" fn rssn_permutations(
    n: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let result = permutations((*n).clone(), (*k).clone());

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_combinations(
    n: *const Expr,
    k: *const Expr,
) -> *mut Expr {

    let result = combinations(&(*n), (*k).clone());

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_catalan_number(n: usize) -> *mut Expr {

    let result = catalan_number(n);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_stirling_number_second_kind(
    n: usize,
    k: usize,
) -> *mut Expr {

    let result = stirling_number_second_kind(n, k);

    Box::into_raw(Box::new(result))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_bell_number(n: usize) -> *mut Expr {

    let result = bell_number(n);

    Box::into_raw(Box::new(result))
}
