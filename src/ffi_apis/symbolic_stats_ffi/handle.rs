//! Handle-based FFI API for symbolic statistics functions.

use crate::symbolic::core::Expr;
use crate::symbolic::stats;

/// Computes the symbolic mean of a set of expressions.
///
/// # Safety
/// The caller must ensure `data` is a valid pointer to an array of `Expr` pointers of size `len`.
#[no_mangle]

pub unsafe extern "C" fn rssn_mean(
    data: *const *const Expr,
    len: usize,
) -> *mut Expr {

    if data.is_null() {

        return std::ptr::null_mut();
    }

    let mut exprs =
        Vec::with_capacity(len);

    for i in 0..len {

        let ptr = *data.add(i);

        if !ptr.is_null() {

            exprs.push((*ptr).clone());
        }
    }

    Box::into_raw(Box::new(
        stats::mean(&exprs),
    ))
}

/// Computes the symbolic variance of a set of expressions.
///
/// # Safety
/// The caller must ensure `data` is a valid pointer to an array of `Expr` pointers of size `len`.
#[no_mangle]

pub unsafe extern "C" fn rssn_variance(
    data: *const *const Expr,
    len: usize,
) -> *mut Expr {

    if data.is_null() {

        return std::ptr::null_mut();
    }

    let mut exprs =
        Vec::with_capacity(len);

    for i in 0..len {

        let ptr = *data.add(i);

        if !ptr.is_null() {

            exprs.push((*ptr).clone());
        }
    }

    Box::into_raw(Box::new(
        stats::variance(&exprs),
    ))
}

/// Computes the symbolic standard deviation of a set of expressions.
///
/// # Safety
/// The caller must ensure `data` is a valid pointer to an array of `Expr` pointers of size `len`.
#[no_mangle]

pub unsafe extern "C" fn rssn_std_dev(
    data: *const *const Expr,
    len: usize,
) -> *mut Expr {

    if data.is_null() {

        return std::ptr::null_mut();
    }

    let mut exprs =
        Vec::with_capacity(len);

    for i in 0..len {

        let ptr = *data.add(i);

        if !ptr.is_null() {

            exprs.push((*ptr).clone());
        }
    }

    Box::into_raw(Box::new(
        stats::std_dev(&exprs),
    ))
}

/// Computes the symbolic covariance of two sets of expressions.
///
/// # Safety
/// The caller must ensure `data1` and `data2` are valid pointers to arrays of `Expr` pointers.
#[no_mangle]

pub unsafe extern "C" fn rssn_covariance(
    data1: *const *const Expr,
    len1: usize,
    data2: *const *const Expr,
    len2: usize,
) -> *mut Expr {

    if data1.is_null()
        || data2.is_null()
    {

        return std::ptr::null_mut();
    }

    let mut exprs1 =
        Vec::with_capacity(len1);

    for i in 0..len1 {

        let ptr = *data1.add(i);

        if !ptr.is_null() {

            exprs1.push((*ptr).clone());
        }
    }

    let mut exprs2 =
        Vec::with_capacity(len2);

    for i in 0..len2 {

        let ptr = *data2.add(i);

        if !ptr.is_null() {

            exprs2.push((*ptr).clone());
        }
    }

    Box::into_raw(Box::new(
        stats::covariance(
            &exprs1, &exprs2,
        ),
    ))
}

/// Computes the symbolic Pearson correlation coefficient.
///
/// # Safety
/// The caller must ensure `data1` and `data2` are valid pointers to arrays of `Expr` pointers.
#[no_mangle]

pub unsafe extern "C" fn rssn_correlation(
    data1: *const *const Expr,
    len1: usize,
    data2: *const *const Expr,
    len2: usize,
) -> *mut Expr {

    if data1.is_null()
        || data2.is_null()
    {

        return std::ptr::null_mut();
    }

    let mut exprs1 =
        Vec::with_capacity(len1);

    for i in 0..len1 {

        let ptr = *data1.add(i);

        if !ptr.is_null() {

            exprs1.push((*ptr).clone());
        }
    }

    let mut exprs2 =
        Vec::with_capacity(len2);

    for i in 0..len2 {

        let ptr = *data2.add(i);

        if !ptr.is_null() {

            exprs2.push((*ptr).clone());
        }
    }

    Box::into_raw(Box::new(
        stats::correlation(
            &exprs1, &exprs2,
        ),
    ))
}
