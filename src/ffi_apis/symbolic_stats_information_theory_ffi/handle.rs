use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_information_theory;
use std::os::raw::c_char;

unsafe fn collect_exprs(data: *const *const Expr, len: usize) -> Vec<Expr> {

    let mut exprs = Vec::with_capacity(len);

    for i in 0..len {

        let ptr = *data.add(i);

        if !ptr.is_null() {

            exprs.push((*ptr).clone());
        }
    }

    exprs
}

#[no_mangle]

pub unsafe extern "C" fn rssn_shannon_entropy(probs: *const *const Expr, len: usize) -> *mut Expr {

    if probs.is_null() {

        return std::ptr::null_mut();
    }

    let p = collect_exprs(probs, len);

    Box::into_raw(Box::new(stats_information_theory::shannon_entropy(&p)))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_kl_divergence(
    p_probs: *const *const Expr,
    p_len: usize,
    q_probs: *const *const Expr,
    q_len: usize,
) -> *mut Expr {

    if p_probs.is_null() || q_probs.is_null() {

        return std::ptr::null_mut();
    }

    let p = collect_exprs(p_probs, p_len);

    let q = collect_exprs(q_probs, q_len);

    match stats_information_theory::kl_divergence(&p, &q) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_cross_entropy(
    p_probs: *const *const Expr,
    p_len: usize,
    q_probs: *const *const Expr,
    q_len: usize,
) -> *mut Expr {

    if p_probs.is_null() || q_probs.is_null() {

        return std::ptr::null_mut();
    }

    let p = collect_exprs(p_probs, p_len);

    let q = collect_exprs(q_probs, q_len);

    match stats_information_theory::cross_entropy(&p, &q) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_gini_impurity(probs: *const *const Expr, len: usize) -> *mut Expr {

    if probs.is_null() {

        return std::ptr::null_mut();
    }

    let p = collect_exprs(probs, len);

    Box::into_raw(Box::new(stats_information_theory::gini_impurity(&p)))
}

#[no_mangle]

pub unsafe extern "C" fn rssn_joint_entropy(joint_probs: *const Expr) -> *mut Expr {

    if joint_probs.is_null() {

        return std::ptr::null_mut();
    }

    match stats_information_theory::joint_entropy(&*joint_probs) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_conditional_entropy(joint_probs: *const Expr) -> *mut Expr {

    if joint_probs.is_null() {

        return std::ptr::null_mut();
    }

    match stats_information_theory::conditional_entropy(&*joint_probs) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(_) => std::ptr::null_mut(),
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_mutual_information(joint_probs: *const Expr) -> *mut Expr {

    if joint_probs.is_null() {

        return std::ptr::null_mut();
    }

    match stats_information_theory::mutual_information(&*joint_probs) {
        Ok(res) => Box::into_raw(Box::new(res)),
        Err(_) => std::ptr::null_mut(),
    }
}
