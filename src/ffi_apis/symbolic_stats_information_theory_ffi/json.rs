use crate::symbolic::core::Expr;
use crate::symbolic::stats_information_theory;
use crate::ffi_apis::common::*;
use std::os::raw::c_char;

#[no_mangle]
pub unsafe extern "C" fn rssn_json_shannon_entropy(probs_json: *const c_char) -> *mut c_char {
    let probs: Option<Vec<Expr>> = from_json_string(probs_json);
    if let Some(p) = probs {
        let res = stats_information_theory::shannon_entropy(&p);
        to_json_string(&res)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_kl_divergence(
    p_probs_json: *const c_char,
    q_probs_json: *const c_char
) -> *mut c_char {
    let p_probs: Option<Vec<Expr>> = from_json_string(p_probs_json);
    let q_probs: Option<Vec<Expr>> = from_json_string(q_probs_json);
    
    if let (Some(p), Some(q)) = (p_probs, q_probs) {
        match stats_information_theory::kl_divergence(&p, &q) {
            Ok(res) => to_json_string(&res),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_cross_entropy(
    p_probs_json: *const c_char,
    q_probs_json: *const c_char
) -> *mut c_char {
    let p_probs: Option<Vec<Expr>> = from_json_string(p_probs_json);
    let q_probs: Option<Vec<Expr>> = from_json_string(q_probs_json);
    
    if let (Some(p), Some(q)) = (p_probs, q_probs) {
        match stats_information_theory::cross_entropy(&p, &q) {
            Ok(res) => to_json_string(&res),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_gini_impurity(probs_json: *const c_char) -> *mut c_char {
    let probs: Option<Vec<Expr>> = from_json_string(probs_json);
    if let Some(p) = probs {
        let res = stats_information_theory::gini_impurity(&p);
        to_json_string(&res)
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_joint_entropy(joint_probs_json: *const c_char) -> *mut c_char {
    let joint: Option<Expr> = from_json_string(joint_probs_json);
    if let Some(j) = joint {
        match stats_information_theory::joint_entropy(&j) {
            Ok(res) => to_json_string(&res),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_conditional_entropy(joint_probs_json: *const c_char) -> *mut c_char {
    let joint: Option<Expr> = from_json_string(joint_probs_json);
    if let Some(j) = joint {
        match stats_information_theory::conditional_entropy(&j) {
            Ok(res) => to_json_string(&res),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}

#[no_mangle]
pub unsafe extern "C" fn rssn_json_mutual_information(joint_probs_json: *const c_char) -> *mut c_char {
    let joint: Option<Expr> = from_json_string(joint_probs_json);
    if let Some(j) = joint {
        match stats_information_theory::mutual_information(&j) {
            Ok(res) => to_json_string(&res),
            Err(_) => std::ptr::null_mut(),
        }
    } else {
        std::ptr::null_mut()
    }
}
