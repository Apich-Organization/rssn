//! JSON-based FFI API for thermodynamics functions.

use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::thermodynamics;

/// Calculates ideal gas Law using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_ideal_gas_law(
    p_json : *const c_char,
    v_json : *const c_char,
    n_json : *const c_char,
    r_json : *const c_char,
    t_json : *const c_char,
) -> *mut c_char {

    let p : Option<Expr> = from_json_string(p_json);

    let v : Option<Expr> = from_json_string(v_json);

    let n : Option<Expr> = from_json_string(n_json);

    let r : Option<Expr> = from_json_string(r_json);

    let t : Option<Expr> = from_json_string(t_json);

    if let (Some(p), Some(v), Some(n), Some(r), Some(t)) = (p, v, n, r, t) {

        to_json_string(&thermodynamics::ideal_gas_law(&p, &v, &n, &r, &t))
    } else {

        std::ptr::null_mut()
    }
}

/// Calculates Gibbs Free Energy using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_gibbs_free_energy(
    h_json : *const c_char,
    t_json : *const c_char,
    s_json : *const c_char,
) -> *mut c_char {

    let h : Option<Expr> = from_json_string(h_json);

    let t : Option<Expr> = from_json_string(t_json);

    let s : Option<Expr> = from_json_string(s_json);

    if let (Some(h), Some(t), Some(s)) = (h, t, s) {

        to_json_string(&thermodynamics::gibbs_free_energy(&h, &t, &s))
    } else {

        std::ptr::null_mut()
    }
}
