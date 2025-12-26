use std::ffi::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::convergence::analyze_convergence;
use crate::symbolic::core::Expr;

#[no_mangle]

pub extern "C" fn rssn_json_analyze_convergence(
    term_json : *const c_char,
    var_json : *const c_char,
) -> *mut c_char {

    let term : Option<Expr> =
        from_json_string(term_json);

    let var : Option<String> =
        from_json_string(var_json);

    if let (Some(t), Some(v)) =
        (term, var)
    {

        let result =
            analyze_convergence(&t, &v);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
