use std::ffi::c_char;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::convergence::analyze_convergence;
use crate::symbolic::core::Expr;

/// Analyzes the convergence of a series using JSON-serialized inputs.

///

/// Takes C-style strings containing JSON-serialized `Expr` for the series term

/// and `String` for the variable.

/// Returns a C-style string containing the JSON-serialized analysis result.

#[no_mangle]

pub extern "C" fn rssn_json_analyze_convergence(
    term_json: *const c_char,

    var_json: *const c_char,
) -> *mut c_char {

    let term: Option<Expr> =
        from_json_string(term_json);

    let var: Option<String> =
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
