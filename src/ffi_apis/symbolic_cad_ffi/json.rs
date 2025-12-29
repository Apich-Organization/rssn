use std::os::raw::c_char;

use serde::Deserialize;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::cad::cad;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct CadInput {
    polys: Vec<Expr>,
    vars: Vec<String>,
}

/// Computes CAD for a set of polynomials via JSON interface.
///
/// Input JSON should be an object: `{"polys": [Expr, ...], "vars": ["x", "y", ...]}`.
#[no_mangle]

pub unsafe extern "C" fn rssn_json_cad(
    input_json: *const c_char
) -> *mut c_char {

    let input: Option<CadInput> =
        from_json_string(input_json);

    if let Some(data) = input {

        let vars_refs: Vec<&str> = data
            .vars
            .iter()
            .map(std::string::String::as_str)
            .collect();

        let mut sparse_polys =
            Vec::new();

        for expr in data.polys {

            let sp = crate::symbolic::polynomial::expr_to_sparse_poly(&expr, &vars_refs);

            sparse_polys.push(sp);
        }

        match cad(
            &sparse_polys,
            &vars_refs,
        ) {
            | Ok(c) => {
                to_json_string(&c)
            },
            | Err(_) => {
                std::ptr::null_mut()
            },
        }
    } else {

        std::ptr::null_mut()
    }
}
