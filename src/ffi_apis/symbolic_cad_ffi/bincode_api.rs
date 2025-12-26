use serde::Deserialize;

use crate::ffi_apis::common::*;
use crate::symbolic::cad::cad;
use crate::symbolic::cad::Cad;
use crate::symbolic::core::Expr;

#[derive(Deserialize)]

struct CadInput {
    polys : Vec<Expr>,
    vars : Vec<String>,
}

/// Computes CAD for a set of polynomials via Bincode interface.
///
/// Input buffer should contain a serialized `CadInput`: `{"polys": [Expr, ...], "vars": ["x", "y", ...]}`.
#[no_mangle]

pub extern "C" fn rssn_bincode_cad(input_buf : BincodeBuffer) -> BincodeBuffer {

    let input : Option<CadInput> = from_bincode_buffer(&input_buf);

    if let Some(data) = input {

        let vars_refs : Vec<&str> = data
            .vars
            .iter()
            .map(|s| s.as_str())
            .collect();

        let mut sparse_polys = Vec::new();

        for expr in data.polys {

            let sp = crate::symbolic::polynomial::expr_to_sparse_poly(&expr, &vars_refs);

            sparse_polys.push(sp);
        }

        match cad(
            &sparse_polys,
            &vars_refs,
        ) {
            | Ok(c) => to_bincode_buffer(&c),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}
