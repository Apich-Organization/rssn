use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::quantum_mechanics::Bra;
use crate::symbolic::quantum_mechanics::Ket;
use crate::symbolic::quantum_mechanics::Operator;
use crate::symbolic::quantum_mechanics::{
    self,
};

/// Computes the expectation value using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_expectation_value(
    op_json : *const c_char,
    psi_json : *const c_char,
) -> *mut c_char {

    let op : Option<Operator> =
        from_json_string(op_json);

    let psi : Option<Ket> =
        from_json_string(psi_json);

    if let (Some(op), Some(psi)) =
        (op, psi)
    {

        to_json_string(&quantum_mechanics::expectation_value(&op, &psi))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the uncertainty using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_uncertainty(
    op_json : *const c_char,
    psi_json : *const c_char,
) -> *mut c_char {

    let op : Option<Operator> =
        from_json_string(op_json);

    let psi : Option<Ket> =
        from_json_string(psi_json);

    if let (Some(op), Some(psi)) =
        (op, psi)
    {

        to_json_string(&quantum_mechanics::uncertainty(&op, &psi))
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the inner product <Bra|Ket> using JSON.
#[no_mangle]

pub extern "C" fn rssn_json_bra_ket(
    bra_json : *const c_char,
    ket_json : *const c_char,
) -> *mut c_char {

    let bra : Option<Bra> =
        from_json_string(bra_json);

    let ket : Option<Ket> =
        from_json_string(ket_json);

    if let (Some(bra), Some(ket)) =
        (bra, ket)
    {

        to_json_string(
            &quantum_mechanics::bra_ket(
                &bra, &ket,
            ),
        )
    } else {

        std::ptr::null_mut()
    }
}
