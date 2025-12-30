use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::convergence::analyze_convergence;
use crate::symbolic::core::Expr;

/// Analyzes the convergence of a series using bincode-serialized inputs.

///

/// Takes a `BincodeBuffer` containing the series term (`Expr`) and another

/// `BincodeBuffer` for the variable (`String`).

/// Returns a `BincodeBuffer` containing the analysis result.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_analyze_convergence(
    term_buf: BincodeBuffer,

    var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let term: Option<Expr> =
        from_bincode_buffer(&term_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    match (term, var)
    { (Some(t), Some(v)) => {

        let result =
            analyze_convergence(&t, &v);

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}
