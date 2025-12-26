use crate::ffi_apis::common::*;
use crate::symbolic::convergence::analyze_convergence;
use crate::symbolic::core::Expr;

#[no_mangle]

pub extern "C" fn rssn_bincode_analyze_convergence(
    term_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
) -> BincodeBuffer {

    let term : Option<Expr> = from_bincode_buffer(&term_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    if let (Some(t), Some(v)) = (term, var) {

        let result = analyze_convergence(&t, &v);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
