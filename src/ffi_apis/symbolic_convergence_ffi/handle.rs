use crate::symbolic::convergence::{
    analyze_convergence,
    ConvergenceResult,
};
use crate::symbolic::core::Expr;

#[no_mangle]

pub extern "C" fn rssn_analyze_convergence_handle(
    term: *const Expr,
    var: *const std::ffi::c_char,
) -> ConvergenceResult {

    let term_ref = unsafe {

        &*term
    };

    let var_str = unsafe {

        if var.is_null() {

            return ConvergenceResult::Inconclusive;
        }

        std::ffi::CStr::from_ptr(var)
            .to_string_lossy()
            .into_owned()
    };

    analyze_convergence(term_ref, &var_str)
}
