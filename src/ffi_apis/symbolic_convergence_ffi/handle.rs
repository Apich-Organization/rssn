use crate::symbolic::convergence::analyze_convergence;
use crate::symbolic::convergence::ConvergenceResult;
use crate::symbolic::core::Expr;

/// Analyzes the convergence of a series using direct pointers.

///

/// Takes a pointer to an `Expr` representing the series term and a C-style string

/// for the variable.

/// Returns a `ConvergenceResult` enum indicating the analysis outcome.

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

    analyze_convergence(
        term_ref,
        &var_str,
    )
}
