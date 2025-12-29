use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::logic::is_satisfiable;
use crate::symbolic::logic::simplify_logic;
use crate::symbolic::logic::to_cnf;
use crate::symbolic::logic::to_dnf;

/// Simplifies a logical expression using bincode-based FFI.
#[no_mangle]

pub extern "C" fn rssn_bincode_simplify_logic(
    expr_buf: BincodeBuffer
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        let result = simplify_logic(&e);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Converts a logical expression to Conjunctive Normal Form (CNF) using bincode-based FFI.
#[no_mangle]

pub extern "C" fn rssn_bincode_to_cnf(
    expr_buf: BincodeBuffer
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        let result = to_cnf(&e);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Converts a logical expression to Disjunctive Normal Form (DNF) using bincode-based FFI.
#[no_mangle]

pub extern "C" fn rssn_bincode_to_dnf(
    expr_buf: BincodeBuffer
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        let result = to_dnf(&e);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Checks if a logical expression is satisfiable using bincode-based FFI.
///
/// Returns a bincode buffer containing:
/// - `Some(true)` if satisfiable
/// - `Some(false)` if unsatisfiable
/// - `None` if the expression contains quantifiers (undecidable)
#[no_mangle]

pub extern "C" fn rssn_bincode_is_satisfiable(
    expr_buf: BincodeBuffer
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    if let Some(e) = expr {

        let result = is_satisfiable(&e);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
