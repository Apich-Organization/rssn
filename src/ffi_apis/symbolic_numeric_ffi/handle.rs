use crate::symbolic::core::Expr;
use crate::symbolic::numeric::evaluate_numerical;

/// Numerically evaluates a symbolic expression.

///

/// Takes a raw pointer to an `Expr` as input,

/// and returns an `f64` representing the numerical evaluation of that expression.

#[no_mangle]

pub extern "C" fn rssn_evaluate_numerical_handle(
    expr: *const Expr
) -> f64 {

    let expr_ref = unsafe {

        &*expr
    };

    evaluate_numerical(expr_ref)
        .unwrap_or(f64::NAN)
}
