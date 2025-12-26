use crate::symbolic::core::Expr;
use crate::symbolic::numeric::evaluate_numerical;

#[no_mangle]

pub extern "C" fn rssn_evaluate_numerical_handle(
    expr : *const Expr
) -> f64 {

    let expr_ref = unsafe {

        &*expr
    };

    evaluate_numerical(expr_ref)
        .unwrap_or(f64::NAN)
}
