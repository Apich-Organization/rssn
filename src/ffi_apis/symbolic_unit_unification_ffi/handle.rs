use crate::symbolic::core::Expr;
use crate::symbolic::unit_unification::unify_expression;

#[no_mangle]
pub extern "C" fn rssn_unify_expression_handle(expr: *const Expr) -> *mut Expr {
    let expr_ref = unsafe { &*expr };
    match unify_expression(expr_ref) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}
