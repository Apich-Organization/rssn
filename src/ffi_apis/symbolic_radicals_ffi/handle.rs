use crate::symbolic::core::Expr;
use crate::symbolic::radicals::denest_sqrt;
use crate::symbolic::radicals::simplify_radicals;

/// Simplifies radical expressions (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_simplify_radicals_handle(
    expr: *const Expr
) -> *mut Expr {

    let expr_ref = unsafe {

        if expr.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*expr
    };

    let result =
        simplify_radicals(expr_ref);

    Box::into_raw(Box::new(result))
}

/// Denests a nested square root (Handle)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_denest_sqrt_handle(
    expr: *const Expr
) -> *mut Expr {

    let expr_ref = unsafe {

        if expr.is_null() {

            return std::ptr::null_mut(
            );
        }

        &*expr
    };

    let result = denest_sqrt(expr_ref);

    Box::into_raw(Box::new(result))
}
