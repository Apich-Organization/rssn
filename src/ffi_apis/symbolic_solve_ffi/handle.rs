use crate::symbolic::core::Expr;
use crate::symbolic::solve::{solve, solve_linear_system, solve_system};

#[no_mangle]
pub extern "C" fn rssn_solve_handle(
    expr: *const Expr,
    var: *const std::ffi::c_char,
) -> *mut Vec<Expr> {
    let expr_ref = unsafe { &*expr };
    let var_str = unsafe {
        if var.is_null() {
            return std::ptr::null_mut();
        }
        std::ffi::CStr::from_ptr(var).to_string_lossy().into_owned()
    };

    let result = solve(expr_ref, &var_str);
    Box::into_raw(Box::new(result))
}

#[no_mangle]
pub extern "C" fn rssn_solve_system_handle(
    equations: *const Vec<Expr>,
    vars: *const Vec<String>,
) -> *mut Vec<(Expr, Expr)> {
    let eqs_ref = unsafe { &*equations };
    let vars_ref = unsafe { &*vars };
    let vars_str: Vec<&str> = vars_ref.iter().map(|s| s.as_str()).collect();

    match solve_system(eqs_ref, &vars_str) {
        Some(result) => Box::into_raw(Box::new(result)),
        None => std::ptr::null_mut(),
    }
}

#[no_mangle]
pub extern "C" fn rssn_solve_linear_system_handle(
    system: *const Expr,
    vars: *const Vec<String>,
) -> *mut Vec<Expr> {
    let sys_ref = unsafe { &*system };
    let vars_ref = unsafe { &*vars };

    match solve_linear_system(sys_ref, vars_ref) {
        Ok(result) => Box::into_raw(Box::new(result)),
        Err(_) => std::ptr::null_mut(),
    }
}
