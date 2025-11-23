use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::solve::{solve, solve_linear_system, solve_system};

#[no_mangle]
pub extern "C" fn rssn_bincode_solve(expr_buf: BincodeBuffer, var_buf: BincodeBuffer) -> BincodeBuffer {
    let expr: Option<Expr> = from_bincode_buffer(&expr_buf);
    let var: Option<String> = from_bincode_buffer(&var_buf);
    
    if let (Some(e), Some(v)) = (expr, var) {
        let result = solve(&e, &v);
        to_bincode_buffer(&result)
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_solve_system(equations_buf: BincodeBuffer, vars_buf: BincodeBuffer) -> BincodeBuffer {
    let equations: Option<Vec<Expr>> = from_bincode_buffer(&equations_buf);
    let vars: Option<Vec<String>> = from_bincode_buffer(&vars_buf);
    
    if let (Some(eqs), Some(vs)) = (equations, vars) {
        let vars_str: Vec<&str> = vs.iter().map(|s| s.as_str()).collect();
        match solve_system(&eqs, &vars_str) {
            Some(result) => to_bincode_buffer(&result),
            None => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}

#[no_mangle]
pub extern "C" fn rssn_bincode_solve_linear_system(system_buf: BincodeBuffer, vars_buf: BincodeBuffer) -> BincodeBuffer {
    let system: Option<Expr> = from_bincode_buffer(&system_buf);
    let vars: Option<Vec<String>> = from_bincode_buffer(&vars_buf);
    
    if let (Some(sys), Some(vs)) = (system, vars) {
        match solve_linear_system(&sys, &vs) {
            Ok(result) => to_bincode_buffer(&result),
            Err(_) => BincodeBuffer::empty(),
        }
    } else {
        BincodeBuffer::empty()
    }
}
