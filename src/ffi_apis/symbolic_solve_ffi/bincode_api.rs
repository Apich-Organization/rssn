use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::solve::solve;
use crate::symbolic::solve::solve_linear_system;
use crate::symbolic::solve::solve_system;

/// Solves an equation for a given variable.

///

/// Takes bincode-serialized `Expr` (equation) and `String` (variable).

/// Returns a bincode-serialized `Expr` representing the solution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_solve(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    match (expr, var)
    { (Some(e), Some(v)) => {

        let result = solve(&e, &v);

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Solves a system of equations for given variables.

///

/// Takes bincode-serialized `Vec<Expr>` (equations) and `Vec<String>` (variables).

/// Returns a bincode-serialized `Expr` representing the solution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_solve_system(
    equations_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let equations: Option<Vec<Expr>> =
        from_bincode_buffer(
            &equations_buf,
        );

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    match (equations, vars)
    { (Some(eqs), Some(vs)) => {

        let vars_str: Vec<&str> = vs
            .iter()
            .map(std::string::String::as_str)
            .collect();

        match solve_system(
            &eqs,
            &vars_str,
        ) {
            | Some(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | None => {
                BincodeBuffer::empty()
            },
        }
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Solves a linear system of equations.

///

/// Takes bincode-serialized `Expr` (system) and `Vec<String>` (variables).

/// Returns a bincode-serialized `Expr` representing the solution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_solve_linear_system(
    system_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let system: Option<Expr> =
        from_bincode_buffer(
            &system_buf,
        );

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    match (system, vars)
    { (Some(sys), Some(vs)) => {

        match solve_linear_system(
            &sys, &vs,
        ) {
            | Ok(result) => {
                to_bincode_buffer(
                    &result,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
    } _ => {

        BincodeBuffer::empty()
    }}
}
