use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::optimize::find_constrained_extrema;
use crate::symbolic::optimize::find_extrema;
use crate::symbolic::optimize::hessian_matrix;

/// Finds extrema of a function (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_find_extrema(
    expr_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    match (expr, vars) {
        | (Some(e), Some(v)) => {

            let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

            match find_extrema(
                &e,
                &vars_refs,
            ) {
                | Ok(points) => {
                    to_bincode_buffer(
                        &points,
                    )
                },
                | Err(_) => {
                    BincodeBuffer::empty(
                    )
                },
            }
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes Hessian matrix (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_hessian_matrix(
    expr_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    match (expr, vars) {
        | (Some(e), Some(v)) => {

            let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

            let hessian =
                hessian_matrix(
                    &e,
                    &vars_refs,
                );

            to_bincode_buffer(&hessian)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Finds constrained extrema (Bincode)
#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_find_constrained_extrema(
    expr_buf: BincodeBuffer,
    constraints_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let constraints: Option<Vec<Expr>> =
        from_bincode_buffer(
            &constraints_buf,
        );

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    match (
        expr,
        constraints,
        vars,
    ) {
        | (
            Some(e),
            Some(c),
            Some(v),
        ) => {

            let vars_refs: Vec<&str> = v
            .iter()
            .map(std::string::String::as_str)
            .collect();

            match find_constrained_extrema(
            &e,
            &c,
            &vars_refs,
        ) {
            | Ok(solutions) => {
                to_bincode_buffer(
                    &solutions,
                )
            },
            | Err(_) => {
                BincodeBuffer::empty()
            },
        }
        },
        | _ => BincodeBuffer::empty(),
    }
}
