use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::ffi_apis::common::BincodeBuffer;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_regression;

/// Performs a simple linear regression.

///

/// Takes a bincode-serialized `Vec<(Expr, Expr)>` representing the data points.

/// Returns a bincode-serialized `Vec<Expr>` containing the intercept and slope coefficients.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_simple_linear_regression(
    data_buf: BincodeBuffer
) -> BincodeBuffer {

    let data: Option<
        Vec<(Expr, Expr)>,
    > = from_bincode_buffer(&data_buf);

    if let Some(data) = data {

        let (b0, b1) = stats_regression::simple_linear_regression_symbolic(&data);

        to_bincode_buffer(&vec![b0, b1])
    } else {

        BincodeBuffer::empty()
    }
}

/// Performs a polynomial regression.

///

/// Takes a bincode-serialized `Vec<(Expr, Expr)>` representing the data points

/// and a `usize` for the degree of the polynomial.

/// Returns a bincode-serialized `Vec<Expr>` containing the coefficients of the polynomial.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_polynomial_regression(
    data_buf: BincodeBuffer,
    degree: usize,
) -> BincodeBuffer {

    let data: Option<
        Vec<(Expr, Expr)>,
    > = from_bincode_buffer(&data_buf);

    if let Some(data) = data {

        match stats_regression::polynomial_regression_symbolic(&data, degree) {
            | Ok(coeffs) => to_bincode_buffer(&coeffs),
            | Err(_) => BincodeBuffer::empty(),
        }
    } else {

        BincodeBuffer::empty()
    }
}

/// Performs a nonlinear regression.

///

/// Takes bincode-serialized `Vec<(Expr, Expr)>` (data points), `Expr` (model),

/// `Vec<String>` (variables), and `Vec<String>` (parameters).

/// Returns a bincode-serialized `Vec<Expr>` representing the optimized parameter values.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_nonlinear_regression(
    data_buf: BincodeBuffer,
    model_buf: BincodeBuffer,
    vars_buf: BincodeBuffer,
    params_buf: BincodeBuffer,
) -> BincodeBuffer {

    let data: Option<
        Vec<(Expr, Expr)>,
    > = from_bincode_buffer(&data_buf);

    let model: Option<Expr> =
        from_bincode_buffer(&model_buf);

    let vars: Option<Vec<String>> =
        from_bincode_buffer(&vars_buf);

    let params: Option<Vec<String>> =
        from_bincode_buffer(
            &params_buf,
        );

    match (
        data,
        model,
        vars,
        params,
    ) { (
        Some(data),
        Some(model),
        Some(vars),
        Some(params),
    ) => {

        let vars_refs: Vec<&str> = vars
            .iter()
            .map(std::string::String::as_str)
            .collect();

        let params_refs: Vec<&str> =
            params
                .iter()
                .map(std::string::String::as_str)
                .collect();

        match stats_regression::nonlinear_regression_symbolic(
            &data,
            &model,
            &vars_refs,
            &params_refs,
        ) {
            | Some(solutions) => to_bincode_buffer(&solutions),
            | None => BincodeBuffer::empty(),
        }
    } _ => {

        BincodeBuffer::empty()
    }}
}
