use crate::ffi_apis::common::BincodeBuffer;
use crate::ffi_apis::common::from_bincode_buffer;
use crate::ffi_apis::common::to_bincode_buffer;
use crate::symbolic::core::Expr;
use crate::symbolic::series::analytic_continuation;
use crate::symbolic::series::analyze_convergence;
use crate::symbolic::series::asymptotic_expansion;
use crate::symbolic::series::fourier_series;
use crate::symbolic::series::laurent_series;
use crate::symbolic::series::product;
use crate::symbolic::series::summation;
use crate::symbolic::series::taylor_series;

/// Computes the Taylor series expansion of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (center), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the Taylor series.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_taylor_series(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
    center_buf: BincodeBuffer,
    order_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    let center: Option<Expr> =
        from_bincode_buffer(
            &center_buf,
        );

    let order: Option<usize> =
        from_bincode_buffer(&order_buf);

    match (
        expr,
        var,
        center,
        order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(c),
            Some(o),
        ) => {

            let result = taylor_series(
                &e, &v, &c, o,
            );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the Laurent series expansion of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (center), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the Laurent series.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_laurent_series(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
    center_buf: BincodeBuffer,
    order_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    let center: Option<Expr> =
        from_bincode_buffer(
            &center_buf,
        );

    let order: Option<usize> =
        from_bincode_buffer(&order_buf);

    match (
        expr,
        var,
        center,
        order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(c),
            Some(o),
        ) => {

            let result = laurent_series(
                &e, &v, &c, o,
            );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the Fourier series expansion of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (period), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the Fourier series.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_fourier_series(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
    period_buf: BincodeBuffer,
    order_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    let period: Option<Expr> =
        from_bincode_buffer(
            &period_buf,
        );

    let order: Option<usize> =
        from_bincode_buffer(&order_buf);

    match (
        expr,
        var,
        period,
        order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(p),
            Some(o),
        ) => {

            let result = fourier_series(
                &e, &v, &p, o,
            );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the summation of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (lower bound), and `Expr` (upper bound).

/// Returns a bincode-serialized `Expr` representing the summation.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_summation(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
    lower_buf: BincodeBuffer,
    upper_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    let lower: Option<Expr> =
        from_bincode_buffer(&lower_buf);

    let upper: Option<Expr> =
        from_bincode_buffer(&upper_buf);

    match (
        expr, var, lower, upper,
    ) {
        | (
            Some(e),
            Some(v),
            Some(l),
            Some(u),
        ) => {

            let result = summation(
                &e, &v, &l, &u,
            );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the product of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (lower bound), and `Expr` (upper bound).

/// Returns a bincode-serialized `Expr` representing the product.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_product(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
    lower_buf: BincodeBuffer,
    upper_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    let lower: Option<Expr> =
        from_bincode_buffer(&lower_buf);

    let upper: Option<Expr> =
        from_bincode_buffer(&upper_buf);

    match (
        expr, var, lower, upper,
    ) {
        | (
            Some(e),
            Some(v),
            Some(l),
            Some(u),
        ) => {

            let result =
                product(&e, &v, &l, &u);

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Analyzes the convergence of a series.

///

/// Takes bincode-serialized `Expr` (series expression) and `String` (variable).

/// Returns a bincode-serialized `Expr` representing the convergence analysis result.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_series_bincode_analyze_convergence(
    series_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
) -> BincodeBuffer {

    let series: Option<Expr> =
        from_bincode_buffer(
            &series_buf,
        );

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    match (series, var) {
        | (Some(s), Some(v)) => {

            let result =
                analyze_convergence(
                    &s, &v,
                );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the asymptotic expansion of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (point), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the asymptotic expansion.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_asymptotic_expansion(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
    point_buf: BincodeBuffer,
    order_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    let point: Option<Expr> =
        from_bincode_buffer(&point_buf);

    let order: Option<usize> =
        from_bincode_buffer(&order_buf);

    match (
        expr, var, point, order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(p),
            Some(o),
        ) => {

            let result =
                asymptotic_expansion(
                    &e, &v, &p, o,
                );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}

/// Computes the analytic continuation of a series.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (original center), `Expr` (new center), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the analytic continuation.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_analytic_continuation(
    expr_buf: BincodeBuffer,
    var_buf: BincodeBuffer,
    orig_center_buf: BincodeBuffer,
    new_center_buf: BincodeBuffer,
    order_buf: BincodeBuffer,
) -> BincodeBuffer {

    let expr: Option<Expr> =
        from_bincode_buffer(&expr_buf);

    let var: Option<String> =
        from_bincode_buffer(&var_buf);

    let orig_center: Option<Expr> =
        from_bincode_buffer(
            &orig_center_buf,
        );

    let new_center: Option<Expr> =
        from_bincode_buffer(
            &new_center_buf,
        );

    let order: Option<usize> =
        from_bincode_buffer(&order_buf);

    match (
        expr,
        var,
        orig_center,
        new_center,
        order,
    ) {
        | (
            Some(e),
            Some(v),
            Some(oc),
            Some(nc),
            Some(o),
        ) => {

            let result =
                analytic_continuation(
                    &e, &v, &oc, &nc, o,
                );

            to_bincode_buffer(&result)
        },
        | _ => BincodeBuffer::empty(),
    }
}
