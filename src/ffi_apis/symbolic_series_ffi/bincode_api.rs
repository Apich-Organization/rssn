use crate::ffi_apis::common::*;
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

#[no_mangle]

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

    if let (
        Some(e),
        Some(v),
        Some(c),
        Some(o),
    ) = (
        expr,
        var,
        center,
        order,
    ) {

        let result = taylor_series(
            &e, &v, &c, o,
        );

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Laurent series expansion of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (center), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the Laurent series.

#[no_mangle]

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

    if let (
        Some(e),
        Some(v),
        Some(c),
        Some(o),
    ) = (
        expr,
        var,
        center,
        order,
    ) {

        let result = laurent_series(
            &e, &v, &c, o,
        );

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the Fourier series expansion of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (period), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the Fourier series.

#[no_mangle]

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

    if let (
        Some(e),
        Some(v),
        Some(p),
        Some(o),
    ) = (
        expr,
        var,
        period,
        order,
    ) {

        let result = fourier_series(
            &e, &v, &p, o,
        );

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the summation of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (lower bound), and `Expr` (upper bound).

/// Returns a bincode-serialized `Expr` representing the summation.

#[no_mangle]

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

    if let (
        Some(e),
        Some(v),
        Some(l),
        Some(u),
    ) = (
        expr, var, lower, upper,
    ) {

        let result =
            summation(&e, &v, &l, &u);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the product of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (lower bound), and `Expr` (upper bound).

/// Returns a bincode-serialized `Expr` representing the product.

#[no_mangle]

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

    if let (
        Some(e),
        Some(v),
        Some(l),
        Some(u),
    ) = (
        expr, var, lower, upper,
    ) {

        let result =
            product(&e, &v, &l, &u);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Analyzes the convergence of a series.

///

/// Takes bincode-serialized `Expr` (series expression) and `String` (variable).

/// Returns a bincode-serialized `Expr` representing the convergence analysis result.

#[no_mangle]

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

    if let (Some(s), Some(v)) =
        (series, var)
    {

        let result =
            analyze_convergence(&s, &v);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the asymptotic expansion of an expression.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (point), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the asymptotic expansion.

#[no_mangle]

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

    if let (
        Some(e),
        Some(v),
        Some(p),
        Some(o),
    ) = (
        expr, var, point, order,
    ) {

        let result =
            asymptotic_expansion(
                &e, &v, &p, o,
            );

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the analytic continuation of a series.

///

/// Takes bincode-serialized `Expr` (expression), `String` (variable),

/// `Expr` (original center), `Expr` (new center), and `usize` (order).

/// Returns a bincode-serialized `Expr` representing the analytic continuation.

#[no_mangle]

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

    if let (
        Some(e),
        Some(v),
        Some(oc),
        Some(nc),
        Some(o),
    ) = (
        expr,
        var,
        orig_center,
        new_center,
        order,
    ) {

        let result =
            analytic_continuation(
                &e, &v, &oc, &nc, o,
            );

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
