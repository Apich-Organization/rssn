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

#[no_mangle]

pub extern "C" fn rssn_bincode_taylor_series(
    expr_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
    center_buf : BincodeBuffer,
    order_buf : BincodeBuffer,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    let center : Option<Expr> = from_bincode_buffer(&center_buf);

    let order : Option<usize> = from_bincode_buffer(&order_buf);

    if let (Some(e), Some(v), Some(c), Some(o)) = (
        expr,
        var,
        center,
        order,
    ) {

        let result = taylor_series(&e, &v, &c, o);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_laurent_series(
    expr_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
    center_buf : BincodeBuffer,
    order_buf : BincodeBuffer,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    let center : Option<Expr> = from_bincode_buffer(&center_buf);

    let order : Option<usize> = from_bincode_buffer(&order_buf);

    if let (Some(e), Some(v), Some(c), Some(o)) = (
        expr,
        var,
        center,
        order,
    ) {

        let result = laurent_series(&e, &v, &c, o);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_fourier_series(
    expr_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
    period_buf : BincodeBuffer,
    order_buf : BincodeBuffer,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    let period : Option<Expr> = from_bincode_buffer(&period_buf);

    let order : Option<usize> = from_bincode_buffer(&order_buf);

    if let (Some(e), Some(v), Some(p), Some(o)) = (
        expr,
        var,
        period,
        order,
    ) {

        let result = fourier_series(&e, &v, &p, o);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_summation(
    expr_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
    lower_buf : BincodeBuffer,
    upper_buf : BincodeBuffer,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    let lower : Option<Expr> = from_bincode_buffer(&lower_buf);

    let upper : Option<Expr> = from_bincode_buffer(&upper_buf);

    if let (Some(e), Some(v), Some(l), Some(u)) = (
        expr, var, lower, upper,
    ) {

        let result = summation(&e, &v, &l, &u);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_product(
    expr_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
    lower_buf : BincodeBuffer,
    upper_buf : BincodeBuffer,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    let lower : Option<Expr> = from_bincode_buffer(&lower_buf);

    let upper : Option<Expr> = from_bincode_buffer(&upper_buf);

    if let (Some(e), Some(v), Some(l), Some(u)) = (
        expr, var, lower, upper,
    ) {

        let result = product(&e, &v, &l, &u);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_series_bincode_analyze_convergence(
    series_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
) -> BincodeBuffer {

    let series : Option<Expr> = from_bincode_buffer(&series_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    if let (Some(s), Some(v)) = (series, var) {

        let result = analyze_convergence(&s, &v);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_asymptotic_expansion(
    expr_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
    point_buf : BincodeBuffer,
    order_buf : BincodeBuffer,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    let point : Option<Expr> = from_bincode_buffer(&point_buf);

    let order : Option<usize> = from_bincode_buffer(&order_buf);

    if let (Some(e), Some(v), Some(p), Some(o)) = (
        expr, var, point, order,
    ) {

        let result = asymptotic_expansion(&e, &v, &p, o);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_analytic_continuation(
    expr_buf : BincodeBuffer,
    var_buf : BincodeBuffer,
    orig_center_buf : BincodeBuffer,
    new_center_buf : BincodeBuffer,
    order_buf : BincodeBuffer,
) -> BincodeBuffer {

    let expr : Option<Expr> = from_bincode_buffer(&expr_buf);

    let var : Option<String> = from_bincode_buffer(&var_buf);

    let orig_center : Option<Expr> = from_bincode_buffer(&orig_center_buf);

    let new_center : Option<Expr> = from_bincode_buffer(&new_center_buf);

    let order : Option<usize> = from_bincode_buffer(&order_buf);

    if let (Some(e), Some(v), Some(oc), Some(nc), Some(o)) = (
        expr,
        var,
        orig_center,
        new_center,
        order,
    ) {

        let result = analytic_continuation(&e, &v, &oc, &nc, o);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
