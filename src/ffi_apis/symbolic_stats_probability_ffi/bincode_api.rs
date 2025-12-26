use crate::ffi_apis::common::*;
use crate::symbolic::core::{
    Distribution,
    Expr,
};
use crate::symbolic::stats_probability::{
    Bernoulli,
    Beta,
    Binomial,
    Exponential,
    Gamma,
    Normal,
    Poisson,
    StudentT,
    Uniform,
};
use std::os::raw::c_char;
use std::sync::Arc; // Technically unused here but consistent with other files

// Helper to safely convert BincodeBuffer to Expr
fn parse_expr(buf: BincodeBuffer) -> Option<Expr> {

    from_bincode_buffer(&buf)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_normal(
    mean_buf: BincodeBuffer,
    std_dev_buf: BincodeBuffer,
) -> BincodeBuffer {

    let mean = parse_expr(mean_buf).unwrap_or(Expr::Constant(0.0));

    let std_dev = parse_expr(std_dev_buf).unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(Arc::new(Normal {
        mean,
        std_dev,
    }));

    to_bincode_buffer(&dist)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_uniform(
    min_buf: BincodeBuffer,
    max_buf: BincodeBuffer,
) -> BincodeBuffer {

    let min = parse_expr(min_buf).unwrap_or(Expr::Constant(0.0));

    let max = parse_expr(max_buf).unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(Arc::new(Uniform {
        min,
        max,
    }));

    to_bincode_buffer(&dist)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_binomial(
    n_buf: BincodeBuffer,
    p_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n = parse_expr(n_buf).unwrap_or(Expr::Constant(1.0));

    let p = parse_expr(p_buf).unwrap_or(Expr::Constant(0.5));

    let dist = Expr::Distribution(Arc::new(Binomial {
        n,
        p,
    }));

    to_bincode_buffer(&dist)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_poisson(rate_buf: BincodeBuffer) -> BincodeBuffer {

    let rate = parse_expr(rate_buf).unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(Arc::new(Poisson {
        rate,
    }));

    to_bincode_buffer(&dist)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_bernoulli(p_buf: BincodeBuffer) -> BincodeBuffer {

    let p = parse_expr(p_buf).unwrap_or(Expr::Constant(0.5));

    let dist = Expr::Distribution(Arc::new(
        Bernoulli { p },
    ));

    to_bincode_buffer(&dist)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_exponential(rate_buf: BincodeBuffer) -> BincodeBuffer {

    let rate = parse_expr(rate_buf).unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(Arc::new(
        Exponential { rate },
    ));

    to_bincode_buffer(&dist)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_gamma(
    shape_buf: BincodeBuffer,
    rate_buf: BincodeBuffer,
) -> BincodeBuffer {

    let shape = parse_expr(shape_buf).unwrap_or(Expr::Constant(1.0));

    let rate = parse_expr(rate_buf).unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(Arc::new(Gamma {
        shape,
        rate,
    }));

    to_bincode_buffer(&dist)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_beta(
    alpha_buf: BincodeBuffer,
    beta_buf: BincodeBuffer,
) -> BincodeBuffer {

    let alpha = parse_expr(alpha_buf).unwrap_or(Expr::Constant(1.0));

    let beta = parse_expr(beta_buf).unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(Arc::new(Beta {
        alpha,
        beta,
    }));

    to_bincode_buffer(&dist)
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_student_t(nu_buf: BincodeBuffer) -> BincodeBuffer {

    let nu = parse_expr(nu_buf).unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(Arc::new(StudentT {
        nu,
    }));

    to_bincode_buffer(&dist)
}

// --- Methods ---

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_pdf(
    dist_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let dist_expr = parse_expr(dist_buf);

    let x_expr = parse_expr(x_buf).unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) = dist_expr {

        let result = d.pdf(&x_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_cdf(
    dist_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let dist_expr = parse_expr(dist_buf);

    let x_expr = parse_expr(x_buf).unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) = dist_expr {

        let result = d.cdf(&x_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_expectation(dist_buf: BincodeBuffer) -> BincodeBuffer {

    let dist_expr = parse_expr(dist_buf);

    if let Some(Expr::Distribution(d)) = dist_expr {

        let result = d.expectation();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_variance(dist_buf: BincodeBuffer) -> BincodeBuffer {

    let dist_expr = parse_expr(dist_buf);

    if let Some(Expr::Distribution(d)) = dist_expr {

        let result = d.variance();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_dist_mgf(
    dist_buf: BincodeBuffer,
    t_buf: BincodeBuffer,
) -> BincodeBuffer {

    let dist_expr = parse_expr(dist_buf);

    let t_expr = parse_expr(t_buf).unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) = dist_expr {

        let result = d.mgf(&t_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
