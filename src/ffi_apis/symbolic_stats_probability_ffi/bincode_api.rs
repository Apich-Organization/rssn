use std::sync::Arc;

use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::core::Expr;
use crate::symbolic::stats_probability::Bernoulli;
use crate::symbolic::stats_probability::Beta;
use crate::symbolic::stats_probability::Binomial;
use crate::symbolic::stats_probability::Exponential;
use crate::symbolic::stats_probability::Gamma;
use crate::symbolic::stats_probability::Normal;
use crate::symbolic::stats_probability::Poisson;
use crate::symbolic::stats_probability::StudentT;
use crate::symbolic::stats_probability::Uniform; // Technically unused here but consistent with other files

// Helper to safely convert BincodeBuffer to Expr
pub(crate) fn parse_expr(
    buf: BincodeBuffer
) -> Option<Expr> {

    from_bincode_buffer(&buf)
}

/// Creates a normal distribution.

///

/// Takes bincode-serialized `Expr` representing the mean and standard deviation.

/// Returns a bincode-serialized `Expr` representing the normal distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_normal(
    mean_buf: BincodeBuffer,
    std_dev_buf: BincodeBuffer,
) -> BincodeBuffer {

    let mean = parse_expr(mean_buf)
        .unwrap_or(Expr::Constant(0.0));

    let std_dev =
        parse_expr(std_dev_buf)
            .unwrap_or(Expr::Constant(
                1.0,
            ));

    let dist = Expr::Distribution(
        Arc::new(Normal {
            mean,
            std_dev,
        }),
    );

    to_bincode_buffer(&dist)
}

/// Creates a uniform distribution.

///

/// Takes bincode-serialized `Expr` representing the minimum and maximum values.

/// Returns a bincode-serialized `Expr` representing the uniform distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_uniform(
    min_buf: BincodeBuffer,
    max_buf: BincodeBuffer,
) -> BincodeBuffer {

    let min = parse_expr(min_buf)
        .unwrap_or(Expr::Constant(0.0));

    let max = parse_expr(max_buf)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Uniform {
            min,
            max,
        }),
    );

    to_bincode_buffer(&dist)
}

/// Creates a binomial distribution.

///

/// Takes bincode-serialized `Expr` representing `n` (number of trials) and `p` (probability of success).

/// Returns a bincode-serialized `Expr` representing the binomial distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_binomial(
    n_buf: BincodeBuffer,
    p_buf: BincodeBuffer,
) -> BincodeBuffer {

    let n = parse_expr(n_buf)
        .unwrap_or(Expr::Constant(1.0));

    let p = parse_expr(p_buf)
        .unwrap_or(Expr::Constant(0.5));

    let dist = Expr::Distribution(
        Arc::new(Binomial {
            n,
            p,
        }),
    );

    to_bincode_buffer(&dist)
}

/// Creates a Poisson distribution.

///

/// Takes a bincode-serialized `Expr` representing the rate parameter (λ).

/// Returns a bincode-serialized `Expr` representing the Poisson distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_poisson(
    rate_buf: BincodeBuffer
) -> BincodeBuffer {

    let rate = parse_expr(rate_buf)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Poisson {
            rate,
        }),
    );

    to_bincode_buffer(&dist)
}

/// Creates a Bernoulli distribution.

///

/// Takes a bincode-serialized `Expr` representing `p` (probability of success).

/// Returns a bincode-serialized `Expr` representing the Bernoulli distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_bernoulli(
    p_buf: BincodeBuffer
) -> BincodeBuffer {

    let p = parse_expr(p_buf)
        .unwrap_or(Expr::Constant(0.5));

    let dist = Expr::Distribution(
        Arc::new(Bernoulli {
            p,
        }),
    );

    to_bincode_buffer(&dist)
}

/// Creates an exponential distribution.

///

/// Takes a bincode-serialized `Expr` representing the rate parameter (λ).

/// Returns a bincode-serialized `Expr` representing the exponential distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_exponential(
    rate_buf: BincodeBuffer
) -> BincodeBuffer {

    let rate = parse_expr(rate_buf)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Exponential {
            rate,
        }),
    );

    to_bincode_buffer(&dist)
}

/// Creates a gamma distribution.

///

/// Takes bincode-serialized `Expr` representing the shape and rate parameters.

/// Returns a bincode-serialized `Expr` representing the gamma distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_gamma(
    shape_buf: BincodeBuffer,
    rate_buf: BincodeBuffer,
) -> BincodeBuffer {

    let shape = parse_expr(shape_buf)
        .unwrap_or(Expr::Constant(1.0));

    let rate = parse_expr(rate_buf)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Gamma {
            shape,
            rate,
        }),
    );

    to_bincode_buffer(&dist)
}

/// Creates a beta distribution.

///

/// Takes bincode-serialized `Expr` representing the alpha and beta parameters.

/// Returns a bincode-serialized `Expr` representing the beta distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_beta(
    alpha_buf: BincodeBuffer,
    beta_buf: BincodeBuffer,
) -> BincodeBuffer {

    let alpha = parse_expr(alpha_buf)
        .unwrap_or(Expr::Constant(1.0));

    let beta = parse_expr(beta_buf)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Beta {
            alpha,
            beta,
        }),
    );

    to_bincode_buffer(&dist)
}

/// Creates a Student's t-distribution.

///

/// Takes a bincode-serialized `Expr` representing the degrees of freedom (ν).

/// Returns a bincode-serialized `Expr` representing the Student's t-distribution.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_student_t(
    nu_buf: BincodeBuffer
) -> BincodeBuffer {

    let nu = parse_expr(nu_buf)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(StudentT {
            nu,
        }),
    );

    to_bincode_buffer(&dist)
}

// --- Methods ---

/// Computes the probability density function (PDF) of a distribution.

///

/// Takes bincode-serialized `Expr` representing the distribution and the value `x`.

/// Returns a bincode-serialized `Expr` representing the PDF at `x`.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_pdf(
    dist_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let dist_expr =
        parse_expr(dist_buf);

    let x_expr = parse_expr(x_buf)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.pdf(&x_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the cumulative distribution function (CDF) of a distribution.

///

/// Takes bincode-serialized `Expr` representing the distribution and the value `x`.

/// Returns a bincode-serialized `Expr` representing the CDF at `x`.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_cdf(
    dist_buf: BincodeBuffer,
    x_buf: BincodeBuffer,
) -> BincodeBuffer {

    let dist_expr =
        parse_expr(dist_buf);

    let x_expr = parse_expr(x_buf)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.cdf(&x_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the expectation (mean) of a distribution.

///

/// Takes a bincode-serialized `Expr` representing the distribution.

/// Returns a bincode-serialized `Expr` representing the expectation.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_expectation(
    dist_buf: BincodeBuffer
) -> BincodeBuffer {

    let dist_expr =
        parse_expr(dist_buf);

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.expectation();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the variance of a distribution.

///

/// Takes a bincode-serialized `Expr` representing the distribution.

/// Returns a bincode-serialized `Expr` representing the variance.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_variance(
    dist_buf: BincodeBuffer
) -> BincodeBuffer {

    let dist_expr =
        parse_expr(dist_buf);

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.variance();

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

/// Computes the moment generating function (MGF) of a distribution.

///

/// Takes bincode-serialized `Expr` representing the distribution and the variable `t`.

/// Returns a bincode-serialized `Expr` representing the MGF.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_dist_mgf(
    dist_buf: BincodeBuffer,
    t_buf: BincodeBuffer,
) -> BincodeBuffer {

    let dist_expr =
        parse_expr(dist_buf);

    let t_expr = parse_expr(t_buf)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.mgf(&t_expr);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
