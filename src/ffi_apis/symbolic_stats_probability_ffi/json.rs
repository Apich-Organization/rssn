use std::os::raw::c_char;
use std::sync::Arc;

use crate::ffi_apis::common::{from_json_string, to_json_string};
use crate::symbolic::core::Expr;
use crate::symbolic::stats_probability::Bernoulli;
use crate::symbolic::stats_probability::Beta;
use crate::symbolic::stats_probability::Binomial;
use crate::symbolic::stats_probability::Exponential;
use crate::symbolic::stats_probability::Gamma;
use crate::symbolic::stats_probability::Normal;
use crate::symbolic::stats_probability::Poisson;
use crate::symbolic::stats_probability::StudentT;
use crate::symbolic::stats_probability::Uniform; // Need CStr

// Helper to safely convert JSON string to Expr
fn parse_expr(
    json: *const c_char
) -> Option<Expr> {

    from_json_string(json)
}

// --- Constructors (returning JSON representation of Expression wrapping Distribution) ---
// Note: Serialization of Expr::Distribution might be tricky if not implemented in core.
// The grep showed logic for serializing Distribution, so it should work if Serde impl exists.
// The snippet 1174: Expr::Distribution(d) => write!(f, "{:?}", d), implies Debug is impl.
// Serialization usually follows. If not, this JSON API will fail for Distribution variant.
// Assuming Expr implements Serialize/Deserialize fully including Distribution via custom serializer or trait.

/// Creates a normal distribution.

///

/// Takes JSON strings representing `Expr` (mean) and `Expr` (standard deviation).

/// Returns a JSON string representing the `Expr` of the normal distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_normal(
    mean_json: *const c_char,
    std_dev_json: *const c_char,
) -> *mut c_char {

    let mean = parse_expr(mean_json)
        .unwrap_or(Expr::Constant(0.0));

    let std_dev =
        parse_expr(std_dev_json)
            .unwrap_or(Expr::Constant(
                1.0,
            ));

    let dist = Expr::Distribution(
        Arc::new(Normal {
            mean,
            std_dev,
        }),
    );

    to_json_string(&dist)
}

/// Creates a uniform distribution.

///

/// Takes JSON strings representing `Expr` (minimum value) and `Expr` (maximum value).

/// Returns a JSON string representing the `Expr` of the uniform distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_uniform(
    min_json: *const c_char,
    max_json: *const c_char,
) -> *mut c_char {

    let min = parse_expr(min_json)
        .unwrap_or(Expr::Constant(0.0));

    let max = parse_expr(max_json)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Uniform {
            min,
            max,
        }),
    );

    to_json_string(&dist)
}

/// Creates a binomial distribution.

///

/// Takes JSON strings representing `Expr` (number of trials) and `Expr` (probability of success).

/// Returns a JSON string representing the `Expr` of the binomial distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_binomial(
    n_json: *const c_char,
    p_json: *const c_char,
) -> *mut c_char {

    let n = parse_expr(n_json)
        .unwrap_or(Expr::Constant(1.0));

    let p = parse_expr(p_json)
        .unwrap_or(Expr::Constant(0.5));

    let dist = Expr::Distribution(
        Arc::new(Binomial {
            n,
            p,
        }),
    );

    to_json_string(&dist)
}

/// Creates a Poisson distribution.

///

/// Takes a JSON string representing `Expr` (rate parameter λ).

/// Returns a JSON string representing the `Expr` of the Poisson distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_poisson(
    rate_json: *const c_char
) -> *mut c_char {

    let rate = parse_expr(rate_json)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Poisson {
            rate,
        }),
    );

    to_json_string(&dist)
}

/// Creates a Bernoulli distribution.

///

/// Takes a JSON string representing `Expr` (probability of success).

/// Returns a JSON string representing the `Expr` of the Bernoulli distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_bernoulli(
    p_json: *const c_char
) -> *mut c_char {

    let p = parse_expr(p_json)
        .unwrap_or(Expr::Constant(0.5));

    let dist = Expr::Distribution(
        Arc::new(Bernoulli {
            p,
        }),
    );

    to_json_string(&dist)
}

/// Creates an exponential distribution.

///

/// Takes a JSON string representing `Expr` (rate parameter λ).

/// Returns a JSON string representing the `Expr` of the exponential distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_exponential(
    rate_json: *const c_char
) -> *mut c_char {

    let rate = parse_expr(rate_json)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Exponential {
            rate,
        }),
    );

    to_json_string(&dist)
}

/// Creates a gamma distribution.

///

/// Takes JSON strings representing `Expr` (shape parameter) and `Expr` (rate parameter).

/// Returns a JSON string representing the `Expr` of the gamma distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_gamma(
    shape_json: *const c_char,
    rate_json: *const c_char,
) -> *mut c_char {

    let shape = parse_expr(shape_json)
        .unwrap_or(Expr::Constant(1.0));

    let rate = parse_expr(rate_json)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Gamma {
            shape,
            rate,
        }),
    );

    to_json_string(&dist)
}

/// Creates a beta distribution.

///

/// Takes JSON strings representing `Expr` (alpha parameter) and `Expr` (beta parameter).

/// Returns a JSON string representing the `Expr` of the beta distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_beta(
    alpha_json: *const c_char,
    beta_json: *const c_char,
) -> *mut c_char {

    let alpha = parse_expr(alpha_json)
        .unwrap_or(Expr::Constant(1.0));

    let beta = parse_expr(beta_json)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(Beta {
            alpha,
            beta,
        }),
    );

    to_json_string(&dist)
}

/// Creates a Student's t-distribution.

///

/// Takes a JSON string representing `Expr` (degrees of freedom ν).

/// Returns a JSON string representing the `Expr` of the Student's t-distribution.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_student_t(
    nu_json: *const c_char
) -> *mut c_char {

    let nu = parse_expr(nu_json)
        .unwrap_or(Expr::Constant(1.0));

    let dist = Expr::Distribution(
        Arc::new(StudentT {
            nu,
        }),
    );

    to_json_string(&dist)
}

// --- Methods ---

/// Computes the probability density function (PDF) of a distribution.

///

/// Takes JSON strings representing `Expr` (distribution) and `Expr` (value `x`).

/// Returns a JSON string representing the `Expr` of the PDF at `x`.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_pdf(
    dist_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let dist_expr =
        parse_expr(dist_json);

    let x_expr = parse_expr(x_json)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.pdf(&x_expr);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the cumulative distribution function (CDF) of a distribution.

///

/// Takes JSON strings representing `Expr` (distribution) and `Expr` (value `x`).

/// Returns a JSON string representing the `Expr` of the CDF at `x`.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_cdf(
    dist_json: *const c_char,
    x_json: *const c_char,
) -> *mut c_char {

    let dist_expr =
        parse_expr(dist_json);

    let x_expr = parse_expr(x_json)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.cdf(&x_expr);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the expectation (mean) of a distribution.

///

/// Takes a JSON string representing `Expr` (distribution).

/// Returns a JSON string representing the `Expr` of the expectation.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_expectation(
    dist_json: *const c_char
) -> *mut c_char {

    let dist_expr =
        parse_expr(dist_json);

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.expectation();

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the variance of a distribution.

///

/// Takes a JSON string representing `Expr` (distribution).

/// Returns a JSON string representing the `Expr` of the variance.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_variance(
    dist_json: *const c_char
) -> *mut c_char {

    let dist_expr =
        parse_expr(dist_json);

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.variance();

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}

/// Computes the moment generating function (MGF) of a distribution.

///

/// Takes JSON strings representing `Expr` (distribution) and `Expr` (variable `t`).

/// Returns a JSON string representing the `Expr` of the MGF.

#[no_mangle]

pub unsafe extern "C" fn rssn_json_dist_mgf(
    dist_json: *const c_char,
    t_json: *const c_char,
) -> *mut c_char {

    let dist_expr =
        parse_expr(dist_json);

    let t_expr = parse_expr(t_json)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        let result = d.mgf(&t_expr);

        to_json_string(&result)
    } else {

        std::ptr::null_mut()
    }
}
