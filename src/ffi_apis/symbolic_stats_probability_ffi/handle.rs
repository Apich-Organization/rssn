use std::ffi::CStr;
use std::ffi::CString;
use std::os::raw::c_char;
use std::sync::Arc;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Distribution;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_probability::Bernoulli;
use crate::symbolic::stats_probability::Beta;
use crate::symbolic::stats_probability::Binomial;
use crate::symbolic::stats_probability::Exponential;
use crate::symbolic::stats_probability::Gamma;
use crate::symbolic::stats_probability::Normal;
use crate::symbolic::stats_probability::Poisson;
use crate::symbolic::stats_probability::StudentT;
use crate::symbolic::stats_probability::Uniform;

// --- Helper to convert a raw pointer to a boxed Expr ---
unsafe fn ptr_to_expr(
    ptr : *const Expr
) -> Option<Expr> {

    if ptr.is_null() {

        None
    } else {

        Some((*ptr).clone())
    }
}

// --- Generic Helper to wrap a Distribution in Expr ---
fn wrap_dist<
    D : Distribution + 'static,
>(
    dist : D
) -> *mut Expr {

    Box::into_raw(Box::new(
        Expr::Distribution(Arc::new(
            dist,
        )),
    ))
}

// --- Constructors for Distributions ---

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_normal(
    mean : *const Expr,
    std_dev : *const Expr,
) -> *mut Expr {

    let mean = ptr_to_expr(mean)
        .unwrap_or(Expr::Constant(0.0));

    let std_dev = ptr_to_expr(std_dev)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Normal {
        mean,
        std_dev,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_uniform(
    min : *const Expr,
    max : *const Expr,
) -> *mut Expr {

    let min = ptr_to_expr(min)
        .unwrap_or(Expr::Constant(0.0));

    let max = ptr_to_expr(max)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Uniform {
        min,
        max,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_binomial(
    n : *const Expr,
    p : *const Expr,
) -> *mut Expr {

    let n = ptr_to_expr(n)
        .unwrap_or(Expr::Constant(1.0));

    let p = ptr_to_expr(p)
        .unwrap_or(Expr::Constant(0.5));

    wrap_dist(Binomial {
        n,
        p,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_poisson(
    rate : *const Expr
) -> *mut Expr {

    let rate = ptr_to_expr(rate)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Poisson {
        rate,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_bernoulli(
    p : *const Expr
) -> *mut Expr {

    let p = ptr_to_expr(p)
        .unwrap_or(Expr::Constant(0.5));

    wrap_dist(Bernoulli {
        p,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_exponential(
    rate : *const Expr
) -> *mut Expr {

    let rate = ptr_to_expr(rate)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Exponential {
        rate,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_gamma(
    shape : *const Expr,
    rate : *const Expr,
) -> *mut Expr {

    let shape = ptr_to_expr(shape)
        .unwrap_or(Expr::Constant(1.0));

    let rate = ptr_to_expr(rate)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Gamma {
        shape,
        rate,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_beta(
    alpha : *const Expr,
    beta : *const Expr,
) -> *mut Expr {

    let alpha = ptr_to_expr(alpha)
        .unwrap_or(Expr::Constant(1.0));

    let beta = ptr_to_expr(beta)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Beta {
        alpha,
        beta,
    })
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_student_t(
    nu : *const Expr
) -> *mut Expr {

    let nu = ptr_to_expr(nu)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(StudentT {
        nu,
    })
}

// --- Methods on Distributions ---

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_pdf(
    dist : *const Expr,
    x : *const Expr,
) -> *mut Expr {

    let dist_expr = ptr_to_expr(dist);

    let x_expr = ptr_to_expr(x)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        Box::into_raw(Box::new(
            d.pdf(&x_expr),
        ))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_cdf(
    dist : *const Expr,
    x : *const Expr,
) -> *mut Expr {

    let dist_expr = ptr_to_expr(dist);

    let x_expr = ptr_to_expr(x)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        Box::into_raw(Box::new(
            d.cdf(&x_expr),
        ))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_expectation(
    dist : *const Expr
) -> *mut Expr {

    let dist_expr = ptr_to_expr(dist);

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        Box::into_raw(Box::new(
            d.expectation(),
        ))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_variance(
    dist : *const Expr
) -> *mut Expr {

    let dist_expr = ptr_to_expr(dist);

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        Box::into_raw(Box::new(
            d.variance(),
        ))
    } else {

        std::ptr::null_mut()
    }
}

#[no_mangle]

pub unsafe extern "C" fn rssn_dist_mgf(
    dist : *const Expr,
    t : *const Expr,
) -> *mut Expr {

    let dist_expr = ptr_to_expr(dist);

    let t_expr = ptr_to_expr(t)
        .unwrap_or(Expr::Constant(0.0));

    if let Some(Expr::Distribution(d)) =
        dist_expr
    {

        Box::into_raw(Box::new(
            d.mgf(&t_expr),
        ))
    } else {

        std::ptr::null_mut()
    }
}
