use std::sync::Arc;

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
/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

unsafe fn ptr_to_expr(
    ptr: *const Expr
) -> Option<Expr> { unsafe {

    if ptr.is_null() {

        None
    } else {

        Some((*ptr).clone())
    }
}}

// --- Generic Helper to wrap a Distribution in Expr ---
fn wrap_dist<
    D: Distribution + 'static,
>(
    dist: D
) -> *mut Expr {

    Box::into_raw(Box::new(
        Expr::Distribution(Arc::new(
            dist,
        )),
    ))
}

// --- Constructors for Distributions ---

/// Creates a normal distribution.

///

/// Takes raw pointers to `Expr` representing the mean and standard deviation.

/// Returns a raw pointer to an `Expr` representing the normal distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_normal(
    mean: *const Expr,
    std_dev: *const Expr,
) -> *mut Expr { unsafe {

    let mean = ptr_to_expr(mean)
        .unwrap_or(Expr::Constant(0.0));

    let std_dev = ptr_to_expr(std_dev)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Normal {
        mean,
        std_dev,
    })
}}

/// Creates a uniform distribution.

///

/// Takes raw pointers to `Expr` representing the minimum and maximum values.

/// Returns a raw pointer to an `Expr` representing the uniform distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_uniform(
    min: *const Expr,
    max: *const Expr,
) -> *mut Expr { unsafe {

    let min = ptr_to_expr(min)
        .unwrap_or(Expr::Constant(0.0));

    let max = ptr_to_expr(max)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Uniform {
        min,
        max,
    })
}}

/// Creates a binomial distribution.

///

/// Takes raw pointers to `Expr` representing `n` (number of trials) and `p` (probability of success).

/// Returns a raw pointer to an `Expr` representing the binomial distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_binomial(
    n: *const Expr,
    p: *const Expr,
) -> *mut Expr { unsafe {

    let n = ptr_to_expr(n)
        .unwrap_or(Expr::Constant(1.0));

    let p = ptr_to_expr(p)
        .unwrap_or(Expr::Constant(0.5));

    wrap_dist(Binomial {
        n,
        p,
    })
}}

/// Creates a Poisson distribution.

///

/// Takes a raw pointer to an `Expr` representing the rate parameter (λ).

/// Returns a raw pointer to an `Expr` representing the Poisson distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_poisson(
    rate: *const Expr
) -> *mut Expr { unsafe {

    let rate = ptr_to_expr(rate)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Poisson {
        rate,
    })
}}

/// Creates a Bernoulli distribution.

///

/// Takes a raw pointer to an `Expr` representing `p` (probability of success).

/// Returns a raw pointer to an `Expr` representing the Bernoulli distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_bernoulli(
    p: *const Expr
) -> *mut Expr { unsafe {

    let p = ptr_to_expr(p)
        .unwrap_or(Expr::Constant(0.5));

    wrap_dist(Bernoulli {
        p,
    })
}}

/// Creates an exponential distribution.

///

/// Takes a raw pointer to an `Expr` representing the rate parameter (λ).

/// Returns a raw pointer to an `Expr` representing the exponential distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_exponential(
    rate: *const Expr
) -> *mut Expr { unsafe {

    let rate = ptr_to_expr(rate)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Exponential {
        rate,
    })
}}

/// Creates a gamma distribution.

///

/// Takes raw pointers to `Expr` representing the shape and rate parameters.

/// Returns a raw pointer to an `Expr` representing the gamma distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_gamma(
    shape: *const Expr,
    rate: *const Expr,
) -> *mut Expr { unsafe {

    let shape = ptr_to_expr(shape)
        .unwrap_or(Expr::Constant(1.0));

    let rate = ptr_to_expr(rate)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Gamma {
        shape,
        rate,
    })
}}

/// Creates a beta distribution.

///

/// Takes raw pointers to `Expr` representing the alpha and beta parameters.

/// Returns a raw pointer to an `Expr` representing the beta distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_beta(
    alpha: *const Expr,
    beta: *const Expr,
) -> *mut Expr { unsafe {

    let alpha = ptr_to_expr(alpha)
        .unwrap_or(Expr::Constant(1.0));

    let beta = ptr_to_expr(beta)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(Beta {
        alpha,
        beta,
    })
}}

/// Creates a Student's t-distribution.

///

/// Takes a raw pointer to an `Expr` representing the degrees of freedom (ν).

/// Returns a raw pointer to an `Expr` representing the Student's t-distribution.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_student_t(
    nu: *const Expr
) -> *mut Expr { unsafe {

    let nu = ptr_to_expr(nu)
        .unwrap_or(Expr::Constant(1.0));

    wrap_dist(StudentT {
        nu,
    })
}}

// --- Methods on Distributions ---

/// Computes the probability density function (PDF) of a distribution.

///

/// Takes raw pointers to `Expr` representing the distribution and the value `x`.

/// Returns a raw pointer to an `Expr` representing the PDF at `x`.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_pdf(
    dist: *const Expr,
    x: *const Expr,
) -> *mut Expr { unsafe {

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
}}

/// Computes the cumulative distribution function (CDF) of a distribution.

///

/// Takes raw pointers to `Expr` representing the distribution and the value `x`.

/// Returns a raw pointer to an `Expr` representing the CDF at `x`.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_cdf(
    dist: *const Expr,
    x: *const Expr,
) -> *mut Expr { unsafe {

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
}}

/// Computes the expectation (mean) of a distribution.

///

/// Takes a raw pointer to an `Expr` representing the distribution.

/// Returns a raw pointer to an `Expr` representing the expectation.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_expectation(
    dist: *const Expr
) -> *mut Expr { unsafe {

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
}}

/// Computes the variance of a distribution.

///

/// Takes a raw pointer to an `Expr` representing the distribution.

/// Returns a raw pointer to an `Expr` representing the variance.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_variance(
    dist: *const Expr
) -> *mut Expr { unsafe {

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
}}

/// Computes the moment generating function (MGF) of a distribution.

///

/// Takes raw pointers to `Expr` representing the distribution and the variable `t`.

/// Returns a raw pointer to an `Expr` representing the MGF.

#[unsafe(no_mangle)]

/// # Safety
///
/// This function is unsafe because it dereferences raw pointers as part of the FFI boundary.
/// The caller must ensure:
/// 1. All pointer arguments are valid and point to initialized memory.
/// 2. The memory layout of passed structures matches the expected C-ABI layout.
/// 3. Any pointers returned by this function are managed according to the API's ownership rules.

pub unsafe extern "C" fn rssn_dist_mgf(
    dist: *const Expr,
    t: *const Expr,
) -> *mut Expr { unsafe {

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
}}
