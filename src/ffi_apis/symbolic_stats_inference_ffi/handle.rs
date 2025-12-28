use std::os::raw::c_char;

use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_inference::HypothesisTest;
use crate::symbolic::stats_inference::{
    self,
};

unsafe fn collect_exprs(
    data: *const *const Expr,
    len: usize,
) -> Vec<Expr> {

    let mut exprs =
        Vec::with_capacity(len);

    for i in 0 .. len {

        let ptr = *data.add(i);

        if !ptr.is_null() {

            exprs.push((*ptr).clone());
        }
    }

    exprs
}

// Convert HypothesisTest to a boxed Expr (representing a struct/map)
// Since HypothesisTest is a struct, we maybe should return it as a serialized string or abstract handle?
// Handle API typically returns *mut Expr. We can wrap HypothesisTest in Expr::Custom or just return test_statistic?
// The user likely wants the full test result.
// A better way is to return *mut HypothesisTest if we want to expose it as an opaque type, or return a standardized Expr representation (e.g. Map).
// However, Expr doesn't have a Map variant.
// Let's assume we return the `test_statistic` for now, OR we add a FFI function to get specific fields.
// Actually, let's just return the test statistic as the primary result, or create separate functions for p-value.

// BUT, to be comprehensive, let's assume the user wants the full object.
// We can serialize it to JSON and return string in JSON API.
// In Handle API, we might return an opaque handle to HypothesisTest, but `Expr` is the standard handle type.
// Let's assume we return `test_statistic` as the Expr handle for simple usage.
// AND provide a specific function to get p-value expr.

// Better yet, let's implement `rssn_one_sample_t_test` returning `*mut HypothesisTest` (opaque).
// But existing pattern uses `*mut Expr`.
// Let's define an opaque wrapper or just return `test_statistic` Expr for now to keep simple consistency with "symbolic math" theme,
// AND provide `rssn_one_sample_t_test_details` to get p-value etc.

// Wait, the prompt implies "implement FFI for stats_inference".
// Let's create proper opaque handles for `HypothesisTest` ?
// Usually, we stick to `Expr` if possible.
// Let's wrap HypothesisTest in `Expr::Tuple` or similar if we strictly return `Expr`.
// Tuple(stat, p_value, df) seems reasonable?

/// Performs a one-sample t-test.

///

/// Takes a raw pointer to an array of `Expr` (data), its length,

/// and a raw pointer to an `Expr` (target mean).

/// Returns a raw pointer to an `Expr` tuple containing the test statistic,

/// p-value formula, and degrees of freedom.

#[no_mangle]

pub unsafe extern "C" fn rssn_one_sample_t_test(
    data: *const *const Expr,
    len: usize,
    target_mean: *const Expr,
) -> *mut Expr {

    if data.is_null() {

        return std::ptr::null_mut();
    }

    let sample =
        collect_exprs(data, len);

    let target =
        if target_mean.is_null() {

            Expr::Constant(0.0)
        } else {

            (*target_mean).clone()
        };

    let result = stats_inference::one_sample_t_test_symbolic(&sample, &target);

    // Return Tuple(statistic, p_value_formula, df)
    let df = result
        .degrees_of_freedom
        .unwrap_or(Expr::Constant(0.0)); // 0 if None
    Box::into_raw(Box::new(
        Expr::Tuple(vec![
            result.test_statistic,
            result.p_value_formula,
            df,
        ]),
    ))
}

/// Performs a two-sample t-test.

///

/// Takes raw pointers to two arrays of `Expr` (data sets), their lengths,

/// and a raw pointer to an `Expr` (hypothesized difference in means).

/// Returns a raw pointer to an `Expr` tuple containing the test statistic,

/// p-value formula, and degrees of freedom.

#[no_mangle]

pub unsafe extern "C" fn rssn_two_sample_t_test(
    data1: *const *const Expr,
    len1: usize,
    data2: *const *const Expr,
    len2: usize,
    mu_diff: *const Expr,
) -> *mut Expr {

    if data1.is_null()
        || data2.is_null()
    {

        return std::ptr::null_mut();
    }

    let sample1 =
        collect_exprs(data1, len1);

    let sample2 =
        collect_exprs(data2, len2);

    let diff = if mu_diff.is_null() {

        Expr::Constant(0.0)
    } else {

        (*mu_diff).clone()
    };

    let result = stats_inference::two_sample_t_test_symbolic(
        &sample1,
        &sample2,
        &diff,
    );

    let df = result
        .degrees_of_freedom
        .unwrap_or(Expr::Constant(0.0));

    Box::into_raw(Box::new(
        Expr::Tuple(vec![
            result.test_statistic,
            result.p_value_formula,
            df,
        ]),
    ))
}

/// Performs a z-test.

///

/// Takes a raw pointer to an array of `Expr` (data), its length,

/// a raw pointer to an `Expr` (target mean), and a raw pointer to an `Expr` (population standard deviation).

/// Returns a raw pointer to an `Expr` tuple containing the test statistic,

/// p-value formula, and a placeholder for degrees of freedom.

#[no_mangle]

pub unsafe extern "C" fn rssn_z_test(
    data: *const *const Expr,
    len: usize,
    target_mean: *const Expr,
    pop_std_dev: *const Expr,
) -> *mut Expr {

    if data.is_null()
        || pop_std_dev.is_null()
    {

        return std::ptr::null_mut();
    }

    let sample =
        collect_exprs(data, len);

    let target =
        if target_mean.is_null() {

            Expr::Constant(0.0)
        } else {

            (*target_mean).clone()
        };

    let sigma = (*pop_std_dev).clone();

    let result = stats_inference::z_test_symbolic(
        &sample,
        &target,
        &sigma,
    );

    // Z-test has no DF, so return Tuple(stat, p_value, null)
    Box::into_raw(Box::new(
        Expr::Tuple(vec![
            result.test_statistic,
            result.p_value_formula,
            Expr::NoSolution, /* Placeholder for None */
        ]),
    ))
}
