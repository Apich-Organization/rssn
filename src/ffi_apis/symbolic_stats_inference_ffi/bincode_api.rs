use crate::ffi_apis::common::{BincodeBuffer, from_bincode_buffer, to_bincode_buffer};
use crate::symbolic::core::Expr;
use crate::symbolic::stats_inference::{
    self,
};

/// Performs a one-sample t-test.

///

/// Takes bincode-serialized `Vec<Expr>` (data) and `Expr` (target mean).

/// Returns a bincode-serialized `HypothesisTest` representing the test result.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_one_sample_t_test(
    data_buf: BincodeBuffer,
    target_mean_buf: BincodeBuffer,
) -> BincodeBuffer {

    let data: Option<Vec<Expr>> =
        from_bincode_buffer(&data_buf);

    let target: Option<Expr> =
        from_bincode_buffer(
            &target_mean_buf,
        );

    match (data, target)
    { (Some(data), Some(target)) => {

        let result = stats_inference::one_sample_t_test_symbolic(&data, &target);

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Performs a two-sample t-test.

///

/// Takes bincode-serialized `Vec<Expr>` (two data sets) and `Expr` (hypothesized difference in means).

/// Returns a bincode-serialized `HypothesisTest` representing the test result.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_two_sample_t_test(
    data1_buf: BincodeBuffer,
    data2_buf: BincodeBuffer,
    mu_diff_buf: BincodeBuffer,
) -> BincodeBuffer {

    let data1: Option<Vec<Expr>> =
        from_bincode_buffer(&data1_buf);

    let data2: Option<Vec<Expr>> =
        from_bincode_buffer(&data2_buf);

    let diff: Option<Expr> =
        from_bincode_buffer(
            &mu_diff_buf,
        );

    match (data1, data2, diff)
    { (
        Some(d1),
        Some(d2),
        Some(diff),
    ) => {

        let result = stats_inference::two_sample_t_test_symbolic(&d1, &d2, &diff);

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}

/// Performs a z-test.

///

/// Takes bincode-serialized `Vec<Expr>` (data), `Expr` (target mean), and `Expr` (population standard deviation).

/// Returns a bincode-serialized `HypothesisTest` representing the test result.

#[unsafe(no_mangle)]

pub extern "C" fn rssn_bincode_z_test(
    data_buf: BincodeBuffer,
    target_mean_buf: BincodeBuffer,
    pop_std_dev_buf: BincodeBuffer,
) -> BincodeBuffer {

    let data: Option<Vec<Expr>> =
        from_bincode_buffer(&data_buf);

    let target: Option<Expr> =
        from_bincode_buffer(
            &target_mean_buf,
        );

    let sigma: Option<Expr> =
        from_bincode_buffer(
            &pop_std_dev_buf,
        );

    match (data, target, sigma)
    { (
        Some(data),
        Some(target),
        Some(sigma),
    ) => {

        let result = stats_inference::z_test_symbolic(
            &data,
            &target,
            &sigma,
        );

        to_bincode_buffer(&result)
    } _ => {

        BincodeBuffer::empty()
    }}
}
