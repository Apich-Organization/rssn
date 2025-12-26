use crate::ffi_apis::common::*;
use crate::symbolic::core::Expr;
use crate::symbolic::stats_inference::HypothesisTest;
use crate::symbolic::stats_inference::{
    self,
};

#[no_mangle]

pub extern "C" fn rssn_bincode_one_sample_t_test(
    data_buf : BincodeBuffer,
    target_mean_buf : BincodeBuffer,
) -> BincodeBuffer {

    let data : Option<Vec<Expr>> = from_bincode_buffer(&data_buf);

    let target : Option<Expr> = from_bincode_buffer(&target_mean_buf);

    if let (Some(data), Some(target)) = (data, target) {

        let result = stats_inference::one_sample_t_test_symbolic(&data, &target);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_two_sample_t_test(
    data1_buf : BincodeBuffer,
    data2_buf : BincodeBuffer,
    mu_diff_buf : BincodeBuffer,
) -> BincodeBuffer {

    let data1 : Option<Vec<Expr>> = from_bincode_buffer(&data1_buf);

    let data2 : Option<Vec<Expr>> = from_bincode_buffer(&data2_buf);

    let diff : Option<Expr> = from_bincode_buffer(&mu_diff_buf);

    if let (Some(d1), Some(d2), Some(diff)) = (data1, data2, diff) {

        let result = stats_inference::two_sample_t_test_symbolic(&d1, &d2, &diff);

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}

#[no_mangle]

pub extern "C" fn rssn_bincode_z_test(
    data_buf : BincodeBuffer,
    target_mean_buf : BincodeBuffer,
    pop_std_dev_buf : BincodeBuffer,
) -> BincodeBuffer {

    let data : Option<Vec<Expr>> = from_bincode_buffer(&data_buf);

    let target : Option<Expr> = from_bincode_buffer(&target_mean_buf);

    let sigma : Option<Expr> = from_bincode_buffer(&pop_std_dev_buf);

    if let (Some(data), Some(target), Some(sigma)) = (data, target, sigma) {

        let result = stats_inference::z_test_symbolic(
            &data,
            &target,
            &sigma,
        );

        to_bincode_buffer(&result)
    } else {

        BincodeBuffer::empty()
    }
}
